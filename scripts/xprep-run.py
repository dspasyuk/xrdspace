#!/usr/bin/env python3
# Copyright (c) 2026 Denis Spasyuk. MIT License.
#
# Minimal pexpect driver for Bruker XPREP used by tests/xrdspace-xprep.js.
#
# Runs XPREP on a SHELX .fcf file, navigates the interactive menu far enough to
# let XPREP determine the space group (accepting XPREP's own recommended option),
# quits, and prints a JSON summary on stdout:
#
#   {"ok": true, "sgNumber": 121, "hm": "I-42m", "chosen": "E",
#    "candidates": [{"opt": "A", "hm": "I-4", "number": 82}, ...]}
#
# Usage: xprep-run.py <file.fcf>

import json
import os
import re
import signal
import sys
import time

import pexpect

XPREP = os.environ.get("XPREP_BIN", "/home/denis/CODE/xdsgo/executables/xprep")

ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")
PROMPT_RE = re.compile(r"\[[^\]\n]*\]\s*:\s*$")


def strip_ansi(s):
    return ANSI_RE.sub("", s)


def parse_prp(path):
    """Parse an XPREP .prp log: return (chosen_letter, candidates list)."""
    if not os.path.exists(path):
        return None, []
    text = open(path, "r", errors="ignore").read()
    text = strip_ansi(text)
    chosen = None
    m = re.search(r"Option\s*\[([A-Za-z0-9])\]\s*chosen", text)
    if m:
        chosen = m.group(1)
    candidates = []
    seen = set()
    # Full candidate row, e.g.:
    #   [A] I-4            # 82  non-cen  1    99  0.012     24  0.0 / 53.4  31.06
    # columns: opt  space-group  #no  type  axes  csd  R(sym)  N(eq)  Syst.Abs  CFOM
    # The "Syst.Abs" column is the only "x / y" pair; it separates the leading
    # columns (R(sym), N(eq) are its two predecessors) from the trailing CFOM.
    for raw in text.splitlines():
        line = raw.strip()
        m = re.match(r"\[([A-Za-z0-9])\]\s+(\S+)\s+#\s*(\d+)\s+(\S+)", line)
        if not m:
            continue
        opt = m.group(1)
        if opt in seen:
            continue
        seen.add(opt)
        toks = line.split()
        def num(idx):
            try:
                return float(toks[idx])
            except (IndexError, ValueError):
                return None
        rsym = cfom = None
        for i in range(2, len(toks) - 3):
            if toks[i + 1] == "/":
                # toks[i-2]=R(sym)  toks[i-1]=N(eq)  toks[i]=x  '/'  toks[i+2]=y  toks[i+3]=CFOM
                rsym = num(i - 2)
                cfom = num(i + 3)
                break
        candidates.append({
            "opt": opt,
            "hm": m.group(2),
            "number": int(m.group(3)),
            "type": m.group(4),
            "rsym": rsym,
            "cfom": cfom,
        })
    return chosen, candidates


def decide(plain, state):
    """Given the recent (ANSI-stripped) output, return the reply string."""
    low = plain[-800:].lower()
    if os.environ.get("XPREP_DEBUG"):
        sys.stderr.write("FLAGS sg=%s no=%s type=%s axes=%s len=%d\n" % (
            "space group" in low, "no." in low, "type" in low, "axes" in low,
            len(plain)))
    if "process dataset" in low:
        return ""
    if "use cell from reflection data file" in low:
        return ""
    if "type <cr> to continue" in low:
        return ""

    # After the space group is chosen we are done; the caller breaks out and
    # parses the .prp. Remaining branches are only reached before that point.
    if state.get("sg_prompted"):
        return "Q"

    # Space-group candidate table: rows look like "[E] I-42m  # 121  non-cen".
    # The prompt that follows picks one of these options; accept XPREP's
    # recommended default. Some datasets offer a single candidate, so also
    # recognize the table by its header.
    tail = plain[-2500:]
    rows = re.findall(r"\[[A-Za-z0-9]\]\s+\S+\s+#\s*\d+", tail)
    header = tail.lower()
    if len(rows) >= 2 or ("space group" in header and "no." in header and
                          "type" in header and "axes" in header):
        state["sg_prompted"] = True
        state["accept_sg"] = True
        return ""

    # Lattice / metric-symmetry prompts: accept the default recommendation.
    if "lattice type" in low:
        return ""
    if re.search(r"select option \[[pabcifor]\]\s*:", low):
        return ""

    # Space-group submenu: choose "Determine SPACE GROUP" (S).
    if "input known space group" in low:
        return ""
    if "determine or input space group" in low:  # main menu
        return "S"
    if "determine space group" in low:
        return ""

    return ""


def interact(child, state):
    """Navigate the XPREP menus until the space group is chosen."""
    buf = ""
    seen = {}
    deadline = time.time() + 12
    while time.time() < deadline:
        try:
            chunk = child.read_nonblocking(size=16384, timeout=0.5)
            if chunk:
                buf += chunk
        except pexpect.TIMEOUT:
            if not child.isalive():
                break
            continue
        except pexpect.EOF:
            break

        plain = strip_ansi(buf)
        if PROMPT_RE.search(plain) or plain.rstrip().endswith("Type <CR> to continue"):
            reply = decide(plain, state)
            if os.environ.get("XPREP_DEBUG"):
                sys.stderr.write("REPLY %r :: %s\n" % (
                    reply, plain.replace("\n", " ")[-140:]))
            # Bail out if the menu interaction starts cycling (e.g. XPREP
            # loops on an F-superlattice warning without ever offering a
            # space-group table).
            sig = (reply + "|" + plain.replace("\n", " ")[-80:])
            seen[sig] = seen.get(sig, 0) + 1
            if seen[sig] > 8:
                break
            try:
                child.sendline(reply)
            except Exception:
                break
            buf = ""
            if state.pop("accept_sg", False):
                # XPREP writes the chosen space group to the .prp as it accepts
                # the option; give it a moment to flush, then stop.
                time.sleep(0.8)
                break


def _on_alarm(signum, frame):
    raise TimeoutError("xprep-run timed out")


def main():
    if len(sys.argv) < 2:
        print(json.dumps({"ok": False, "error": "usage: xprep-run.py <file.fcf>"}))
        return 2
    fcf = os.path.abspath(sys.argv[1])
    prp = os.path.splitext(fcf)[0] + ".prp"
    if os.path.exists(prp):
        try:
            os.remove(prp)
        except OSError:
            pass

    state = {"sg_prompted": False, "accept_sg": False}
    child = pexpect.spawn(XPREP, [fcf], encoding="utf-8", timeout=20,
                          dimensions=(200, 400))
    old = signal.signal(signal.SIGALRM, _on_alarm)
    signal.alarm(15)
    try:
        interact(child, state)
    except Exception:
        pass
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
        try:
            child.kill(signal.SIGKILL)
        except Exception:
            pass
        try:
            child.close(force=True)
        except Exception:
            pass

    chosen, candidates = parse_prp(prp)
    if chosen is None:
        result = {"ok": False, "error": "no chosen option in .prp",
                  "candidates": candidates}
    else:
        sel = next((c for c in candidates if c["opt"] == chosen), None)
        result = {"ok": sel is not None, "chosen": chosen,
                  "sgNumber": sel["number"] if sel else None,
                  "hm": sel["hm"] if sel else None,
                  "candidates": candidates}
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
