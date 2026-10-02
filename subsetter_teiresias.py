"""
Step 2: run teiresias_char (pattern discovery).

teiresias_char.exe is used as is. It cannot open files with long paths or
paths containing non-ASCII characters, so it runs inside a temporary folder
with short file names (in.dat, equiv.txt, out.thr). The result is then
copied to the requested location.

If teiresias_char.unix is given (the Linux program the Perl used), it runs
through WSL on Windows. The .exe and the .unix usually find the same
patterns, but not always (see the README).

Parameters as in subsetter.pl:
  -l3   at least 3 defined positions (not dots) in each elementary pattern
  -w6   within a span of at most 6 positions
  -c2   overlap when joining elementary patterns (= l - 1)
  -k2   the pattern must occur in at least 2 sequences
  -v    -k counts sequences, not occurrences
  -n1   at most one equivalence class, e.g. [KRH], per pattern
  -b    the amino acid equivalence file (imgt11.aaequiv)

Usage:  python subsetter_teiresias.py <in.dat> <out.thr> [equivalence file] [executable]
"""

import os
import shutil
import subprocess
import sys
import tempfile

from perl_compat import utf8_output

HERE = os.path.dirname(os.path.abspath(__file__))

L, W, K, N = 3, 6, 2, 1

# In WSL the program runs in a Linux folder: writing there is much faster
# than writing to the Windows disk.
WSL_SCRIPT = ('d=$(mktemp -d) && cp in.dat equiv.txt teiresias_char.unix "$d"/ && '
              'here=$(pwd) && cd "$d" && chmod +x teiresias_char.unix && '
              './teiresias_char.unix "$@"; rc=$?; '
              '[ -f out.thr ] && cp out.thr "$here"/; cd "$here"; rm -rf "$d"; exit $rc')


def default_executable():
    if os.name == "nt":
        return os.path.join(HERE, "teiresias_char.exe")
    return os.path.join(HERE, "teiresias_char.unix")


def run(dat_path, thr_path, equiv_path=None, executable=None, l=L, w=W, k=K, n=N):
    executable = os.path.abspath(executable or default_executable())
    if equiv_path is None:
        equiv_path = os.path.join(HERE, "imgt11.aaequiv")
    if not os.path.exists(executable):
        raise FileNotFoundError(f"{executable} not found")

    options = ["-iin.dat", "-oout.thr", f"-l{l}", f"-c{l - 1}", f"-k{k}",
               f"-w{w}", f"-n{n}", "-v", "-bequiv.txt"]
    with tempfile.TemporaryDirectory(prefix="teir") as work:
        shutil.copyfile(dat_path, os.path.join(work, "in.dat"))
        shutil.copyfile(equiv_path, os.path.join(work, "equiv.txt"))
        if os.name == "nt" and executable.endswith(".unix"):
            shutil.copyfile(executable, os.path.join(work, "teiresias_char.unix"))
            args = ["wsl", "-e", "sh", "-c", WSL_SCRIPT, "sh"] + options
        else:
            args = [executable] + options
        proc = subprocess.run(args, cwd=work, capture_output=True, text=True, errors="replace")
        out_thr = os.path.join(work, "out.thr")
        if proc.returncode != 0 or not os.path.exists(out_thr):
            raise RuntimeError(f"teiresias_char failed (exit code {proc.returncode})\n"
                               f"{proc.stdout}\n{proc.stderr}")
        shutil.copyfile(out_thr, thr_path)
    return proc.stdout


if __name__ == "__main__":
    utf8_output()
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    run(sys.argv[1], sys.argv[2],
        sys.argv[3] if len(sys.argv) > 3 else None,
        sys.argv[4] if len(sys.argv) > 4 else None)
