"""
Test: the Python port must write the same files as the Perl tool.

The reference files in test/expected/ are the ones shipped with the Perl
tool (test.fasta.e.vq.table.*).

Test A: the whole pipeline, with teiresias_char.exe.
  On this test data the Windows .exe finds the same patterns as the Linux
  program that produced the reference files, but in a slightly different
  order. So .thr, .hits and .pairs.layout are compared ignoring order. All
  other files, including the final table, must be identical byte for byte.

Test B: steps 3-5 starting from the reference .thr.
  ALL files must be identical byte for byte.

Test C: the whole pipeline with teiresias_char.unix through WSL (only if WSL exists).
  ALL files must be identical byte for byte.

Usage:  python test/run_test.py
Results are written to test/output/.
"""

import os
import re
import shutil
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(TEST_DIR))

import subsetter  # noqa: E402
import subsetter_clustering  # noqa: E402
import subsetter_input  # noqa: E402
import subsetter_nrgrep  # noqa: E402
import subsetter_pairs  # noqa: E402
from perl_compat import utf8_output  # noqa: E402

INPUT = "test.fasta.e.vq.table"
EXPECTED = os.path.join(TEST_DIR, "expected")
ORDER_FREE = {".ign2.thr", ".ign2.thr.hits", ".ign2.thr.hits.ide0.5sim0.7.pairs.layout"}


def read_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


def without_order(suffix, data):
    """The content of a file in a form where order does not matter."""
    text = data.decode("latin-1").replace("\r\n", "\n")
    if suffix == ".ign2.thr.hits":
        blocks = re.split(r"(?m)^(?=#)", text)
        return sorted(b for b in blocks if b)
    if suffix.endswith(".layout"):
        text = re.sub(r" \[\d+\]", "", text)   # the pattern number depends on the order
    return sorted(text.splitlines())


def compare(out_dir, order_free):
    failures = 0
    for name in sorted(os.listdir(EXPECTED)):
        suffix = name[len(INPUT):]
        produced = os.path.join(out_dir, name)
        if not os.path.exists(produced):
            print(f"  MISSING    {name}")
            failures += 1
            continue
        expected, actual = read_bytes(os.path.join(EXPECTED, name)), read_bytes(produced)
        if expected == actual:
            print(f"  same       {name}")
        elif suffix in order_free and without_order(suffix, expected) == without_order(suffix, actual):
            print(f"  same*      {name}   (* same content, different pattern order)")
        else:
            print(f"  DIFFERENT  {name}")
            failures += 1
    return failures


def fresh_dir(name):
    path = os.path.join(TEST_DIR, "output", name)
    shutil.rmtree(path, ignore_errors=True)
    os.makedirs(path)
    shutil.copyfile(os.path.join(TEST_DIR, INPUT), os.path.join(path, INPUT))
    return path


def test_full_pipeline():
    print("Test A: whole pipeline with teiresias_char.exe")
    out = fresh_dir("full")
    subsetter.run(os.path.join(out, INPUT), "IG", log=lambda _msg: None)
    return compare(out, ORDER_FREE)


def test_from_reference_thr():
    print("Test B: steps 3-5 from the reference .thr")
    out = fresh_dir("from_thr")
    infile = os.path.join(out, INPUT)
    subsetter_input.run(infile)
    shutil.copyfile(os.path.join(EXPECTED, INPUT + ".ign2.thr"), infile + ".ign2.thr")
    subsetter_nrgrep.run(infile + ".ign2.thr", infile + ".ign2.board4nrgrep", "0.5", "0.7")
    layout = subsetter_pairs.run(infile + ".ign2.thr.hits", infile + ".orig.board", "0.5", "0.7",
                                 2, 0, 0, subsetter.PHYLO_EQUIVALENCES)
    table = subsetter_clustering.run(layout, infile + ".orig.board", "0", "0")
    shutil.copyfile(table, infile + ".ss.table")
    return compare(out, set())


def test_linux_teiresias():
    print("Test C: whole pipeline with teiresias_char.unix through WSL")
    unix = os.path.join(os.path.dirname(TEST_DIR), "teiresias_char.unix")
    if os.name != "nt" or not shutil.which("wsl") or not os.path.exists(unix):
        print("  skipped (no WSL or no teiresias_char.unix)")
        return 0
    out = fresh_dir("linux_teiresias")
    subsetter.run(os.path.join(out, INPUT), "IG", teiresias=unix, log=lambda _msg: None)
    return compare(out, set())


if __name__ == "__main__":
    utf8_output()
    failures = test_full_pipeline()
    print()
    failures += test_from_reference_thr()
    print()
    failures += test_linux_teiresias()
    print()
    print("OK: everything matches the Perl" if failures == 0 else f"FAILED: {failures} files differ")
    sys.exit(1 if failures else 0)
