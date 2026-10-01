"""
subsetter - antigen receptor sequence clustering (ARResT/Subsetter)
Agathangelidis and Darzentas, http://bat.infspire.org
Python port of subsetter.pl and the subsetter--*.pl scripts.

Usage:
  python subsetter.py --receptor IG --in <V-QUEST table> [options]

  --receptor     IG (Ig / BcR / B-cell) or TR (TcR / T-cell)
  --in           the IMGT/V-QUEST output table
  --ide          minimum identity for a link (default 0.5)
  --sim          minimum similarity for a link (default 0.7)
  --ignorephylo  ignore the IGHV phylogenetic groups
  --annotation   file with "seq id <TAB> annotation" lines
  --out_table    where to write the final table (default <in>.ss.table)
  --teiresias    another teiresias_char executable (default: teiresias_char.exe in this folder)

Intermediate files are written next to the input file, with the same names
the Perl used.
"""

import argparse
import os
import shutil
import sys

import subsetter_clustering
import subsetter_input
import subsetter_nrgrep
import subsetter_pairs
import subsetter_teiresias
from perl_compat import utf8_output

HERE = os.path.dirname(os.path.abspath(__file__))

# The 13 fixed parameters of subsetter.pl ($args[0] to $args[12])
AA2IGNORE = 2                 # arg 0: amino acids to ignore (1 from each end)
TEIRESIAS_L = 3               # arg 1: l
TEIRESIAS_W = 6               # arg 2: w
TEIRESIAS_K = 2               # arg 3: k
TEIRESIAS_N = 1               # arg 4: n
EQUIVALENCES = os.path.join(HERE, "imgt11.aaequiv")         # arg 5
DEFAULT_IDE = "0.5"           # arg 6: identity threshold
DEFAULT_SIM = "0.7"           # arg 7: similarity threshold
LENGTH_DIFF = 0               # arg 8: allowed CDR3 length difference
OFFSET_DIFF = 0               # arg 9: allowed difference in pattern position
SCORE_THRESHOLD = "0"         # arg 10
CLUSTER_LINK_RATIO = "0"      # arg 11
PHYLO_EQUIVALENCES = os.path.join(HERE, "subsetter.phyloequivs")  # arg 12


def run(in_path, receptor, ide=None, sim=None, ignorephylo=False, annotation=None,
        out_table=None, teiresias=None, log=print):
    """Runs all steps. Returns the path of the final table, or None if no
    patterns were found."""
    ide = DEFAULT_IDE if ide is None else ide
    sim = DEFAULT_SIM if sim is None else sim
    phylo = None if (receptor == "TR" or ignorephylo) else PHYLO_EQUIVALENCES
    out_table = out_table or f"{in_path}.ss.table"
    ign = f"{in_path}.ign{AA2IGNORE}"

    log("1/5 reading the V-QUEST table")
    subsetter_input.run(in_path, annotation or None)

    log("2/5 TEIRESIAS: pattern discovery")
    subsetter_teiresias.run(f"{ign}.dat", f"{ign}.thr", EQUIVALENCES, teiresias,
                            l=TEIRESIAS_L, w=TEIRESIAS_W, k=TEIRESIAS_K, n=TEIRESIAS_N)

    log("3/5 searching the patterns in the sequences")
    flt_patterns = subsetter_nrgrep.run(f"{ign}.thr", f"{ign}.board4nrgrep", ide, sim)
    if flt_patterns == 0:
        log("(!) no patterns to look at... maybe try more relaxed parameters - exiting")
        return None

    log("4/5 sequence pairs")
    layout = subsetter_pairs.run(f"{ign}.thr.hits", f"{in_path}.orig.board", ide, sim,
                                 AA2IGNORE, LENGTH_DIFF, OFFSET_DIFF, phylo)

    log("5/5 clustering")
    table = subsetter_clustering.run(layout, f"{in_path}.orig.board",
                                     SCORE_THRESHOLD, CLUSTER_LINK_RATIO)
    shutil.copyfile(table, out_table)
    log(f"done: {out_table}")
    return out_table


def main(argv=None):
    utf8_output()
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--receptor", required=True, help="IG or TR")
    parser.add_argument("--in", dest="in_path", required=True, help="V-QUEST table")
    parser.add_argument("--ide", help="minimum identity (e.g. 0.5)")
    parser.add_argument("--sim", help="minimum similarity (e.g. 0.7)")
    parser.add_argument("--ignorephylo", action="store_true", help="ignore the phylogenetic groups")
    parser.add_argument("--annotation", help="annotation file")
    parser.add_argument("--out_table", help="output file")
    parser.add_argument("--teiresias", help="path to teiresias_char")
    parser.add_argument("--pass", dest="passphrase", help=argparse.SUPPRESS)  # ignored, as in the Perl
    args = parser.parse_args(argv)

    try:
        run(args.in_path, args.receptor, args.ide, args.sim, args.ignorephylo,
            args.annotation, args.out_table, args.teiresias)
    except Exception as error:  # noqa: BLE001 - Perl-style message and a failing exit code
        print(f"\n(!) error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
