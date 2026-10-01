"""
Step 1: subsetter--input.pl in Python.

Reads the IMGT/V-QUEST table and writes four files next to it:

  <in>.orig.board         label + full junction (with C and W)
  <in>.ign2.board         label + junction without its first and last amino acid
  <in>.ign2.dat           the same, in the input format of teiresias_char
  <in>.ign2.board4nrgrep  one line per sequence: label@sequence

The label of each sequence is:
  >seq id <TAB> annotation <TAB> V <TAB> D <TAB> J <TAB> V identity % <TAB> frame <TAB> CDR3 length

Usage:  python subsetter_input.py <V-QUEST table> [annotation file]
"""

import os
import re
import sys

from perl_compat import WORD, at, perl_num, perl_split, read_lines, s, uc, utf8_output, write_open

AA2IGNORE = 2  # amino acids trimmed in total (1 from each end), fixed as in the Perl


def read_annotation(path):
    """Annotation file: lines "seq id <TAB> value". Lines without exactly
    one TAB are ignored, as in the Perl."""
    annotation = {}
    if not path:
        return annotation
    if not os.path.exists(path):
        print(f"(!) annotation file {path} not found, continuing without it", file=sys.stderr)
        return annotation
    for line in read_lines(path):
        if WORD.search(line) and line.count("\t") == 1:
            fields = perl_split("\t", line)
            annotation[s(at(fields, 0))] = at(fields, 1)
    return annotation


def write_boards(in_path, annotation):
    """Writes <in>.ign2.board and <in>.orig.board."""
    cut = AA2IGNORE // 2
    header = []
    annotated = {}
    with write_open(f"{in_path}.ign{AA2IGNORE}.board") as cdr3aaseq_out, \
            write_open(f"{in_path}.orig.board") as orig_out:
        for line in read_lines(in_path):
            if re.match(r"sequence\.seq id", line):
                header = perl_split("\t", line)
            elif WORD.search(line):
                data = perl_split("\t", line)
                for i, name in enumerate(header):
                    annotated[name] = at(data, i)

                # the higher of the two V-region identities
                id_plain = annotated.get("sequence.V-REGION.id %")
                id_indel = annotated.get("sequence.V-REGION.id % (with ins/del events)")
                idpc = id_plain if perl_num(id_plain) >= perl_num(id_indel) else id_indel

                seq_id = s(annotated.get("sequence.seq id"))
                if annotation.get(seq_id) is None:
                    annotation[seq_id] = "UNK"

                label = ">" + "\t".join([
                    seq_id,
                    annotation[seq_id],
                    s(annotated.get("sequence.V-GENE and allele")),
                    s(annotated.get("sequence.D-GENE and allele")),
                    s(annotated.get("sequence.J-GENE and allele")),
                    s(idpc),
                    s(annotated.get("sequence.JUNCTION.frame")),
                    s(annotated.get("sequence.CDR3-IMGT.len")),
                ])

                cdr3aaseq = s(annotated.get("sequence.JUNCTION.aa seq")).replace(" ", "")
                if WORD.search(cdr3aaseq):
                    orig_out.write(f"{label}\n{cdr3aaseq}\n")

                subcdr3aaseq = cdr3aaseq[cut:-cut]
                if WORD.search(subcdr3aaseq):
                    cdr3aaseq_out.write(f"{label}\n{subcdr3aaseq}\n")


def write_dat(in_path):
    """<in>.ign2.board -> <in>.ign2.dat (input of teiresias_char).

    In the labels every character other than a letter or digit becomes _,
    and the sequence number is appended at the end.
    """
    board = f"{in_path}.ign{AA2IGNORE}.board"
    index = 0
    with write_open(f"{in_path}.ign{AA2IGNORE}.dat") as dat:
        for line in read_lines(board):
            if line.startswith(">"):
                line = re.sub(r"^>+", "", line)
                line = re.sub(r" +", "_", line)
                line = re.sub(r"[^A-Za-z0-9]", "_", line)
                if index != 0:
                    dat.write("\n")
                dat.write(f">{line} {index}\n")
                index += 1
            else:
                dat.write(line)
        dat.write("\n")


def write_board4nrgrep(in_path):
    """<in>.ign2.board -> <in>.ign2.board4nrgrep, one line "label@sequence".

    Only the 20 amino acids, X and the comma are kept in the sequence
    (the comma because it is inside the character class of the Perl).
    """
    board = f"{in_path}.ign{AA2IGNORE}.board"
    check = 0
    with write_open(f"{in_path}.ign{AA2IGNORE}.board4nrgrep") as out:
        for line in read_lines(board):
            if line.startswith(">"):
                line = re.sub(r"[\r\f]", "", line)
                if check != 0:
                    out.write("\n")
                out.write(line + "@")
                check += 1
            else:
                line = uc(line)
                line = re.sub(r"[\r\f]", "", line)
                line = re.sub(r"[^A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X,x]", "", line)
                out.write(line)
        out.write("\n")


def run(in_path, annotation_path=None):
    annotation = read_annotation(annotation_path)
    write_boards(in_path, annotation)
    write_dat(in_path)
    write_board4nrgrep(in_path)


if __name__ == "__main__":
    utf8_output()
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
