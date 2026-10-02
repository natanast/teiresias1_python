"""
Step 4: subsetter--pairs.pl in Python.

For each pattern in <hits>:
  1. Finds exactly where the pattern falls in each CDR3 and highlights it
     in upper case (cARGPNISGGIAaLIKGFRAGDWYFDLw).
  2. Builds groups: each sequence is joined with those that have the same
     CDR3 length, the pattern at the same position, and enough coverage
     (ide, sim).
  3. Gives each group a 12-digit composite score: 6 digits for identity and
     6 for similarity, e.g. 96.1538 % and 100 % -> 961538999999
     (100 % is written as 999999).

Finally, for each pair of sequences found together, keeps the highest score
and the patterns that gave it, and writes them to
<hits>.ide<ide>sim<sim>.pairs.layout:

  "1099  cARG...DLw  999999999999  IGHV1-2*04 F ...  UNK"  "2099  ..."  999999999999  ARGP...DL [0]

Perl quirks kept on purpose (for identical results):
  * The last pattern of the <hits> file is never processed, because the
    Perl only processes a pattern when it reaches the next "#".
  * If a pattern is not found exactly in the CDR3, the sequence keeps the
    position of the previous sequence (variable cdr3offset).
  * Links (patternblast) are shared across all patterns.

Usage:  python subsetter_pairs.py <hits> <orig.board> <ide> <sim> <aa2ignore>
                                  <length diff> <offset diff> [phyloequivs]
"""

import re
import sys

from perl_compat import (at, lc, num_asc, num_desc, perl_num, perl_split, perl_split_ws,
                         read_lines, s, uc, utf8_output, write_open)
from subsetter_nrgrep import pattern_counts


def read_board(board_path):
    """The orig.board: label, genes and CDR3 (lower case) of each sequence."""
    mem = {}
    sample = ""
    for line in read_lines(board_path):
        if line.startswith(">"):
            label = line[1:]
            labelarray = perl_split("\t", label)
            sample = re.sub(r" +", "", s(at(labelarray, 0)))
            entry = mem.setdefault(sample, {})
            entry["label"] = label
            entry["vdj"] = f"{s(at(labelarray, 2))} {s(at(labelarray, 3))} {s(at(labelarray, 4))}"
            entry["vsubgroup"] = at(perl_split("-", at(labelarray, 2)), 0)   # IGHV1
            entry["vgene"] = at(perl_split(r"\*", at(labelarray, 2)), 0)     # IGHV1-2
        else:
            entry = mem.setdefault(sample, {})
            entry["cdr3aaseq"] = lc(line)
            entry["cdr3len"] = len(entry["cdr3aaseq"])
    return mem


def read_phylo(phylo_path):
    """subsetter.phyloequivs: each line is one phylogenetic group.

    Returns {name: group number} and the level (vsubgroup, vgene or
    vallele), which, as in the Perl, is set by the last word of the file.
    """
    phyloequivs = {}
    phylolevel = None
    phylogroup = 1
    for line in read_lines(phylo_path):
        for bit in perl_split(r"[ \t\n\r\f\v]+", line):
            phyloequivs[bit] = phylogroup
            if "-" in bit:
                phylolevel = "vgene"
            elif "*" in bit:
                phylolevel = "vallele"
            else:
                phylolevel = "vsubgroup"
        phylogroup += 1
    return phyloequivs, phylolevel


def length_key(value):
    """Perl hash keys are strings (undef -> "")."""
    return "" if value is None else str(value)


def run(hits_path, board_path, ide, sim, aa2ignore=2, lendiff=0, offsetdiff=0, phylo_path=None):
    ide_n, sim_n = perl_num(ide), perl_num(sim)
    aa2ignore, lendiff, offsetdiff = int(aa2ignore), int(lendiff), int(offsetdiff)
    half = aa2ignore / 2

    mem = read_board(board_path)

    def M(sample, key):
        return mem.get(sample, {}).get(key)

    phyloflag = phylo_path is not None
    phyloequivs, phylolevel = read_phylo(phylo_path) if phyloflag else ({}, None)

    # State kept from pattern to pattern (global variables in the Perl)
    patternblast = {}       # sample -> set of samples it was linked with
    samples_only_blo = {}   # score -> pattern -> group -> [entries]
    patternindices = {}     # pattern -> sequence number
    annotation = {}
    patternindex = 0
    cdr3offset = None       # NOT reset between sequences, as in the Perl

    # State of the current pattern
    flag = False
    patternflag = False
    count = 0
    pattern = None
    per_pattern = {}        # CDR3 length -> sample -> line
    cdr3offsetmem = {}
    highlighted_cdr3 = {}
    identical = similar = patternlen = 0
    patternarray = []
    lc_regex = None

    def flush_pattern():
        """Grouping for the pattern that just ended (lines 70-148 of the Perl)."""
        nonlocal patternindex, patternflag
        finalcount, toprint, grouplensum, samples_for_blo = {}, {}, {}, {}
        printed = set()
        group = 0
        for length in num_asc(per_pattern):
            for sample in sorted(per_pattern[length]):
                finalcount[group] = finalcount.get(group, 0) + 1
                grouplensum[group] = grouplensum.get(group, 0) + perl_num(M(sample, "cdr3len")) - aa2ignore
                toadd2print = {per_pattern[length_key(M(sample, "cdr3len"))][sample]}
                samples_for_blo.setdefault(group, set()).add(sample)
                basesample = sample
                base_len = perl_num(M(basesample, "cdr3len"))
                for i in range(base_len, base_len + lendiff + 1):
                    for other in sorted(per_pattern.get(str(i), {})):
                        if other == basesample:
                            continue
                        if abs(perl_num(cdr3offsetmem.get(basesample))
                               - perl_num(cdr3offsetmem.get(other))) > offsetdiff:
                            continue
                        base_trim = base_len - aa2ignore
                        other_trim = perl_num(M(other, "cdr3len")) - aa2ignore
                        if (similar / base_trim >= sim_n and similar / other_trim >= sim_n
                                and identical / base_trim >= ide_n and identical / other_trim >= ide_n):
                            patternblast.setdefault(basesample, set()).add(other)
                            finalcount[group] += 1
                            grouplensum[group] += other_trim
                            toadd2print.add(per_pattern[length_key(M(other, "cdr3len"))][other])
                            samples_for_blo[group].add(other)
                toprint[group] = "".join(sorted(toadd2print))
                group += 1

        for group in sorted(samples_for_blo):
            if finalcount[group] <= 1:
                continue
            # groups with exactly the same members count once
            if toprint[group] in printed:
                continue
            printed.add(toprint[group])
            cdr3lenave = grouplensum[group] / finalcount[group]
            ide2print = "%.4f" % (identical / cdr3lenave * 100)
            sim2print = "%.4f" % (similar / cdr3lenave * 100)
            ide4score = "999999" if float(ide2print) == 100 else ide2print
            sim4score = "999999" if float(sim2print) == 100 else sim2print
            patternscore = (ide4score + sim4score).replace(".", "")
            patternflag = True
            for sample in sorted(samples_for_blo[group]):
                entry = (f"{sample}  {s(highlighted_cdr3.get(sample))}  {patternscore}  "
                         f"{s(M(sample, 'vdj'))}  {s(annotation.get(sample))}")
                samples_only_blo.setdefault(patternscore, {}).setdefault(pattern, {}) \
                    .setdefault(group, []).append(entry)
        if patternflag:
            patternindices[pattern] = patternindex
            patternindex += 1

    for line in read_lines(hits_path):
        if line.startswith("#"):
            patternflag = False
            if flag and count > 1:
                flush_pattern()
            patternflag = False
            pattern = at(perl_split_ws(line), 1)
            lc_regex = re.compile(lc(pattern))
            pattern4len, patternlen, identical, similar = pattern_counts(pattern)
            patternarray = list(pattern4len)
            if identical / patternlen >= ide_n and similar / patternlen >= sim_n:
                patternflag = True
            flag = True
            count = 0
            per_pattern = {}
            cdr3offsetmem = {}
            highlighted_cdr3 = {}

        if line.startswith(">") and patternflag:
            sample = re.sub(r" +", "", line[1:])
            hit = M(sample, "label")
            cdr3 = s(M(sample, "cdr3aaseq"))
            cdr3array = list(cdr3)
            skip = False
            for m in lc_regex.finditer(cdr3):
                pos = m.end()
                cdr3offset = pos - patternlen + 1      # start position, counting from 1
                cdr3endcheck = len(cdr3) - pos + 1
                # the pattern must not touch the C (start) or the W/F (end)
                if not (cdr3offset > half and cdr3endcheck > half):
                    skip = True
                    break
                for i, ch in enumerate(patternarray):
                    if ("A" <= ch <= "Z" or "a" <= ch <= "z") and cdr3offset + i - 1 < len(cdr3array):
                        cdr3array[cdr3offset + i - 1] = uc(cdr3array[cdr3offset + i - 1])
            if skip:
                continue
            cdr3 = "".join(cdr3array)
            highlighted_cdr3[sample] = cdr3
            count += 1
            cdr3offsetmem[sample] = cdr3offset
            annotation[sample] = at(perl_split("\t", hit), 1)
            per_pattern.setdefault(length_key(M(sample, "cdr3len")), {})[sample] = \
                f"{patternindex}\t@\t%\t{pattern}\t{s(hit)}\t{cdr3}\n"
    # NOTE: the Perl does not flush the last pattern here.

    def phylo_group(sample):
        value = M(sample, phylolevel) if phylolevel is not None else None
        return phyloequivs.get(s(value))

    def linked(a, b):
        return b in patternblast.get(a, ()) or a in patternblast.get(b, ())

    # For each pair: the highest score and the patterns that gave it
    pairwise = {}   # entry a -> entry b -> [score, set of patterns]
    for patternscore in num_desc(samples_only_blo):
        score_n = perl_num(patternscore)
        for pattern in sorted(samples_only_blo[patternscore]):
            for group in sorted(samples_only_blo[patternscore][pattern]):
                entries = samples_only_blo[patternscore][pattern][group]
                for entry_a in entries:
                    sample_a = entry_a.split("  ")[0]
                    for entry_b in entries:
                        sample_b = entry_b.split("  ")[0]
                        if phyloflag:
                            group_a, group_b = phylo_group(sample_a), phylo_group(sample_b)
                            same_phylo = (group_a is not None and group_b is not None and group_a == group_b) \
                                or (group_a is None and group_b is None)
                            if not same_phylo:
                                continue
                        if not linked(sample_a, sample_b) or entry_a == entry_b:
                            continue
                        record = pairwise.setdefault(entry_a, {}).get(entry_b)
                        if record is None:
                            record = pairwise[entry_a][entry_b] = [-1, set()]
                        if score_n >= perl_num(record[0]):
                            record[0] = patternscore
                            record[1].add(pattern)

    layout_path = f"{hits_path}.ide{ide}sim{sim}.pairs.layout"
    written = set()
    with write_open(layout_path) as out:
        for entry_a in sorted(pairwise):
            for entry_b in sorted(pairwise[entry_a]):
                if (entry_b, entry_a) in written:
                    continue
                score, patterns = pairwise[entry_a][entry_b]
                max_patterns = ",".join(f"{p} [{patternindices.get(p, '')}]" for p in sorted(patterns))
                out.write(f'"{entry_a}"\t"{entry_b}"\t{score}\t{max_patterns}\n')
                written.add((entry_a, entry_b))
    return layout_path


if __name__ == "__main__":
    utf8_output()
    if len(sys.argv) < 8:
        sys.exit(__doc__)
    run(*sys.argv[1:8], phylo_path=sys.argv[8] if len(sys.argv) > 8 else None)
