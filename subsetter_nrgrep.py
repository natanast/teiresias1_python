"""
Step 3: subsetter--nrgrep.pl in Python, without the nrgrep program.

For each pattern in <thr> finds which sequences contain it and keeps only
those that the pattern covers enough:

  defined amino acids of the pattern / sequence length                 >= ide
  (defined + equivalence-class positions) / sequence length            >= sim

If at least 2 sequences remain, writes to <thr>.hits:

  #  ARGPNISGGIA.LI[KRH]GFRAGDWYFDL
  >1099
  >2099
  ...

Difference from the Perl: the Perl called "nrgrep -k 1i", which allows one
extra letter and, because of internal optimisations in nrgrep, sometimes
misses exact occurrences at the end of the sequence. Here the search is
exact: a sequence is found if and only if it contains the pattern. That is
also all the next step (pairs) can use. On the test data the results are
identical. See the README for details.

Usage:  python subsetter_nrgrep.py <thr> <board4nrgrep> <ide> <sim>
"""

import bisect
import re
import sys

from perl_compat import at, perl_num, perl_split, perl_split_ws, read_lines, utf8_output, write_open

TOKEN = re.compile(r">\S+", re.ASCII)


def pattern_counts(pattern):
    """Counts the positions of a TEIRESIAS pattern.

    "RG..[AVLI]" -> length 5, 2 identical (R, G), 3 similar (R, G, [AVLI]).
    Returns (pattern4len, patternlen, identical, similar), where every
    [..] class in pattern4len has been replaced by x.
    """
    pattern4len = re.sub(r"(\[[A-Z.]*\])", "x", pattern)
    patternlen = len(pattern4len)
    # Perl's tr/[A-Z]/i/: counts upper-case letters and brackets
    identical = sum(1 for ch in pattern4len if "A" <= ch <= "Z" or ch in "[]")
    equivalent = pattern4len.count("x")
    return pattern4len, patternlen, identical, identical + equivalent


def pattern_regex(pattern):
    """TEIRESIAS pattern -> Python regular expression.

    The syntax is the same (. any character, [..] one of these). Only
    characters that are special in Python but not in nrgrep are escaped.
    """
    out = []
    i = 0
    while i < len(pattern):
        ch = pattern[i]
        if ch == "[":
            end = pattern.find("]", i + 1)
            if end < 0:
                raise ValueError(f"pattern with an unclosed bracket: {pattern}")
            out.append("[" + "".join(re.escape(c) for c in pattern[i + 1:end]) + "]")
            i = end + 1
            continue
        out.append(ch if ch in ".*+?|()" or ch.isalnum() else re.escape(ch))
        i += 1
    return re.compile("".join(out))


class Board:
    """The board4nrgrep file, prepared for fast searching.

    Sequences are joined into one text separated by \\n, so that each
    pattern is searched with a single regular expression call.
    """

    def __init__(self, path):
        self.lengths = {}   # >id -> sequence length (the Perl %lengths)
        self.tokens = []    # the >... words of each record, as the Perl picked them up
        searchable = []
        for line in read_lines(path):
            if line.startswith(">"):
                seq = at(perl_split("@", line), 1)
                self.lengths[at(perl_split_ws(line), 0)] = None if seq is None else len(seq)
            # nrgrep searched "@.*pattern": the pattern must come after the first @
            at_sign = line.find("@")
            searchable.append(line[at_sign + 1:] if at_sign >= 0 else "")
            self.tokens.append(TOKEN.findall(line))

        # Records are grouped by sequence length, so that a pattern is only
        # searched in sequences short enough to pass the ide/sim filter.
        # This is faster and does not change the result.
        by_length = {}
        for record, tokens in enumerate(self.tokens):
            if tokens:  # a record without >... never produces output
                shortest = min(self.lengths.get(token) or 0 for token in tokens)
                by_length.setdefault(shortest, []).append(record)
        self.groups = []    # (length, text, start of each record in the text, records)
        for length in sorted(by_length):
            records = by_length[length]
            starts, offset = [], 0
            for record in records:
                starts.append(offset)
                offset += len(searchable[record]) + 1
            text = "\n".join(searchable[record] for record in records)
            self.groups.append((length, text, starts, records))

    def records_matching(self, regex, max_length=None):
        """The records (in file order) that contain the pattern, among those
        whose sequence is at most max_length long."""
        found = []
        for length, text, starts, records in self.groups:
            if max_length is not None and length > max_length:
                break
            last = -1
            for m in regex.finditer(text):
                i = bisect.bisect_right(starts, m.start()) - 1
                if i != last:
                    found.append(records[i])
                    last = i
        found.sort()
        return found


def run(thr_path, board_path, ide, sim, hits_path=None):
    """Writes <thr>.hits and returns the number of patterns written."""
    ide_n, sim_n = perl_num(ide), perl_num(sim)
    hits_path = hits_path or f"{thr_path}.hits"
    board = Board(board_path)
    written = 0
    with write_open(hits_path) as hits:
        for line in read_lines(thr_path):
            if line.startswith("#"):
                continue
            fields = perl_split_ws(line)
            if len(fields) < 3:
                continue  # empty line
            grep_seqlet = fields[2]
            _, patternlen, identical, similar = pattern_counts(grep_seqlet)
            if not (identical / patternlen >= ide_n and similar / patternlen >= sim_n):
                continue

            # Longest sequence that can pass the filter below, plus 1 as a
            # margin for rounding (the exact filter is applied below).
            bounds = ([similar / sim_n] if sim_n > 0 else []) + ([identical / ide_n] if ide_n > 0 else [])
            max_length = min(bounds) + 1 if bounds else None

            fltgrep = []
            for record in board.records_matching(pattern_regex(grep_seqlet), max_length):
                for token in board.tokens[record]:
                    length = board.lengths.get(token)
                    if not length:
                        raise ValueError(f"unknown sequence length for {token}")
                    if similar / length >= sim_n and identical / length >= ide_n:
                        fltgrep.append(token)
            if len(fltgrep) > 1:
                hits.write(f"#  {grep_seqlet}\n")
                hits.write("".join(f"{token}\n" for token in fltgrep))
                written += 1
    return written


if __name__ == "__main__":
    utf8_output()
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
