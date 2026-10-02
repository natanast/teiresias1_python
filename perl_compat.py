"""
Helpers that reproduce Perl behaviour.

The original scripts rely on Perl details (how split works, how a string
becomes a number, that chomp only removes \n). They are collected here so
that the Python port writes exactly the same files.
"""

import re

# Files are read and written as latin-1: one byte is one character. Bytes go
# out exactly as they came in, and strings sort in byte order, as in Perl.
ENCODING = "latin-1"

# Perl's \w (no locale) matches only ASCII letters, digits and _
WORD = re.compile(r"\w", re.ASCII)

# Perl's \s
_PERL_SPACE = " \t\n\r\f\v"

_UPPER = str.maketrans("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
_LOWER = str.maketrans("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz")


def utf8_output():
    """Print non-ASCII text (e.g. file paths) even when output is redirected
    (on Windows a redirected stream defaults to cp1252)."""
    import sys
    for stream in (sys.stdout, sys.stderr):
        encoding = (getattr(stream, "encoding", None) or "").lower().replace("-", "")
        if hasattr(stream, "reconfigure") and encoding != "utf8":
            stream.reconfigure(encoding="utf-8", errors="replace")


def read_lines(path):
    """Lines of a file without the final \\n (like chomp).

    A trailing \\r is kept, as with Perl on Linux.
    """
    with open(path, encoding=ENCODING, newline="\n") as fh:
        for line in fh:
            yield line[:-1] if line.endswith("\n") else line


def write_open(path):
    """Open a file for writing with \\n line endings (not \\r\\n), also on Windows."""
    return open(path, "w", encoding=ENCODING, newline="\n")


def uc(text):
    """Perl uc: changes ASCII letters only."""
    return text.translate(_UPPER)


def lc(text):
    """Perl lc: changes ASCII letters only."""
    return text.translate(_LOWER)


def s(value):
    """A Perl undef inside a string becomes the empty string."""
    return "" if value is None else value


def at(items, index):
    """$array[index]: None (undef) if the element does not exist."""
    return items[index] if 0 <= index < len(items) else None


def perl_split(pattern, text):
    """split /pattern/, text

    As in Perl, empty trailing fields are removed.
    """
    if text is None:
        return []
    parts = re.split(pattern, text)
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def perl_split_ws(text):
    """split " ", text: split into words, ignoring leading whitespace."""
    if text is None:
        return []
    return perl_split("[" + re.escape(_PERL_SPACE) + "]+", text.lstrip(_PERL_SPACE))


_NUMBER = re.compile(r"[ \t\n\r\f\v]*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)")
_INTEGER = re.compile(r"[+-]?\d+")


def perl_num(value):
    """The numeric value of a string, as Perl sees it.

    "97.13" -> 97.13, "12abc" -> 12, "" or undef -> 0.
    Integers stay integers (int) so that sums are exact.
    """
    if value is None:
        return 0
    if isinstance(value, (int, float)):
        return value
    m = _NUMBER.match(value)
    if not m:
        return 0
    number = m.group(1)
    if _INTEGER.fullmatch(number):
        return int(number)
    return float(number)


def perl_divide(total, count):
    """$total / $count as 64-bit Perl computes it for integer arguments."""
    if isinstance(total, int) and isinstance(count, int) and count and total % count == 0:
        return float(total // count)
    return float(total) / float(count)


def num_desc(keys):
    """sort {$b <=> $a} keys"""
    return sorted(keys, key=perl_num, reverse=True)


def num_asc(keys):
    """sort {$a <=> $b} keys"""
    return sorted(keys, key=perl_num)


def nsort(keys):
    """nsort from Sort::Naturally: numbers inside the text compare as numbers,
    e.g. CLUSTER-0-9999 before CLUSTER-0-10000."""
    def natural_key(text):
        return [(0, int(chunk), "") if chunk[0] in "0123456789" else (1, 0, lc(chunk))
                for chunk in re.findall(r"[0-9]+|[^0-9]+", text)]
    return sorted(keys, key=natural_key)
