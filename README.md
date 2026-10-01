# teiresias1_python

The ARResT/Subsetter tool (folder `teiresias1/`) rewritten in Python.
The files of the original folder are unchanged.

## How to run

Requires Python 3.8 or newer, with no extra packages.

```
py subsetter.py --receptor IG --in test\test.fasta.e.vq.table
```

The options are the same as in `subsetter.pl`:

| Option | What it does | Default |
|---|---|---|
| `--receptor` | `IG` or `TR`. With `TR` the IGHV phylogenetic groups are not applied | required |
| `--in` | the IMGT/V-QUEST table | required |
| `--ide` | minimum identity | 0.5 |
| `--sim` | minimum similarity | 0.7 |
| `--ignorephylo` | ignore the phylogenetic groups | no |
| `--annotation` | file with `seq id <TAB> annotation` lines | none (`UNK` is written) |
| `--out_table` | where the final table is written | `<in>.ss.table` |
| `--teiresias` | another teiresias_char executable, e.g. `teiresias_char.unix` (runs through WSL) | `teiresias_char.exe` in this folder |

As in the Perl, all intermediate files are written next to the input file,
with the same names.

## Which file corresponds to which

| Perl (`teiresias1/`) | Python (here) | Step |
|---|---|---|
| `subsetter.pl` | `subsetter.py` | orchestrator, the 13 fixed parameters |
| `subsetter--input.pl` | `subsetter_input.py` | 1. V-QUEST table → `.orig.board`, `.ign2.board`, `.ign2.dat`, `.ign2.board4nrgrep` |
| `teiresias_char.exe` (unchanged) | `subsetter_teiresias.py` | 2. pattern discovery → `.ign2.thr` |
| `subsetter--nrgrep.pl` + `nrgrep` | `subsetter_nrgrep.py` | 3. which sequences contain each pattern → `.thr.hits` |
| `subsetter--pairs.pl` | `subsetter_pairs.py` | 4. pairs and scores → `.pairs.layout` |
| `subsetter--clustering.pl` | `subsetter_clustering.py` | 5. clusters L0, L1, ... → `.lite.clusters`, `.lite.clusters.table` |
| (Perl behaviour) | `perl_compat.py` | split, numbers from strings etc., so that the files come out identical |

Variable names are the same as in the Perl (`per_pattern`, `patternblast`,
`cluster2member`...), so the two can be read side by side. Each step can
also be run on its own, e.g. `py subsetter_input.py table.txt`.

## Test

```
py test\run_test.py
```

Runs the tool on `test/test.fasta.e.vq.table` and compares every file with
the Perl output in `test/expected/`. The results are kept in `test/output/`.

## Differences from the Perl

### 1. nrgrep was replaced by an exact search

The Perl searched each pattern with `nrgrep -k 1i '@.*PATTERN'`, i.e. "find
the pattern, allowing one extra letter". The original nrgrep (the Linux
build, run through WSL) was tested and does not do this consistently.
Depending on the pattern, it internally picks a different algorithm:

| Pattern | Sequence | Found? |
|---|---|---|
| `ABCDEFGH` | `QQABCDEFGXHQQ` (extra X in the second half) | yes |
| `ABCDEFGH` | `QQABXCDEFGHQQ` (extra X in the first half) | **no** |
| `ABCD` | `QQABCDQ` | yes |
| `ABCD` | `QQABCD` (the pattern ends at the end of the sequence) | **no**, although it is an exact occurrence |

In addition, the next step (pairs) only works with exact occurrences. When a
sequence comes from nrgrep with an extra letter, pairs cannot find where the
pattern falls and gives it the position of the previous sequence.

The Python port therefore uses an exact search: a sequence is found if and
only if it contains the pattern. On the test data the result is identical.
On larger data it can differ (see "Comparison on 3000 sequences" below).

### 2. teiresias_char.exe (Windows) and teiresias_char.unix (Linux)

The Perl used `teiresias_char.unix`. The Windows `.exe` does not produce
exactly the same file:

* On the test data it finds the same 23 patterns, in a slightly different
  order. Only the pattern number `[n]` in `.pairs.layout` changes. The final
  table is identical.
* On 3000 CLL sequences they find 484,598 (Linux) and 484,455 (Windows)
  patterns. About 0.2% differ, all of them with an equivalence class such as
  `[AVLI]`. All of them were checked to be valid patterns in both cases; the
  two programs simply choose differently between equivalent options. The
  final table has 800 lines with the Linux program and 807 with the Windows
  one.

The default is the `.exe`, which runs directly on Windows. For results
identical to the Perl, pass the Linux program, which runs through WSL:

```
py subsetter.py --receptor IG --in table.txt --teiresias teiresias_char.unix
```

With it, all 10 files of the test data are identical to the Perl output,
byte for byte.

### 3. Fixed order in the clustering

In ten places `subsetter--clustering.pl` iterates over a hash without
sorting it. Since version 5.18, Perl shuffles hash order on every run, so
when two links had the same score, the same connectivity and the same
length, which one was processed first was random. The Python port always
processes them in alphabetical order (numerical for clusters), so every run
gives the same result.

## Perl quirks kept on purpose

They were kept so that the results are identical. They are the first
candidates for fixing in the improvement phase. Each one has a comment in
the code.

1. **The last pattern of `.hits` is never used** (`subsetter_pairs.py`).
   The Perl processes a pattern when it reaches the next `#`, and there is
   none after the last one.
2. **A "sticky" pattern position** (`subsetter_pairs.py`, `cdr3offset`). If
   the pattern is not found exactly in a CDR3, the position of the previous
   sequence is kept. With the exact search of step 3 this practically no
   longer happens.
3. **At level L1 the "linked with all members?" check almost always fails**
   (`subsetter_clustering.py`, `joins_all`). It looks at the links between
   sequences instead of the links between clusters, so a new cluster is
   created almost every time.
4. **`connectivity` is shared between sequences and clusters and is not
   reset between levels.** If sequence IDs are plain numbers (e.g. `5`),
   they get mixed up with cluster numbers.
5. **In the final table, if a sequence belongs to two clusters at one level,
   only the last one is followed to the next levels** (`write_table`).

## Comparison on 3000 sequences

The test data has only 9 sequences and one cluster. To also check levels
L1-L3 and the phylogenetic constraint, a table was built from the 3000 CLL
sequences in `teiresias2/data` (junction from position 104 to the W of the
WGxG motif). The original Perl was run in WSL with the Linux programs.

| Check | Result |
|---|---|
| Step 1 (input), Perl vs Python | identical files, byte for byte |
| Steps 4-5 with the same `.hits` (pairs, clustering) | identical files, byte for byte: 17,355 pairs, 810 table lines* |
| Step 3: exact search instead of nrgrep | 5 of 246 L0 clusters change |
| Time of step 3 | Perl 886 seconds, Python 22 |

\* The comparison used a copy of the Perl in which the loops iterate in
alphabetical order (see "Fixed order in the clustering"). The original Perl
on the same data gave 3 different tables in 5 runs (814, 805 and 810 lines).
Any two orders, Perl or Python, agree on about 93-98% of the L0 clusters.

### Example from the 5 clusters that changed

Pattern `ARD.[AVLI]I`:

| Sequence | Junction | Contains the pattern? |
|---|---|---|
| DE-02-3501-H1 | `CARDFLIGW` | yes: `ARD` `F` `L` `I` |
| GR-01-0348-H1 | `CARDPITIF` | no: `ARD` `P` `I` **`T`** `I` |

nrgrep found both, because it allows one extra letter (the `T`). pairs did
not find the pattern in the second one, gave it the position of the first,
and put them in the same cluster. In the Perl `.pairs.layout` the second
sequence has no highlighted positions:

```
"DE-02-3501-H1  cARDfLIgw ..."   "GR-01-0348-H1  cardpitif ..."   571429714286   ARD.[AVLI]I
```

With the exact search of the Python port this link is not made.
