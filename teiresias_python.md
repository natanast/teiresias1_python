# Notes: teiresias1_python

The Perl tool in `teiresias1/` (ARResT/Subsetter) rewritten in Python, in the
folder `teiresias1_python/`. The original folder was not changed.

The algorithm is the same, step by step. Given the same input, steps 1 and 4
write the same files as the Perl, byte for byte, and so does step 5 when the
Perl's random order (2.3) is fixed the same way. Three things were changed on
purpose or are outside the port's control; they are explained in section 2.

---

## 1. How to run it

### The command

The Python command takes exactly the same arguments as the Perl one. Only
the start changes.

Perl (old):

```
/apps/teiresias1/subsetter.pl --receptor 'IG' --in 'Total_sequences_MCL.txt' --sim 0.70 --ide 0.50 --annotation
```

Python (new), from a terminal opened in the folder that contains the data:

```
py "C:\Users\natan\Υπολογιστής\teiresias1\teiresias1_python\subsetter.py" --receptor 'IG' --in 'Total_sequences_MCL.txt' --sim 0.70 --ide 0.50 --annotation
```

It prints the five steps and the name of the final table:

```
1/5 reading the V-QUEST table
2/5 TEIRESIAS: pattern discovery
3/5 searching the patterns in the sequences
4/5 sequence pairs
5/5 clustering
done: Total_sequences_MCL.txt.ss.table
```

Requirements: Python 3.8 or newer, nothing else. No packages and no virtual
environment are needed, because only the standard library is used.

### Options

| Option | Meaning | Default |
|---|---|---|
| `--receptor` | `IG` or `TR`. With `TR` the IGHV phylogenetic groups are not used | required |
| `--in` | the IMGT/V-QUEST table | required |
| `--ide` | minimum identity | 0.5 |
| `--sim` | minimum similarity | 0.7 |
| `--ignorephylo` | ignore the IGHV phylogenetic groups | off |
| `--annotation` | file with `seq id <TAB> value` lines. Without a value it means "no annotation" and the column says `UNK` | no annotation |
| `--out_table` | where to write the final table | `<in>.ss.table` |
| `--teiresias` | which teiresias_char to use, e.g. `teiresias_char.unix` (new, see 2.2) | `teiresias_char.exe` on Windows, `teiresias_char.unix` on Linux |

### Try a small part first

In PowerShell, in the folder of the data:

```
# number of lines (sequences = lines - 1, for the header)
(Get-Content 'Total_sequences_MCL.txt' | Measure-Object -Line).Lines

# the header and the first 1000 sequences
Get-Content 'Total_sequences_MCL.txt' -TotalCount 1001 | Set-Content 'MCL_1000.txt'
```

Then run the same command with `--in 'MCL_1000.txt'`.

### Run times on this PC

Measured on real CLL sequences (Ryzen 7 5700U, 16 GB RAM), whole pipeline:

| Sequences | Time |
|---|---|
| 1,000 | 15 seconds |
| 3,000 | 1.5 minutes |
| 10,000 | 15.5 minutes |

Memory stayed under 0.5 GB even for 10,000 sequences. Most of the time is
TEIRESIAS, and its time roughly quadruples when the number of sequences
doubles (rough estimate: 20,000 sequences about 1 hour, 40,000 about 4
hours). For bigger files a server is a better fit.

### In Positron

Open the data folder (File → Open Folder), choose Python 3.11 if asked, and
paste in the Python console:

```python
import sys
sys.path.insert(0, r"C:\Users\natan\Υπολογιστής\teiresias1\teiresias1_python")
import subsetter

table = subsetter.run("Total_sequences_MCL.txt", "IG", ide="0.50", sim="0.70")

import pandas as pd
df = pd.read_csv(table, sep="\t")
```

### In WSL or on a Linux server

Same command with `python3` instead of `py`. The Windows disk is under
`/mnt/c`:

```
python3 "/mnt/c/Users/natan/Υπολογιστής/teiresias1/teiresias1_python/subsetter.py" --receptor 'IG' --in 'Total_sequences_MCL.txt' --sim 0.70 --ide 0.50 --annotation
```

On Linux the tool uses `teiresias_char.unix` automatically, the same program
the Perl used. For large files inside WSL, keep the data in the WSL home
folder (e.g. `~/data`) rather than under `/mnt/c`, because WSL writes slowly
to the Windows disk.

### Check that everything works

```
py "C:\Users\natan\Υπολογιστής\teiresias1\teiresias1_python\test\run_test.py"
```

It runs the tool on the Perl test data and compares every file with the Perl
output. The last line must be `OK: everything matches the Perl`.

### What it writes

As in the Perl, all files are written next to the input file, with the same
names:

| File | Step | Content |
|---|---|---|
| `.orig.board` | 1 | label + full junction |
| `.ign2.board`, `.ign2.dat`, `.ign2.board4nrgrep` | 1 | junction without the first and last amino acid, in three formats |
| `.ign2.thr` | 2 | the patterns found by TEIRESIAS |
| `.ign2.thr.hits` | 3 | for each pattern, the sequences that contain it |
| `...pairs.layout` | 4 | linked pairs of sequences with their best score |
| `...lite.clusters` | 5 | clusters at each level |
| `.ss.table` | 5 | **the final table** |

One line of the final table, from `Final_27.txt`:

| Column | Value | Meaning |
|---|---|---|
| sequence.seq id | `CLLGreek-E21358_42` | |
| V, D, J genes | `IGHV3-13*05 F`, `IGHD3-3*01 F`, `IGHJ6*02 F` | |
| junction | `cARDSGNQYYDFWSGYYRDYYYYYGMDVw` | upper case = positions fixed by the shared pattern, lower case = the rest |
| composite score | `999999999999` | 6 digits identity + 6 digits similarity (100 % is written 999999) |
| identity, similarity | `100.00`, `100.00` | the same two numbers as percentages |
| annotation | `UNK` | from `--annotation`, `UNK` if none |
| L0, L1, L2, ... | `CLUSTER-0-0000`, `CLUSTER-1-0003`, ... | the cluster at each level |

A sequence can appear on more than one line if it belongs to more than one
cluster.

### What L0, L1, L2, ... mean

They are the levels of a hierarchy. L0 holds the tightest clusters; each
higher level joins clusters of the level below into bigger, looser ones.

Example: the first line of `Final_27.txt.ss.table`, sequence
`CLLGreek-E21358_42`, has `CLUSTER-0-0000`, `CLUSTER-1-0003`,
`CLUSTER-2-0000`, `CLUSTER-3-0001`, `CLUSTER-4-0000`. That is:

```
L4  CLUSTER-4-0000 ............................ 103 sequences
 └ L3  CLUSTER-3-0001 ......................... 78
    └ L2  CLUSTER-2-0000 ...................... 18
       └ L1  CLUSTER-1-0003 ................... 12
          ├ L0  CLUSTER-0-0000 ................ 10   CARDSGNQYYDFWSGYYRDYYYYYGMDVW
          └ L0  CLUSTER-0-0096 ................ 3    CARDGGGRYYDFWSGYYRFPYYYYGMDVW
```

**L0: the strict cluster.** Every sequence is linked directly with every
other one: a shared pattern at the same position, the same CDR3 length,
identity ≥ 50 %, similarity ≥ 70 %, and the same IGHV phylogenetic group.
This corresponds to a stereotyped subset. `CLUSTER-0-0000` has 10 almost
identical sequences:

```
CARDSGNQYYDFWSGYYRDYYYYYGMDVW   (7 identical)
CARDSALQYYDFWSGYYRDYYYYYGMDVW
CARDSAFQYYDFWSGYYPDYYYYYGMDVW
CARDSGNQYYDFWSGYSSLYYYYYGMDVW
```

**L1: a cluster of L0 clusters.** Two L0 clusters are joined when some of
their sequences are linked, even if not all of them are. `CLUSTER-1-0003`
joins `CLUSTER-0-0000` with `CLUSTER-0-0096`, which share
`YYDFWSGYY...YYYYGMDV` but differ at the start (`SGNQ` vs `GGGR`, `LGNT`).
Related, not identical. Because of bug 3 in section 3, an L1 cluster
practically always contains exactly two L0 clusters: in `Final_27.txt` all
199 of them do.

**L2, L3, L4:** the same again, clusters of clusters, each level looser. At
L4 the whole neighbourhood of 103 sequences with this motif is one cluster.
The process stops when nothing else can be joined.

Across the whole file:

| Level | Clusters |
|---|---|
| L0 | 354 |
| L1 | 199 |
| L2 | 36 |
| L3 | 7 |
| L4 | 1 |

How to read it:

* For subsets, look at L0. That is the main result.
* L1 and above show which subsets resemble each other.
* The numbers (`0000`, `0096`) are just names, not a ranking.
* If a line stops before a level (no cluster in that column), its cluster
  was not joined with any other at that level.
* A sequence can be in two L0 clusters: `CLLGreek-E21358_42` is in both
  `CLUSTER-0-0000` and `CLUSTER-0-0096`, so it appears on more than one line.

---

## 2. What changed compared to the Perl tool

| Change | Why | Effect on results |
|---|---|---|
| 2.1 nrgrep replaced by an exact search | nrgrep.exe does not run here, and nrgrep made false links | small: 97% of L0 clusters identical on `Final_27.txt` |
| 2.2 teiresias_char.exe by default instead of .unix | the .unix only runs on Linux | about 0.2% of patterns differ (on 3000 sequences); avoidable with `--teiresias teiresias_char.unix` |
| 2.3 fixed order in the clustering | the Perl order was random | the Perl itself gives different results from run to run |
| 2.4 smaller changes | usability, speed | none |

### 2.1 nrgrep replaced by an exact search

The Perl searched each pattern with `nrgrep -k 1i`, which means "allow one
extra letter". Testing the original nrgrep showed that it does not do this
consistently, because it picks a different internal algorithm depending on
the pattern:

| Pattern | Sequence | nrgrep finds it? |
|---|---|---|
| `ABCDEFGH` | `QQABCDEFGXHQQ` (extra X in the second half) | yes |
| `ABCDEFGH` | `QQABXCDEFGHQQ` (extra X in the first half) | **no** |
| `ABCD` | `QQABCD` (exact match at the end of the sequence) | **no** |

Also, the next step (pairs) only understands exact matches. When nrgrep
returns a sequence that contains the pattern only with an extra letter,
pairs cannot find the pattern position and silently reuses the position of
the previous sequence. That creates false links. Real example, pattern
`ARD.[AVLI]I`:

| Sequence | Junction | Contains the pattern? |
|---|---|---|
| DE-02-3501-H1 | `CARDFLIGW` | yes: `ARD` `F` `L` `I` |
| GR-01-0348-H1 | `CARDPITIF` | no: `ARD` `P` `I` **`T`** `I` |

The Perl put these two in the same cluster. In its `.pairs.layout` the second
one has no upper-case positions (`cardpitif`), which is the sign that the
pattern was never found in it.

The Python port finds a sequence if and only if it contains the pattern.

### 2.2 teiresias_char.exe and teiresias_char.unix

The Perl used `teiresias_char.unix` (Linux). On Windows the Python port uses
`teiresias_char.exe` by default. Both find valid patterns, but they do not
always choose the same ones:

* Test data: the same 23 patterns, in a slightly different order. The final
  table is identical.
* `Final_27.txt`: 881,431 patterns (Linux) vs 881,334 (Windows).

For results identical to the Perl, add `--teiresias teiresias_char.unix` on
Windows (it runs through WSL). With it, all 10 files of the test data are
identical to the Perl output.

The .exe cannot open long paths or paths with Greek characters
(`Υπολογιστής`). The Python port runs it in a temporary folder with short
names, so this no longer matters.

### 2.3 Fixed order in the clustering

In ten places `subsetter--clustering.pl` goes through a list without sorting
it. Perl shuffles that order on every run, so when two links tie (same
score, same connectivity, same length) the winner was random. On 3000 CLL
sequences, 5 runs of the same Perl gave 3 different final tables.

The Python port always uses alphabetical order (numerical for clusters), so
every run gives the same table. Its result is one of the results the Perl
could have given.

### 2.4 Smaller changes

* Step 3 is much faster: on 3000 sequences, 22 seconds instead of 886. The
  sequences are grouped by length and each pattern is only searched in
  sequences short enough to pass the identity/similarity filter. The output
  is identical.
* `--annotation` without a value works, as in the Perl.
* New option `--teiresias` (see 2.2).
* If the input file does not exist, the tool stops with
  `(!) error: input file not found: ...` and writes nothing. (Before this
  fix it left two empty files, like the `Total_sequences_MCL.txt.*` files in
  `test_run_t1_py/`, which can be deleted.)
* If the annotation file does not exist, a warning is printed. The Perl
  ignored it silently.
* On Linux, `teiresias_char.unix` works even if it lost its execute
  permission when copied from Windows.
* Each step can run on its own, like the Perl scripts, e.g.
  `py subsetter_input.py table.txt`.
* README, comments and messages are in English.

---

## 3. What did not change: Perl bugs kept on purpose

They were kept so that the results match the Perl. They are the first
candidates for the improvement phase. Each one is marked in the code.

1. **The last pattern of the `.hits` file is never used** (`subsetter_pairs.py`).
   The Perl processes a pattern when it reaches the next `#` line, and there
   is none after the last pattern.
2. **Sticky pattern position** (`subsetter_pairs.py`, `cdr3offset`). If a
   pattern is not found exactly, the position of the previous sequence is
   reused. With the exact search of 2.1 this practically no longer happens.
3. **At level L1 the "linked with all members?" check almost always fails**
   (`subsetter_clustering.py`, `joins_all`). It looks up links between
   sequences instead of links between clusters, so a new cluster is created
   almost every time.
4. **`connectivity` is shared between sequences and clusters and never reset
   between levels.** Sequence IDs that are plain numbers (e.g. `5`) get mixed
   up with cluster numbers.
5. **In the final table, if a sequence is in two clusters at one level, only
   the last one is followed to the next levels** (`write_table`).

---

## 4. How the port was checked

| Data | Check | Result |
|---|---|---|
| Perl test data (9 sequences) | whole pipeline with `teiresias_char.unix` | all 10 files identical, byte for byte |
| Perl test data | whole pipeline with `teiresias_char.exe` | final table identical |
| 3000 CLL sequences from `teiresias2/data` | steps 4-5 from the same `.hits` as the Perl | identical, byte for byte (17,355 pairs, 810 lines)* |
| 3000 CLL sequences | effect of the exact search alone | 5 of 246 L0 clusters change |
| `Final_27.txt` (2,473 sequences) | step 4 from the same `.hits` as the Perl | identical, byte for byte |
| `Final_27.txt` | final table vs the Perl run in `big_run/` | 348 of 357 L0 clusters identical (97%) |

\* Compared with a copy of the Perl in which the ten loops of 2.3 are sorted
the same way, since the original Perl gives a different answer each run.

On `Final_27.txt` the remaining 3% splits like this: with the same `.hits`
as the Perl, 352 of 357 L0 clusters are identical (only the order of 2.3
differs); the exact search of 2.1 accounts for the rest. Three sequences
that the Perl clustered are not clustered by the Python port (642 vs 639).
All three were clustered by the Perl only through the false links of 2.1.
In fact all 5 L0 clusters that change because of the exact search are
pairs in which one of the two sequences does not contain the pattern (no
upper-case positions in the Perl output, e.g. `carsrtmgyydssgyyygpgyyyygmdvw`
for the pattern `[KRH]..GYYDSSGYYY..YYYYGMDV`). In two cases the false link
even replaced a real one: the Perl paired `BLOOD-HC80_782` with
`CLLGreek-IG8228_893` (false, score 666667703704) instead of
`HAEM-naB-2_13391` (real, score 629630703704); the Python port keeps the real
one.

With `--teiresias teiresias_char.unix` the `.thr` is identical to the Perl
one, and the 96 sequence-pattern hits that nrgrep added (54 patterns) are the
only difference in the `.hits`. The .exe and the .unix give the same L0
clusters on this file; they differ in 1 L1 and 2 L2 clusters. Five runs of
the original Perl on the same `.pairs.layout` gave five different final
tables.

---

## 5. Open items

* **Memory in step 4 (pairs) on large data.** The port keeps the Perl data
  structures, so the problem described in `notes.md` section 7 (labels used
  as keys) is still there. It is the first thing to fix before running the
  very large datasets.
* **The five Perl bugs of section 3.**
* **A deliberate tie-break rule for the clustering** (2.3), instead of
  alphabetical order, if a biologically meaningful one is preferred.
