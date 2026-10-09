# Notes: how teiresias1 works

teiresias1 is the Perl tool (ARResT/Subsetter) in the folder `teiresias1/`. It
groups sequences whose CDR3 regions share a pattern.

It works in **two phases**:

1. It builds a **list of links** between sequences. This is done once.
2. It builds **groups** by reading that list.

After phase 1 it **never compares sequences again**. Phase 2 only looks things
up in the list.

---

## Phase 1: build the list of links (done once)

| Step | File | What it does |
|---|---|---|
| 1 | `subsetter--input.pl` | keeps only the CDR3 of each sequence (the junction, minus one amino acid at each end) |
| 2 | `teiresias_char.exe` | finds **patterns**: pieces that occur in at least two sequences |
| 3 | `subsetter--nrgrep.pl` | for each pattern, finds which sequences contain it |
| 4 | `subsetter--pairs.pl` | turns shared patterns into **linked pairs**, each with a score |

Two sequences become a linked pair only if all of these hold:

- they share a pattern **at the same position** in the CDR3,
- their CDR3 has **the same length**,
- their V genes are in **the same phylogenetic group** (IGHV1/5/7, IGHV2/4/6 or
  IGHV3; this check is skipped for TR or with `--ignorephylo`).

The result is a list like this (a higher score means more similar):

```
A–B  95
D–E  95
A–C  90
B–C  90
C–D  60
```

Here sequence F shares no pattern with anyone, so it has no link and does not
appear in the list at all.

**What "match" means here:** two sequences match when they share a pattern at
the same position of the CDR3. It does not mean "similar overall", as in
teiresias2.

---

## Phase 2: build the groups from the list

`subsetter--clustering.pl` reads the list from the strongest score to the
weakest. It does not look at the sequences any more, only at the list.

**The rule:** a sequence joins a group only if the list contains a link between
it and **every** member of that group.

| Line of the list | What it checks | Result |
|---|---|---|
| A–B | does either one already have a group? No | new group {A, B} |
| D–E | same | new group {D, E} |
| A–C | is C–A in the list? Yes. Is C–B in the list? Yes | C joins: {A, B, C} |
| B–C | they are already together | nothing |
| C–D | they are in different groups | records that the two groups are linked (used for level 1) |

### The three cases when a line is read

Every line of the list names two sequences. What happens depends on whether
they already belong to a group.

**Case 1: neither has a group.** They form a new group together.

```
reading A–B
before:   (no groups yet)
after:    {A, B}
```

**Case 2: one has a group, the other does not.** The one without a group asks
to join. It may join only if the list has a link between it and **every**
member of that group.

When it can join (reading A–C in the example above):

```
reading A–C:   C has no group, A is in {A, B}
check:         C–A in the list?   yes  ✓
               C–B in the list?   yes  ✓
before:        {A, B}   {D, E}
after:         {A, B, C}   {D, E}
```

When it cannot join (the same list, but without the line B–C):

```
reading A–C:   C has no group, A is in {A, B}
check:         C–A in the list?   yes  ✓
               C–B in the list?   no   ✗
before:        {A, B}   {D, E}
after:         {A, B}   {D, E}   {A, C}        ← A is now in two groups
```

One link is missing, so C does **not** join {A, B}. But the link A–C is real,
so it is not thrown away: the tool makes a **new, separate group** with just
those two. It also keeps a note that {A, B} and {A, C} are linked (through A).
This is why one sequence can appear on more than one row of the output table.

**Case 3: both already have groups.** The tool does **not** merge the two
groups here, and does not make a new one. It only keeps a **note** that the
two groups are linked. The notes are used one level up, where linked groups
join.

```
reading C–D:   C is in {A, B, C}, D is in {D, E}
before:        {A, B, C}   {D, E}
after:         {A, B, C}   {D, E}        note: {A, B, C} ↔ {D, E}
```

If both are already in the same group (B–C in the example), nothing changes.

Both versions of the example (with and without B–C) were run through the real
`subsetter--clustering.pl`, and these are the groups it produced.

### A common misreading

It is tempting to think the tool puts two sequences together and then compares
that pair with all the remaining sequences again. **It does not.** When C's
turn comes, the tool does not compare C with A and B. It looks up whether the
lines C–A and C–B exist in the list. If even one is missing, C does not join.

### Levels: L0, L1, L2 ...

When the whole list has been read, the groups above form **level 0**:

```
L0:   [A B C]   [D E]
```

Then each level-0 group becomes a single unit, and groups that have a note
between them (case 3) join together. Because of the note left by C–D:

```
L1:   [ [A B C] + [D E] ]
```

So at level 0 a link between two groups is only written down; the groups are
joined one level higher.

The tool keeps going up a level while new groups are still being formed, and
stops when there is nothing left to join. It never forces everything into one
top group, so there is **no root**.

F has no links, so it **never appears in the output**, not even as "alone".

These levels are the columns `L0`, `L1`, `L2` ... of the output table. Level 0
is the most detailed; each higher level is coarser.

### How strict each level is

The rule gets looser as you go up:

| Level | The units are | To join a group, a unit needs a link to ... |
|---|---|---|
| L0 | sequences | **every** sequence of the group |
| L1 | L0 groups | **every** L0 group of the group |
| L2 and up | groups of the level below | **one** member of the group is enough |

So:

- **At L0** only sequences that are all linked to each other are together.
  This is where the groups are separated best.
- **At L1** it is enough that two L0 groups **touch**: at least one sequence of
  one group is linked to at least one sequence of the other. The sequences of
  an L1 group do not all need to be linked to each other.
- **From L2 up** a single link is enough, so the groups become broad families.

The "every member" check is applied only when the level is 0 or 1 (lines 124
and 186 of `subsetter--clustering.pl`).

---

## Reading the output table: a real example

### One row per path

Each row of the output table is **one path** from L0 upwards. If a sequence's
L0 group belongs to three L1 groups, the sequence is written on three rows:
same sequence, same L0, only L1 (and above) differs. The three rows are the
same information repeated, not three different results.

**When counting** how many sequences a group has from this table, count
**unique** sequence IDs, otherwise each sequence is counted once per path.

### An L0 group that is in three L1 groups

From `test_run_t1_py/Final_27.txt.ss.table`. The L0 group 0274 belongs to
three L1 groups, and all three go up to the same L2 group:

```
L2                              CLUSTER-2-0024
                    ┌─────────────────┼─────────────────┐
L1           CLUSTER-1-0124     CLUSTER-1-0132     CLUSTER-1-0138
              ┌────┴────┐        ┌────┴────┐        ┌────┴────┐
L0          0274       0005    0274       0279    0274       0260
```

| L0 group | Sequences |
|---|---|
| 0274 | DE-02-3594, ES-01-0512 |
| 0005 | DE-02-4492, DE-02-7182, DE-03-7692, BE-02-0039 |
| 0279 | DE-03-10569, DE-03-9900 |
| 0260 | BLOOD-HC66_1651, DE-03-9900 |

All of them share the piece `YDFWSGY…YYGMDV` of the CDR3.

The other L0 groups (0005, 0279, 0260) are **separate groups**, with their own
sequences. They are **not** inside 0274; they meet 0274 only at L1.

Each L1 group here is **exactly** two L0 groups and nothing else:

```
L0    0274 = {DE-02-3594, ES-01-0512}
      0005 = {DE-02-4492, DE-02-7182, DE-03-7692, BE-02-0039}

L1    0124 = 0274 + 0005                        → 6 sequences, nothing else

L2    0024 = 0124 + 0132 + 0138 + any others    → other sequences join here
```

Inside 0005 all four sequences are linked to each other (L0 rule). Inside 0124
the six sequences are **not** all linked to each other: it is enough that 0274
and 0005 touch (L1 rule).

0274 ends up in three L1 groups by case 2, one level up: 0274 is linked to
0005, 0279 and 0260, but those three are not linked to each other. So each time
one of them asks to join, it is not linked to every member, and a new separate
L1 group is made. The same thing happens at L0: DE-03-9900 is in both 0279 and
0260.

**To find the other L0 groups of an L1 group** in your own table: filter the
`L1` column for that L1 group (for example `CLUSTER-1-0048`) and look at the
`L0` column. Every value there other than your own L0 group is another group.

### The structure depends on the run

In this run the four sequences DE-02-3594, IT-01-0361, DE-03-9206 and
NY-01-0413 form the L0 group 0004, which belongs to **one** L1 group only:

```
L1    0125 = 0004 + 0184
```

In another run (a Perl run, not in this folder) the same four sequences formed
`CLUSTER-0-0000`, which belonged to **three** L1 groups (0048, 0049, 0050), so
each of them was written on three rows. How many paths a group has depends on
the run: on false links from nrgrep in the Perl, and on the order of tied
links.

### A suspicious L0 link (from the Perl run)

In that Perl run, IT-01-0361 was in two L0 groups: `CLUSTER-0-0000` (with the
three sequences above) and `CLUSTER-0-0088`, together with
CLLGreek-IG8220_224. This is case 2:

```
reading:   IT-0361 – CLLGreek               (score 70.37, read after the
                                              stronger 100 and 88.89 links)
           CLLGreek has no group; IT-0361 is already in CLUSTER-0-0000

check:     CLLGreek needs a link to ALL 4 members of CLUSTER-0-0000.
           At least one is missing                                    ✗

after:     CLUSTER-0-0000 = {DE-3594, IT-0361, DE-9206, NY-0413}
           CLUSTER-0-0088 = {IT-0361, CLLGreek}                ← new group
```

But this link itself is very likely **false**:

```
shared piece:   DYDFWSGYYPN YYYY  GMDV
IT-0361:        DYDFWSGYYPN YYYY  GMDV     ✓ contains it exactly, position 10
CLLGreek:       dydfwsgyypn yyyyy gmdv     ✗ one extra Y
```

- In the table, upper-case letters mark where the pattern was found. The CDR3
  of CLLGreek is **all lower-case**, so the pattern was never found in it.
- It has one extra Y. nrgrep allows "one extra letter", so it let it through;
  the next step could not find the position and reused the position of the
  previous sequence (see `teiresias_python.md`, section 2.1).
- Position by position the two CDR3s are identical in only 11 of 27
  positions. A real link at 50% identity needs at least 14.

With an exact search, as in the Python port, this link would not exist.

---

## Ties: why two runs can give different results

The tool does **not** put sequences into groups at random. The rule is fixed:
the strongest link is read first. Chance comes in only when links are
**tied**: same score, same number of links, and same CDR3 length. Perl reads
tied links in an order that it shuffles on every run. Because groups are built
on a first-come basis (the three cases above), a different order can give
different groups.

### Example: four sequences in a square

All four links have score 90, and there are no diagonals (no A–C, no B–D):

```
A ──── B
│      │
D ──── C
```

The real `subsetter--clustering.pl` was run 12 times on exactly this input. It
gave three different level-0 results:

| Level 0 | Runs |
|---|---|
| {A B} {B C} {C D} | 5 |
| {A B} {C D} | 4 |
| {A D} {B C} | 3 |

To see why, follow two different reading orders.

**Order 1: A–B, C–D, B–C, D–A**

| Line | What happens | Groups after |
|---|---|---|
| A–B | case 1: neither has a group, new group | {A, B} |
| C–D | case 1: neither has a group, new group | {A, B}   {C, D} |
| B–C | case 3: both have groups, note only | {A, B}   {C, D}   note: {A, B} ↔ {C, D} |
| D–A | case 3: both have groups, note only | {A, B}   {C, D}   note: {A, B} ↔ {C, D} |

Result: **{A B} {C D}**.

**Order 2: A–B, B–C, C–D, D–A**

| Line | What happens | Groups after |
|---|---|---|
| A–B | case 1: new group | {A, B} |
| B–C | case 2: C asks to join {A, B}. C–B ✓, C–A ✗, so a new group | {A, B}   {B, C} |
| C–D | case 2: D asks to join {B, C}. D–C ✓, D–B ✗, so a new group | {A, B}   {B, C}   {C, D} |
| D–A | case 3: both have groups, note only | {A, B}   {B, C}   {C, D} |

Result: **{A B} {B C} {C D}**.

Both results follow the rules. Only the reading order was different.

### What this means for the output files

A difference at level 0 carries on into L1, L2 and the final table. That is why
the output files can differ slightly between two runs on the same input. On
3,000 CLL sequences, five runs of the same Perl gave three different final
tables.

The order is left unsorted in two kinds of places in
`subsetter--clustering.pl`:

- lines 102 and 105 decide **which tied link is read first**;
- lines 126 and 188 decide **which group is tried first**, when a sequence
  already belongs to more than one.

The Python port (`teiresias1_python/`, see `teiresias_python.md`, section 2.3)
always reads in alphabetical order, so every run gives the same table. That
table is one of the tables the Perl could have produced.

---

## Where it happens in the code

All in `teiresias1/subsetter--clustering.pl`:

| What | Line |
|---|---|
| sort the links from the highest score | 99 |
| read tied links, without sorting | 102 and 105 |
| case 1: neither sequence has a group, so a new group | 109 |
| case 2: one has a group, try to join it | 122 and 184 |
| "every member" check only at L0 and L1 (`if ($level <= 1)`) | 124 and 186 |
| which group is tried first, without sorting | 126 and 188 |
| check "is the candidate linked to every member?" | 131 and 193 |
| case 2: could not join, so a new separate group | 146 and 208 |
| case 3: both have groups, so only a note | 246 |
| print only groups with more than one member | 279 |
| "if at least two groups were formed, go up one level" | 326 |
