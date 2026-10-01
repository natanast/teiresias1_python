"""
Step 5: subsetter--clustering.pl in Python.

Takes the pairs in <layout> and builds clusters in levels:

  Level 0: links are processed from the highest score to the lowest. A
  sequence joins an existing cluster only if it is linked with ALL of its
  members, otherwise it starts a new cluster.
  Levels 1, 2, ...: the clusters of the previous level become the new units
  and are joined the same way. Stops when a level creates fewer than two
  clusters.

Writes:
  <layout>.<threshold>.ratio<ratio>.lite.clusters        cluster -> member -> score
  <layout>.<threshold>.ratio<ratio>.lite.clusters.table  the final table (L0, L1, ...)

Perl quirks kept on purpose (for identical results):
  * The connectivity dictionary is shared between sequences and clusters and
    is not reset between levels.
  * At level 1 the "linked with all members?" check looks at the links of
    the sequences, not of the clusters, so it almost always fails and a new
    cluster is created.
  * In a few places the Perl iterates over a hash without sorting it, so the
    order was random. Here the order is always alphabetical (numerical for
    clusters), so that every run gives the same result.

Usage:  python subsetter_clustering.py <layout> <orig.board> [threshold] [ratio]
"""

import re
import sys

from perl_compat import (WORD, at, nsort, num_asc, num_desc, perl_divide, perl_num, perl_split,
                         perl_split_ws, read_lines, s, utf8_output, write_open)


def read_board(board_path):
    """From the orig.board: V, D, J genes and junction length of each sequence."""
    genes, lengths = {}, {}
    sample = ""
    for line in read_lines(board_path):
        if line.startswith(">"):
            fields = perl_split("\t", line)
            sample = s(at(fields, 0)).replace(">", "", 1).replace(" ", "")
            genes[sample] = "\t".join(s(at(fields, i)) for i in (2, 3, 4))
        elif WORD.search(line):
            lengths[sample] = len(line)
    return genes, lengths


def key(value):
    """Perl hash keys are strings (undef -> "")."""
    return "" if value is None else str(value)


def nested_add(tree, *path):
    """$tree{a}{b}...{z} = 1"""
    for step in path[:-1]:
        tree = tree.setdefault(step, {})
    tree[path[-1]] = 1


class Pair:
    """The link between two clusters (the Perl $cluster2cluster{one}{two})."""
    __slots__ = ("sum", "count", "members", "members_count")

    def __init__(self):
        self.sum = 0
        self.count = 0
        self.members = set()
        self.members_count = 0

    def add(self, score, sample_one, sample_two):
        self.sum += perl_num(score)
        self.count += 1
        if sample_one not in self.members and sample_two not in self.members:
            self.members_count += 1
            self.members.add(sample_one)
            self.members.add(sample_two)


def cluster_levels(per_score, per_score_extra, connectivity, pairs, threshold, ratio, lite):
    """The level loop (lines 93-333 of the Perl). Returns the last level."""
    threshold_n, ratio_n = perl_num(threshold), perl_num(ratio)
    level = 0
    while True:
        clusterindex = -1
        cluster2member = {}   # cluster -> score -> member -> entry
        member2cluster = {}   # member -> score -> {cluster: 1}
        cluster2cluster = {}  # cluster -> cluster -> Pair
        clustersize = {}

        def value(score, one, two, which):
            # the full entry at level 0; not needed at later levels (undef in the Perl)
            item = per_score[score][one][two]
            return item[which] if isinstance(item, dict) else None

        def new_cluster(score, one, two):
            nonlocal clusterindex
            clusterindex += 1
            ci = str(clusterindex)
            cluster2member.setdefault(ci, {}).setdefault(score, {})[one] = value(score, one, two, "one")
            cluster2member[ci][score][two] = value(score, one, two, "two")
            clustersize[ci] = clustersize.get(ci, 0) + 2
            nested_add(member2cluster, one, score, ci)
            nested_add(member2cluster, two, score, ci)
            return ci

        def clusters_of(sample):
            """The clusters of a member, in the order the Perl visits them."""
            for score1 in num_desc(member2cluster[sample]):
                for cluster in num_asc(member2cluster[sample][score1]):
                    yield cluster

        def link_to_new(member, ci, score, one, two):
            """The new cluster ci is linked with every cluster of member (including itself)."""
            for cluster in list(clusters_of(member)):
                cluster2cluster.setdefault(cluster, {}).setdefault(ci, Pair()).add(score, one, two)
                connectivity[cluster] = connectivity.get(cluster, 0) + 1
                connectivity[ci] = connectivity.get(ci, 0) + 1

        def joins_all(cluster, candidate):
            """Is candidate linked with every member of the cluster?"""
            for members in cluster2member[cluster].values():
                for member in members:
                    if candidate not in pairs.get(member, ()) and member not in pairs.get(candidate, ()):
                        return False
            return True

        def attach(member, joining, score, which, one, two):
            """Sequence joining enters a cluster of member (case "only one of the two has a cluster")."""
            if level <= 1:
                for cluster in list(clusters_of(member)):
                    if joins_all(cluster, joining):
                        cluster2member[cluster].setdefault(score, {})[joining] = value(score, one, two, which)
                        clustersize[cluster] += 1
                        nested_add(member2cluster, joining, score, cluster)
                        return
                ci = new_cluster(score, one, two)
                link_to_new(member, ci, score, one, two)
            else:
                for cluster in list(clusters_of(member)):
                    cluster2member[cluster].setdefault(score, {})[joining] = value(score, one, two, which)
                    clustersize[cluster] += 1
                    nested_add(member2cluster, joining, score, cluster)

        for score in num_desc(per_score_extra):
            by_conne_one = per_score_extra[score]
            for conne_one in num_desc(by_conne_one):
                for len_one in num_desc(by_conne_one[conne_one]):
                    for sample_one in sorted(by_conne_one[conne_one][len_one]):
                        by_conne_two = by_conne_one[conne_one][len_one][sample_one]
                        for conne_two in num_desc(by_conne_two):
                            for len_two in num_desc(by_conne_two[conne_two]):
                                for sample_two in sorted(by_conne_two[conne_two][len_two]):
                                    if sample_one == sample_two:
                                        continue
                                    has_one = sample_one in member2cluster
                                    has_two = sample_two in member2cluster
                                    if not has_one and not has_two:
                                        new_cluster(score, sample_one, sample_two)
                                    elif has_one and not has_two:
                                        attach(sample_one, sample_two, score, "two", sample_one, sample_two)
                                    elif not has_one and has_two:
                                        attach(sample_two, sample_one, score, "one", sample_one, sample_two)
                                    else:
                                        for onecluster in list(clusters_of(sample_one)):
                                            for twocluster in list(clusters_of(sample_two)):
                                                cluster2cluster.setdefault(onecluster, {}) \
                                                    .setdefault(twocluster, Pair()) \
                                                    .add(score, sample_one, sample_two)
                                                connectivity[onecluster] = connectivity.get(onecluster, 0) + 1
                                                connectivity[twocluster] = connectivity.get(twocluster, 0) + 1

        # print the clusters of this level
        for cluster in num_asc(cluster2member):
            if clustersize.get(cluster, 0) <= 1:
                continue
            cluster2print = "CLUSTER-%d-%04d" % (level, int(cluster))
            for score in num_desc(cluster2member[cluster]):
                for member in sorted(cluster2member[cluster][score]):
                    if level > 0:
                        member2print = "CLUSTER-%d-%04d" % (level - 1, int(member))
                    else:
                        member2print = '"%s"' % s(cluster2member[cluster][score][member])
                    lite.write("%s\t%s\t%.0f\n" % (cluster2print, member2print, perl_num(score)))

        # the clusters become the units of the next level
        per_score, per_score_extra = {}, {}
        for onecluster, partners in cluster2cluster.items():
            for twocluster, pair in partners.items():
                linked_enough = (pair.members_count / clustersize[onecluster] >= ratio_n
                                 and pair.members_count / clustersize[twocluster] >= ratio_n)
                if linked_enough or level >= 1:
                    score = "%.0f" % perl_divide(pair.sum, pair.count)
                    if perl_num(score) > threshold_n:
                        nested_add(per_score, score, onecluster, twocluster)
                        nested_add(per_score_extra, score,
                                   key(connectivity.get(onecluster)), key(clustersize.get(onecluster)), onecluster,
                                   key(connectivity.get(twocluster)), key(clustersize.get(twocluster)), twocluster)

        if clusterindex >= 1:
            level += 1
            continue
        return level


def write_table(lite_path, table_path, genes, annotation, maxlevel):
    """The final table: one line per sequence and cluster, with columns L0, L1, ..."""
    excelify = {}   # level -> node -> cluster -> score
    for line in read_lines(lite_path):
        line = line.replace('"', "")
        fields = perl_split("\t", line)
        if line.startswith("CLUSTER-0"):
            node_fields = perl_split("  ", s(at(fields, 1)))
            sample = s(at(node_fields, 0))
            mod_node = "\t".join([sample, s(at(node_fields, 1)), s(at(node_fields, 2)),
                                  s(annotation.get(sample))])
            cluster = s(at(fields, 0))
            maxlevel = s(at(perl_split("-", cluster), 1))
            excelify.setdefault(maxlevel, {}).setdefault(mod_node, {})[cluster] = at(fields, 2)
        elif re.search(r"CLUSTER-[1-9]+-", line):
            cluster1, cluster2 = s(at(fields, 0)), s(at(fields, 1))
            maxlevel = s(at(perl_split("-", cluster1), 1))
            excelify.setdefault(maxlevel, {}).setdefault(cluster2, {})[cluster1] = at(fields, 2)
    maxlevel = perl_num(maxlevel)

    alter2order = {}   # L0 -> composite score -> lines
    for node in sorted(excelify.get("0", {})):
        current = node
        alters = {node + "\t"}
        level = 0
        while level <= maxlevel:
            parents = excelify.get(str(level), {}).get(current)
            if parents is not None:
                new_alters = set()
                for parent in sorted(parents):
                    for alter in sorted(alters):
                        new_alters.add(alter + parent + "\t")
                    current = parent   # as in the Perl, only the last cluster is followed
                alters = new_alters
            level += 1

        for alter in alters:
            alter = alter[:-1]
            sample = s(at(perl_split("\t", alter), 0)).replace("(", r"\(").replace(")", r"\)")
            replacement = sample + "\t" + s(genes.get(sample))
            alter = re.sub(sample, lambda _m: replacement, alter, count=1)
            composcore = s(at(perl_split("\t", alter), 5))
            idepc = "%.2f" % (perl_num(composcore[0:6]) / 10000)
            simpc = "%.2f" % (perl_num(composcore[6:12]) / 10000)
            alter = re.sub(re.escape(composcore), lambda _m: f"{composcore}\t{idepc}\t{simpc}", alter, count=1)
            l0 = s(at(perl_split("\t", alter), 9))
            alter2order.setdefault(l0, {}).setdefault(composcore, set()).add(alter)

    with write_open(table_path) as table:
        table.write("# sequence.seq id\tsequence.V-GENE and allele\tsequence.D-GENE and allele\t"
                    "sequence.J-GENE and allele\tsequence.JUNCTION.aa seq\tcomposite score\t"
                    "identity\tsimilarity\tannotation\tL0\tL1\tL2\tL3\tL4\tL5\tetc\n")
        for l0 in nsort(alter2order):
            for composcore in num_desc(alter2order[l0]):
                for alter in sorted(alter2order[l0][composcore]):
                    table.write(alter + "\n")


def run(layout_path, board_path, threshold="0", cluster_link_ratio="0"):
    genes, lengths = read_board(board_path)

    per_score = {}       # score -> sample1 -> sample2 -> {"one": entry1, "two": entry2}
    connectivity = {}    # number of links of each sequence (and later of each cluster)
    pairs = {}           # sample1 -> set of sample2
    annotation = {}
    threshold_n = perl_num(threshold)
    for line in read_lines(layout_path):
        line = line.replace('"', "")
        fields = perl_split("\t", line)
        one, two = s(at(fields, 0)), s(at(fields, 1))
        sample_one = s(at(perl_split_ws(one), 0))
        sample_two = s(at(perl_split_ws(two), 0))
        array_one, array_two = perl_split("  ", one), perl_split("  ", two)
        annotation[sample_one] = array_one[-1] if array_one else None
        annotation[sample_two] = array_two[-1] if array_two else None
        score = s(at(fields, 2))
        if perl_num(score) >= threshold_n:
            per_score.setdefault(score, {}).setdefault(sample_one, {})[sample_two] = {"one": one, "two": two}
            if sample_two not in pairs.get(sample_one, ()):
                connectivity[sample_one] = connectivity.get(sample_one, 0) + 1
                connectivity[sample_two] = connectivity.get(sample_two, 0) + 1
            pairs.setdefault(sample_one, set()).add(sample_two)

    per_score_extra = {}
    for score in per_score:
        for sample_one in per_score[score]:
            for sample_two in per_score[score][sample_one]:
                nested_add(per_score_extra, score,
                           key(connectivity.get(sample_one)), key(lengths.get(sample_one)), sample_one,
                           key(connectivity.get(sample_two)), key(lengths.get(sample_two)), sample_two)

    base = f"{layout_path}.{threshold}.ratio{cluster_link_ratio}.lite.clusters"
    with write_open(base) as lite:
        maxlevel = cluster_levels(per_score, per_score_extra, connectivity, pairs,
                                  threshold, cluster_link_ratio, lite)
    write_table(base, base + ".table", genes, annotation, maxlevel)
    return base + ".table"


if __name__ == "__main__":
    utf8_output()
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    run(*sys.argv[1:5])
