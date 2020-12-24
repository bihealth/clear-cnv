#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, numpy as np


# In[2]:


overlap = 0.8


# In[3]:


def set_overlap(val):
    if val <= 1 and val >= 0:
        global overlap
        overlap = val


# In[4]:


class Target:
    def __init__(self, chr, start, end, gene):
        self.chr = str(chr)
        self.start = int(start)
        self.end = int(end)
        self.gene = gene

    def __str__(self):
        return "\t".join(map(str, [self.chr, self.start, self.end, self.gene]))

    def __hash__(self):
        return hash(str(self))

    def __ne__(self, other):
        return not self == other

    def __eq__(self, other):
        global overlap
        alpha = (self.end - self.start) * (1.0 - overlap)
        al, ar = self.start - alpha, self.start + alpha
        bl, br = self.end - alpha, self.end + alpha
        return (
            self.chr == other.chr
            and other.start > al
            and other.start < ar
            and other.end > bl
            and other.end < br
        )

    def __lt__(self, other):
        if self.chr == other.chr:
            return (
                self.start + (self.end - self.start) / 2
                < other.start + (other.end - other.start) / 2
            )
        elif self.chr.isdigit() and other.chr.isdigit():
            return int(self.chr) < int(other.chr)
        elif not self.chr.isdigit() and not other.chr.isdigit():
            return self.chr < other.chr
        elif self.chr.isdigit() and not other.chr.isdigit():
            return True
        elif not self.chr.isdigit() and other.chr.isdigit():
            return False

    def __gt__(self, other):
        return self == other or not (self < other)

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    def __add__(self, other):
        if self.chr != other.chr:
            raise ValueError("Chromosomes have to be the same of both targets.")
        newgene = ""
        if self.gene in other.gene:
            newgene = other.gene
        elif other.gene in self.gene:
            newgene = self.gene
        else:
            ",".join(set([self.gene, other.gene]))
        return Target(self.chr, min([self.start, other.start]), max(self.end, other.end), newgene)

    def to_list(self):
        return list(map(str, [self.chr, self.start, self.end, self.gene]))


# In[5]:


class CNV:
    def __init__(self, chr, start, end, gene, abb, score=1.0, size=None):
        Target.__init__(self, chr, start, end, gene)
        self.abb = abb
        self.size = size
        self.score = float(score)

    def __str__(self):
        if self.size:
            return "\t".join(
                map(
                    str,
                    [self.chr, self.start, self.end, self.gene, self.abb, self.score, self.size],
                )
            )
        return "\t".join(
            map(str, [self.chr, self.start, self.end, self.gene, self.abb, self.score])
        )

    def __hash__(self):
        return hash("\t".join(str(self).split("\t")[:-1]))

    def __eq__(self, other):
        global overlap
        alpha = (self.end - self.start) * (1.0 - overlap)
        al, ar = self.start - alpha, self.start + alpha
        bl, br = self.end - alpha, self.end + alpha
        return (
            self.chr == other.chr
            and other.start > al
            and other.start < ar
            and other.end > bl
            and other.end < br
            and self.abb == other.abb
        )

    def __ne__(self, other):
        return not self == other

    def to_list(self):
        if self.size:
            return list(
                map(
                    str,
                    [self.chr, self.start, self.end, self.gene, self.abb, self.score, self.size],
                )
            )
        return list(map(str, [self.chr, self.start, self.end, self.gene, self.abb, self.score]))

    def get_size(self):
        return self.end - self.start


# In[6]:


class Hit:
    def __init__(self, sample, CNV):
        self.sample = sample
        self.cnv = CNV

    def __str__(self):
        return str(self.sample) + "\t" + str(self.cnv)

    def __hash__(self):
        return hash("\t".join(str(self).split("\t")[:-1]))

    def __eq__(self, other):
        return self.sample == other.sample and self.cnv == other.cnv

    def __ne__(self, other):
        return not self == other

    def to_list(self):
        return [str(self.sample), *self.cnv.to_list()]

    def get_size(self):
        return self.cnv.get_size()


# In[7]:


class Sample:
    def __init__(self, name, CNVs):
        self.name = str(name)
        self.cnvs = CNVs

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return self.name

    def __hash__(self):
        return hash(str(self))

    def has_cnv(self, cnv):
        return cnv in self.cnvs

    def get_hits(self):
        return [Hit(self, cnv) for cnv in self.cnvs]


# In[8]:


def load_cnvs_from_df(df):
    samples = []
    header = ["sample", "chr", "start", "end", "gene", "aberration", "score", "size"]
    df.columns = header
    for sample in set(df["sample"]):
        sub_df = df[df["sample"] == sample]
        cnvs = []
        for target in sub_df.index:
            cnvs.append(CNV(*sub_df.loc[target, header[1:]]))
        samples.append(Sample(sample, cnvs))
    return samples


# In[9]:


def union_hits(hits):
    if len(hits) <= 1:
        return hits[0]
    return set(hits[0]).union(union_hits(hits[1:]))


def difference_left(hits, hits_other):
    return set(hits).difference(set(difference_left(hits_other)))


def difference_hits(hits):
    if len(hits) <= 1:
        return hits[0]
    return set(hits[0]).difference(difference_hits(hits[1:]))


def intersect_hits(hits):
    if len(hits) <= 1:
        return hits[0]
    return set(hits[0]).intersection(intersect_hits(hits[1:]))


# In[11]:


def find_samples(hit, hits):
    X = []
    for a in hits:
        if hit.cnv == a.cnv and hit.sample != a.sample:
            X.append(str(a.sample))
    return X


# In[12]:


def read_sample(file, sep="\t"):
    cnvs = []
    with open(file) as f:
        for line in f:
            if sum([s.isdigit() for s in line.split(sep)]) == 0:
                continue
            cnvs.append(CNV(*(line.split(sep))))
    return Sample(os.path.basename(file).split(".")[0], cnvs)


def read_bedfile(file, sep="\t"):
    targets = []
    with open(file) as f:
        for line in f:
            if sum([s.isdigit() for s in line.split(sep)]) == 0:
                continue
            targets.append(Target(*(line.split(sep))))
    return targets


# In[13]:


def overlapping(c1, c2):
    return c1.chr == c2.chr and (min([c1.end, c2.end]) - (max([c2.start, c1.start])) > 0)


def is_hit_unique(hit, hits):
    for h in hits:
        if overlapping(hit.cnv, h.cnv) and hit.sample.name == h.sample.name:
            return False
    return True


def hitsA_not_in_hitsB(hitsA, hitsB):
    R = []
    for h in hitsA:
        if is_hit_unique(h, hitsB):
            R.append(h)
    return R


# bed needs to be formatted: chr, start, end, gene
def num_targets(bed, target, tolerance=10):
    return bed[
        (bed["chr"].astype(str) == str(target.chr))
        & (bed["start"] + tolerance >= int(target.start))
        & (bed["end"] - tolerance <= int(target.end))
    ].shape[0]


def distance(t0, t1):
    if t0.chr == t1.chr:
        return abs((t0.end + t1.start) - (t1.end - t1.start))
    else:
        return 0


# In[10]:


def flatten_list(L):
    return [val for l in L for val in l]
