import subprocess
import pandas as pd
import os

# =============================================================================
#  BED file merging

def merge_bedfile(bedfile,merged,cores=4):
    f = open(merged, "wt")
    p1 = subprocess.Popen(["bedtools","merge","-i",bedfile,"-c",str(cores),"-o","collapse"], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["sort", "-V", "-k1,1", "-k2,2"], stdin=p1.stdout, stdout=f)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0]
    f.close()

# =============================================================================
#  ANNOTATIONS

def gc_content(reference, bedfile_merged, GC_annotation):
    with open(GC_annotation, "wt") as f:
        subprocess.call(["bedtools", "nuc", "-fi", reference, "-bed", bedfile_merged], stdout=f)

def mappability(bedfile_merged, kmer_align, mappability_annotation):
    with open(mappability_annotation, "wt") as f:
        subprocess.call(["bedtools", "coverage", "-a", bedfile_merged, "-b", kmer_align], stdout=f)

def merge_annotations(mappability_annotation, GC_annotation, annotations):
    pd.read_csv(mappability_annotation,sep='\t',header=None).iloc[:,[0,1,2,3,6,7]].join(pd.read_csv(GC_annotation,sep='\t').iloc[:,5]).to_csv(annotations,sep='\t',header=None,index=False)

# =============================================================================
#  COVERAGES

# WRITE TEST!!!
def count_pairs(bamfiles, bedfile, outdir):
    for bam in bamfiles:
        count_pair(bam,bedfile,outdir+os.path.basename(bam).split('.')[0]+'.rtbed')

def count_pair(bamfile, bedfile, rtbed):
    with open(rtbed, "wt") as f:
        subprocess.call(["bedtools", "multicov", "-p", "-q", "20", "-bams", bamfile, "-bed", bedfile], stdout=f)


def merge_rtbeds(bedfile, rtbed_paths, coverages):
    D = pd.read_csv(bedfile, sep='\t',header=None)
    for i in range(len(rtbed_paths)):
        D[5+i] = pd.read_csv(rtbed_paths[i], sep='\t', header=None)[4]
    D.columns=['chr','start','end','gene']+[os.path.basename(s).split('.')[0] for s in rtbed_paths]
    D.to_csv(coverages, sep='\t',index=False)



# =============================================================================
#  tests

def test():
    P = "/fast/work/users/vmay_m/development/clearCNV/test/testdata/"
    print("merge bed file test")
    merge_bedfile(
        bedfile=P+"bedfile.bed",
        merged=P+"merged.bed",
        cores=4)

    print("gc content test")
    gc_content(
        reference=P+"reference.fa",
        bedfile_merged=P+"merged.bed",
        GC_annotation=P+"GC_annotations.bed")

    print("mappability test")
    mappability(
        kmer_align=P+"mappability.eq.1.bed",
        bedfile_merged=P+"merged.bed",
        mappability_annotation=P+"mappability.bed")

    print("merge annotations test")
    merge_annotations(
        mappability_annotation=P+"mappability.bed",
        GC_annotation=P+"GC_annotations.bed",
        annotations=P+"annotations.bed")

    print("merge .rtbeds test")
    merge_rtbeds(
        bedfile=P+"merged.bed",
        rtbed_paths=[P+"sample01.rtbed",P+"sample02.rtbed",P+"sample03.rtbed"],
        coverages=P+"annotations.bed")

    print("count coverage test")
    count_pair(
        bedfile=P+"targets.bed",
        bamfile=P+"small.bam",
        rtbed=P+"small.rtbed")


# =============================================================================
# main
test()
