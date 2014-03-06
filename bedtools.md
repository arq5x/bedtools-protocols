% bedtools Tutorial
% Aaron Quinlan
% November 22, 2013

BEDTools: arithmetic on the genome.

Abstract
========
Technological advances have enabled the use of DNA sequencing as a flexible tool to characterize genetic variation and to measure the activity of diverse cellular phenomena such as gene isoform expression and transcription factor binding. Extracting biological insight from the experiments enabled by these advances demands the analysis of large, multi-dimensional datasets. This unit describes the use of the BEDTools toolkit for the exploration of high-throughput genomics datasets. I present several protocols for common genomic analyses and demonstrate how multiple BEDTools operations may be combined to quickly address complex analysis tasks.


Introduction
============
Modern genomics research combines high throughput DNA sequencing technologies with computational and statistical analysis to gain insight into genome biology. Interpreting genomic datasets is complicated by the spectrum of experimental assays that are now possible (*cite Pachter*) and the scale of the data generated. Additional complexity comes from the fact that sever genomics data formats, such as BED, GFF, VCF, BAM, and BigWig are used by the research community. While these data formats differ in their intent and structure, they each describe the attributes of one or more genome intervals. A genome interval represents a consecutive stretch of nucleotides on a chromosome or assembly scaffold (**Figure**), and most genomic analyses can be distilled to "genome arithmetic": that is, the comparison of sets of many genome intervals. For example, quantifying transcript expression in RNA-seq experiments is fundamentally a process of counting the number of cDNA alignments (i.e., intervals in BAM format) that overlap transcript annotations (i.e., intervals in GFF or BED format).

BEDTools is an open source software package comprised of multiple tools for comparing and exploring genomic datasets via fundamental "genome arithmetic" tasks. While the individual tools in the BEDTools suite are each focused on a relatively simple task (e.g., interval intersection), complex analyses are possible by properly combining multiple tools. The goals of this unit are to introduce the basic concepts of genome arithmetic with BEDTools and to demonstrate, via biologically relevant examples, how analytical power is conferred through clever combinations of individual tools. I emphasize that this unit is intended to give new users a sense of what is possible with the BEDTools suite. Interested users are encouraged to read the BEDTools documentation (bedtools.readthedocs.org), as there are currently 39 subcommands, yet only a subset are covered in this unit.


Strategic Planning
==================
The protocols described in this unit require a computer with a UNIX, Linux, or Apple OSX operating system. Microsoft Windows users may also complete the unit if they first install Cygwin, but Windows usage is not directly supported. In the following sections, I will describe how to install BEDTools and other required software, as well as the basic usage concepts.

Conventions
-------------------
Throughout this unit, I will demonstrate BEDTools usage via 
commands issued on the UNIX command line. Such commands will use a different font as follows. Also, the "$" character is intended
to represent the command prompt and should not be typed.

~~~~ {#mycode .bash}
    $ bedtools --help
~~~~

Additionally, this unit will include commands in the R programming language, primarily as a means for creating plots
describing the results of BEDTools analyses. Such commands will use the same font, yet will be preceded by a ">" to denote the
R command prompt. The ">" should likewise not be typed.

~~~~ {#mycode .R}
    > x <- c(1,2,3)
~~~~


Background knowledge
--------------------
This unit assumes that the reader has previous experience working on the UNIX command line, as well as a basic understanding of common genomics file formats such as BED, VCF, GFF, and BAM. If not, I encourage you to first read the bedtools documentation (bedtools.readthedocs.org), as well as the manuscripts describing the above formats. There are also many freely available tutorials on the Internet that describe the basics of working on the UNIX command line.


Installing BEDTools
-------------------
BEDTools is freely available software that is archived and maintained on GitHub. The latest version of BEDTools can always be found at the following URL: https://github.com/arq5x/bedtools2/releases. At the time this unit as written, the latest release was 2.19.1. Future readers should check for subsequent releases and adjust the installation commands below accordingly.

**Necessary Resources**
- a C/C++ compiler such as GCC. For OSX users, this typically requires the installation of the Xcode developer tools.
- the zlib and zlib-devel compression libraries (installed by default on many systems)

1. Download and compile the latest version of bedtools. In this example, we are downloading version 2.19.1.

~~~~ {#mycode .bash}
    curl -OL https://github.com/arq5x/bedtools2/releases/download/v2.19.0/bedtools-2.19.0.tar.gz
    tar -zxvf bedtools-2.19.0.tar.gz
    cd bedtools2-2.19.0
    make
    
    # At this point, you need to make the bedtools executable 
    # accessible on your system. If you have administrator 
    # privileges, you can do this by copying the bedtools
    # executable to /usr/local/bin
    cp bin/bedtools /usr/local/bin/

    # if you do not have administrator privileges, you should
    # copy the executable to a "bin" directory located within
    # your home directory and then update your PATH to know
    # where bedtools can be found
    mkdir ~/bin
    cp bin/bedtools ~/bin
    PATH=PATH:~/bin
~~~~

2. Alternatively, bedtools is also available via package management software available on most UNIX systems.

~~~~ {#mycode .bash}
    # Fedora / Centos
    $ yum install bedtools

    # Debian / Ubuntu
    $ apt-get install bedtools

    # OSX (via HomeBrew)
    $ brew install bedtools
~~~~


Downloading datasets for this unit
----------------------------------
To do.
cpg.bed
exons.bed


The bedtools help
-----------------
The BEDTools software package is comprised of many subtools. One can be reminded of the tools available and a brief
summary of their functionality by typing the following on the command line:

~~~~ {#mycode .bash}
    $ bedtools --help
~~~~

If bedtools has been installed correctly, on your system, you should see several examples of the bedtools "subcommands:. If not, please refer back to step 1 above.

~~~~ {#mycode .bash}
    $ bedtools --help
    bedtools: flexible tools for genome arithmetic and DNA sequence analysis.
    usage:    bedtools <subcommand> [options]

    The bedtools sub-commands include:

    [ Genome arithmetic ]
        intersect     Find overlapping intervals in various ways.
        window        Find overlapping intervals within a window around an interval.
        closest       Find the closest, potentially non-overlapping interval.
        coverage      Compute the coverage over defined intervals.
        map           Apply a function to a column for each overlapping interval.
        genomecov     Compute the coverage over an entire genome.
        merge         Combine overlapping/nearby intervals into a single interval.
        cluster       Cluster (but do not merge) overlapping/nearby intervals.
        complement    Extract intervals _not_ represented by an interval file.
        subtract      Remove intervals based on overlaps b/w two files.
        slop          Adjust the size of intervals.
        flank         Create new intervals from the flanks of existing intervals.
        sort          Order the intervals in a file.
        random        Generate random intervals in a genome.
        shuffle       Randomly redistrubute intervals in a genome.
        sample        Sample random records from file using reservoir sampling.
        annotate      Annotate coverage of features from multiple files.
~~~~

In order to conduct an analysis of genomic intervals with BEDTools, one must employ one of the BEDTools subcommands. I will illustrate this with the `intersect` subcommand (more details will be provided in the first protocol). For example, to intersect BED files representing Alu elements and CpG islands, one would use the following command.


~~~~ {#mycode .bash}
    $ bedtools intersect -a alu.bed -b cpg.bed
~~~~

One may also request a detailed help menu regarding the specifics of each tool as follows.

~~~~ {#mycode .bash}
    $ bedtools intersect -h
~~~~

All other subcommands follow the same basic convention.

~~~~ {#mycode .bash}
    $ bedtools [SUBCOMMAND] [OPTIONS]
~~~~


Working with "genome-sorted" datasets
-------------------------------------
The default algorithm that BEDTools employes for detecting overlapping genomic intervals loads one file into an R-tree data structure. While fast, it can consume substantial memory, especially for very large files. For this reason, we provide an alternative, yet very fast and memory efficient algorithm that requires one's input files to be "genome-sorted": that is, sorted first by chromosome and then by start position. When both input files are genome-sorted, the algorithm can “sweep” through the data and detect overlaps on the fly in a manner much like the way database systems join two tables. This algorithm is invoked via the `-sorted` option and its use is demonstrated through the manuscript. The performance gains conferred through the use of "genome-sorted" datasets are illustrated in Figure 1.

![](https://bedtools.readthedocs.org/en/latest/_images/memory-comparo.png)


Genome files
------------
Some of the BEDTools subcommand need to know the size of each of the chromosomes from the organism with which you are working. These "genome" files must be tab delimited; the first column must be the chromosome label and the second column must be the length of the chromosome. For example, below is an example "genome" file for build 37 (a.k.a "hg19") of the human genome.

~~~~ {#mycode .bash}
    $ head -24 human.hg19.genome
    chr1    249250621
    chr2    243199373
    chr3    198022430
    chr4    191154276
    chr5    180915260
    chr6    171115067
    chr7    159138663
    chrX    155270560
    chr8    146364022
    chr9    141213431
    chr10   135534747
    chr11   135006516
    chr12   133851895
    chr13   115169878
    chr14   107349540
    chr15   102531392
    chr16   90354753
    chr17   81195210
    chr18   78077248
    chr20   63025520
    chrY    59373566
    chr19   59128983
    chr22   51304566
    chr21   48129895
~~~~



BP1: Intersecting genome feature files 
=====================================================

The `intersect` command is the workhorse of the BEDTools suite. It compares two BED/VCF/GFF/BAM files and reports all of the genomic where the features in the two files overlap (that is, share at least one base pair in common).

By default, `intersect` reports the subset of intervals that are common between your two files. The "A" file is considered the "query" file, whereas the "B" file is considered the "database" file. The BEDTools convention is to, for the most part, report results with respect the the "query" (A) file.

To demonstrate, let's identify all of the CpG islands that overlap exons.

~~~~ {#mycode .bash}
    bedtools intersect -a cpg.bed -b exons.bed | head -5
    chr1    29320   29370   CpG:_116
    chr1    135124  135563  CpG:_30
    chr1    327790  328229  CpG:_29
    chr1    327790  328229  CpG:_29
    chr1    327790  328229  CpG:_29
~~~~


AP1a: Reporting the original feature in each file.
--------------------------------------------------
The `-wa` (write A) and `-wb` (write B) options allow one to see the original records from the A and B files that overlapped.  As such, instead of solely showing you *where* the intersections occurred, these options show you the original intervals that intersected from the two files.

~~~~ {#mycode .bash}
    bedtools intersect -a cpg.bed -b exons.bed -wa -wb \
    | head -5
    chr1    28735   29810   CpG:_116    chr1    29320   29370   NR  _024540_exon_10_0_chr1_29321_r    0   -
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +
~~~~


AP1b: How many base pairs of overlap were there?
------------------------------------------------
The `-wo` (write overlap) option allows one to also report the *number* of base pairs of overlap between the features that overlap between each of the files.

~~~~ {#mycode .bash}
    bedtools intersect -a cpg.bed -b exons.bed -wo \
    | head -5
    chr1    28735   29810   CpG:_116    chr1    29320   29370   NR  _024540_exon_10_0_chr1_29321_r    0   -   50
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +   439
~~~~


AP1c: Counting the number of overlapping features.
--------------------------------------------------
We can also count, for each feature in the "A" file, the number of overlapping features in the "B" file. This is handled with the `-c` option.

~~~~ {#mycode .bash}
    bedtools intersect -a cpg.bed -b exons.bed -c \
    | head -5
    chr1    28735   29810   CpG:_116    1
    chr1    135124  135563  CpG:_30 1
    chr1    327790  328229  CpG:_29 3
    chr1    437151  438164  CpG:_84 0
    chr1    449273  450544  CpG:_99 0
~~~~

AP1d: Find features that DO NOT overlap
--------------------------------------------
Often we want to identify those features in our "A" file that **do not** overlap features in the B file. The `-v` option is your friend in this case.

~~~~ {#mycode .bash}
    bedtools intersect -a cpg.bed -b exons.bed -v \
    | head -5
    chr1    437151  438164  CpG:_84
    chr1    449273  450544  CpG:_99
    chr1    533219  534114  CpG:_94
    chr1    544738  546649  CpG:_171
    chr1    801975  802338  CpG:_24
~~~~


AP1e: Require a minimal fraction of overlap.
--------------------------------------------
Recall that the default is to report overlaps between features in A and B so long as *at least one basepair* of overlap exists. However, the `-f` option allows you to specify what fraction of each feature in A should be overlapped by a feature in B before it is reported.

Let's be more strict and require 50% of overlap.

~~~~ {#mycode .bash}
    bedtools intersect \
             -a cpg.bed \
             -b exons.bed \
             -wo \
             -f 0.50 \
    | head -5
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +   439
    chr1    788863  789211  CpG:_28 chr1    788770  794826  NR_ 047525_exon_4_0_chr1_788771_f    0   +   348
~~~~


BP2: Assessing coverage in DNA sequencing experiments 
=====================================================

http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html

~~~~ {#mycode .bash}
    # Download alignments for sample NA19146 in BAM format
    # from the 1000 Genomes Project.
    # Filesize: 
    # Estimate: 15-30 minute download
    KGFTP=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/NA19146/alignment
    
    curl -O $KGFTP/NA19146.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam

    # create a symbolic link for brevity
    ln -s NA19146.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam NA19146.bam
    
    # Calculate a histogram of coverage for each chromosome
    # as well as genome-wide.
    # Estimate: 30 minutes
    bedtools genomecov \
             -ibam NA19146.bam \
    > NA19146.coverage.hist.txt
~~~~

At this point, we can use R to plot the genome-wide distribution of coverage. 

~~~~ {#mycode .R}
    > cov = read.table('NA19146.coverage.hist.txt')
    > head(cov)
      V1 V2       V3        V4        V5
    1  1  0 29172804 249250621 0.1170420
    2  1  1  7196069 249250621 0.0288708
    3  1  2  9698769 249250621 0.0389117
    4  1  3 11608275 249250621 0.0465727
    5  1  4 12960686 249250621 0.0519986
    6  1  5 13769135 249250621 0.0552421

    # extract the genome-wide (i.e., no the per-chromosome)   histogram entries
    > gcov = cov[cov[,1] == 'genome',]
    > head(gcov)
              V1 V2        V3         V4        V5
    91947 genome  0 321856083 3137454505 0.1025850
    91948 genome  1 117285058 3137454505 0.0373822
    91949 genome  2 152949464 3137454505 0.0487495
    91950 genome  3 176290526 3137454505 0.0561890
    91951 genome  4 189358028 3137454505 0.0603540
    91952 genome  5 194473449 3137454505 0.0619845

    # plot a density function for the genome-wide coverage
    > plot(gcov[1:51,2], gcov[1:51,5], 
    +      type='h', col='darkgreen', lwd=3,
    +      xlab="Depth", ylab="Fraction of genome at depth",)
    > axis(1,at=c(1,5,10,15,20,25,30,35,40,45,50))
~~~~

This should result in a plot such as the following.

![](protocols/coverage/genomecov-pdf.png)

We can also turn this into a cumulative density function to assess what fraction of the genome is covered by more than 1, 10, 20, 50, etc. sequences.

~~~~ {#mycode .R}
    # Create a cumulative distribution from the "raw" hist 
    # (truncate at depth >=50)
    > gcov_cumul = 1 - cumsum(gcov[,5])

    # Create a plot of the CDF
    > plot(gcov[2:51,2], gcov_cumul[1:50], 
    +      col='darkgreen', type='l', lwd=3, 
    +      xlab="Depth", ylab="Fraction of genome >= depth", 
    +      ylim=c(0,1.0)
    + )
    > axis(1,at=c(1,5,10,15,20,25,30,35,40,45,50))
    > axis(2,at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
~~~~

This should result in a plot such as the following.

![](protocols/coverage/genomecov-cdf.png)


AP2a: Which genome regions had reduced or excessive coverage? 
====================================================================

~~~~ {#mycode .bash}
    # Calculate a histogram of coverage for each chromosome
    # as well as genome-wide.
    # Estimate: 30 minutes
    bedtools genomecov \
             -ibam NA19146.bam \
             -bga \
    > NA19146.coverage.bedg
~~~~

Since we now have a BEDGRAPH representing the coverage through the entire genome of NA19146, we can use simple
awk filters to extract intervals with insufficent coverage.
Based on the distribution shown above (Figure/Panel), we might define insufficent coverage as those regions with less than five aligned sequences. In the BEDGRAPH format,
the observed coverage for a given interval is represented in the fourth column.

~~~~ {#mycode .bash}
    head -10 NA19146.coverage.bedg
    1   0   9994    0
    1   9994    9996    1
    1   9996    9999    2
    1   9999    10000   4
    1   10000   10001   42
    1   10001   10002   82
    1   10002   10003   112
    1   10003   10004   145
    1   10004   10005   184
    1   10005   10006   233

    awk '$4 < 5' NA19146.coverage.bedg \
    | head -10
    1   0   9994    0
    1   9994    9996    1
    1   9996    9999    2
    1   9999    10000   4
    1   10525   10526   4
    1   10526   10534   3
    1   10534   10550   2
    1   10550   10576   1
    1   10576   10589   2
    1   10589   10617   1
~~~~

Similarly, we can identify regions with excessive coverage, which based on above distribution, we will
define as intervals with more than 25 aligned sequences.

~~~~ {#mycode .bash}
    awk '$4 > 25' NA19146.coverage.bedg \
    | head -10
    1   10000   10001   42
    1   10001   10002   82
    1   10002   10003   112
    1   10003   10004   145
    1   10004   10005   184
    1   10005   10006   233
    1   10006   10007   250
    1   10007   10008   258
    1   10008   10009   306
    1   10009   10010   313
~~~~

Notice that many of the intervals with excessive coverage
are adjacent to one another. In such cases, we may wish to
merge these adjacent bases into single, high-coverage intervals. This can easily be accomplished with the "merge" tool:

~~~~ {#mycode .bash}
    awk '$4 > 25' NA19146.coverage.bedg \
    | bedtools merge -i - \
    | head -10
    1   10000   10287
    1   10337   10465
    1   11740   11741
    1   11742   11822
    1   11834   13924
    1   13931   15266
    1   15305   15308
    1   15310   15311
    1   15317   15528
    1   15568   15576
~~~~

AP2b: Assessing coverage in exome capture experiments 
=====================================================

~~~~ {#mycode .bash}
    # Download exome capture alignments for sample
    # NA12891 in BAM format from the 1000 Genomes Project.
    KGFTP=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/phase3_EX_or_LC_only_alignment/data/NA12891/exome_alignment/

    curl -O $KGFTP/NA12891.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam

    # create a symbolic link for brevity
    ln -s NA12891.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam NA12891.exome.bam

    # Download the exome capture targets
    KGFTP=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference
    
    curl -O $KGFTP/exome_pull_down_targets/20130108.exome.targets.bed

    # Remove the "chr" from the chromosome labels.
    sed -e 's/chr//' 20130108.exome.targets.bed \
    > targets.numeric.chroms.bed

    # Compute the histogram of coverage for each target, as 
    # well as across all targets.
    bedtools coverage \
             -hist \
             -abam NA12891.exome.bam \
             -b 20130108.exome.targets.numeric.chroms.bed \
    > NA12891.exome.coverage.hist.txt
~~~~

We can now use R to assess the fraction of the exome
capture targets that were covered by 0,1,2...N sequences.

~~~~ {#mycode .R}
    > cov = read.table('NA12891.exome.coverage.hist.txt');

    > gcov = cov[cov[,1] == 'all',]
    
    # Create a cumulative distribution from the "raw" hist 
    # (truncate at depth >=1000)
    > gcov_cumul = 1 - cumsum(gcov[,5])

    # Create a plot of the CDF
    > plot(gcov[2:401,2], gcov_cumul[1:400],
        + col='darkred', type='l', lwd=2,
        + xlab="Depth", 
        + ylab="Fraction of capture target bases >= depth", 
        + ylim=c(0,1.0),
    + )
    # add gridlines
    > abline(v = 20, col = "gray60")
    > abline(v = 50, col = "gray60")
    > abline(v = 80, col = "gray60")
    > abline(v = 100, col = "gray60")
    > abline(h = 0.50, col = "gray60")
    > abline(h = 0.90, col = "gray60")

    > axis(1, at=c(20,50,80), labels=c(20,50,80))
    > axis(2, at=c(0.90), labels=c(0.90))
    > axis(2, at=c(0.50), labels=c(0.50))
~~~~

The resulting plot should look something like:

![](protocols/coverage/exomecov-cdf.png)


AP2c: Which target regions were not covered at all? 
=====================================================

Blah blah blah

~~~~ {#mycode .bash}
    bedtools genomecov \
             -ibam NA12891.exome.bam \
             -bga \
    | awk '$4 == 0' \
    > NA12891.uncovered.bedg

    bedtools intersect \
             -a targets.numeric.chroms.bed \
             -b NA12891.uncovered.bedg \
             -sorted \
    > unsequenced.exome.target.intervals.bed

    head unsequenced.exome.target.intervals.bed
    1   15903   15930
    1   16714   16719
    1   69551   69560
    1   129150  129152
    1   877795  877813
    1   878015  878017
    1   878068  878123
    1   896053  896180
    1   896634  896659
    1   899929  899939
~~~~

At this point, we can check these results by inspecting the aligned sequence coverage directly in IGV.
**FIGURE!!!**


AP2c: Coverage in genomic windows normalized by GC content.
===========================================================



BP3: Plot transcription factor occupancy surrounding the transcription start site.
==================================================================================

Goal: plot TF binding occupancy around transcription start sites (TSS) for Sp1 and a control.
Uses: bedtools, R, and free ENCODE data

First, we must create a BED file of the TSSs. To do this, we will query the UCSC Genome Browser and choose the correct TSS based on 
the transcript's strand.

~~~~ {.bash}
    # -N : no headers
    # -B : tab-delimted output
    # uniq to remove duplicate TSSs across tmultiple transcripts
    # grep -v "_" to remove unplaced contigs
    mysql --user genome \
          --host genome-mysql.cse.ucsc.edu \
          -N \
          -B \
          -D hg19 \
          -e  "SELECT chrom, txStart, txEnd, \
                      X.geneSymbol, 1, strand \
               FROM knownGene as K, kgXref as X \
               WHERE txStart != txEnd \
               AND X.kgID = K.name" \
    | awk 'BEGIN{OFS=FS="\t"} \
           { if ($6 == "+") \
             { print $1,$2,$2+1,$4,$5,$6 } \
             else if ($6 == "-") \
             { print $1,$3-1,$3,$4,$5,$6 } \
           }' \
    | sort -k1,1 -k2,2n \
    | uniq \
    | grep -v "_" \
    > tss.bed

    head -5 tss.bed
    chr1    11873   11874   DDX11L1 1   +
    chr1    16764   16765   WASH7P  1   -
    chr1    17750   17751   WASH7P  1   -
    chr1    18060   18061   WASH7P  1   -
    chr1    19758   19759   WASH7P  1   -
~~~~

Now, let's add 1000 bp upstream and downstream of each TSS. To do this, we use the bedtools "slop" command.

~~~~ {.bash}
    bedtools slop \
             -b 1000 \
             -i tss.bed \
             -g hg19.chromsizes \
    > tss.plusminus.1000bp.bed

    head -5 
    chr1    10873   12874   DDX11L1 1   +
    chr1    15764   17765   WASH7P  1   -
    chr1    16750   18751   WASH7P  1   -
    chr1    17060   19061   WASH7P  1   -
    chr1    18758   20759   WASH7P  1   -
~~~~

To provide greater resolution to the plot we will produce, let's
break up each 2000bp interval flanking each TSS into 5bp sub-wondows.
We can easily do this with the `makewindows` command.

~~~~ {.bash}
    # the tr statement makes the window number the 5th column
    # this will be used to summarize the coverage
    # observed at each of the 2000 bases flanking the TSS
    # across all TSSs.
    bedtools makewindows \
             -b tss.plusminus.1000bp.bed \
             -w 5 \
             -i srcwinnum \
    | sort -k1,1 -k2,2n \
    | tr "_" "\t" \
    > tss.plusminus.1000bp.5bp.windows.bed
~~~~

Download BigWig files for the Sp1 transcription factor and
a negative control.

~~~~ {.bash}
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsH1hescSp1Pcr1xRawRep1.bigWig
    
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsH1hescRxlchPcr1xRawRep1.bigWig
~~~~

Map the Sp1 transcription factor.

~~~~ {.bash}
    # -c 4 -o mean: get the mean of the coverage
    # -null 0: if no overlap with bigwig, set to zero
    bedtools map \
             -a tss.plusminus.1000bp.5bp.windows.bed \
             -b <(bigWigToBedGraph  wgEncodeHaibTfbsH1hescSp1Pcr1xRawRep1.bigWig stdout) \
             -c 4 \
             -o mean \
             -null 0 \
    > sp1.tss.window.coverage.bedg
~~~~

Map the negative control.

~~~~ {.bash}
    # -c 4 -o mean: get the mean of the coverage
    # -null 0: if no overlap with bigwig, set to zero
    bedtools map \
             -a tss.plusminus.1000bp.5bp.windows.bed \
             -b <(bigWigToBedGraph  wgEncodeHaibTfbsH1hescRxlchPcr1xRawRep1.bigWig stdout) \
             -c 4 \
             -o mean \
             -null 0 \
    > Rxl.tss.window.coverage.bedg
~~~~

Summarize the coverage observed around the TSS 
from both Sp1 and the negative control.

~~~~ {.bash}
    # sort by the window number
    # -t$'\t' to specify that TABS should 
    # be used as the delimiter
    sort -t$'\t' -k5,5n sp1.tss.window.coverage.bedg \
    | bedtools groupby \
               -i - \
               -g 5 \
               -c 6 \
               -o sum \
    > sp1.tss.window.counts.txt

    sort -t$'\t' -k5,5n Rxl.tss.window.coverage.bedg \
    | bedtools groupby \
               -i - \
               -g 5 \
               -c 6 \
               -o sum \
    > Rxl.tss.window.counts.txt
~~~~

Plot the coverage around the TSS for SP1 and the control.

~~~~ {.R}
    > sp1 <- read.table('sp1.tss.window.counts.txt')
    > rxl <- read.table('Rxl.tss.window.counts.txt')
    
    # Plot SP1 in red
    > plot(sp1[,1], sp1[,2], 
      +   col='darkred', 
      +   xaxt = "n", 
      +   xlab="Distance from TSS", 
      +   ylab="Depth")
    # Add the control in grey
    > points(rxl[,1], rxl[,2], col='darkgrey')
    
    # adjust labels based on distance to TSS.
    # recall that the window size is 5 base pairs
    > axis(1, at=seq(0,400,40), labels=seq(-1000,1000,200))
    
    # add a vertical line at TSS
    > abline(v = 200, col = "gray60", lwd=3, lty=3)
    
    # add a legend.
    > legend('topright', 
             c("SP1 Trans. factor","Reverse cross-link Control"),
             lty=1, 
             col=c('darkred', 'darkgrey'), 
             bty='n')
~~~~

The resulting figure should look like this:

![](protocols/map/tsscov.png)



BP4 : Comparing and exploring the relationships among many datasets.
===================================================================

Download the sample BED files I have provided.

    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/maurano.dnaseI.tgz

Now, we need to extract all of the 20 Dnase I hypersensitivity BED files from the "tarball" named
`maurano.dnaseI.tgz`.

    tar -zxvf maurano.dnaseI.tgz
    rm maurano.dnaseI.tgz


Your directory should now contain 20 BED files, which reflect Dnase I hypersensitivity sites measured in twenty different fetal tissue samples from the brain, heart, intestine, kidney, lung, muscle, skin, and stomach.

Example finding the regions in the fetal intestine samples.

~~~~ {.bash}
    bedtools multiinter -i fIntestine*.bed \
                        -header \
                        -names DS16559 DS16712 DS16822 DS17808 DS18495 \
    | head
    chrom   start   end num list    DS16559 DS16712 DS16822 DS17808 DS18495
    chr1    10148   10150   2   DS16712,DS17808 0   1   0   1   0
    chr1    10150   10151   3   DS16712,DS16822,DS17808 0   1   1   1   0
    chr1    10151   10284   4   DS16559,DS16712,DS16822,DS17808 1   1   1   1   0
    chr1    10284   10315   3   DS16559,DS16712,DS17808 1   1   0   1   0
    chr1    10315   10353   2   DS16712,DS17808 0   1   0   1   0
    chr1    237719  237721  1   DS16822 0   0   1   0   0
    chr1    237721  237728  3   DS16712,DS16822,DS17808 0   1   1   1   0
    chr1    237728  237783  4   DS16559,DS16712,DS16822,DS17808 1   1   1   1   0
    chr1    237783  237784  3   DS16559,DS16712,DS16822 1   1   1   0   0
~~~~

What portion of the DnaseI hypersensitivity sites are found in 1,2,3,...20 of the cell types assayed?

~~~~ {.bash}
    bedtools multiinter -i *.bed \
        | awk '{print $4"\t"$3-$2}' \
        | sort -k1,1n \
        | bedtools groupby -g 1 -c 2 -o sum \
    > dnase.occupancy.dist.txt

    cat dnase.occupancy.dist.txt
    1   172639699
    2   70626095
    3   51945770
    4   35992709
    5   27751090
    6   18118487
    7   13375483
    8   10759033
    9   8887535
    10  7421260
    11  6229434
    12  5512967
    13  4964581
    14  4499490
    15  4252603
    16  4181704
    17  4430412
    18  4892327
    19  6148254
        20  10996900
~~~~


Plot the distribution of occupancy.

~~~~ {.R}
    > dnase_occ <- read.table('dnase.occupancy.dist.txt')

    # Plot the fraction of bases found to
    # be DnaseI heypersensitive in 1, 2, 3, ... 20 cells assayed.
    > plot(dnase_occ[,1], dnase_occ[,2] / sum(dnase_occ[,2]), 'h', 
      +   col="darkred", 
      +   lwd=4, 
      +   xlab="No. of assayed cells", 
      +   ylab="Fraction of bases")
~~~~

The resulting figure should look like this:

![](protocols/multiple-datasets/dnaseI-occupancy.png)

This demonstrates that among these 20 assayed cells the majority of hypersensitive bases are exclusive to a single cell.


AP4a: What is the DnaseI hypersensitivity signal from each cell type at each interval? 
======================================================================================

Example finding the regions in the fetal intestine samples.

~~~~ {.bash}
    bedtools unionbed -i fIntestine*.bed \
                        -header \
                        -names DS16559 DS16712 DS16822 DS17808 DS18495 \
    | head
    chrom   start   end DS16559 DS16712 DS16822 DS17808 DS18495
    chr1    10148   10150   0   7.76326 0   12.6573 0
    chr1    10150   10151   0   7.76326 9.704   12.6573 0
    chr1    10151   10284   9.71568 7.76326 9.704   12.6573 0
    chr1    10284   10315   9.71568 7.76326 0   12.6573 0
    chr1    10315   10353   0   7.76326 0   12.6573 0
    chr1    237719  237721  0   0   7.36415 0   0
    chr1    237721  237728  0   11.4351 7.36415 7.88268 0
    chr1    237728  237783  8.57969 11.4351 7.36415 7.88268 0
    chr1    237783  237784  8.57969 11.4351 7.36415 0   0
~~~~


BP5 : Distance measures.
========================
To do.


\


BP6 : Measuring dataset similarity.
====================================

We will use the bedtools implementation of a Jaccard statistic to meaure the 
similarity of two datasets. Briefly, the Jaccard statistic measures the ratio 
of the number of *intersecting* base pairs to the *total* number of base 
pairs in the two sets.  As such, the score ranges from 0.0 to 1.
0; lower values reflect lower similarity, whereas higher values reflect 
higher similarity.

Let's walk through an example: we would expect the Dnase hypersensivity 
sites to be rather similar between two samples of the **same** fetal tissue 
type.  Let's test:

~~~~ {.bash}
    bedtools jaccard \
        -a fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed \
        -b fHeart-DS15839.hotspot.twopass.fdr0.05.merge.bed
    intersection    union   jaccard
    81269248    160493950   0.50637
~~~~

But what about the similarity of two **different** tissue types?

~~~~ {.bash}
    bedtools jaccard \
        -a fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed \
        -b fSkin_fibro_bicep_R-DS19745.hg19.hotspot.twopass.fdr0.05.merge.bed
    intersection    union   jaccard
    28076951    164197278   0.170995
~~~~

Hopefully this demonstrates how the Jaccard statistic can be used as a simple 
statistic to reduce the dimensionality of the comparison between two large (e.g., 
often containing thousands or millions of intervals) feature sets.

We are going to take this a bit further and use the Jaccard statistic to 
measure the similarity of all 20 tissue samples against all other 20 samples. We will 
use BASH script loops to compute a Jaccard statistic for the 400 (20*20) pairwise 
comparisons among the fetal tissue samples.

~~~~ {.bash}
    file_labels=`ls *.bed | sed -e 's/.hotspot.twopass.fdr0.05.merge.bed//g' \
                                       -e 's/.hg19//g'`
    echo name" "$file_labels >> pairwise_jaccard.txt
    for file1 in `ls *.bed`
    do
        # make reasonable file labels
        file1_short=`echo $file1 \
                    | sed -e 's/.hotspot.twopass.fdr0.05.merge.bed//g' \
                    -e 's/.hg19//g'`
        echo -n $file1_short >> pairwise_jaccard.txt

        for file2 in `ls *.bed`;
        do
            # compute the jaccard stat for these two files.
            jaccard=`bedtools jaccard \
                       -a $file1 \
                       -b $file2 \
                       -valueOnly`
            
            # report the jaccard stat for these two files
            echo -n " "$jaccard >> pairwise_jaccard.txt
        done
        echo >> pairwise_jaccard.txt
    done
~~~~
     
Since the Jaccard statistic serves as a measure of the similarity of two datasets,
we can use a simple heatmap to graphically convey the overall similarity of all 20 
DnaseI hypersensitivity patterns.

~~~~ {.R}
    # install RColorBrewer if missing
    if (!require("RColorBrewer")) {
        install.packages("RColorBrewer")
        library(RColorBrewer)
    }

    jaccard_table <- read.table('pairwise_jaccard.txt', header=TRUE)
    row.names(jaccard_table) <- jaccard_table$name
    jaccard_table <- jaccard_table[, -1]
    jaccard_matrix <- as.matrix(jaccard_table)
    heatmap.2(jaccard_matrix, 
              col=brewer.pal(9,"Blues"), 
              margins = c(14, 14),
              density.info = "none",
              lhei = c(2, 8),
              trace="none")
~~~~

![](protocols/multiple-datasets/heatmap.png)

