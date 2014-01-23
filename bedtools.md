% bedtools Tutorial
% Aaron Quinlan
% November 22, 2013

Abstract
========


Introduction
============

- put everything in a Github GIST

We therefore encourage you to read the bedtools [documentation](http://bedtools.readthedocs.org/en/latest/).

The bedtools help
-----------------
To bring up the help, just type

    bedtools

As you can see, there are multiple "subcommands" and for bedtools to
work you must tell it which subcommand you want to use. Examples:

    bedtools intersect
    bedtools merge
    bedtools subtract



Strategic Planning
==================
- setup
- basic usage
- conventions:
    black font: command line
    blue font: R console
- dependencies
    - bedtools
    - R
    - R packages?
    - samtools?
    - mysql
- sorting 
    -show chromsweep comparison plot
- genome files



BP1: Intersecting genome feature files 
=====================================================

The `intersect` command is the workhorse of the `bedtools` suite. It compares two BED/VCF/GFF files (or a BAM file and one of the aforementioned files) and identifies all the regions in the gemome where the features in the two files overlap (that is, share at least one base pair in common).

By default, `intersect` reports the intervals that represent overlaps between your two files.  To demonstrate, let's identify all of the CpG islands that overlap exons.

    bedtools intersect -a cpg.bed -b exons.bed | head -5
    chr1    29320   29370   CpG:_116
    chr1    135124  135563  CpG:_30
    chr1    327790  328229  CpG:_29
    chr1    327790  328229  CpG:_29
    chr1    327790  328229  CpG:_29

AP1a: Reporting the original feature in each file.
--------------------------------------------------
The `-wa` (write A) and `-wb` (write B) options allow one to see the original records from the A and B files that overlapped.  As such, instead of not only showing you *where* the intersections occurred, it shows you *what* intersected.

    bedtools intersect -a cpg.bed -b exons.bed -wa -wb \
    | head -5
    chr1    28735   29810   CpG:_116    chr1    29320   29370   NR  _024540_exon_10_0_chr1_29321_r    0   -
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +

AP1b: How many base pairs of overlap were there?
------------------------------------------------
The `-wo` (write overlap) option allows one to also report the *number* of base pairs of overlap between the features that overlap between each of the files.

    bedtools intersect -a cpg.bed -b exons.bed -wo \
    | head -5
    chr1    28735   29810   CpG:_116    chr1    29320   29370   NR  _024540_exon_10_0_chr1_29321_r    0   -   50
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +   439

AP1c: Counting the number of overlapping features.
--------------------------------------------------
We can also count, for each feature in the "A" file, the number of overlapping features in the "B" file. This is handled with the `-c` option.

    bedtools intersect -a cpg.bed -b exons.bed -c \
    | head -5
    chr1    28735   29810   CpG:_116    1
    chr1    135124  135563  CpG:_30 1
    chr1    327790  328229  CpG:_29 3
    chr1    437151  438164  CpG:_84 0
    chr1    449273  450544  CpG:_99 0

\

AP1d: Find features that DO NOT overlap
--------------------------------------------
Often we want to identify those features in our A file that **do not** overlap features in the B file. The `-v` option is your friend in this case.

    bedtools intersect -a cpg.bed -b exons.bed -v \
    | head -5
    chr1    437151  438164  CpG:_84
    chr1    449273  450544  CpG:_99
    chr1    533219  534114  CpG:_94
    chr1    544738  546649  CpG:_171
    chr1    801975  802338  CpG:_24

AP1e: Require a minimal fraction of overlap.
--------------------------------------------
Recall that the default is to report overlaps between features in A and B so long as *at least one basepair* of overlap exists. However, the `-f` option allows you to specify what fraction of each feature in A should be overlapped by a feature in B before it is reported.

Let's be more strict and require 50% of overlap.

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

AP1f: Require intersections with the same strand
------------------------------------------------

AP1g: Require intersections with opposite strands
-------------------------------------------------

AP1h: Working with sorted data.
-------------------------------------------------
Figure comparing bedtools to bedops.

\


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


BP3: Something with the map tool.
=====================================================

Goal: plot TF binding occupancy around TSS for Sp1 and a control.
Uses: bedtools, R, and free ENCODE data

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

Make windows.

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



BP3: Analysis of ChIP datasets.
=====================================================
To do.



BP4 : Comparing and exploring the relationships among many datasets.
===================================================================
To do.
multintersect
unionbedgraph
jaccard



BP5 : Distance measures.
========================
To do.



Sophistication through chaining multiple bedtools
=================================================
Analytical power in `bedtools` comes from the ability to "chain" together multiple tools in order to construct rather sophisicated analyses with very little programming - you just need **genome arithmetic**!  Have a look at the examples [here](http://bedtools.readthedocs.org/en/latest/content/advanced-usage.html).


\


Principal component analysis
=============================

We will use the bedtools implementation of a Jaccard statistic to meaure the similarity of two 
datasets. Briefly, the Jaccard statistic measures the ratio of the number of *intersecting* base 
pairs to the *total* number of base pairs in the two sets.  As such, the score ranges from 0.0 to 1.
0; lower values reflect lower similarity, whereas higher values reflect higher similarity.

Let's walk through an example: we would expect the Dnase hypersensivity sites to be rather similar 
between two samples of the **same** fetal tissue type.  Let's test:

    bedtools jaccard \
        -a fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed \
        -b fHeart-DS15839.hotspot.twopass.fdr0.05.merge.bed
    intersection    union   jaccard
    81269248    160493950   0.50637

But what about the similarity of two **different** tissue types?

    bedtools jaccard \
        -a fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed \
        -b fSkin_fibro_bicep_R-DS19745.hg19.hotspot.twopass.fdr0.05.merge.bed
    intersection    union   jaccard
    28076951    164197278   0.170995

Hopefully this demonstrates how the Jaccard statistic can be used as a simple statistic to reduce 
the dimensionality of the comparison between two large (e.g., often containing thousands or 
millions of intervals) feature sets.


\


A Jaccard statistic for all 400 pairwise comparisons.
------------------------------------------------------


We are going to take this a bit further and use the Jaccard statistic to measure the similarity of 
all 20 tissue samples against all other 20 samples.  Once we have a 20x20 matrix of similarities, 
we can use dimensionality reduction techniques such as hierarchical clustering or principal 
component analysis to detect higher order similarities among **all** of the datasets.


We will use GNU parallel to compute a Jaccard statistic for the 400 (20*20) pairwise comparisons 
among the fetal tissue samples.

But first, we need to install [GNU parallel](http://www.gnu.org/software/parallel/).

    brew install parallel

Next, we need to install a tiny script I wrote for this analysis.

    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/make-matrix.py


Now, we can use `parallel` to, you guessed it, compute the 400 pairwise Jaccard statistics in parallel using as many processors as you have available.

    parallel "bedtools jaccard -a {1} -b {2} \
             | awk 'NR>1' \
             | cut -f 3 \
             > {1}.{2}.jaccard" \
             ::: `ls *.merge.bed` ::: `ls *.merge.bed`

This command will create a single file containing the pairwise Jaccard measurements from all 400 tests.

    find . \
        | grep jaccard \
        | xargs grep "" \
        | sed -e s"/\.\///" \
        | perl -pi -e "s/.bed./.bed\t/" \
        | perl -pi -e "s/.jaccard:/\t/" \
        > pairwise.dnase.txt

A bit of cleanup to use more intelligible names for each of the samples.

    cat pairwise.dnase.txt \
    | sed -e 's/.hotspot.twopass.fdr0.05.merge.bed//g' \
    | sed -e 's/.hg19//g' \
    > pairwise.dnase.shortnames.txt
 
Now let's make a 20x20 matrix of the Jaccard statistic. This will allow the data to play nicely with R.

    awk 'NF==3' pairwise.dnase.shortnames.txt \
    | awk '$1 ~ /^f/ && $2 ~ /^f/' \
    | python make-matrix.py \
    > dnase.shortnames.distance.matrix
 
Let's also make a file of labels for each dataset so that we can label each dataset in our R plot.

    cut -f 1 dnase.shortnames.distance.matrix | cut -f 1 -d "-" | cut -f 1 -d "_" > labels.txt
 
Now start up R. (This assumes you have installed the `ggplot2` package).

    R

You should see something very similar to this:


    R version 2.15.1 (2012-06-22) -- "Roasted Marshmallows"
    Copyright (C) 2012 The R Foundation for Statistical Computing
    ISBN 3-900051-07-0
    Platform: x86_64-apple-darwin12.0.0 (64-bit)
    
    R is free software and comes with ABSOLUTELY NO WARRANTY.
    You are welcome to redistribute it under certain conditions.
    Type 'license()' or 'licence()' for distribution details.
    
      Natural language support but running in an English locale
    
    R is a collaborative project with many contributors.
    Type 'contributors()' for more information and
    'citation()' on how to cite R or R packages in publications.
    
    Type 'demo()' for some demos, 'help()' for on-line help, or
    'help.start()' for an HTML browser interface to help.
    Type 'q()' to quit R.
    
    >

No paste these commands into the R console:

    library(ggplot2)
    library(RColorBrewer)
    blues <- colorRampPalette(c('dark blue', 'light blue'))
    greens <- colorRampPalette(c('dark green', 'light green'))
    reds <- colorRampPalette(c('pink', 'dark red'))
     
    setwd("~/Desktop/bedtools-demo")
    x <- read.table('dnase.shortnames.distance.matrix')
    labels <- read.table('labels.txt')
    ngroups <- length(unique(labels))
    pca <- princomp(x)
    qplot(pca$scores[,1], pca$scores[,2], color=labels[,1],     geom="point", size=1) +
      scale_color_manual(values = c(blues(4), greens(5), reds(5))) 

You should see this:


Et voila. Note that PCA was used in this case as an example of what PCA does for the CSHL Adv. Seq. course. Heatmaps are a more informative visualization in this case.



