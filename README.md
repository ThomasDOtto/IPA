# IPA
Script to improve long read (pacbio) assemblies

This wrapper script was written to cleanup PacBio assemblies produced with HGAP v3. 
We used in to assemble the P. cynomolgi assembly, many P. falciparum genomes. 

It was also test for Leishmania, where Karen (thanks!) compared her merges with the output 
of the script, confirming all the edits.

Following step are done:
Ignore contigs smaller than 5kb.
Delete contigs if they are contained with more than 90% of there overlap (IPA.getOverlap.pl). 
Merge contigs if the Illumina reads coverage is 50% of the median coverage. (findoverlaps_ver3.pl - Author Susanne Reimering)
Abacas2 to order the contigs (at least rename them)
icorn2 to correct homopolymer track from the PacBio assembly


## Requirements
You need following package installed and in the path
abacas2 - https://github.com/sanger-pathogens/ABACAS2
icorn2 - icorn.sourceforge.net/ : set ICORN2_THREADS=8; export ICORN2_THREADS;
circlator - https://github.com/sanger-pathogens/circlator
fastaq - https://github.com/sanger-pathogens/Fastaq
smalt - http://www.sanger.ac.uk/science/tools/smalt-0
samtools - you know where...
bwa - last version with the bwa mem
blastall - ncbi distribution that has megablast...

## Comments
This scripts works well on our computer farm. I adapted it for general use. Certain steps,
like the circularization are programmed for Plasmodium. 
If you have problems to run the script please let me know.
