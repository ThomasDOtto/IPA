#!/bin/bash
# (c) Copyright: Genome Research Ltd. 2017
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
#
# For question please see the readme file.

set -e

##########
### Improve / finish up PacBio Plasmodium assemblies

assembly=$1
dir=$2
correctedReads=$3
name=$4
reference=$5
illuminaReads=$6

## path must be set 
#PATH=/nfs/pathogen003/tdo/Tools/ABACAS2/working:$PATH; export PATH

if [ -z "$illuminaReads" ] ; then 
echo "usage: IPA <assembly>  <Result directory> <Pacbio corrected_reads> <name of result> <reference for abacas> <root of illumina reads>";
exit 1;
fi;

##########
# 1. Do an abacas with overlap to split - idea: resolve trivial gaps, get the chromosome names
# 2. megablast -F F 
# 2a: delete contained contigs
# 2b: find overlaps 
# 1. Do an abacas with overlap to split - idea: resolve trivial gaps, get the chromosome names
# 3 map corrected reads
# 4 circulator if MT and API are in set
# 5 icorn

doAbacas1=0;
doAbacas2=1;

mkdir $dir
cd $dir
ln -s ../$assembly 00.assembly.fa 
#### 1. abacas
if [ "$doAbacas1" -eq 1 ] ; then
#reference=/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Feb2012/NoSubtelomers/Pf3D7_v3.nonSubtelomers.fasta
ABA_CHECK_OVERLAP=0; export ABA_CHECK_OVERLAP
abacas2.nonparallel.sh $reference 00.assembly.fa 1000 98
### need to break  and delete N's
fastaq scaffolds_to_contigs Genome.abacas.fasta 01.assembly.fa
else 
ln -s 00.assembly.fa 01.assembly.fa
fi

#### 2. megablast
## we want to exclude overlapping and 
formatdb -p F -i 01.assembly.fa
megablast -W 40 -F F -a 8 -m 8 -e 1e-80 -d 01.assembly.fa -i 01.assembly.fa | awk '$3>98 && $4>500 && $1 != $2' > comp.self1.blast
perl IPA.addLengthBlast.pl 01.assembly.fa 01.assembly.fa comp.self1.blast &> /dev/null

# 2a: delete contained contigs
# we want the query to be always the smaller one
awk '$3>99 && $4>500 && $13 < $14' comp.self1.blast.length | perl IPA.getOverlap.pl > 02.ListContained.txt

cat 02.ListContained.txt | awk '$4>90' | cut -f 1 > List.Contained.fofn

perl IPA.deleteContigs.pl List.Contained.fofn 01.assembly.fa 02.assembly.fa
### now delete the contigs from the blast
# also filter for just hits > 2kb > 99%
cat comp.self1.blast.length | awk '$3> 99 && $4 > 2000' | perl IPA.deleteEntryinBlast.pl List.Contained.fofn > Blast.merge.blast.length

# 2b: find overlaps 
### Get a bam files to do the merge
SMALT_PARAMETER=" -n 8 "; export SMALT_PARAMETER
cp -s ../$illuminaReads*fastq .
IPA.runSMALT.sh 02.assembly.fa 20 3 $illuminaReads\_1.fastq $illuminaReads\_2.fastq first 800 

perl findoverlaps_ver3.pl Blast.merge.blast.length first.bam 02.assembly.fa OUT
mv mergedseq_OUT.fasta 03.assembly.fa

#### 3. abacas mostly to rename the sequences
if [ "$doAbacas2" -eq 1 ] ; then
#reference=/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Feb2012/NoSubtelomers/Pf3D7_v3.nonSubtelomers.fasta
PATH=/nfs/pathogen003/tdo/Tools/ABACAS2/working:$PATH; export PATH
ABA_CHECK_OVERLAP=1; export ABA_CHECK_OVERLAP
abacas2.nonparallel.sh $reference 03.assembly.fa 1000 98
### need to break  and delete N's
fastaq trim_Ns_at_end Genome.abacas.fasta 03b.assembly.fa
else 
ln -s 03.assembly.fa 03b.assembly.fa
fi

### 4. Do icorn now
icorn2.serial.sh $illuminaReads 500 03b.assembly.fa 1 5
ln -s ICORN2.03b.assembly.fa.6  04.assembly.fa

### 5. map the corrected reads 
bwa index  04.assembly.fa
bwa mem -t 8 -x pacbio 04.assembly.fa ../$correctedReads  > Mapped.corrected.04.sam
awk '$3=="Pf_MT" || $3=="Pf_API"'  Mapped.corrected.04.sam | awk '{ print ">"$1"\n"$10  }' > ForCirc.reads.fasta
samtools faidx 04.assembly.fa Pf_MT Pf_API > ForCirc.Ref.fasta
circlator all ForCirc.Ref.fasta ForCirc.reads.fasta Out.Circ
# delete the plastids from the current version (04.assembly.fa)
echo "Pf_MT" > List.MT_API.fofn
echo "Pf_API" >> List.MT_API.fofn
perl IPA.deleteContigs.pl List.MT_API.fofn 04.assembly.fa 05.assembly.fa
### merge in good trust of mh12 and his master piece
cat Out.Circ/05.*a >> 05.assembly.fa
samtools faidx 04.assembly.fa
#doBAM.sh Mapped.corrected.04 04.assembly.fa.fai

#### 5. Rename the sequences
## need to adapted this on your genome
export name
cat 05.assembly.fa | perl -nle 'if (/>Pf(\S+)$/){ $n=$1; $n =~ s/_v3//g; $n=~ s/3D7//g; print ">Pf".$ENV{name}.$n } else { print }' |  sed 's/|quiver//g' > 06.assembly.fa
geecee 06.assembly.fa GC.txt
fasta2singleLine.pl  06.assembly.fa | perl -nle 'if (/>/){ print } else { $l="-"; $r="+"; $left=substr($_,0,1000); $right=substr($_,(length($_)-1000),1000);if ($left =~ /^.*CCCTAAACCCTAAACCCTAAA/i) {$l="L";} if ($right =~ /TTTAGGGTTTAGGGTTTAGGG.*$/i){ $r="R"} print "$l\t$r"}' > 06.TelomerCounts.txt
fasta2singleLine.pl  06.assembly.fa |  perl -nle 'if (/>/){ print } elsif (length($_)>400000) { $l="-"; $r="+"; $left=substr($_,0,1000); $right=substr($_,(length($_)-1000),1000);if ($left =~ /^.*CCCTAAACCCTAAACCCTAAA/i) {$l="L";} if ($right =~ /TTTAGGGTTTAGGGTTTAGGG.*$/i){ $r="R"} print "$l\t$r"}' > 06.TelomerCounts.large.txt
num1=$(egrep -c L 06.TelomerCounts.txt)
num2=$(egrep -c R 06.TelomerCounts.txt)
num=$(($num1+$num2))
echo "Amount of telomer repeats $num"
num1=$(egrep -c L 06.TelomerCounts.large.txt)
num2=$(egrep -c R 06.TelomerCounts.large.txt)
num=$(($num1+$num2))
echo "Amount of telomer repeats attached to 400kb $num"
num=$(egrep -c "L.*R"  06.TelomerCounts.large.txt)
echo "$num chromosome have both telomers attached"
#### 7. Run icorn


#Countcountigsize.pl 06.assembly.fa 5000 07.assembly.fa

cat 07.assembly.fa | perl IPA.renameAssembly.pl $name > $name.V1.fasta 
