#!/bin/bash
# Authors: José Luis Ruiz & Thomas Dan Otto  (APPROPRIATE ORDER)
# License: GNU General Public License v3.0 (APPROPRIATE? I REMOVED GENOME RESEARCH STATEMENTS? ALSO DO THE SAME IN THE HEADER OF THE REST OF FILES?)
# https://github.com/ThomasDOtto/IPA



########## IPA
# Pipeline for  Improving / finish up PacBio Assemblies

##### Arguments / Variables
assembly=$1
dir=$2
correctedReads=$3
name=$4
reference=$5
illuminaReads=$6
cores=$7
seqs_circl_1=$8
seqs_circl_2=$9
number_iterations_icorn=${10}
contigs_threshold_size=${11}
InsertsizeRange=${12}

if [ "$#" -ne 12 ]; then
echo -e "Not all parameters received as input. Proceed with care and check the 'Arguments / Variables' section in the pipeline if needed. All parameters are positional. If you want to skip a parameter provide an empty string. If not provided, some parameters will be set to default (check Arguments / Variables' section in the pipeline)\n"
echo -e "Example usage: IPA <assembly> <Result directory> <Pacbio corrected_reads> <name of result> <reference for ABACAS2> <root of illumina reads> <number of cores> <first sequence name to circularize> <second sequence name to circularize> <number of iterations for iCORN2> <Size threshold for discarding contigs> <Insert size range for Illumina short reads>\n"
echo "Be aware many perl scripts wihin IPA folder are used, and you may need to manually change the interpreter in the shebang. Many software and dependencies are also required and will be automatically checked. The pipeline will exit if any required software is not found in the PATH. Smalt has to be updated enough to accept gzipped reads"
# exit 1;
fi

if [ -z "$number_iterations_icorn" ] ; then
number_iterations_icorn=3
echo "Number of icorn2 iterations: "$number_iterations_icorn
fi

if [ -z "$cores" ] ; then
cores=40
echo "Number of threads: "$cores
fi

if [ -z "$seqs_circl_1" ] ; then
seqs_circl_1="MT|M76611"
echo "Seqs to circularize: "$seqs_circl_1
fi

if [ -z "$seqs_circl_2" ] ; then
seqs_circl_2="API"
echo "Seqs to circularize: "$seqs_circl_2
fi

if [ -z "$reference" ] ; then
echo "No reference genome provided. ABACAS2 will be skipped ..."
doAbacas2=0
else
doAbacas2=1
fi

if [ -z "$contigs_threshold_size" ] ; then
contigs_threshold_size=5000
echo "Length threshold to discard contigs (bp): "$contigs_threshold_size
fi

if [ -z "$InsertsizeRange" ] ; then
InsertsizeRange=800
echo "Insert size range Illumina short reads (bp): "$InsertsizeRange
fi

if [[ $assembly == /* ]] ; then
echo -e "\n"
else
assembly=$PWD/$assembly
fi

if [[ $dir == /* ]] ; then
echo -e "\n"
else
dir=$PWD/$dir
fi

if [[ $correctedReads == /* ]] ; then
echo -e "\n"
else
correctedReads=$PWD/$correctedReads
fi

if [[ $illuminaReads == /* ]] ; then
echo -e "\n"
else
  if [ -z "$illuminaReads" ]; then
    echo "Illumina reads not detected. Some steps of the pipeline will not be executed"
  else
    illuminaReads=$PWD/$illuminaReads
  fi
fi

if [ -f $illuminaReads\_1.fastq ]; then
    echo "Be aware of the naming required for the Illumina reads: _1.fastq.gz and _2.fastq.gz"
else
    echo -e "Be aware of the naming required for the Illumina reads: _1.fastq.gz and _2.fastq.gz\n"
    echo "IPA IS NOT DETECTING THE ILLUMINA READS FILES. If you want to use Illumina reads for polishing, please check naming, paths and that the files do exist. Otherwise ignore this if you don't want to use Illumina reads"
fi

# Check installed/available programs:
# PATH with paths to all tools must be properly set: export PATH=$PATH:...
# Check IPA .sh files to make sure which is the required software, but quick checking:
type formatdb >/dev/null 2>&1 || { echo >&2 "I require formatdb but it's not installed/available. Aborting."; exit 1; }
type megablast >/dev/null 2>&1 || { echo >&2 "I require megablast but it's not installed/available. Aborting."; exit 1; }
type samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed/available. Aborting."; exit 1; }
type smalt >/dev/null 2>&1 || { echo >&2 "I require smalt but it's not installed/available. Aborting."; exit 1; }
type abacas2.nonparallel.sh >/dev/null 2>&1 || { echo >&2 "I require ABACAS2 and dependencies but something is not installed/available. Aborting."; exit 1; }
type fastaq >/dev/null 2>&1 || { echo >&2 "I require Fastaq but it's not installed/available. Aborting."; exit 1; }
type icorn2.serial_bowtie2.sh >/dev/null 2>&1 || { echo >&2 "I require iCORN2 and dependencies but something is not installed/available. Aborting."; exit 1; }
type bowtie2 >/dev/null 2>&1 || { echo >&2 "I require Bowtie2 but it's not installed/available. Aborting."; exit 1; }
type java >/dev/null 2>&1 || { echo >&2 "I require java, particularly v1.7.0 for iCORN2, but it's not installed/available. Aborting."; exit 1; }
type bwa >/dev/null 2>&1 || { echo >&2 "I require BWA but it's not installed/available. Aborting."; exit 1; }
type circlator >/dev/null 2>&1 || { echo >&2 "I require Circlator but it's not installed/available. Aborting."; exit 1; }
type gtool.py >/dev/null 2>&1 || { echo >&2 "I require gtool.py for getting GC content but it's not installed/available. Aborting."; exit 1; }
type assembly-stats >/dev/null 2>&1 || { echo >&2 "I require assembly-stats but it's not installed/available. Aborting."; exit 1; }
type fasta2singleLine.pl >/dev/null 2>&1 || { echo >&2 "I require fasta2singleLine.pl but it's not installed/available. Aborting."; exit 1; }
# type fastq_info.sh >/dev/null 2>&1 || { echo >&2 "I require fastq_info.sh but it's not installed/available. Aborting."; exit 1; }
# type fasta_to_fastq.pl >/dev/null 2>&1 || { echo >&2 "I require fasta_to_fastq.pl but it's not installed/available. Aborting."; exit 1; }



##### STEPS:
# 1.  Discard smaller contigs
# 2.  megablast -F F
# 2a. Delete contained contigs
# 2b. Find overlaps
# 3.  ABACAS2 with overlap checking (resolve trivial gaps, reorder contigs, get the chromosome names, rename the sequences...)
# 4.  iCORN2 - error correction
# 5.  Circularization or organelles or any sequence reqired
# 6.  Rename sequences, evaluate the assemblies, get telomere counts, stats...
#### 1. Discard contigs smaller than a threshold:
mkdir -p $dir/1.Filtering; cd $dir/1.Filtering
echo -e "\nSTEP 1: Starting...\n"
echo -e "### Excluded contigs based on length threshold: (removesmalls.pl)" > Excluded.contigs.fofn
perl -S removesmalls.pl $contigs_threshold_size $assembly | sed 's/|/_/g' > 01.assembly.fa
mv Excluded.contigs.fofn ../Excluded.contigs.fofn
echo -e "\nSTEP 1: DONE\n"

#### 2. megablast
echo -e "\nSTEP 2: Starting...\n"
mkdir -p $dir/2.MegaBLAST; cd $dir/2.MegaBLAST
formatdb -p F -i $dir/1.Filtering/01.assembly.fa
megablast -W 40 -F F -a $cores -m 8 -e 1e-80 -d $dir/1.Filtering/01.assembly.fa -i $dir/1.Filtering/01.assembly.fa | awk '$3>98 && $4>500 && $1 != $2' > comp.self1.blast
IPA.addLengthBlast.pl $dir/1.Filtering/01.assembly.fa $dir/1.Filtering/01.assembly.fa comp.self1.blast &> /dev/null

# 2a. delete contained contigs
# We want the query to be always the smaller one
echo -e "\nSTEP 2a: Starting...\n"
awk '$3>99 && $4>500 && $13 < $14' comp.self1.blast.length | IPA.getOverlap.pl > 02.ListContained.txt
cat 02.ListContained.txt | awk '$4>90' | cut -f 1 > List.Contained.fofn
echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Excluded contigs that are contained in others: (MegaBLAST/IPA.getOverlap.pl/IPA.deleteContigs.pl)" >> ../Excluded.contigs.fofn
cat List.Contained.fofn >> ../Excluded.contigs.fofn
# Now delete the contigs from the blast, also filter for just hits > 2kb > 99%
IPA.deleteContigs.pl List.Contained.fofn $dir/1.Filtering/01.assembly.fa 02.assembly.fa
cat comp.self1.blast.length | awk '$3> 99 && $4 > 2000' | IPA.deleteEntryinBlast.pl List.Contained.fofn > Blast.merge.blast.length

# 2b. find overlaps
# Get a bam files to do the merge
echo -e "\nSTEP 2b: Starting...\n"
SMALT_PARAMETER="-n "$cores; export SMALT_PARAMETER
if [ -z "$illuminaReads" ]; then
echo -e "Illumina reads not provided and no filtering findoverlaps_ver3.pl\n"
cp 02.assembly.fa 03.assembly.fa
else
IPA.runSMALT_V2.sh 02.assembly.fa 20 3 $illuminaReads\_1.fastq.gz $illuminaReads\_2.fastq.gz first $InsertsizeRange $cores
findoverlaps_ver3.pl Blast.merge.blast.length first.bam 02.assembly.fa OUT
echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Contigs not covered: (findoverlaps_ver3.pl)" >> ../Excluded.contigs.fofn
cat notcovered_OUT.fasta | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' | grep ">" >> ../Excluded.contigs.fofn
echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### New contigs merging overlapping contigs(findoverlaps_ver3.pl):" >> ../Excluded.contigs.fofn
if [  $(cat log_OUT.txt | grep ^sequences: | awk '{ gsub("sequences: " , "") ; print }') -eq 0 ] ; then
echo "No filtering findoverlaps_ver3.pl"; else
cat results_OUT.txt >> ../Excluded.contigs.fofn; fi
mv mergedseq_OUT.fasta 03.assembly.fa
fi
echo -e "\nSTEP 2: DONE\n"

#### 3. ABACAS2
echo -e "\nSTEP 3: Starting...\n"
mkdir -p $dir/3.ABACAS2; cd $dir/3.ABACAS2
if [ "$doAbacas2" -eq 1 ] ; then
echo "ABACAS2 parameters are: ABA_CHECK_OVERLAP=0, Min_Alignment_Length=1000, Identity_Cutoff=98. Please check ABACAS2 help and change manually within the pipeline (section 3) these parameters if needed"
abacas2.nonparallel.sh
ABA_CHECK_OVERLAP=0; export ABA_CHECK_OVERLAP; Min_Alignment_Length=1000; Identity_Cutoff=98
abacas2.nonparallel.sh $reference $dir/2.MegaBLAST/03.assembly.fa $Min_Alignment_Length $Identity_Cutoff
# Break  and delete N's
fastaq trim_Ns_at_end Genome.abacas.fasta 03b.assembly.fa
# Extract info
# echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Overlapping contigs_2 (ABACAS2)" >> ../Excluded.contigs.fofn # Uncomment if ABA_CHECK_OVERLAP=1
# zcat *.overlap_report.gz >> ../Excluded.contigs.fofn # Uncomment if ABA_CHECK_OVERLAP=1
echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Final sequences not mapped to the reference: (ABACAS2)" >> ../Excluded.contigs.fofn
cat Res.abacasBin.fna | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' | grep ">" >> ../Excluded.contigs.fofn
echo -e "\n" >> ../Excluded.contigs.fofn; echo -e "### Final sequences mapped to the reference: (ABACAS2)" >> ../Excluded.contigs.fofn
for i in $(ls | grep .contigs.gff); do
echo "# "$(echo $i | sed 's/^[^.]*.\([^.]*\)..*/\1/')":" >> ../Excluded.contigs.fofn
cat $i | grep contig | awk '{ print $9 }' | sed 's/^[^"]*"\([^"]*\)".*/\1/' | sed 's/.*=//' >> ../Excluded.contigs.fofn
done; else
ln -s 03.assembly.fa 03b.assembly.fa; fi
echo -e "\nSTEP 3: DONE\n"

#### 4. iCORN2
echo -e "\nSTEP 4: Starting...\n"
mkdir -p $dir/4.iCORN2; cd $dir/4.iCORN2
if [ -z "$illuminaReads" ]; then
echo -e "\n"
else
ln -s $dir/3.ABACAS2/03b.assembly.fa 04.assembly.fa
echo -e "Be aware a particular version of java, 1.7.0, is required for iCORN2 and any error may be due to incorrect java version\n"
icorn2.serial_bowtie2.sh $illuminaReads 500 $dir/3.ABACAS2/03b.assembly.fa 1 $number_iterations_icorn
ln -s ICORN2.03b.assembly.fa.[$(expr $number_iterations_icorn + 1)] 04.assembly.fa
mkdir -p $dir/6.stats; icorn2.collectResults.pl $PWD > ../6.stats/iCORN2.final_corrections.results.txt
echo -e "### Total SNP: " >> ../6.stats/iCORN2.final_corrections.results.txt; cat ../6.stats/iCORN2.final_corrections.results.txt | awk '{sum+=$6;} END{print sum;}' >> ../6.stats/iCORN2.final_corrections.results.txt
echo -e "### Total INS: " >> ../6.stats/iCORN2.final_corrections.results.txt; cat ../6.stats/iCORN2.final_corrections.results.txt | awk '{sum+=$7;} END{print sum;}' >> ../6.stats/iCORN2.final_corrections.results.txt
echo -e "### Total DEL: " >> ../6.stats/iCORN2.final_corrections.results.txt; cat ../6.stats/iCORN2.final_corrections.results.txt | awk '{sum+=$8;} END{print sum;}' >> ../6.stats/iCORN2.final_corrections.results.txt
echo -e "### Total HETERO: " >> ../6.stats/iCORN2.final_corrections.results.txt; cat ../6.stats/iCORN2.final_corrections.results.txt | awk '{sum+=$9;} END{print sum;}' >> ../6.stats/iCORN2.final_corrections.results.txt
fi
echo -e "\nSTEP 4: DONE\n"

#### 5. Circlator for organelles
echo -e "\nSTEP 5: Starting...\n"
mkdir -p $dir/5.Circlator; cd $dir/5.Circlator
# Map the corrected reads
bwa index $dir/4.iCORN2/04.assembly.fa; bwa mem -t $cores -x pacbio $dir/4.iCORN2/04.assembly.fa $correctedReads > Mapped.corrected.04.sam
# Circlator:
if grep -q -E "$seqs_circl_1|$seqs_circl_2" Mapped.corrected.04.sam; then
seq_ids=$(awk '{print $3}' Mapped.corrected.04.sam | grep -E "$seqs_circl_1|$seqs_circl_2" | tail -n +2 | sort | uniq)
for i in $seq_ids; do cat Mapped.corrected.04.sam | grep $i | awk '{ print ">"$1"\n"$10 }' >> ForCirc.reads.fasta; done
for i in $seq_ids; do samtools faidx $dir/4.iCORN2/04.assembly.fa $i >> ForCirc.Ref.fasta; done
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' ForCirc.reads.fasta > ForCirc.reads_2.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' ForCirc.Ref.fasta > ForCirc.Ref_2.fasta
circlator all ForCirc.Ref_2.fasta ForCirc.reads_2.fasta Out.Circ
# Delete the plastids/organelles/circular sequences from the current version (04.assembly.fa)
for i in $seq_ids; do echo $i >> List.circular_sequences.fofn; done
IPA.deleteContigs.pl List.circular_sequences.fofn $dir/4.iCORN2/04.assembly.fa 05.assembly.fa
cat Out.Circ/06.fixstart.fasta >> 05.assembly.fa; else
ln -s $dir/4.iCORN2/04.assembly.fa 05.assembly.fa; fi
echo -e "\nSTEP 5: DONE\n"

#### 6. Rename sequences, evaluate the assemblies, get telomere counts, gc stats...
echo -e "\nSTEP 6: Starting...\n"
mkdir -p $dir/6.stats; cd $dir/6.stats; export name
{ IFS= read -r telomere_seq_1 && IFS= read -r telomere_seq_2; } < $(echo $PATH | awk '{gsub(/:/,"\n",$0)}1' | egrep IPA | uniq)/telomere_sequences_IPA.txt
if [ -z "$telomere_seq_1" ] ; then
{ IFS= read -r telomere_seq_1 && IFS= read -r telomere_seq_2; } < $(dirname $0)/telomere_sequences_IPA.txt
fi
if [ -z "$telomere_seq_2" ] ; then
telomere_seq_1="left: CCCTAAACCCTAAACCCTAAA"
telomere_seq_2="right: TTTAGGGTTTAGGGTTTAGGG"
echo "Check out the file telomere_sequences_IPA.txt is located in the same folder than the IPA script or in a folder in the PATH containing the string 'IPA'"
fi
printf "The telomere sequences used are:\n"$telomere_seq_1"\n"$telomere_seq_2"\nChange the file 'telomere_sequences_IPA.txt' if you want to use other sequences\n"
# cat $dir/5.Circlator/05.assembly.fa | IPA.renameAssembly.pl $name > 06.assembly.fa
if [ -z "$reference" ] ; then
cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_".$n } else { print }' > 06.assembly.fa; else
cat $dir/5.Circlator/05.assembly.fa | perl -nle 'if (/>(\S+)$/){ $n=$1; print ">".$ENV{name}."_ref_".$n } else { print }' > 06.assembly.fa
fi
gtool.py -sg 06.assembly.fa | tee contigs_GC_size.txt
for i in $(grep "^>" 06.assembly.fa | awk 'sub(/^>/, "")'); do cat 06.assembly.fa | gtool.py -sg -c $i - >> GC_length.txt; done
fasta2singleLine.pl 06.assembly.fa - | perl -nle 'if (/>/){ print } else { $l="-"; $r="+"; $left=substr($_,0,1000); $right=substr($_,(length($_)-1000),1000);if ($left =~ /^.*$ENV{'telomere_seq_1'}/i) {$l="L";} if ($right =~ /$ENV{'telomere_seq_2'}.*$/i){ $r="R"} print "$l\t$r"}' > 06.TelomerCounts.txt
fasta2singleLine.pl 06.assembly.fa - | perl -nle 'if (/>/){ print } elsif (length($_)>400000) { $l="-"; $r="+"; $left=substr($_,0,1000); $right=substr($_,(length($_)-1000),1000);if ($left =~ /^.*$ENV{'telomere_seq_1'}/i) {$l="L";} if ($right =~ /$ENV{'telomere_seq_2'}.*$/i){ $r="R"} print "$l\t$r"}' > 06.TelomerCounts.large.txt
num1=$(egrep -c L 06.TelomerCounts.txt); num2=$(egrep -c R 06.TelomerCounts.txt); num=$(($num1+$num2))
echo "Amount of telomer repeats $num" | tee TelomersSummary.txt
num1=$(egrep -c L 06.TelomerCounts.large.txt); num2=$(egrep -c R 06.TelomerCounts.large.txt); num=$(($num1+$num2))
echo "Amount of telomer repeats attached to 400kb $num" | tee -a TelomersSummary.txt
num=$(egrep -c "L.*R"  06.TelomerCounts.large.txt);
echo "$num chromosomes have both telomers attached" | tee -a TelomersSummary.txt
assembly-stats $assembly | head -n 3
assembly-stats $assembly > assembly_stats_original_correction_IPA.txt
cp 06.assembly.fa ../$name.IPA.fasta
assembly-stats ../$name.IPA.fasta | head -n 3
assembly-stats ../$name.IPA.fasta >> assembly_stats_original_correction_IPA.txt
echo -e "\nSTEP 6: DONE\n"





######## BETA:

# 6b. Get sequencing depth: (reads need to be uncompressed)
# https://github.com/wudustan/fastq-info
# fastq_info.sh $correctedReads &> fastq_info_correctedReads.txt
# if [ -z "$illuminaReads" ]; then
# echo -e "\n"
# else
# fastq_info_2.sh $illuminaReads\_1.fastq $illuminaReads\_2.fastq &> fastq_info_IlluminaReads.txt
# illumina_read_length=$(head -n 2 R1.fastq |tail -n 2|wc -c)
# echo "Illumina_read_length $illumina_read_length"
# fastq-info.sh -r $illumina_read_length $illuminaReads\_1.fastq $illuminaReads\_2.fastq ../$name.IPA.fasta &> fastq_info_IlluminaReads_assembly.txt
# fi





#### 7. GATK4 comparison with the reference:
# mkdir -p $dir/7.GATK; cd $dir/7.GATK
# Mapping to reference:
# perl -S fasta_to_fastq.pl ../$name.IPA.fasta ? > $name.IPA.fasta.fq # assuming default quality 30 (? symbol, https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)
# cp $reference ref.fa
# bwa index ref.fa; bwa mem -t $cores -x intractg ref.fa $name.IPA.fasta.fq | samtools sort -@ $cores -o Mapped.assembly.ref.bam -
# samtools view -@ $cores -F 0x800 Mapped.assembly.ref.bam -o Mapped.assembly.prim.ref.bam # only primary alignments
# Preparing for Variant calling
# samtools faidx ref.fa
# java -XX:ParallelGCThreads=$cores -jar picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict
# java -XX:ParallelGCThreads=$cores -jar picard.jar AddOrReplaceReadGroups I=Mapped.assembly.prim.ref.bam O=Mapped.assembly.ref.with.RG.bam RGLB=IPA_1 RGPL=IPA_2 RGPU=barcode RGSM=INFO_2
# gatk MarkDuplicatesSpark -I Mapped.assembly.ref.with.RG.bam -O Mapped.assembly.ref.marked.duplicates.bam -M marked.dup.metrics.txt --spark-master local[$cores] --conf 'spark.executor.cores='$cores
# GATK Best Practice Data Pre-processing
# https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
# HaplotypeCaller
# https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
# The version that can be multithreaded is HaplotypeCallerSpark, still in beta version:
# https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_HaplotypeCallerSpark.php
# https://javadoc.io/doc/org.broadinstitute/gatk/latest/org/broadinstitute/hellbender/tools/HaplotypeCallerSpark.html
# Variant Calling
# gatk HaplotypeCaller -R ref.fa -I Mapped.assembly.ref.marked.duplicates.bam -O output_HaplotypeCaller.vcf -ERC GVCF --native-pair-hmm-threads $cores
# gatk CountVariants -V output_HaplotypeCaller.vcf | tail -n 2 > GATK_CountVariants_output_number.txt
# gatk VariantsToTable -V output_HaplotypeCaller.vcf -F CHROM -F POS -F END -F REF -F ALT -O GATK_Variants_final_table.txt
# gatk HaplotypeCallerSpark -R ref.fa -I Mapped.assembly.ref.marked.duplicates.bam -O output_HaplotypeCallerSpark.vcf -ERC GVCF --native-pair-hmm-threads $cores --spark-master local[$cores] --conf 'spark.executor.cores='$cores

#### 8. Busco
# mkdir -p $dir/8.busco; cd $dir/8.busco
# https://busco.ezlab.org/busco_userguide.html
# busco -m genome -i ../$name.IPA.fasta -o OUTPUT_busco -l plasmodium_odb10 --cpu $cores
# busco -m genome -i ../$name.IPA.fasta -o OUTPUT_busco --auto-lineage-euk --cpu $cores
# busco -m genome -i ../$name.IPA.fasta -o OUTPUT_busco --auto-lineage-prok --cpu $cores
# busco -m genome -i ../$name.IPA.fasta -o OUTPUT_busco --auto-lineage --cpu $cores
