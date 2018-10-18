# Identify spliced leaders in *H.glycines* transcripts, associate with  genes and genomic strata


## Direct BLAST of spliced leaders to transcripts
### BLAST of full length SL's to transcripts using same BLAST params as below
```
#/work/GIF/remkv6/SplicedLeaders/22_FLblastTrinity

blastn -db ../21_ReadMappingToTranscripts/01_BLAST/nematode_transcripts_Trinity.blastdb -query query22bp.fasta -evalue 10000 -num_alignments 10000 -word_size 5 -dust no -task blastn-short -num_threads 16 -outfmt 6 -out 22bpSLTranscripts.blastout

#number of transcripts that meet standards
 awk '$7<5 && $9<27 && $4>20' 22bpSLTranscripts.blastout |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|wc
     15      15     240
```
### BLAST of SLs to transcripts with 11bp 3' ends of SL sequences and no dust
```
#/work/GIF/remkv6/SplicedLeaders/17_SL2Transcripts

blastn -db ../21_ReadMappingToTranscripts/01_BLAST/nematode_transcripts_Trinity.blastdb -query query.fasta -evalue 10000 -num_alignments 10000 -word_size 5 -dust no -task blastn-short -num_threads 16 -outfmt 6 -out SL2Transcriptsv2Filtered.blast.out

less SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sort|uniq|grep  -f - ../../Baum/CamTechGenomeComparison/30_transcriptMappingUpdated/nematode_transcripts_Trinity.fasta.gff >TranscriptsMapped2Genome.SLpositive.gff &

#number of transcripts that meet standards
less SL2Transcriptsv2Filtered.blast.out|awk '$7<13 && $9<22 {print $2}' |sed 's/ID=//g'|sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|wc
   2076    2076   33120
less SL2Transcriptsv2Filtered.blast.out|awk '$7<13 && $9<22 {print $2}' |sed 's/ID=//g'|sed 's/_/\t/2' |awk '{print $1}' |sort|uniq>SLReadsInTranscripts.list
```


# Identify spliced leaders in RNA seq paired-end reads
Got a new list of SL's from stacey, ones that are high confidence. Repeat a blast of these SLs to the R1 and R2 reads.\\

### Blast spliced leaders to paired reads
```
#/work/GIF/remkv6/SplicedLeaders/02_SLAlignment

#blast with 11bp putative SL query to second mate of read
blastn -db ../01_OldAnalyses/RNAseq2.blastdb/Heterodera_glycines_reallyraw_rnaseq.blastdb -query 11bpQuery.fasta  -evalue 10000 -num_alignments 10000 -word_size 5 -dust no -task blastn-short -num_threads 16 -outfmt 6 -out SlList11bpR2.blastout

#blast with 11bp putative SL query first mate of read
blastn -db all_R1.fastq.blastdb -query 11bpQuery.fasta  -evalue 10000 -task blastn-short -num_alignments 10000 -word_size 5 -dust no -num_threads 16 -outfmt 6  -out SlList11bpR1.blastout
```

# How many reads do we hit with 11bp 3' spliced leaders/Extract reads
```
#Methods were proven in 20_SLtest
#How many reads did the
less SlList11bpR1.blastout | awk '{if($10 > $9 && $9 < 13 && $4>10) {print $0} else if($9 > $10 && $9 > 88 && $4>10) {print $0}}' |awk '{print $2}' |sort|uniq|wc
  27158   27158 1174923
less SlList11bpR2.blastout | awk '{if($10 > $9 && $9 < 13 && $4>10) {print $0} else if($9 > $10 && $9 > 88 && $4>10) {print $0}}' |awk '{print $2}' |sort|uniq|wc
  58712   58712 2565525

less SlList11bpR1.blastout | awk '{if($10 > $9 && $9 < 13 && $4>10) {print $0} else if($9 > $10 && $9 > 88 && $4>10) {print $0}}' |awk '{print $2}' |sort|uniq>R1.READS.list
less SlList11bpR2.blastout | awk '{if($10 > $9 && $9 < 13 && $4>10) {print $0} else if($9 > $10 && $9 > 88 && $4>10) {print $0}}' |awk '{print $2}' |sort|uniq>R2.READS.list

seqtk subseq  all.R1.1col.fasta R1.READS.list >11bp.extracted.reads.R1
grep -c ">" 11bp.extracted.reads.R1
27158

seqtk subseq  ../01_OldAnalyses/RNAseq2.blastdb/all_R2.fasta R2.READS.list >11bp.extracted.reads.R2
grep -c ">" 11bp.extracted.reads.R2
58712
```

# BLAST of 11bp SL-reads to transcripts to see if lower counts are being caused by stringency in aligner
```
#/work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST

#truncated the Trinity fasta assembly to only the first 150 bp of each transcript, in order to reduce the blast time and subsequent file size.
grep ">" nematode_transcripts_Trinity.fasta |sed 's/>//g' |awk '{print "samtools faidx nematode_transcripts_Trinity.fasta "$1":1-150 &"}' |sed '0~15 s/$/\nwait/g'>faidxTruncateTrinity.sh
sh faidxTruncateTrinity.sh >nematode_transcripts_TrinityTruncation.fasta &
#did it finish? yes
grep -c ">" nematode_transcripts_TrinityTruncation.fasta
180003

#blast the reads to the truncated trinity assembly
makeblastdb -in nematode_transcripts_TrinityTruncation.fasta -dbtype nucl -out nematode_transcripts_TrinityTruncation.blastdb
cat ../../02_SLAlignment/11bp.extracted.reads.R1 ../../02_SLAlignment/11bp.extracted.reads.R2 >extracted.reads.combine.fasta
blastn -db nematode_transcripts_TrinityTruncation.blastdb -query extracted.reads.combine.fasta -evalue 10000 -num_alignments 100 -word_size 5 -dust no -task blastn-short -num_threads 16 -outfmt 6 -out reads2Transcripts.blastout &


#How many transcripts had spliced leaders based on the read analysis using blast
less reads2Transcripts.blastout |awk '$4>80' |awk '$9<12 ||$10<12' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|wc
   1635    1635   26111

#Tests of other parameters
####less reads2Transcripts.blastout |awk '$4>80' |awk '$9<14 ||$10<14' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|wc
####   1661    1661   26526
####less reads2Transcripts.blastout |awk '$4>80' |awk '$9<19 ||$10<19' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|wc
####   1721    1721   27485


#How many overlapped with the direct blast analysis? --eventually used these parameters
less reads2Transcripts.blastout |awk '$4>80' |awk '$9<12 ||$10<12' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|cat - <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq)|sort|uniq -c|awk '$1==2' |wc
   1281    2562   30709

#Tests of other parameters
####[remkv6@condo041 01_BLAST]$ less reads2Transcripts.blastout |awk '$4>80' |awk '$9<14 ||$10<14' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|cat - <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq)|sort|uniq -c|awk '$1==2' |wc
####   1287    2574   30853
####[remkv6@condo041 01_BLAST]$ less reads2Transcripts.blastout |awk '$4>80' |awk '$9<19 ||$10<19' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|cat - <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq)|sort|uniq -c|awk '$1==2' |wc
####   1303    2606   31237

It looks like as I low the stringency, the error rate goes up significantly.  
I will stick with these parameters  "awk '$4>80' |awk '$9<12 ||$10<12'
less reads2Transcripts.blastout |awk '$4>80' |awk '$9<12 ||$10<12' |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq >transcriptCt.list


#What is the new total for spliced leader transcripts?
cat ../../05_EST/EST2TrinityTranscripts.list transcriptCt.list <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq)|sort|uniq|wc
   2532    2532   40402

#How many did not overlap with the direct blast analysis?
 less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|grep -w -v -f - transcriptCt.list |wc
    354     354    5650
```


### Use mapping gffs to convert transcripts and ESTs to genes
```
#getting Transcripts that the EST's aligned to.
#/work/GIF/remkv6/SplicedLeaders/05_EST
less Filtered.blastout |sort -u -k1,1 |awk '{print $2}' |grep -f - ../../Baum/CamTechGenomeComparison/30_transcriptMappingUpdated/nematode_transcripts_Trinity.fasta.gff >estConvertTranscript.gff


#/work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST  
#getting gff of the transcripts that were found in both reads,the direct blast, and ESTs
 awk '$3=="gene" ||$3=="exon"' ../../../Baum/CamTechGenomeComparison/30_transcriptMappingUpdated/nematode_transcripts_Trinity.fasta.gff |tr " " "\t"|sed 's/_/ /2' |sed 's/ID=/ID= /g' >nematode_transcripts_Trinity.fastaMod.gff
less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|cat -  transcriptCt.list |cat - ../../05_EST/EST2TrinityTranscripts.list |sort|uniq |while read line; do echo "awk '\$10==\""$line"\"' nematode_transcripts_Trinity.fastaMod.gff >>SLnematode_transcripts_Trinity.fastaMod.gff ";done |sed '0~15 s/$/\nwait/g' >getgff.sh
sh getgff.sh
less SLnematode_transcripts_Trinity.fastaMod.gff |sed 's/ //1' |sed 's/ /_/g' >SLnematode_transcripts_Trinity.gff

#Grabs only those genes that have a transcript overlapping with the first exon of the gene.
bedtools intersect -wo -a  <(awk '$3=="exon"' SLnematode_transcripts_Trinity.gff)  -b <(awk '$3=="exon"' ../../../Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 ) |cut -f 21| sed 's/\./\t/2' |sed 's/;/\t/g' |grep -w "exon1" |awk '{print $1}' |sed 's/\./\t/g' |sed 's/ID=//g' |cut -f 1 |sort|uniq >SLGenesExon1.list

#This is a representation of the below file in bedtools intersect format.
bedtools intersect -wo -a  <(awk '$3=="exon"' SLnematode_transcripts_Trinity.gff)  -b <(awk '$3=="exon"' ../../../Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 ) |sed 's/ID=/ID=\t/2' |sed 's/\./\t/8' |grep -w -f SLGenesExon1.list - >SLContainingGenes.gff

less SLGenesExon1.list |while read line; do echo "awk '\$10==\""$line"\"' augustusMod.gff3 >> SLGenes.gff3";done |sed '0~15 s/$/\nwait/g' >getgff2exon1.sh
less SLGenes.gff3 |sed 's/\t//9'|sed 's/\t//9'|sed 's/\t//9'|sed 's/\t//9'|cut -f 1-9 >SLGenesFixed.gff3

#How many genes had transcripts overlapping with the first exon?
awk '$3=="gene"' SLGenesFixed.gff3 |wc
   9042   81378  476232
```

### How many of the genes have exons overlapping with repeats?
```
bedtools intersect -wo -a <(awk '$3=="exon"' SLContainingGenes_GenesOnlyTrans.gff) -b ../../Baum/CamTechGenomeComparison/58_Renamatorium/2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |cut -f 9,18- |awk '($4-$3)>100 ' |sed 's/\./\t/2'|awk '{print $4}' |uniq|sort|uniq -c |sort -k1,1nr |awk '{print $1}' |summary.sh
Total:  7,827
Count:  396
Mean:   19
Median: 6
Min:    1
Max:    407

module load GIF/perl/5.24.1
bedtools intersect -wo -a <(awk '$3=="exon"' SLContainingGenes_GenesOnlyTrans.gff) -b ../../Baum/CamTechGenomeComparison/58_Renamatorium/2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |cut -f 9,18- |awk '($4-$3)>100 ' |sed 's/\./\t/2'|awk '{print $4}' |uniq|sort|uniq -c |sort -k1,1nr |awk '{print $1}' |st --sd |less
43.2181 == standard deviation
43.22*3 +19 == 149
bedtools intersect -wo -a <(awk '$3=="exon"' SLContainingGenes_GenesOnlyTrans.gff) -b ../../Baum/CamTechGenomeComparison/58_Renamatorium/2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |cut -f 9,18- |awk '($4-$3)>100 ' |sed 's/\./\t/2'|awk '{print $4}' |uniq|sort|uniq -c |sort -k1,1nr |awk '$1>1391 {print $1/69182,$2}'  |awk '{print $2}' |sed 's/"//g' |sed 's/Motif://g' |grep -f - <(grep ">" ../../Baum/CamTechGenomeComparison/22_RepeatModeler/RM_3002.WedNov301615082016/consensi.fa.classified)|cat <(bedtools intersect -wo -a <(awk '$3=="exon"' SLContainingGenes_GenesOnlyTrans.gff) -b ../../Baum/CamTechGenomeComparison/58_Renamatorium/2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |cut -f 9,18- |awk '($4-$3)>100 ' |sed 's/\./\t/2'|awk '{print $4}' |uniq|sort|uniq -c |sort -k1,1nr |awk '$1>149 {print $1/8435,$2}' ) - >3SDRepeatsCtPctNameType.list
[remkv6@condo189 19_SLTranscriptOverlapAnalyses]$ less <(bedtools intersect -wo -a <(awk '$3=="exon"' SLContainingGenes_GenesOnlyTrans.gff) -b ../../Baum/CamTechGenomeComparison/58_Renamatorium/2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |cut -f 9,18- |awk '($4-$3)>100 ' |sed 's/\./\t/2'|awk '{print $4}' |uniq|sort|uniq -c |sort -k1,1nr |awk '$1>149 {print $1/8435,$2}' ) - >3SDRepeatsCtPctNameType.list
Those that  are above 3SD(171).
0.0482513 "Motif:rnd-4_family-65" Unknown
0.0368702 "Motif:rnd-4_family-352" Unknown
0.0336692 "Motif:rnd-4_family-1273" LINE/CR1
0.0309425 "Motif:rnd-5_family-2862" Unknown
0.0248963 "Motif:rnd-3_family-228" Unknown
0.0237107 "Motif:rnd-4_family-1152" DNA/TcMar-Tc2
0.0203912 "Motif:rnd-4_family-299" Unknown
0.0203912 "Motif:rnd-4_family-71" DNA/Mule-MuDR
0.0197985 "Motif:rnd-3_family-757" DNA/MuLE-MuDR
0.0182573 "Motif:rnd-3_family-42" Unknown
```

### functional analysis
```
#/work/GIF/remkv6/SplicedLeaders/19_SLTranscriptOverlapAnalyses/01_Ontologizer
wget http://purl.obolibrary.org/obo/go.obo
ln -s ../../21_ReadMappingToTranscripts/01_BLAST/SLGenesExon1.list
simpleformat.ids -> ../../../Baum/CamTechGenomeComparison/57_secretome/ontologenizer/simpleformat.ids
population -> ../../../Baum/CamTechGenomeComparison/57_secretome/ontologenizer/population
Ontologizer.jar -> ../../../Baum/CamTechGenomeComparison/57_secretome/ontologenizer/Ontologizer.jar

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s SLGenesExon1.list
```

```
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0016773      16030   447     9042    180     856     240     1       false   2.129109577352885E-17   4.7223650425686995E-14  1.7731639413095818E-256 "phosphotransferase activity, alcohol group as acceptor"
GO:0016301      16030   480     9042    189     863     243     2       false   3.9829775046754863E-17  8.834244105370228E-14   1.402943450201154E-256  "kinase activity"
GO:0044267      16030   895     9042    379     2393    796     2       false   2.9624121947224525E-13  6.5706302478944E-10     0.0     "cellular protein metabolic process"
GO:0017076      16030   1085    9042    424     1566    532     1       false   4.522455783858459E-11   1.0030806928598062E-7   0.0     "purine nucleotide binding"
GO:0032553      16030   1079    9042    421     1568    532     2       false   8.347452792627276E-11   1.8514650294047297E-7   0.0     "ribonucleotide binding"
GO:0035639      16030   1024    9042    399     1602    547     2       false   2.956651868176443E-8    6.557853843615351E-5    0.0     "purine ribonucleoside triphosphate binding"
GO:0006261      16030   17      9042    10      259     23      1       false   5.025360003269504E-8    1.114624848725176E-4    5.7282504081053664E-27  "DNA-dependent DNA replication"
GO:0006414      16030   127     9042    64      917     270     2       false   6.583056075757906E-8    1.4601218376031035E-4   1.7167452646656717E-159 "translational elongation"
GO:0043412      16030   810     9042    328     2428    808     1       false   7.312438757494017E-8    1.621898916412173E-4    0.0     "macromolecule modification"
GO:0003723      16030   191     9042    81      1805    470     1       false   1.3646149469977655E-7   3.026715952441044E-4    6.424296199953532E-264  "RNA binding"
GO:0005575      16030   2016    9042    775     7387    2495    1       false   1.3907423608678334E-7   3.0846665564048545E-4   0.0     "cellular_component"
GO:0005198      16030   197     9042    100     6526    2163    1       false   1.632676340562515E-7    3.6212761233676583E-4   0.0     "structural molecule activity"
GO:0097747      16030   46      9042    19      370     49      1       false   2.2890516454945808E-7   5.07711654970698E-4     7.434414868994774E-60   "RNA polymerase activity"
GO:0032991      16030   508     9042    242     2016    775     1       false   6.513714142728797E-7    0.001444741796857247    0.0     "macromolecular complex"
GO:0005515      16030   964     9042    381     4267    1417    1       false   1.6630056246693103E-6   0.0036885464755165302   0.0     "protein binding"
GO:0043604      16030   213     9042    90      1069    310     2       false   2.4642144726819073E-6   0.005465627700408471    4.7911230994740565E-231 "amide biosynthetic process"
GO:0097367      16030   1081    9042    421     4267    1417    1       false   2.5356851445634256E-6   0.005624149650641678    0.0     "carbohydrate derivative binding"
GO:0042578      16030   107     9042    39      348     72      1       false   2.595419441078604E-6    0.005756640320312344    1.1884467014208498E-92  "phosphoric ester hydrolase activity"
GO:0071840      16030   463     9042    201     4821    1642    1       false   6.860307033627184E-6    0.015216161000585093    0.0     "cellular component organization or biogenesis"
GO:0030529      16030   150     9042    86      1165    483     2       false   1.9996988708968514E-5   0.044353320956492165    1.4603944775325897E-193 "intracellular ribonucleoprotein complex"
```

### EST found to have spliced leaders attributed to transcripts#
```
#/work/GIF/remkv6/SplicedLeaders/05_EST
#getting Transcripts that the EST's aligned to.


makeblastdb -in ../17_SL2Transcripts/nematode_transcripts_Trinity.fasta -dbtype nucl -out nematode_transcripts_Trinity.blastdb
blastn -db nematode_transcripts_Trinity.blastdb -query SLESTs.fa -outfmt 6 -dust no -out SLESTs2TrinityTranscripts.blastout

#How many unique hits are there?
less SLESTs2TrinityTranscripts.blastout |sed 's/_/ /2' |sort -u -k2,2V |wc
    342    5130   29360

#How many meet the filtering criteria of alignment length vs total query length?
less SLESTs2TrinityTranscripts.blastout |sed 's/_/ /2' |sort -u -k2,2V |awk '($5/$14)> .5 && ($5/$14)<1.5' |wc
    187    2805   15889

less SLESTs2TrinityTranscripts.blastout |sed 's/_/ /2' |sort -u -k2,2V |awk '($5/$14)> .5 && ($5/$14)<1.5' |awk '{print $2}' |sort|uniq>EST2TrinityTranscripts.list

less SLESTs2TrinityTranscripts.blastout |sed 's/_/ /2' |sort -u -k2,2V |awk '($5/$14)> .5 && ($5/$14)<1.5' >Filtered.blastout

#Create gff for further comparison if needed
less Filtered.blastout |sort -u -k1,1 |awk '{print $2}' |grep -f - ../../Baum/CamTechGenomeComparison/30_transcriptMappingUpdated/nematode_transcripts_Trinity.fasta.gff >estConvertTranscript.gff

#how many are not found already in the direct blast and the read oriented analysis.
#/work/GIF/remkv6/SplicedLeaders/19_SLTranscriptOverlapAnalyses
cat SLReadsInTranscripts.list <(awk '{print $2}' transcriptCt.list) |sort|uniq |grep -w -v -f - EST2TrinityTranscripts.list |wc                           
104
```


### Extract transcript fasta sequences for Blast to Go transcript Comparison
```
#/work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST
#total number of transcripts with EST's, direct blasts, and read analysis.
 wc AllSLTranscriptsCombine.list
 2532  2532 40402 AllSLTranscriptsCombine.list

grep -w -f AllSLTranscriptsCombine.list <(grep ">" nematode_transcripts_Trinity.fasta|sed 's/_/\t/2'|sed 's/>//g') |sort -u -k1,1V|sed 's/\t/_/1' |awk '{print $1}' |sort|uniq|cdbyank nematode_transcripts_Trinity.fasta.cidx  >B2GAllSLTranscripts.fasta &
```

### Effector overlap
```
#/work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST

#how many spliced leader genes are known effectors?
less SLGenesExon1.list |grep -w -f - ../../../Baum/CamTechGenomeComparison/58_Renamatorium/18_effectorRedo/121EffectorOldNewGeneNames.list |wc
     29      58     713

#what is the distribution of these effectors
less SLGenesExon1.list |grep -w -f - ../../../Baum/CamTechGenomeComparison/58_Renamatorium/18_effectorRedo/121EffectorOldNewGeneNames.list |awk '{print $2}' |grep -w -f - ../../../Baum/CamTechGenomeComparison/58_Renamatorium/18_effectorRedo/4SebastianEffvsGenes.list|sort -u -k1,1V |awk '{print $2}' |sort|uniq -c |sort -k1,1nr |less
# # # # # # # # # # # # # # #
5 11A06
      3 16A01
      3 GLAND5
      2 10C02
      2 4F01
      1 10A06
      1 17G06
      1 19B10
      1 21E12
      1 2D01
      1 33A09
      1 33E05
      1 5A08(missing5’)
      1 5D08
      1 8A07
      1 flGSB3
      1 GLAND11
      1 GLAND12
      1 GLAND4
# # # # # # # #

#What is the total of effectors in the genome
less ../../../Baum/CamTechGenomeComparison/58_Renamatorium/18_effectorRedo/4SebastianEffvsGenes.list|sort -u -k1,1V |awk '{print $2}' |sort|uniq -c |sort -k1,1nr |less

################
6 GLAND5
5 11A06
5 16A01
5 32E03
5 4D06
4 15A10
4 GLAND13
4 GLAND6
4 GLAND9
3 16B09
3 5D08
3 8C06
3 GLAND14
2 10A07
2 10C02
2 17G06
2 25A01
2 29D09
2 2A05
2 2B10
2 30C02
2 30D08
2 33E05
2 3H07
2 4F01
2 6F06
2 GLAND1
2 GLAND10
1 10A06
1 12H04
1 13A06
1 13C08
1 16H02
1 17C07
1 17G01
1 18H08
1 19B10
1 19C07
1 20E03
1 21E12
1 22C12
1 26D05
1 2D01
1 2E04
1 31A08
1 33A09
1 3B05
1 3D11
1 45D07
1 4D09
1 5A08(missing5’)
1 5D06
1 7E05
1 8A07
1 flGSB3
1 GLAND11
1 GLAND12
1 GLAND15
1 GLAND16
1 GLAND17
1 GLAND2
1 GLAND3
1 GLAND4
1 GLAND8
################
```

### Locational genomic clustering of genes giving rise to trans-spliced transcripts
```
bedtools intersect -wo -a SLGenesFixed.gff3 -b <(awk '$3=="gene"' ../../../Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3) |cut -f 10,11,12,13,14,15,16,17,18 |uniq|sort|uniq|awk '{print $1,$4/50000}' |sed 's/\./\t/1' |awk '{print $1,$2*50000}' |sort |uniq -c |sort -k1,1nr |less
# # # # # # # # # #
    18 000001 800000
     16 000002 700000
     16 000118 200000
     16 000221 0
     16 000613K 0
     16 002293 100000
     15 000118 150000
     15 24syntmer 100000
     14 000171 200000
     14 000309 50000
     13 000001 250000
     13 000058 600000
     13 000066 500000
     13 000108 50000
     13 000118 50000
     13 000133 50000
     13 000168 0
     13 000203 100000
     13 000220 0
     13 000220 100000
     13 000252 100000
     13 000286 50000
     13 000502 50000
     13 000831K 0
     13 23syntmer 0
     13 24syntmer 150000
     13 35syntmer 450000
     13 39syntmer 100000
     13 46syntmer 100000
     12 000014 50000
     12 000049 200000
     12 000058 900000
     12 000108 200000
     12 000138 100000
     12 000177 250000
     12 000188 100000
     12 000194 100000
     12 000380 100000
     12 000467K 50000
     12 000502 250000
     12 000640K 0
     12 35syntmer 100000
     12 35syntmer 400000
     12 35syntmer 50000
     12 41syntmer 50000
     11 000001 750000
     11 000012 450000
     11 000015 50000
     11 000028 50000
     11 000029 350000
     11 000038 100000
     11 000066 600000
     11 000091 100000
     11 000131 150000
     11 000138 200000
     11 000139 50000
     11 000171 1000000
     11 000217 50000
     11 000220 200000
     11 000252 0
     11 000252 50000


# # # # # # # # # #
counts for 50kb bins of spliced leader genes.
bedtools intersect -wo -a SLGenesFixed.gff3 -b <(awk '$3=="gene"' ../../../Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3) |cut -f 10,11,12,13,14,15,16,17,18 |uniq|sort|uniq|awk '{print $1,$4/50000}' |sed 's/\./\t/1' |awk '{print $1,$2*50000}' |sort |uniq -c |sort -k1,1nr |awk '{print $1}' |summary.sh
Total:  9,085
Count:  2,381
Mean:   3
Median: 3
Min:    1
Max:    18


#standard deviation
bedtools intersect -wo -a SLGenesFixed.gff3 -b <(awk '$3=="gene"' ../../../Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3) |cut -f 10,11,12,13,14,15,16,17,18 |uniq|sort|uniq|awk '{print $1,$4/50000}' |sed 's/\./\t/1' |awk '{print $1,$2*50000}' |sort |uniq -c |sort -k1,1nr |awk '{print $1}' |st --sd

2.76625

3 standard deviations (3+2.8+2.8+2.8)= 11.4

#which bins are left
 bedtools intersect -wo -a SLGenesFixed.gff3 -b <(awk '$3=="gene"' ../../../Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3) |cut -f 10,11,12,13,14,15,16,17,18 |uniq|sort|uniq|awk '{print $1,$4/50000}' |sed 's/\./\t/1' |awk '{print $1,$2*50000}' |sort |uniq -c |sort -k1,1nr |awk '$1>11.4' |less


       18 000001 800000
     16 000002 700000
     16 000118 200000
     16 000221 0
     16 000613K 0
     16 002293 100000
     15 000118 150000
     15 24syntmer 100000
     14 000171 200000
     14 000309 50000
     13 000001 250000
     13 000058 600000
     13 000066 500000
     13 000108 50000
     13 000118 50000
     13 000133 50000
     13 000168 0
     13 000203 100000
     13 000220 0
     13 000220 100000
     13 000252 100000
     13 000286 50000
     13 000502 50000
     13 000831K 0
     13 23syntmer 0
     13 24syntmer 150000
     13 35syntmer 450000
     13 39syntmer 100000
     13 46syntmer 100000
     12 000014 50000
     12 000049 200000
     12 000058 900000
     12 000108 200000
     12 000138 100000
     12 000177 250000
     12 000188 100000
     12 000194 100000
     12 000380 100000
     12 000467K 50000
     12 000502 250000
     12 000640K 0
     12 35syntmer 100000
     12 35syntmer 400000
     12 35syntmer 50000
     12 41syntmer 50000

# # # # # # # # #

#20kb bins may be more appropriate, as the ratio of sl gene/distance is increased.  (12 scaffold_16 700000)
```
### Venn Diagram for promiscuously spliced trans-spliced leaders
I need to find the count for the each SL that had a read with a 5' spliced leader, mapping to the 5' end of transcripts -- 20 ESTs with spliced leaders did not have a corresponding transcript, so they were added to the venn diagram subsequently.
```
/work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST
cat ../../05_EST/EST2TrinityTranscripts.list transcriptCt.list <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq)|sort|uniq|wc
   2532    2532   40402


#Total to direct blast
less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq|wc
   2076    2076   33120

#Total to the read oriented blast
wc transcriptCt.list
 1635  1635 26111 transcriptCt.list

#Total to the ESTs
   wc EST2TrinityTranscripts.list
 187  187 2992 EST2TrinityTranscripts.list


#Unique to the direct blast
cat ../../05_EST/EST2TrinityTranscripts.list transcriptCt.list |grep -w -v -f - <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq) |wc
    785     785   12499

#Unique to the read oriented analysis
cat ../../05_EST/EST2TrinityTranscripts.list  <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq)  |grep -w -v -f - transcriptCt.list |wc
    343     343    5474

#what is unique to the EST'samples
cat <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq) transcriptCt.list |grep -w -v -f - ../../05_EST/EST2TrinityTranscripts.list|wc
    102     102    1632



#what is shared between reads and direct blast
cat <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq) transcriptCt.list |sort|uniq -c|awk '$1==2' |wc
   1281    2562   30709

#what is shared between reads and ESTs
cat transcriptCt.list ../../05_EST/EST2TrinityTranscripts.list |sort|uniq -c|awk '$1==2' |wc
     75     150    1800


#What is shared between the direct blast and ESTs
cat <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq) ../../05_EST/EST2TrinityTranscripts.list |sort|uniq -c|awk '$1==2' |wc
     74     148    1776



#found in all three samples   
cat ../../05_EST/EST2TrinityTranscripts.list <(less ../../17_SL2Transcripts/SL2Transcriptsv2Filtered.blast.out |awk '{print $2}' |sed 's/_/\t/2' |awk '{print $1}' |sort|uniq) transcriptCt.list |sort|uniq -c|awk '$1==3' |wc
     64     128    1536


```

### Promiscuity analysis to see how often each transcript receives a different spliced leader

I need to find the count for the each SL that had a read with a 5' spliced leader, mapping to the 5' end of transcripts

```
#/work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST
###less reads2Transcripts.blastout |awk '$9<20'  |awk '$10>$9 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' |awk '{print $1}' |sort|uniq|grep -f - ../02_SLAlignment/SlList11bpR1.blastout >SLContainingReadsR1.blastout &
###less SLContainingReadsR1.blastout|awk '{print $1}' |sort |uniq -c |sort -k1,1 |less

less reads2TranscriptsFiltered.blastout |sed 's/_/\t/2' |awk '{print $1,$2,$3}' >reads2TranscriptsForGrep.list
less transcriptCt.list |awk '{print $2}' |while read line; do echo "awk '\$2==\"$line\"' reads2TranscriptsForGrep.list >>reads2transcriptsExtracted1539.list";done |sed '0~15 s/$/\nwait/g' >Extractreads2Transcripts.sh

#/work/GIF/remkv6/SplicedLeaders/02_SLAlignment

#Create the grep database
less SlList11bpR1.blastout | awk '{if($10 > $9 && $9 < 13 && $4>10) {print $0} else if($9 > $10 && $9 > 88 && $4>10) {print $0}}' |sed 's/_/\t/2' |awk '{print $1,$2}' |sort|uniq >grepsubjectPromiscuity
less SlList11bpR2.blastout | awk '{if($10 > $9 && $9 < 13 && $4>10) {print $0} else if($9 > $10 && $9 > 88 && $4>10) {print $0}}' |sed 's/_/\t/2'|awk '{print $1,$2}' |sort|uniq >> grepsubjectPromiscuity
sort grepsubjectPromiscuity|uniq >grepsubjectPromiscuity2

#align so that spliced leader match read names and read names match transcripts
less /work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST/reads2Transcripts.blastout |awk '$4>80' |awk '$9<12 ||$10<12' |sed 's/_/\t/2' |awk '{print $1,$2}' |sort|uniq|awk '{print $1}' |while read line; do grep -w  $line grepsubjectPromiscuity2;done|sort -k2,2V|paste - <( less /work/GIF/remkv6/SplicedLeaders/21_ReadMappingToTranscripts/01_BLAST/reads2Transcripts.blastout |awk '$4>80' |awk '$9<12 ||$10<12' |sed 's/_/\t/2' |awk '{print $1,$2}' |sort|uniq|sort -k1,1V) >ImperfectPromiscuityOutput
```
