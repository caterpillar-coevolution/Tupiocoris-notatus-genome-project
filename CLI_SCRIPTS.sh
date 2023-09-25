# Tupiocoris notatus genome project 
# command line scripts used in carrying outthe project
# These are not meant to be run as written, but to describe the steps we used

# step 1: assembly with hifiasm
hifiasm -o Tnot_hifiasm_assembly.asm -t 32 m64041_220407_192334.hifi_reads.fasta.gz

# step 2: purging duplicate contigs
# this step was performed twice
module load python/3.9/3.9.10
/purge_dups/scripts/run_purge_dups.py -p bash mirid_hifiasm_config.json /home/u15/jaykgold/purge_dups/src Tnotatus_round2

# step 3: blob tools contamination check
# create coverage files
module load samtools
/home/u15/jaykgold/minimap2-2.24_x64-linux/minimap2 -ax map-hifi mirid_hifi_assembly.purged.purged.fa /Tupiocoris_notatus/raw_HiFi_reads/m64041_220407_192334.hifi_reads.fasta.gz > Tnot_coverage.sam
samtools view -bS Tnot_coverage.sam > Tnot_coverage.bam
# create "hits" file
module load blast/2.13.0  
blastn \
-task 'dc‑megablast' \
-query 	mirid_hifi_assembly.purged.purged.fa \
-db nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads 32 \
-evalue 1e-25 \
-out Tnot.vs.nt.redownloaded.dcmegablast.out
# acutally make plots
module load python/3.9/3.9.10
source /blob_tools/blobtools_env/bin/activate
/blob_tools/blobtools/blobtools create \
 -i ./mirid_hifi_assembly.purged.purged.fa \
 -b ./Tnot_coverage_sorted.bam \
 -t ./assembly.vs.uniprot_ref.mts1.1e25.taxified.out \
 -o ./my_first_blobplot \
 --db /blob_tools/blobtools/data/nodesDB.txt \
 -x bestsumorder
 

 
/blob_tools/blobtools/blobtools view \
 -i ./my_first_blobplot.blobDB.json \
 -x bestsumorder \
 -o ./
 
 grep '^##' example/my_first_blobplot.blobDB.table.txt ; \
 grep -v '^##' example/my_first_blobplot.blobDB.table.txt | \
 column -t -s $'\t'
 
/xdisk/judieb/jaykgold/blob_tools/blobtools/blobtools plot \
 -i ./my_first_blobplot.blobDB.json \
 -x bestsumorder \
 -o ./
#step 4: inspector analysis and correction

python3 inspector-correct.py -i Tnot_inspector_run2_out/ --datatype pacbio-hifi -o Tnot_inspector_run2_out/ 

# step 5 repeatmasking
singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif BuildDatabase -name Tnot_twice_purged_wLTRs.DB -engine rmblast mirid_hifi_assembly.purged.purged.fa
singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif RepeatModeler -pa 32 -database Tnot_twice_purged_wLTRs.DB -LTRStruct
singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif RepeatMasker -lib Tnot_twice_purged_wLTRs.DB-families.fa -xsmall -pa 32 -gff -e ncbi mirid_hifi_assembly.purged.purged.fa

# step 6: structural annotation
singularity run --nv helixer-docker_helixer_v0.3.1_cuda_11.2.0-cudnn8.sif Helixer.py --fasta-path Tnot_hifiasm_contig_corrected.fa \
--lineage invertebrate --gff-output-path Tnot_corrected_contigs.gff3 --subsequence-length 213840 --overlap-offset 106920 --overlap-core-length 160380
gff3_to_fasta -g Tnot_corrected_contigs.gff3 -f Tnot_hifiasm_contig_corrected.fa -st all -d simple -o Tnot_corrected_contigs

#step 7: functional annotation
# Interproscan
singularity run \
-B Tnot_filtered_analysis/CORRECTED_ASSEMBLY:/data \
-B interproscan-5.45-80.0/data:/opt/interproscan/data \
interproscan_5.45-80.0_1.sif \
-i Tnot_corrected_contigs_pep.fa \
-d outdir_Tnot_corrected_InterPro \
-f tsv,json,xml,html,gff3,svg \
-g \
-p \
-c \
-n \
-D Tnotatus \
-l

#blastp
module load blast/2.13.0 
blastp -db databases/NCBI_nr/nr -query Tnot_corrected_contigs_pep.fa -outfmt 6 -out Tnot_blastp_nr_for_annotation.out -num_threads 32

# AGAT update gff
singularity run \
    agat_1.1.0--pl5321hdfd78af_1.sif \
    agat_sp_manage_functional_annotation.pl \
    -f Tnot_corrected_contigs.gff3 \
    -b ./Tnot_blastp_for_annotation.out \
    --db databases/uniprot_sprot.fasta \
    -i CORRECTED_ASSEMBLY/outdir_Tnot_corrected_InterPro/Tnot_corrected_contigs_pep.tsv \
    --output Tnot_helixer_w_interpro_blastp

# step 7: BUSCO analysis of asembly and annotation
/busco_5.1.3.sif busco -i ./Tnot_hifiasm_2Xpurged.proteins.fa -l hemiptera_odb10 -o mirid_transcriptome_busco_analysis -m prot -f
/busco_5.1.3.sif busco -i ./mirid_hifi_assembly.fa -l hemiptera_odb10 -o mirid_purged_once_busco_analysis -m genome -f
singularity exec /busco_5.1.3.sif python3 /generate_plot.py --working_directory /path/to/wd/

#step 8: alignment with STAR for diff expression
module load star
module load cufflinks
gffread Tnot_hifiasm_2Xpurged.gff3 -T -o Tnot_hifiasm_2Xpurged.gtf
STAR --runThreadN 15 \
--runMode genomeGenerate \
--genomeDir Tnot/ \
--genomeFastaFiles Tnot_hifiasm_2Xpurged.fa.masked  \
--sjdbGTFfile Tnot_hifiasm_2Xpurged.gtf \
--genomeSAindexNbases 13 \
--sjdbOverhang 74 
module load bcftools
module load samtools
module load perl
module load sratoolkit
#Run trimmomatic
java -jar /xdisk/judieb/jaykgold/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
fasta/SRR4289599_1.fastq.gz fasta/SRR4289599_2.fastq.gz \
fasta/SRR4289599_1_trimmed.fa fasta/SRR4289599_1_unpair_trimmed.fa \
fasta/SRR4289599_2_trimmed.fa fasta/SRR4289599_2_unpair_trimmed.fa \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
STAR --genomeDir Tnot/ \
--runThreadN 15 \
--readFilesIn fasta/SRR4289599_1_trimmed.fa fasta/SRR4289599_2_trimmed.fa \
--outFileNamePrefix results/SRR4289599 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--limitSjdbInsertNsj=3000000 \
--limitOutSJcollapsed=3000000
#index the .bam file
samtools index results/SRR4289599Aligned.sortedByCoord.out.bam




