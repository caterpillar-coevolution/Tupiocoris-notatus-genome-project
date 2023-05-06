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

# step 4 repeatmasking
singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif BuildDatabase -name Tnot_twice_purged_wLTRs.DB -engine rmblast mirid_hifi_assembly.purged.purged.fa
singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif RepeatModeler -pa 32 -database Tnot_twice_purged_wLTRs.DB -LTRStruct
singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif RepeatMasker -lib Tnot_twice_purged_wLTRs.DB-families.fa -xsmall -pa 32 -gff -e ncbi mirid_hifi_assembly.purged.purged.fa

$ step 5: structural annotation
singularity exec /Funannotate/funannotate.sif funannotate train -i Tnot_hifiasm_2Xpurged.fa.masked -o Tnot_funannotate_v2 -l Tn1_1.fq.gz Tn2_1.fq.gz Tnot3_1.fq.gz -r Tn1_2.fq.gz Tn2_2.fq.gz Tnot3_2.fq.gz --cpus 32
singularity exec /Funannotate/funannotate.sif funannotate predict -i Tnot_hifiasm_2Xpurged.fa.masked -o Tnot_funannotate_v2 -s "Tnot_hifiasm_2Xpurged" --cpus 32 --max_intronlen 100000 --organism other --busco_db insecta -d /Funannotate/funannotate_databases --repeats2evm --protein_evidence Apolygus_lucorum.anno.pep.fa

#step 6: functional annotation
nohup /my_interproscan/interproscan-5.57-90.0/interproscan.sh -i 		Tnot_hifiasm_2Xpurged.proteins.fa &
singularity exec /funannotate.sif funannotate annotate -i ./Tnot_funannotate_v2/ --cpus 32 --iprscan Tnot_hifiasm_2Xpurged_updated.proteins.fa.xml

# step 7: BUSCO analysis of asembly and annotation
/busco_5.1.3.sif busco -i ./Tnot_hifiasm_2Xpurged.proteins.fa -l hemiptera_odb10 -o mirid_transcriptome_busco_analysis -m prot -f
/busco_5.1.3.sif busco -i ./mirid_hifi_assembly.fa -l hemiptera_odb10 -o mirid_purged_once_busco_analysis -m genome -f
singularity exec /busco_5.1.3.sif python3 /generate_plot.py --working_directory /path/to/wd/

