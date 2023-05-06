# Tupiocoris notatus genome project 
# command line scripts used in carrying outthe project
# These are not meant to be run as written, but to describe the steps we used

# step 1: assembly with hifiasm
hifiasm -o Tnot_hifiasm_assembly.asm -t 32 m64041_220407_192334.hifi_reads.fasta.gz

# step 2: purging duplicate contigs
# this step was performed twice
module load python/3.9/3.9.10
/home/u15/jaykgold/purge_dups/scripts/run_purge_dups.py -p bash mirid_hifiasm_config.json /home/u15/jaykgold/purge_dups/src Tnotatus_round2
