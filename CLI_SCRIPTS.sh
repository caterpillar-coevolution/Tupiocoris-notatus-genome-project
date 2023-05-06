# Tupiocoris notatus genome project 
# command line scripts used in carrying outthe project
# These are not meant to be run as written, but to describe the steps we used

# step 1: assembly with hifiasm
hifiasm -o Tnot_hifiasm_assembly.asm -t 32 m64041_220407_192334.hifi_reads.fasta.gz
