import random
import os


WORKDIR = config['workdir']
UPDATE_EC_NUMS = config['update_ec_nums']
ANNOTATION_COLS = config['annotation_cols']


# BRENDA_PARSER = BrendaParser()

# # # Temporary code for testing to sample from the big BRENDA dict
# keys = random.sample(BRENDA_PARSER.keys(), 2)

# print (len(BRENDA_PARSER.keys()))
# with open('EC_numbers.txt') as ec_txt:
# 	ec_txt.write(BRENDA_PARSER.keys())

# print (ec_nums)

# sample_d = {k: d[k] for k in keys}
# ec_nums = [x.replace(".", "_") for x in BRENDA_PARSER.keys()]


# If there is no EC numbers txtfile or we want to update it
if UPDATE_EC_NUMS or not os.path.isfile("config/ec_nums.py"):

	print ('Updating the EC numbers')
	from brendapy import BrendaParser
	BRENDA_PARSER = BrendaParser()

	keys = [x for x in BRENDA_PARSER.keys() if x != None]


	with open ('config/ec_nums.py', 'w+') as ec_nums_file:
		ec_nums_file.write (f'ec_nums = {keys}')


from config.ec_nums import ec_nums

# We're going to use these as folder names so get rid of full stops
ec_nums = [x.replace(".", "_") for x in ec_nums]

# Just grab a sample of EC numbers
ec_nums = ['3_5_2_6']


print (ec_nums)

rule all:
        input:
        	expand(WORKDIR + "/{ec_num}/files/{ec_num}.nwk", ec_num = ec_nums),
        	expand(WORKDIR + "/{ec_num}/files/{ec_num}_annotations.txt", ec_num = ec_nums), expand(WORKDIR + "/{ec_num}/files/{ec_num}_filt.nwk", ec_num = ec_nums)

# Get the BRENDA annotations for a given EC number
rule get_brenda_annotations:
    output:
    	WORKDIR + "/{ec_num}/csv/{ec_num}_brenda.csv"
    script:
        "scripts/get_brenda_annotations.py"

# Get the UniProt annotations for the entries in BRENDA
rule get_uniprot_annotations:
    input:
    	WORKDIR + "/{ec_num}/csv/{ec_num}_brenda.csv",
    	workdir = WORKDIR
    output:
    	WORKDIR + "/{ec_num}/csv/{ec_num}_uniprot.csv"
    script:
        "scripts/get_uniprot_annotations.py"

# Write out the sequences to file
rule create_fasta:
	input:
		WORKDIR + "/{ec_num}/csv/{ec_num}_uniprot.csv"
	output:
		WORKDIR + "/{ec_num}/files/{ec_num}.fasta"
	script:
	    "scripts/create_fasta.py"

# Align sequences for a given EC number
rule align_sequences:
	input:
		WORKDIR + "/{ec_num}/files/{ec_num}.fasta"
	output:
		WORKDIR + "/{ec_num}/files/{ec_num}.aln"
	shell:
		"mafft --reorder {input} > {output}"

# Make a tree for a given EC number
rule infer_tree:
    input:
    	WORKDIR + "/{ec_num}/files/{ec_num}.aln"

    output:
    	WORKDIR + "/{ec_num}/files/{ec_num}.nwk"

    shell:
        "FastTree {input} > {output}"

rule create_annotation_file:
	input:
		csv = WORKDIR + "/{ec_num}/csv/{ec_num}_uniprot.csv",
		annotation_cols = config['annotation_cols']
	output:
		tsv = WORKDIR + "/{ec_num}/files/{ec_num}_annotations.txt"
	script:
		"scripts/create_annotation_file.py"

rule create_filtered_file:
	input:
		fasta = WORKDIR + "/{ec_num}/files/{ec_num}.fasta",
		csv = WORKDIR + "/{ec_num}/csv/{ec_num}_uniprot.csv",	
	params:
		filter_parameters = config['filter_parameters'],
		row_num = config['row_num']
	output:
		filtered = WORKDIR + "/{ec_num}/files/{ec_num}_filt.fasta"
	script:
		"scripts/auto_filtering.py"

rule create_filtered_tree:
	input:
		WORKDIR + "/{ec_num}/files/{ec_num}_filt.fasta"
	output:
		WORKDIR + "/{ec_num}/files/{ec_num}_filt.aln"
	shell:
		"mafft --reorder {input} > {output}"

rule infer_filtered_tree:
    input:
    	WORKDIR + "/{ec_num}/files/{ec_num}_filt.aln"

    output:
    	WORKDIR + "/{ec_num}/files/{ec_num}_filt.nwk"

    shell:
        "FastTree {input} > {output}"