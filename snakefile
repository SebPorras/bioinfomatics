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
# ec_nums = ['4_6_1_1','6_3_1_2', '2_7_4_1', '2_5_1_78', '3_1_1_3', '1_10_3_2', '2_4_1_16',
# '5_2_1_8', '2_7_12_2', '1_11_1_15', '1_2_1_12', '1_9_3_1', '3_1_1_73', '5_3_4_1', '3_6_4_12', '2_4_1_17',
# '3_1_3_48', '4_2_2_2', '4_3_1_24', '6_2_1_12', '2_7_1_1', 

ec_nums = ['2_7_7_7', '3_2_1_55', '3_1_1_1', '2_7_7_27', 
'1_5_1_20', '4_2_1_11', '3_2_1_18', '1_1_1_219', '2_3_2_27', '2_7_11_17', '3_1_3_2', '3_2_1_39', '2_7_11_13',
'2_7_7_48', '3_1_3_8', '3_1_4_11', '3_2_1_51', '3_4_21_92','2_7_4_6']
print (ec_nums)

rule all:
        input:
        	expand(WORKDIR + "/{ec_num}/files/{ec_num}_filt.nwk", ec_num = ec_nums),
        	expand(WORKDIR + "/{ec_num}/files/{ec_num}_annotations.txt", ec_num = ec_nums),
        	expand(WORKDIR + "/{ec_num}/files/{ec_num}_filt.nwk", ec_num = ec_nums),

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
# rule create_fasta:
# 	input:
# 		WORKDIR + "/{ec_num}/csv/{ec_num}_uniprot.csv"
# 	output:
# 		WORKDIR + "/{ec_num}/files/{ec_num}.fasta"
# 	script:
# 	    "scripts/create_fasta.py"

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