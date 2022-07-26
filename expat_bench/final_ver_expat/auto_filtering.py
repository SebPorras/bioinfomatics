import pandas as pd
from collections import defaultdict
import numpy as np
from sequence import * 

#store parameters from the snakefile 
THRESHOLD = float(snakemake.params.filter_parameters[0])
MIN_SEQS = int(snakemake.params.filter_parameters[1])
ROW_NUM = int(snakemake.params.row_num)


def display_col_tags(col_name: str) -> array:

    df = pd.read_csv(snakemake.input.csv)

    counts = df[col_name].dropna().unique()
    
    return counts 
    
#enter ec num and threshold to define a cutoff for similarity when grouping tags 
def remove_outliers(data_col:str, entry_limit = 0) -> dict:
    """
    Reads in all seqs from an ec_group, 
    
    rm_outliers: condition that will remove any groups 
    that only have one entry in them 
    
    data_col (str): The particular data column being compared for the sequences. 
    e.g. Cross_reference_InterPro or Cross_reference_SMART 
    
    entry_limit: the minimum number of entries a tag must have to be included 
    """

    enzyme_cols = ['Entry', data_col]

    df = pd.read_csv(snakemake.input.csv)

    #will count each appearance of particular interpro class 
    counters = defaultdict(int)
    df = df[enzyme_cols]
    
    #draft before removing na values 
    #interpro_counts = df[data_col].value_counts()
    
    interpro_counts = df[data_col].dropna().value_counts()
    
    filtered_interpro = dict()
    
    #filter out distant entries 

    for key, value in interpro_counts.items():
        if value > entry_limit:
            filtered_interpro[key] = value
         
    return filtered_interpro 
    
    
def grab_seqs_jaccard(similar_groups: list, col_name: str, index: int) -> list:
    """
    Grabs all sequences that have been put together in a group
    that is deemed to be similar to each other based on the calc_jaccard_matrix 
    method 
    
    similar_groups (list): uses the output from calc_jaccard_matrix to grab all entries 
    that have tags in that group 

    col_name(str): this specifies the data column being compared. 
    e.g. Protein_families 

    index(int): specifies which row of Jaccard matrix to grab
    """
    
    df = pd.read_csv(snakemake.input.csv)
    
    df = df.loc[df[col_name].isin(similar_groups[index])]

    df = df['Entry']
        
    df = df.values.tolist() 

    return df

def grab_seqs(data_col:str, pfam_tag: str):
    
    df = pd.read_csv(snakemake.input.csv)
    
    df = df.loc[df[data_col] == pfam_tag]

    df = df['Entry']
        
    df = df.values.tolist()     
    
    return df
    
def compare_overlap(group_1, group_2):
    
    g1 = set(group_1)

    g2 = set(group_2)

    overlap = g2.intersection(g1)

    percent_overalp = len(overlap)/len(g2)

    return percent_overalp


def calc_jaccard_matrix_auto(threshold: float, filtered_counts:dict, data_col:str) -> list:
    """
    Creates a list containing all InterPro tags and all other tags that they
    are similar to.
    
    threshold (int): Cutoff between 0 and 1 for how similar SMART & InterPro tags must 
    be to each other to be classed as similar
    
    filtered_counts(dict): 

    data_col (str): The particular data column being compared for the sequences. 
    e.g. Cross_reference_InterPro or Cross_reference_SMART 

    """
  
    enzyme_cols = ['Entry', data_col]

    df = pd.read_csv(snakemake.input.csv)

    #will count each appearance of particular interpro class 
    counters = defaultdict(int)

    #splits interpro tags and creates a set to allow comparison of elements   
    tag_list = [set(tag[:-1].split(';')) for tag in filtered_counts.keys()]
    
    #initialise the similarity matrix
    distmat = np.zeros((len(tag_list), len(tag_list)))

    #cacluate JS at for each tag pair 
    for i in range(len(tag_list)):
        for j in range(i+1, len(tag_list)):
            
            A = tag_list[i]
            B = tag_list[j]

            intersection = A.intersection(B)

            union = A.union(B)
            
            #calculate Jaccard similarity and store
            similarity = len(intersection)/len(union)
            distmat[i, j] = distmat[j, i] = similarity 
    
    #IP names stored here so that loop below can reference name with matrix indices 
    interpro_names = [key for key in filtered_counts.keys()]
    
    #groupings of IP tags that are simililar will be stored here 
    similar_groups = []
    
    #iterate through each row in the matrix 
    for i in range(len(interpro_names)):
      
        #current row being compared to columns 
        key_group = str(interpro_names[i])
        
        #include group being compared as diagnonals will be zero 
        tag_group = [key_group]
        
        #compare every row to each column 
        for j in range(len(interpro_names)):
            
            #every column above or equal to threshold gets recorded 
            if distmat[i, j] >= threshold:
                tag_group.append(interpro_names[j])  
        
        similar_groups.append(tag_group)
        
    return similar_groups


def auto_generate_filtered_entries(threshold: int, entry_limit, row_num):
    """
    threshold (int): Cutoff between 0 and 1 for how similar SMART & InterPro tags must 
    be to each other to be classed as similar
    
    entry_limit: the minimum number of entries a interpro or SMART tag must have to be included 
    in the Jaccard Similarity matrix. Used to help remove outliers
    """
    
    #will store a list of entries for each type of filtering used 
   
    #same as above except with InterPro 
    filtered_ip = remove_outliers('Cross_reference_InterPro', entry_limit)

    ip_matrix = calc_jaccard_matrix_auto(threshold, filtered_ip, 'Cross_reference_InterPro')  

    #list of seqs that will be written to the fasta file 
    final_seqs = set(grab_seqs_jaccard(ip_matrix, 'Cross_reference_InterPro', row_num))

    #reference of sequences for other data base entries to be compared to 
    target_ip = grab_seqs_jaccard(ip_matrix, 'Cross_reference_InterPro', row_num)

    #record what protein families are present after filtering 
    df = pd.read_csv(snakemake.input.csv)
    inital_filter = df.loc[df['Entry'].isin(target_ip)]
    families = inital_filter['Protein_families'].dropna().unique()    
    print(f'Families present: {families}')

    #grab any entry annotated with observed protein families 
    unfiltered = pd.read_csv(snakemake.input.csv)
    chosen_seqs = unfiltered.loc[df['Protein_families'].isin(families)]
    chosen_seq_names = set(chosen_seqs['Entry'])
   
    #update list of entries with any entries that were missed
    final_seqs.update(chosen_seq_names)

    return list(final_seqs)


def create_filtered_ec(filtered_entries: list):

    original_fa = readFastaFile(snakemake.input.fasta)

    filtered = [seq for seq in original_fa if seq.name in filtered_entries]

    writeFastaFile(snakemake.output.filtered, filtered)


print(snakemake.params.filter_parameters)

test = auto_generate_filtered_entries(THRESHOLD, MIN_SEQS, ROW_NUM)
create_filtered_ec(test)    

