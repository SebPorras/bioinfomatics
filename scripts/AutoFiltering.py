import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import numpy as np
from sequence import * 



def display_col_tags(ec_num:str, col_name: str):
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    counts = df[col_name].dropna().unique()
    
    return counts 

    
def filtered_annots(ec_num:str, filtered_entries: list, enzyme_prop: str):

    for key in filtered_entries.keys():
      
        df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")
    
        df = df.loc[df['Entry'].isin(filtered_entries.get(key))]

        num_rows = df.shape[0]

        if enzyme_prop in df.columns:
                annot_count = df[enzyme_prop].dropna().shape[0]
                percent_annot = round((annot_count/num_rows) * 100, 2)
                print(f"{ec_num}: {enzyme_prop}: {percent_annot}% of entries annotated in {key}")


#enter ec num and threshold to define a cutoff for similarity when grouping tags 
def remove_outliers(ec_num: str, data_col:str, entry_limit = 0):
    """
    Reads in all seqs from an ec_group, 
    
    rm_outliers: condition that will remove any groups 
    that only have one entry in them 
    
    data_col (str): The particular data column being compared for the sequences. 
    e.g. Cross_reference_InterPro or Cross_reference_SMART 
    
    entry_limit: the minimum number of entries a tag must have to be included 
    """

    enzyme_cols = ['Entry', data_col]

    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

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
    
    
def grab_seqs_jaccard(ec_num: str, similar_groups: list, col_name: str, index: int):
    
    """
    Grabs all sequences that have been put together in a group
    that is deemed to be similar to each other based on the calc_jaccard_matrix 
    method 
    
    similar_groups (list): uses the output from calc_jaccard_matrix to grab all entries 
    that have tags in that group 
    """
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")
    
    df = df.loc[df[col_name].isin(similar_groups[index])]

    df = df['Entry']
        
    df = df.values.tolist() 

    return df


    
def grab_seqs(ec_num: str, data_col:str, pfam_tag: str):
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")
    
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


def calc_jaccard_matrix_auto(ec_num: str, threshold: float, filtered_counts, data_col):
  
    enzyme_cols = ['Entry', data_col]

    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

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
            
            #calculate Jaccard similarity 
            similarity = len(intersection)/len(union)
            
            #updates matrix 
            distmat[i, j] = distmat[j, i] = similarity 
    
    #record how many sequences have each type of IP tag, use for axis ticks 
    interpro_tag_counts = [str(value) + ' - '  \
                           + str(key)for key, value in filtered_counts.items()]
    
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

def display_single_families(ec_num: str, gene_tag:str, data_col: str):
    """
    data col (str): particular filtering group Gene3D, SMART etc
    """
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    enzyme_cols = ['Entry', 'Protein_families', data_col]
    
    #will count each appearance of particular family  
    counters = defaultdict(int)

    df = df[enzyme_cols]
    
    #filters data frame to only include entries with tags matching the HOGENOM tag 
    filtered_df = df.loc[(df[data_col] == gene_tag)]
  
    #record appearances of each family #draft
    
    family_counts = filtered_df['Protein_families'].value_counts()
    
    #place these values in a dict so that it's easier to plot below 
    for index, value in family_counts.items():
        counters[index] = value 

    return counters

def create_filtered_ec(ec_num: str, filtered_entries: list):

    original_fa = readFastaFile(f"/home/seb-porras/expat_bench/workflows/{ec_num}/files/{ec_num}.fasta")

    filtered = [seq for seq in original_fa if seq.name in filtered_entries]

    writeFastaFile(f"/home/seb-porras/expat_bench/workflows/{ec_num}/files/{ec_num}_filt.fasta", filtered)


def auto_generate_filtered_entries(ec_num: str, threshold: int, entry_limit = 0):
    """
    ec_num (str): string version of enzyme number 
    
    threshold (int): Cutoff between 0 and 1 for how similar SMART & InterPro tags must 
    be to each other to be classed as similar
    
    enzyme_prop (str): enzyme property that user wants to have annotated, e.g. BRENDA_KM
    
    entry_limit: the minimum number of entries a interpro or SMART tag must have to be included 
    in the Jaccard Similarity matrix. Used to help remove outliers
    """
    
    #will store a list of entries for each type of filtering used 
   
    #same as above except with InterPro 
    filtered_ip = remove_outliers(ec_num, 'Cross_reference_InterPro', entry_limit)

    ip_matrix = calc_jaccard_matrix_auto(ec_num, threshold, filtered_ip, 'Cross_reference_InterPro')  
    

    #Once user is happy with selections, record which entries are included 
    target_ip = grab_seqs_jaccard(ec_num, ip_matrix, 'Cross_reference_InterPro', 0)
    
  
    #remove groups below the minimum number of entries within the EC group 
    filtered_smart = remove_outliers(ec_num, 'Cross_reference_SMART', entry_limit)

    #create jaccard matrix from filtered list 
    smart_matrix = calc_jaccard_matrix_auto(ec_num, threshold, filtered_smart, 'Cross_reference_SMART')
    
    for i in range(len(smart_matrix)):
     
        target_smart = grab_seqs_jaccard(ec_num, smart_matrix, 'Cross_reference_SMART', i)
        
        print(compare_overlap(target_ip, target_smart))
        print(target_ip)

        if compare_overlap(target_ip, target_smart) > 0.8:
            
            target_ip.extend(target_smart)


    tags = display_col_tags(ec_num, 'Cross_reference_Gene3D')

    for i in tags:

      
        target_g3d = grab_seqs(ec_num, "Cross_reference_Gene3D", i)
        print(target_ip)


        print(compare_overlap(target_ip, target_g3d))

        if compare_overlap(target_ip, target_g3d) > 0.8:
            target_ip.extend(target_g3d)

    tags = display_col_tags(ec_num, 'Cross_reference_Pfam')

    for i in tags:
       
        target_pfam = grab_seqs(ec_num, "Cross_reference_Pfam", i)

        print(compare_overlap(target_ip, target_pfam))
        print(target_ip)


        if compare_overlap(target_ip, target_pfam) > 0.8:
            target_ip.extend(target_pfam)

    
    tags = display_col_tags(ec_num, 'Cross_reference_OrthoDB')

    for i in tags:

        print(target_ip)

        target_ortho = grab_seqs(ec_num, "Cross_reference_OrthoDB", i)
    
        print(compare_overlap(target_ip, target_ortho))

        if compare_overlap(target_ip, target_ortho) > 0.8:
            target_ip.extend(target_ortho)
    
   
    return target_ip

if __name__ == "__main__":
    #print(snakemake.input[0]
    #test = auto_generate_filtered_entries("3_5_2_6", 0.15, 0)
    #create_filtered_ec(3_5_2_6, test)
    test = auto_generate_filtered_entries("3_2_1_14", 0.5, 1)
    create_filtered_ec("3_2_1_14", test)    

