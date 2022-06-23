import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import numpy as np


def display_col_tags(ec_num:str, col_name: str):
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    counts = df[[col_name]].value_counts()
    print(type(counts))
    print(counts)

def calculate_annots(ec_num ,col_name):
 
    print(col_name)

    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    num_rows = df.shape[0]

    if col_name in df.columns:
            annot_count = df[col_name].dropna().shape[0]
            percent_annot = round((annot_count/num_rows) * 100, 2)
            print(f"{ec_num} {col_name}: {percent_annot}% of entries annotated")
        
    else:
        print(f"{ec_num} does not have {col_name} column")

    print("\n")

    
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
    
    
def calc_jaccard_matrix(ec_num: str, threshold: float, filtered_counts, data_col):
  
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
    
    #plot heatmat for reference 
    fig, ax = plt.subplots()
    plt.imshow(distmat, plt.cm.inferno, interpolation='nearest')
    plt.colorbar()
    plt.yticks(np.arange(len(tag_list)), interpro_tag_counts)
    plt.xticks(np.arange(len(tag_list)), [x for x in filtered_counts.keys()], rotation=90)
    plt.subplots_adjust(bottom=0.35)
    plt.title(f'{ec_num}: {data_col} Jaccard Similarity')
    plt.show()
    
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



def display_grouped_families(ec_num, similar_groups, index, data_col: str):
    """
    similar_groups: list containing groupings of tags that are similar to each other 
    
    index: user specifies which group they would like to see on a graph 
    
    """
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    enzyme_cols = ['Entry', 'Protein_families', data_col]
    
    #stores which group of tags user wants to look at 
    ip_ref = similar_groups[index]

    #will count each appearance of particular family from the IP groups  
    counters = defaultdict(int)

    df = df[enzyme_cols]
 
    #filters data frame to only include entries with tags in the similarity group 
    df = df.loc[df[data_col].isin(ip_ref)]

    #record appearances of each family 
    family_counts = df['Protein_families'].value_counts()
   
    #place these values in a dict so that it's easier to plot below 
    for index, value in family_counts.items():
        counters[index] = value 
    
    #plot a histogram of family counts
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.bar(counters.keys(), counters.values())
    plt.xticks(rotation='vertical')
    ax.set_xlabel('Family')
    ax.set_ylabel('# of entries')
    ax.set_title(f'Different families in {ec_num} using {data_col}')
    plt.show()
    
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
    
    #plot a histogram of family counts
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.bar(counters.keys(), counters.values())
    plt.xticks(rotation='vertical')
    ax.set_xlabel('Family')
    ax.set_ylabel('# of entries')
    ax.set_title(f'Different families in {ec_num} using {gene_tag}')
    plt.show()
    
def grab_seqs(ec_num: str, data_col:str, pfam_tag: str):
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")
    
    df = df.loc[df[data_col] == pfam_tag]

    df = df['Entry']
        
    df = df.values.tolist()     
    
    return df
    

def calc_total_overlap(filtered_entries: list):

    sets = {key: set(value) for key, value in filtered_entries.items()}

    for key, value in sets.items():
        print(f'{key}: entries = {len(value)}')

    common_vals = [value for value in sets.values()]

    total_overlap = set.intersection(*common_vals)

    print(f'Total overlap between all groups {len(total_overlap)} \n')
    
    return total_overlap 


def compare_overlap(filtered_entries: list, group_1, group_2):
    
    g1 = set(filtered_entries[group_1])
    
    g2 = set(filtered_entries[group_2])
    
    overlap = g1.intersection(g2)
    
    print(f'{group_1}: entries = {len(filtered_entries[group_1])}')
    
    print(f'{group_2}: entries = {len(filtered_entries[group_2])}')
    
    
    print(f'Total overlap between {group_1} & {group_2} is {len(overlap)}\n')

    return overlap


def generate_filtered_entries(ec_num: str, threshold: int, enzyme_prop: str ,entry_limit = 0):
    """
    ec_num (str): string version of enzyme number 
    
    threshold (int): Cutoff between 0 and 1 for how similar SMART & InterPro tags must 
    be to each other to be classed as similar
    
    enzyme_prop (str): enzyme property that user wants to have annotated, e.g. BRENDA_KM
    
    entry_limit: the minimum number of entries a interpro or SMART tag must have to be included 
    in the Jaccard Similarity matrix. Used to help remove outliers
    
    """
    
    #will store a list of entries for each type of filtering used 
    filtering_stages = {}
    
    #show how many entries contain the chosen enzyme property within the EC group 
    calculate_annots(ec_num, enzyme_prop)  
    
    #remove groups below the minimum number of entries within the EC group 
    filtered_smart = remove_outliers(ec_num, 'Cross_reference_SMART', entry_limit)

    #create jaccard matrix from filtered list 
    smart_matrix = calc_jaccard_matrix(ec_num, threshold, filtered_smart, 'Cross_reference_SMART')

    #same as above except with InterPro 
    filtered_ip = remove_outliers(ec_num, 'Cross_reference_InterPro', entry_limit)

    ip_matrix = calc_jaccard_matrix(ec_num, threshold, filtered_ip, 'Cross_reference_InterPro')     
        
    
    filtering_loop = True 
    
    #next, compare the distributions of families until user deems them to be similar enough 
    while filtering_loop:
        
        choose_smart_tag = int(input('Choose a SMART row to grab entries from: '))
        
        choose_ip_tag = int(input('Choose a InterPro row to grab entries from: '))
        
        #display SMART entries of chosen tags and groups that are similar 
        if len(smart_matrix) > 0:
            display_grouped_families(ec_num, smart_matrix, choose_smart_tag, 'Cross_reference_SMART')
        
        #display interpro entries of chosen tags and groups that are similar 
        display_grouped_families(ec_num, ip_matrix, choose_ip_tag, 'Cross_reference_InterPro')
        
        end_comp = input('Happy with selection(Y/N): ')
        
        if end_comp.upper() == "Y":
            filtering_loop = False 
        
    #Once user is happy with selections, record which entries are included 
    target_ip = list(grab_seqs_jaccard(ec_num, ip_matrix, 'Cross_reference_InterPro', choose_ip_tag))
    
    filtering_stages['IP'] = target_ip
    
    print(f'Number of sequences after IP filtering = {len(target_ip)}')
    if len(smart_matrix) > 0:
        target_smart = list(grab_seqs_jaccard(ec_num, smart_matrix, 'Cross_reference_SMART', choose_smart_tag))

        filtering_stages['SMART'] = target_smart

    
    if len(smart_matrix) > 0:
        print(f'Number of sequences after SMART filtering = {len(target_smart)}')
    

    compare_pfam = False 
    
    #display how many entries have Gene3D annot in the EC group 
    calculate_annots(ec_num, 'Cross_reference_Pfam')  
    
    check_pfam = input('Compare entries to Pfam filtering? (Y/N): ')
    
    if check_pfam.upper() == "Y":
        compare_pfam = True 
        
    while compare_pfam: 
        
        #show user what groups are within Gene3D
        display_col_tags(ec_num, 'Cross_reference_Pfam')

        gene_tag = input('Enter a Gene3D tag: ')
        
        #Show protein families within the chosen G3D tag 
        display_single_families(ec_num, gene_tag, 'Cross_reference_Pfam')

        check_again = input('Would you like to choose another group? (Y/N): ')
        
        if check_again.upper() == 'N':

            compare_pfam = False
            
            #record entries once user is happy with selection based on protein families 
            target_pfam = list(grab_seqs(ec_num, "Cross_reference_Pfam", gene_tag))
            
            filtering_stages['Pfam'] = target_pfam
            
            print(f'The number of entries in the Gene3D group are {len(target_pfam)} \n')
    
    #next create a list of relevent entries based on Gene3D tags 
    
    compare_gene3D = False 
    
    #display how many entries have Gene3D annot in the EC group 
    calculate_annots(ec_num, 'Cross_reference_Gene3D')  
    
    check_3d = input('Compare entries to gene3D filtering? (Y/N): ')
    
    if check_3d.upper() == "Y":
        compare_gene3D = True 
        
    while compare_gene3D: 
        
        #show user what groups are within Gene3D
        display_col_tags(ec_num, 'Cross_reference_Gene3D')

        gene_tag = input('Enter a Gene3D tag: ')
        
        #Show protein families within the chosen G3D tag 
        display_single_families(ec_num, gene_tag, 'Cross_reference_Gene3D')

        check_again = input('Would you like to choose another group? (Y/N): ')
        
        if check_again.upper() == 'N':

            compare_gene3D = False
            
            #record entries once user is happy with selection based on protein families 
            target_g3d = list(grab_seqs(ec_num, "Cross_reference_Gene3D", gene_tag))
            
            filtering_stages['Gene3D'] = target_g3d
            
            print(f'The number of entries in the Gene3D group are {len(target_g3d)} \n')

    #create list of entries based off OrthoDB tags 
    
    compare_ortho_DB = False 
    
    calculate_annots(ec_num, 'Cross_reference_OrthoDB')
    
    check_ortho = input('Compare entries to OrthoDB filtering? (Y/N): ')
    
    if check_ortho.upper() == "Y":
        compare_ortho_DB = True 
    
    while compare_ortho_DB: 

        display_col_tags(ec_num, 'Cross_reference_OrthoDB')

        ortho_tag = input('Enter a OrthoDB tag: ')

        display_single_families(ec_num, ortho_tag, 'Cross_reference_OrthoDB')

        check_again = input('Would you like to choose another group? (Y/N): ')

        if check_again.upper() == 'N':

            compare_ortho_DB = False    
    
            target_ortho = list(grab_seqs(ec_num, "Cross_reference_OrthoDB", ortho_tag))
        
            filtering_stages['Ortho_DB'] = target_ortho

    #display how many entries still contain information for chosen enzyme property 
    filtered_annots(ec_num, filtering_stages, enzyme_prop)
    
    return filtering_stages



    
