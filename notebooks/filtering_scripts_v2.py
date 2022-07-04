import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import defaultdict
import numpy as np


def remove_outliers(ec_num: str, data_col:str, entry_limit) -> dict:
    """
    Reads in all seqs from an ec_group, 
    
    ec_num (str): enzyme classification number to access files 

    rm_outliers: condition that will remove any groups 
    that only have one entry in them 
    
    data_col (str): The particular data column being compared for the sequences. 
    e.g. Cross_reference_InterPro 
    
    entry_limit: the minimum number of entries a interpro must have 
    to be included in the Jaccard Similarity matrix. 
    
    Returns:
        a dictionary mapping InterPro tags to number of seqs with this tag
    """

    enzyme_cols = ['Entry', data_col]

    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    #will count each appearance of particular interpro class 
    df = df[enzyme_cols]
    
    interpro_counts = df[data_col].dropna().value_counts()
    
    filtered_interpro = dict()
    
    #filter out distant entries 
    for key, value in interpro_counts.items():
        if value > entry_limit:
            filtered_interpro[key] = value
       
    return filtered_interpro 

def calc_jaccard_matrix(ec_num: str, threshold: float, filtered_counts: dict) -> list:
    
    """
    Sorts entries based on their Interpro tag and creates a heat map to visualise.

    ec_num (str): enzyme classification number to access files 

    threshold (float): Set by user to determine how similar two tags must be 

    filtered_counts (dict): filtered list of tags that has removed outliers 

    Returns:
        A list containing each row with IP tags that are similar to each other
    
    """  
  
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    #splits interpro tags and creates a set to allow comparison of elements   
    tag_list = [set(tag[:-1].split(';')) for tag in filtered_counts.keys()]
    
    #initialise the similarity matrix
    distmat = np.zeros((len(tag_list), len(tag_list)))

    #cacluate jaccard similarity at for each tag pair 
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
    plt.title(f'{ec_num}: InterPro tags - Jaccard Similarity')
    plt.show()
    
    #IP names stored so loop below can reference name with matrix indices 
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


def display_annotations(ec_num: str, similar_groups: list,\
    index: int, data_col: str) -> None:
    """
    Provides a breakdown of what protein families are present 
    for a particular InterPro tag and all tags that are similar to it 

    ec_num (str): enzyme classification number to access files 

    similar_groups(list): list containing groupings of tags that are
    similar to each other 
    
    index(int): user specifies which group (i.e. row on the graphical heatmap)
    they would like to see a breakdown of

    data_col (str): name of tag list used, e.g. Interpro tags     
    """
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")

    enzyme_cols = ['Entry', data_col, 'Cross_reference_InterPro']
    
    #stores which group of tags user wants to look at 
    ip_ref = similar_groups[index]

    #will count each appearance of particular family from the IP groups  
    counters = defaultdict(int)

    df = df[enzyme_cols]
 
    #filters data frame to only include entries with tags in the similarity group 
    df = df.loc[df['Cross_reference_InterPro'].isin(ip_ref)]

    #record appearances particular property 
    family_counts = df[data_col].value_counts()
   
    #place these values in a dict so that it's easier to plot below 
    for index, value in family_counts.items():
        counters[index] = value 
    
    #plot a histogram of family counts
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.bar(counters.keys(), counters.values())
    plt.xticks(rotation='vertical')
    ax.set_xlabel(data_col)
    ax.set_ylabel('# of entries')
    ax.set_title(f'Different vluaes in {ec_num} using {data_col}')
    plt.show()

def grab_seqs_jaccard(ec_num: str, similar_groups: list, col_name: str, index: int):
    
    """
    Grabs all sequences that have been put together in a group that is deemed 
    to be similar to each other based on calc_jaccard_matrix() 

    ec_num (str): enzyme classification number to access files ec_num(str): number used to reference groups in the data set

    similar_groups(list): list containing groupings of tags that are
    similar to each other 

    col_name (str): The name of the column that sequences will be 
    selected by. e.g. InterPro tag 

    index(int): user specifies which group (i.e. row on the 
    graphical heatmap) they would like to see a breakdown of

    """
    
    df = pd.read_csv(f"../workflows/{ec_num}/csv/{ec_num}_uniprot.csv")
    
    df = df.loc[df[col_name].isin(similar_groups[index])]

    df = df['Entry']
        
    df = df.values.tolist() 

    return df

def calculate_annots(ec_num: str ,col_name: str) -> None:

    """Determines how many entries are annotated with 
    a chosen property and prints this to the user.
    
    ec_num (str): enzyme classification number to access files 

    col_name (str): The name of the column that sequences will be 
    selected by. e.g. InterPro tag 
    """
 
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

def generate_filtered_entries(ec_num: str, threshold: float, index:int, \
    enzyme_prop: str, entry_limit = 0) -> None:
    """
    ec_num (str): enzyme classification number to access files  
    
    threshold (int): Cutoff between 0 and 1 for how similar
    InterPro tags must be to each other to be classed as similar

    index(int): user specifies which group (i.e. row on the 
    graphical heatmap) they would like to see a breakdown of
    
    enzyme_prop (str): enzyme property of interest e.g. BRENDA_KM
    
    entry_limit: the minimum number of entries a interpro must have 
    to be included in the Jaccard Similarity matrix. 
    
    """

    #how many entries contain the enzyme property within EC group 
    calculate_annots(ec_num, enzyme_prop)  
    
    #creates filtered list of EC groups 
    filtered_ip = remove_outliers(ec_num, 'Cross_reference_InterPro', entry_limit)

    #constructs the IP heatmap and returns list of tags that are similiar 
    ip_matrix = calc_jaccard_matrix(ec_num, threshold, filtered_ip)     
    
    #display filtered protein families 
    display_annotations(ec_num, ip_matrix, index, 'Protein_families')

    display_annotations(ec_num, ip_matrix, index, 'BRENDA_KM')

    display_annotations(ec_num, ip_matrix, index, 'BRENDA_IC50')

    display_annotations(ec_num, ip_matrix, index, 'BRENDA_KI')