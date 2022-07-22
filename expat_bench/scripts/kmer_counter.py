from collections import defaultdict
from sequence import *
from scipy.stats import fisher_exact

#oddsratio, pvalue = fisher_exact([[8, 2], [1, 5]])



def kmer_counter(filtered_seqs, k):
    '''
    Return a list of the number of times each possible k-mer appears
    in seq, including overlapping occurrences.
    '''
    counts = {}
    for seq in filtered_seqs:
        for i in range(0, len(seq)-k+1):
            
            kmer = seq[i:i+k]
            current_count = counts.get(kmer, 0)
            counts[kmer] = current_count + 1

    return counts

def positve_kmer_counts(ec_num, k):
    
    filtered_seqs = readFastaFile(f"/home/seb-porras/expat_bench/workflows/{ec_num}/files/{ec_num}_filt.fasta")
    
    num_seqs = len(filtered_seqs)

    counts = kmer_counter(filtered_seqs, k)

    return counts, num_seqs

def negative_kmer_counts(ec_num, k):
    
    unfiltered_seqs = readFastaFile(f"/home/seb-porras/expat_bench/workflows/{ec_num}/files/{ec_num}.fasta")
    
    num_seqs = len(unfiltered_seqs)

    counts = kmer_counter(unfiltered_seqs, k)

    return counts, num_seqs

positive, pos_length = positve_kmer_counts('3_5_2_6', 5)
negative, neg_length = negative_kmer_counts('3_5_2_6', 5)


kmer_pvals = defaultdict()

for key, value in positive.items():
    
    a = value
    b = pos_length

    c = negative[key]
    d = neg_length - c

    oddsratio, pvalue = fisher_exact([[a, b], [c, d]])

    kmer_pvals[key] = pvalue

