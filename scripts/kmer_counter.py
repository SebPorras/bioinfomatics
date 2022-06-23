from collections import defaultdict
from scipy.stats import fisher_exact

oddsratio, pvalue = fisher_exact([[8, 2], [1, 5]])

print(pvalue)


def create_kmers(k_len: int, seq: str) -> dict:
    """
    k_len(int): specify a Kmer length 

    seq (str): can be amino acids or proteins 
    """

    kmers = defaultdict()

    for i in range(len(seq) - k_len + 1): 
        
        kmers[seq[i:i+k_len]] = 0

    return kmers

def count_kmers(k_len: int, seq: str):

    kmers = create_kmers(k_len, seq)

    for k in kmers.keys():
        print(k)


seq = "ATATGTCA"
test = create_kmers(2, seq)

#count_kmers(2, seq)




def kmer_counter(k, seq):
    '''
    Return a list of the number of times each possible k-mer appears
    in seq, including overlapping occurrences.
    '''

    counts = {}
    for i in range(0, len(seq)-k+1):
        
        kmer = seq[i:i+k]
        current_count = counts.get(kmer, 0)
        counts[kmer] = current_count + 1


    return list(counts.items())

thing = kmer_counter(3, seq)

print(thing)