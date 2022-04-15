import sequence # import the binfpy sequence package
# here we load all reviewed sequences and create a dictionary which will allow us to recall each by its name
findmyseq = {}
seqs = sequence.readFastaFile('tpswb_reviewed.fa', sequence.Protein_Alphabet, parse_defline=False)
print('Loaded', len(seqs), 'sequences')
for s in seqs:
    findmyseq[s.name] = s
# next we load each of the sub-family specific FASTA files to collect the names of each sub-family
subfamily = ['a','b','c','d','e','f'] # all sub-families
subfnames = {} # a dictionary with names of members of each sub-family
for f in subfamily:
    fseqs = sequence.readFastaFile('tps' + f + '.fa', sequence.Protein_Alphabet, parse_defline=False)
    subfnames[f] = []
    for s in fseqs:
        subfnames[f].append(s.name)
    print('Loaded', len(fseqs), 'sequences for sub-family tps' + f)
# finally we revisit the original sequence data to amend the name of those entries which belong to a sub-family
for s in seqs:
    for f in subfamily:
        if s.name in subfnames[f]:
            s.info = s.info + ' SF=' + f # add another field to the header-line
            s.name = s.name + '|Tps' + f # extend the sequence name (which is used later)
            break
# save the sequences as a new FASTA file
sequence.writeFastaFile('tpswb_v1.fa', seqs)
