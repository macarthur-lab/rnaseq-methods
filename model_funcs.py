import numpy as np

def detect_and_count(splice_file):
    """Given a file with splice sites, calculate the proportion in which they occur in each gene
        
    :param splice_file (tab separated file):        splice file

    :r splices (numpy array):       location of splice sites in genome (chrom:start-end)
    :r prop (numpy array):          proportion of splice occurence in gene
    """
    splices = np.asarray([], dtype='str')
    counts = np.asarray([])

    current_gene = ''
    current_gene_start = 0
    i = -1
    with open(splice_file, 'r') as f:
        for splice_event in f:
            i += 1
            splice_event = splice_event.rstrip('\n').split()
            gene = splice_event[0]
            chrom = splice_event[3]
            start = splice_event[4]
            end = splice_event[5]
            key = chrom + ':' + start + '-' + end

            if i == 0:
                current_gene = gene

            if key not in splices:
                splices = np.append(splices, key)
                counts = np.append(counts, 0)

            counts[np.where(splices == key)[0][0]] += 1

            if gene != current_gene and i != 0:
                gene_total = np.sum(counts[current_gene_start:i])
                counts[current_gene_start:i] = np.divide(counts[current_gene_start:i], gene_total)
                current_gene_start = i
    # account for very last gene
    gene_total = np.sum(counts[current_gene_start:i+1])
    # print('Current gene start: ' + str(current_gene_start) + 'i: ' + str(i))
    counts[current_gene_start:i] = np.divide(counts[current_gene_start:i+1], gene_total)
    print(counts)
    print(gene)
    print('')

    return splices, counts
