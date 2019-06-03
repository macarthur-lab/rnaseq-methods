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
    current_gene_start = -1

    with open(splice_file, 'r') as f:
        for splice_event in f:
            splice_event = splice_event.rstrip('\n').split()
            gene = splice_event[0]
            chrom = splice_event[3]
            start = splice_event[4]
            end = splice_event[5]
            key = chrom + ':' + start + '-' + end

            if current_gene_start == 1:
                current_gene = gene
                current_gene_start = 0

            if key not in splices:
                splices = np.append(splices, key)
                counts = np.append(counts, 0)

            counts[np.where(splices == key)[0][0]] += 1

            if gene != current_gene:
                gene_total = np.sum(counts[current_gene_start:(counts.size - 1)])
                counts[current_gene_start:(counts.size - 1)] = np.divide(counts[current_gene_start:(counts.size - 1)], gene_total)
                current_gene_start = counts.size - 1
                current_gene = gene

    # account for very last gene
    gene_total = np.sum(counts[current_gene_start:counts.size])
    counts[current_gene_start:counts.size] = np.divide(counts[current_gene_start:counts.size], gene_total)

    return splices, counts
