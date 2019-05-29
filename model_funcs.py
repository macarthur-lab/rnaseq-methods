import numpy as np

def detect_and_count(splice_file):
	"""Given a file with splice sites, calculate the proportion in which they occur in each gene
	
	:param splice_file (tab separated file):	splice file

	:r splices (numpy array):	location of splice sites in genome (chrom:start-end)
	:r prop (numpy array):		proportion of splice occurence in gene
	"""
	splices = np.asarray([], dtype='str')
	counts = np.asarray([])

	with open(splice_file, 'r') as f:
		for splice_event in f:
			splice_event = splice_event.rstrip('\n').split()
			chrom = splice_event[3]
			start = splice_event[4]
			end = splice_event[5]
			key = chrom + ':' + start + '-' + end

			if key not in splices:
				splices = np.append(splices, key)
				counts = np.append(counts, 0)
			
			counts[np.where(splices == key)[0][0]] += 1

	total = np.sum(counts)
	prop = np.divide(counts, total)

	return splices, prop
