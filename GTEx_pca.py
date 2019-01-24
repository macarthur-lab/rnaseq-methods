import os.path
import sys
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler, Normalizer
from sklearn.decomposition import PCA
import pickle 
import datetime

def run_pca(file, transpose=True, tissues = None, tissues_of_interest=None, log_transform = True, scaler = 'standard', n_components = None, save_pca = True, pca_name = "GTEx_pca"):
	""" Given a set of gene expression values, run an exploratory Principal Component Analysis (PCA) and store the loadings.

	:param file (str):			file path of tab separated gene expression values
	:param transpose (bool):		transpose gene expression values matrix (default True)
	:param tissues (str):			file path of tissues corresponding to each sample [id, tissue, subtissue] (optional)
	:param tissues_of_interest (str array):	tissues in which to perform PCA (optional)
	:param log_transform (bool):		perform log2(gene_expression + 1) transformation (default True)
	:param scaler (str):			data scaling method. If 'None', data is not scaled. Accepted scalers 'standard', 'robust', 'normalizer' (default 'standard')
	:param n_components (int):		number of principal components to store. If 'None', all are stored (default None)
	:param save_pca (bool):			store sklearn.decomposition.PCA object fitted to data as pickle object (default True) 			
	:param pca_name (str):			basename of files to be stored (default 'GTEx_pca')

	"""

	if not os.path.exists(file):	
		print("Please provide a valid file")
		sys.exit()

	print(now() + ": Loading " + file)
	tpms = np.loadtxt(file, delimiter='\t')

	print(now() + ": Transposing")
	if transpose:
		tpms = tpms.transpose()
	
	if tissues_of_interest is not None:
		if not os.path.exists(tissues):
			print("Please provide a valid tissue file")
			sys.exit()

		tissues = np.loadtxt(tissues, delimiter='\t', dtype='str')
		used = [i for i in range(0, len(tissues)) if tissues[i,1].strip() in tissues_of_interest]
		tpms = tpms[used,:]		
		if tpms.size == 0:
			print("No tissues of interest were found on " + tissues)
			sys.exit()

		np.savetxt(pca_name + "_used_idx", used)
		np.savetxt(pca_name + "_used_tissues.tsv", tissues[used, :], delimiter='\t', fmt='%s')
	
	
	if log_transform:
		print(now() + ": Performing log2 transformation")
		tpms = np.log2(tpms + 1)

	if scaler is not None:
		if scaler == 'standard':
			# (x_i - mean(X))/stdev(X)
			scaler = StandardScaler()
		elif scaler == 'robust':
			# (x_i - Q_1(X))/(Q_3(X) - Q_1(X))
			scaler = RobustScaler()
		elif scaler == 'normalizer':
			mean_scaler = StandardScaler(with_mean = True, with_std = False)
			mean_scaler.fit(tpms)
			tpms = mean_scaler.transform(tpms)

			# x_i / sqrt(x_i, y_i, z_i, ...)
			scaler = Normalizer()
		else:
			print(scaler + " not recognized. Using standard scaler.")
			scaler = StandardScaler()

		print(now() + ": Scaling data")
		scaler.fit(tpms)
		tpms = scaler.transform(tpms)
	
	print(now() + ": Performing PCA")
	pca = PCA(n_components=n_components, random_state = 42, svd_solver = 'full')
	pca.fit(tpms)
	
	if save_pca:
		print(now() + ": Saving PCA object")
		file_name = str(pca_name) + ".object"
		open(file_name, 'wb')
		pickle.dump(pca, file_name, protocol=4)
		file_name.close()
	
	print(now() + ": Transforming principal components")
	pcs = pca.transform(tpms)
	
	print(now() + ": Saving loadings")
	pcs_file_name = pca_name + "_pcs.tsv"
	np.savetxt(pcs_file_name, pcs, delimiter='\t')
	
def now():
	return(str(datetime.datetime.now()))
