import os.path
import sys
import functools
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler, Normalizer
from sklearn.decomposition import PCA
import pickle 
import datetime

print = functools.partial(print, flush=True)

def run_pca(data=None, file_path=None, transpose=True, tissues = None, tissues_of_interest=None, log_transform = True, scaling = 'standard', n_components = None, save_pca = True, pca_name = "GTEx_pca", save_scaler = True):
	""" Given a set of gene expression values, run an exploratory Principal Component Analysis (PCA) and store the loadings.

	:param data (Numpy arr):		data to be analyzed as numpy array (default None)
	:param file_path (str):			file path of tab separated gene expression values (default None)
	:param transpose (bool):		transpose gene expression values matrix (default True)
	:param tissues (str):			file path of tissues corresponding to each sample [id, tissue, subtissue] (optional)
	:param tissues_of_interest (str array):	tissues in which to perform PCA (optional)
	:param log_transform (bool):		perform log2(gene_expression + 1) transformation (default True)
	:param scaling (str):			data scaling method. If 'None', data is not scaled. Accepted scalers 'standard', 'robust', 'normalizer' (default 'standard')
	:param n_components (int):		number of principal components to store. If 'None', all are stored (default None)
	:param save_pca (bool):			store sklearn.decomposition.PCA object fitted to data as pickle object (default True) 			
	:param pca_name (str):			basename of files to be stored (default 'GTEx_pca')
	:param save_scaler (bool):		store sklearn.preprocessing object fitted to data as pickle object (default True)

	:r pca (obj):				sklearn.decomposition.PCA object fitted to data
	:r scaler (obj):			sklearn.preprocessing object with selected scaler fitted to data. If 'scaling' == 'None', then 'None' is returned
	"""

	A = data is None
	B = file_path is None

	if (A and B) or (not A and not B):
		print("Only one of the following parameters can be set: data or file_path")
		sys.exit()
		
	if not A:
		tpms = np.asarray(data)

	if not B:
		if not os.path.exists(file_path):	
			print("Please provide a valid file")
			sys.exit()

		print(now() + ": Loading " + file_path)
		tpms = np.loadtxt(file_path, delimiter='\t')

	if transpose:
		print(now() + ": Transposing")
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

	if scaling is not None:
		if scaling == 'standard':
			# (x_i - mean(X))/stdev(X)
			scaler = StandardScaler()
		elif scaling == 'robust':
			# (x_i - Q_1(X))/(Q_3(X) - Q_1(X))
			scaler = RobustScaler()
		elif scaling == 'normalizer':

			# x_i / sqrt(x_i, y_i, z_i, ...)
			scaler = Normalizer()
		else:
			print(scaling + " not recognized. Using standard scaler.")
			scaler = StandardScaler()

		print(now() + ": Scaling data")
		scaler.fit(tpms)
		tpms = scaler.transform(tpms)
	else:
		scaler = None

	print(now() + ": Performing PCA")
	pca = PCA(n_components=n_components, random_state = 42, svd_solver = 'full')
	pca.fit(tpms)
	
	if save_pca:
		print(now() + ": Saving PCA object")
		file_name = str(pca_name) + ".object"
		#open(file_name, 'wb')
		pickle.dump(pca, open(file_name, 'wb'), protocol=4)
		# file_name.close()
	
	if save_scaler and scaling is not None:
		print(now() + ": Saving scaler object")
		file_name = str(pca_name) + "." + scaling + "_scaler.object"
		pickle.dump(pca, open(file_name, 'wb'), protocol=4)

	print(now() + ": Finishing pca and scaling execution")
	return pca, scaler	

def now():
	return(str(datetime.datetime.now()))
