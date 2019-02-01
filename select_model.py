import sys
import functools
from GTEx_pca import *
from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import StandardScaler, RobustScaler, Normalizer
# from sklearn.decomposition import PCA
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

print = functools.partial(print, flush=True)


# Naive bayes
# https://scikit-learn.org/stable/modules/generated/sklearn.naive_bayes.GaussianNB.html#sklearn.naive_bayes.GaussianNB
# K-Nearest neighbors
# https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html#sklearn.neighbors.KNeighborsClassifier
# Random forest
# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html#sklearn.ensemble.RandomForestClassifier
# Support vector machine classification
# https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html#sklearn.svm.SVC

# Load data
print(now() + ": Loading data")
X = np.loadtxt('data/genes_tpm_values.gct', delimiter='\t')

# Perform log2 transformation and transpose
print(now() + ": Transposing and performing log2 transformation")
X = X.transpose()
X = np.log2(X + 1)

# Train test split
print(now() + ": Splitting into training and test sets")
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = .20, random_state = 42, stratify = y)

# Get pca and scaler fitted to training set
pca, scaler = run_pca(data = X_train, transpose = False, log_transform = False, scaling = "normalizer", n_components = 500, save_pca = True, save_scaler = True, pca_name = "training_pca")

print(now() + ": Normalizing sets")
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

print(now() + ": Transforming into principal components")
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)

gnb_acc = []
knn_acc = []
rf_acc = []
svc_rbf_acc = []
svc_sig_acc = []

for i in range(5, 501, 5):
	print(now() + ": Obtaining accuracies with " + str(i) + " principal components")
	X_train_slice = X_train[:,:i]
	X_test_slice = X_test[:,:i]

	# Naive Bayes
	print(now() + ":    Naive bayes")
	gnb = GaussianNB()
	gnb.fit(X_train_slice, y_train)
	y_pred = gnb.predict(X_test_slice)
	gnb_acc.append(accuracy_score(y_test, y_pred))

	# K-Nearest Neighbors
	print(now() + ":    K-Nearest neighbors")
	knn = KNeighborsClassifier(n_neighbors=6)
	knn.fit(X_train_slice, y_train)
	y_pred = knn.predict(X_test_slice)
	knn_acc.append(accuracy_score(y_test, y_pred))

	# Random Forest
	print(now() + ":    Random forest")
	rf = RandomForestClassifier(max_depth=10, random_state=42, bootstrap=100)
	rf.fit(X_train_slice, y_train)
	y_pred = rf.predict(X_test_slice)
	rf_acc.append(accuracy_score(y_test, y_pred))

	# SVM classification with radial basis function kernel
	print(now() + ":    SVMC with rbf kernel")
	svc = SVC(kernel = 'rbf', random_state = 42)
	svc.fit(X_train_slice, y_train)
	y_pred = svc.predict(X_test_slice)
	svc_rbf_acc.append(accuracy_score(y_test, y_pred))
	
	# SVM classification with sigmoid kernel
	print(now() + ":    SVMC with sigmoid kernel")
	svc = SVC(kernel = 'sigmoid', random_state = 42)
	svc.fit(X_train_slice, y_train)
	y_pred = svc.predict(X_test_slice)
	svc_sig_acc.append(accuracy_score(y_test, y_pred))
	
print(now() + ": Saving accuracy scores")
np.savetxt('gnb_acc',np.asarray(gnb_acc))
np.savetxt('knn_acc',np.asarray(knn_acc))
np.savetxt('rf_acc',np.asarray(rf_acc))
np.savetxt('svc_rbf_acc',np.asarray(svc_rbf_acc))
np.savetxt('svc_sig_acc',np.asarray(svc_sig_acc))



