import sys
import functools
from GTEx_pca import *
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import accuracy_score

print = functools.partial(print, flush=True)

def get_accuracy(classifier, X_train, y_train, X_test, y_test):
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    return accuracy_score(y_test, y_pred)

classifiers_to_use = {
    'gnb': {'classifier': GaussianNB(), 'acc': []},
    'knn': {'classifier': KNeighborsClassifier(n_neighbors=7), 'acc': []},
    'lr': {'classifier': LogisticRegression(random_state = 42, solver = 'lbfgs', multi_class='ovr', n_jobs = -1), 'acc': []},
    'lin_svc': {'classifier': LinearSVC(multi_class='ovr', random_state = 42, max_iter = 2500), 'acc': []}
}


# Load data
print(now() + ": Loading data")
X = np.loadtxt('data/genes_tpm_values.gct', delimiter='\t')
y = np.loadtxt('data/tissue_types_matched.tsv', delimiter='\t', dtype='str')
y = y[:,1]
y = np.char.strip(y)

# Perform log2 transformation and transpose
print(now() + ": Transposing and performing log2 transformation")
X = X.transpose()
X = np.log2(X + 1)

# Train test split
print(now() + ": Splitting into training and test sets")
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = .20, random_state = 42, stratify = y)

# Get pca and scaler fitted to training set
pca, scaler = run_pca(data = X_train, transpose = False, log_transform = False, scaling = "normalizer", n_components = 40, save_pca = True, save_scaler = True, pca_name = "training_pca_2")

print(now() + ": Normalizing sets")
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

print(now() + ": Transforming into principal components")
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)

pcs = [i for i in range(5,41) if (i % 5 == 0 or i % 5 == 3)]

for i in pcs:
    print(now() + ": Obtaining accuracies with " + str(i) + " principal components")
    X_train_slice = X_train[:,:i]
    X_test_slice = X_test[:,:i]

    for method in classifiers_to_use:
        print(now() + ": " + method)
        classifiers_to_use[method]['acc'].append(get_accuracy(classifiers_to_use[method]['classifier'], X_train_slice, y_train, X_test_slice, y_test))
	
print(now() + ": Saving accuracy scores")

for method in classifiers_to_use:
    np.savetxt(method, np.asarray(classifiers_to_use[method]['acc']))



