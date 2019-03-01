import numpy as np
import pickle

PCA="pca.object"
SCALER="scaler.object"
CLASSIFIER="classifier.object"

# Load PCA and prediction model
pca_file = open(PCA, "rb")
pca = pickle.load(pca_file)

scaler_file = open(SCALING, "rb")
scaler = pickle.load(scaler_file)

classifier_file = open(CLASSIFIER, "rb")
classifier = pickle.load(classifier_file)

# Run RNA-seqc
# subprocess.run(["rnaseqc", REF, BAM, "--coverage", "."])

# Load tpms
tpms = np.loadtxt("gene_tpm.gct"))

# Log2 and PCA transformation
tpms = np.log2(tpms+1)
tpms = tpms.reshape(1,-1)
tpms = scaler.transform(tpms)
pcs = pca.transform(tpms)

# Predict
pred_tissue = classifier.predict(pcs[:,55])

# Return output
print(pred_tissue)

