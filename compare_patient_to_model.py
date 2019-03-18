from model import *

splice_file = '/home/jperezde/work/NEB_prop_model/153BR_JB_M1.splicing.txt'
model_file = '/home/jperezde/work/NEB_prop_model/NEB_model_no_annotation_filtered'

with open(model_file) as f:
	model_splices = f.readline()

model_splices = np.fromstring(model_splices, dtype='str', sep=' ')

model_prop = np.loadtxt(model_file, delimiter=' ', skiprows=1)
num_controls = len(model[0]) - 1

patient_splices, patient_prop = detect_and_count(splice_file)
joined_splices = np.concatenate((model_splices, patient_splices), axis=None)
joined_splices = np.unique(joined_splices)

len_comparison_model = len(joined_splices)
comparison_model = np.zeros((num_controls, len_comparison_model))

for i in range(0, len(joined_splices)):
	key = joined_splices[i]
	patient_idx = np.where(patient_splices == key)
	model_idx = np.where(model_splices == key)
	# splice site in patient but not in model
	if len(model_idx) == 0:
		comparison_model = np.append(comparison_model, zeros)
	else:
		# splice is not found in patient
		if len(model_idx) == 0:
			patient_splices = np.append(patient_splices, key)
			patient_prop = np.append(patient_prop, '0')
		#splice in both patient and model
		comparison_model = np.append(comparison_model, model_prop[model_idx])
	
comparison_model = np.insert(comparison_model, 0, joined_splices)

np.savetxt('NEB_comparison_model', comparison_model, delimiter=' ', fmt='%s')
np.savetxt('NEB_patient_model', np.append(patient_splices, patient_prop), delimiter=' ', fmt='%s')
