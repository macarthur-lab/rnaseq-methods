import argparse
from model_funcs import *

splice_file = '/home/jperezde/work/NEB_prop_model/26I_SK_M1.splicing.txt'
model_file = '/home/jperezde/work/NEB_prop_model/NEB_model_no_annotation_filtered'

def main(args):
    with open(args.model) as f:
        model_splices = f.readline().rstrip('\n').split(' ')

    model_splices = np.array(model_splices) 

    model_prop = np.loadtxt(args.model, delimiter=' ', skiprows=1)
    num_controls = len(model_prop[:,0])

    patient_splices, patient_prop = detect_and_count(args.patient_splices)
    joined_splices = np.concatenate((model_splices, patient_splices), axis=None)
    joined_splices = np.unique(joined_splices)

    len_comparison_model = len(joined_splices)
    comparison_model = np.zeros((num_controls, len_comparison_model))

    for i in range(0, len(joined_splices)):
        key = joined_splices[i]
        patient_idx = np.where(patient_splices == key)[0]
        model_idx = np.where(model_splices == key)[0]

        if len(patient_idx) == 0:
            # splice is not found in patient
            comparison_model[:,i] = model_prop[:,model_idx[0]]
            patient_splices = np.append(patient_splices, key)
            patient_prop = np.append(patient_prop, 0)

        elif len(model_idx) != 0:
        # splice site in patient and model
            comparison_model[:,i] = model_prop[:,model_idx[0]]

    patient_splices = ' '.join(patient_splices)
    joined_splices = ' '.join(joined_splices)

    np.savetxt('NEB_26I_SK_M1_comparison_model', comparison_model, delimiter=' ', header=patient_splices)
    np.savetxt('NEB_26I_SK_M1_patient_model', patient_prop.reshape(1, len(patient_prop)), delimiter=' ', header=patient_splices)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Compare patient splicing against model")
    parser.add_argument('--model', help = "File path of gene model", required = True)
    parser.add_argument('--patient_splices', help = "File with patient splices", required = True)
    parser.add_arguemnt('--base_name', help = "Base name of files to be saved", default = '')
    args=parser.parse_args()
    main(args)
