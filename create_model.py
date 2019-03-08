import functools
import numpy as np

print = functools.partial(print, flush=True)

splices = {}

with open('splices.txt', 'r') as splice_list:
    num_samples = 0
    for splice_file in splice_list:
        splice_file = splice_file.rstrip('\n')
        splice_file = 'clean_splice_lists/' + splice_file
        current_splices = np.asarray([], dtype='str')
        current_counts = np.asarray([])

        print("About to process " + splice_file)
        with open(splice_file, 'r') as f:
            for splice_event in f: 
                splice_event = splice_event.rstrip('\n').split()
                start = splice_event[3]
                end = splice_event[4]
                key = '2:' + start + '-' + end

                if key not in current_splices:
                    current_splices = np.append(current_splices, key)
                    current_counts = np.append(current_counts, 0)

                current_counts[np.where(current_splices == key)[0][0]] += 1

        print("Finished processing " + splice_file)
        current_total = np.sum(current_counts)
        current_prop = np.divide(current_counts, current_total)
	
        joined_splices = np.concatenate((list(splices.keys()), current_splices), axis=None)
        joined_splices = np.unique(joined_splices)

        for key in joined_splices:
            idx = np.where(current_splices == key)[0]

            if len(idx) == 0:
                splices[key].append('0')
            else:
                if key not in splices:
                    splices[key] = ['0'] * num_samples

                idx = idx[0]
                splices[key].append(str(current_prop[idx]))

        num_samples += 1
        print("Finished appending to dictionary")

with open('NEB_model_no_annotation', 'a') as model:
     for key in splices:
        mystr = ' '.join(splices[key])
        model.write(key + ' ' + mystr + '\n')
