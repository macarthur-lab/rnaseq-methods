import functools
import numpy as np
from model_funcs import *

print = functools.partial(print, flush=True)

# filters
min_samples = .05

splices = {}

with open('splices.txt', 'r') as splice_list:
    num_samples = 0
    for splice_file in splice_list:
        splice_file = splice_file.rstrip('\n')
        splice_file = 'clean_splice_lists/' + splice_file

        print("About to process " + splice_file)
        current_splices, current_prop = detect_and_count(splice_file)
        print("Finished processing " + splice_file)
	
        joined_splices = np.concatenate((list(splices.keys()), current_splices), axis=None)
        joined_splices = np.unique(joined_splices)

        for key in joined_splices:
            idx = np.where(current_splices == key)[0]

            if len(idx) == 0:
                splices[key].append(0)
            else:
                if key not in splices:
                    splices[key] = [0] * num_samples

                idx = idx[0]
                splices[key].append(current_prop[idx])

        num_samples += 1
        print("Finished appending to dictionary")

# convert splices dictionary to numpy array
model = np.array(list(splices.values()))

# filter
to_keep = []
for i in range(0, len(model)):
    props = model[i]
    if len(props[props > 0]) / num_samples > min_samples:
        to_keep.append(i)

splice_sites = np.array(list(splices.keys()))
splice_sites = ' '.join(splice_sites[to_keep])

print("Saving finished model")
np.savetxt('NEB_model_no_annotation_filtered', model[to_keep].transpose(), delimiter=' ', header=splice_sites)
