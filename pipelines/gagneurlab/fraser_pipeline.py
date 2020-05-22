
# https://rpy2.github.io/doc/latest/html/introduction.html
import rpy2
from rpy2.robjects.packages import importr

#%%

# https://rpy2.github.io/doc/latest/html/notebooks.html
#import rpy2.ipython.html
#rpy2.ipython.html.init_printing()

#%%

from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

#%%

# import R's "base" package
base = importr('base')

# import R's "utils" package
utils = importr('utils')


fraser = importr('FRASER')

#%%
print("test")

#%%

import os
os.chdir(os.path.expanduser("~/project__rnaseq/data/samples/expression/bams"))


#%%

#[p for p in os.listdir(".") if p.endswith(".bam")]

fraser.FraserDataSet()
fraser.countRNAData
