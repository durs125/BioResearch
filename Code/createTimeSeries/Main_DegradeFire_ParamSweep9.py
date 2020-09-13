#!/usr/bin/env python
#yr

import multiprocessing as mp
safeProcessors = max(1, mp.cpu_count() - 1)
pool2 = mp.Pool(safeProcessors)
# Degrade and fire
from  Functions_DegradeFire4 import *
from  Classes_DegradeFire4 import *

import pandas as pd
import os
import math
import time
from pathlib import Path
import numpy as np

#from numpy import random# this is not used

#from  Functions_DegradeFire4 import gillespie_sim #this is not necessary, already imported
mean_range = np.linspace(5, 10, 16)
cv_range = np.linspace(0, .5, 16)
alpha = 300
beta = .1
R0 = 1
C0 = 10
yr =80
par_range = np.linspace(150, 600, 2)  # alpha
param = 'alpha' #THE parameter


for par in par_range:#This loops over THE parameter in question
    alpha = par #This shoud be the parameter
    path1 = 'PostProcessing/Simulations/{}{}'.format(param,par)
    #path2 = 'Simulations/'+ param + str(par)
    Path(path1).mkdir(parents=True, exist_ok=True)
    pd.DataFrame([mean_range]).to_csv(path1 + '/0metadata.csv', header=False, index=False)
    pd.DataFrame([cv_range]).to_csv(path1 + '/0metadata.csv', mode='a', header=False, index=False)
    row_of_file_names = []
    for mu in mean_range:
        for cv in cv_range:
            file_name = path1+ '/mean=' + str(mu) + '_CV=' + str(cv)  + '.csv'
            row_of_file_names.append(file_name)
        paths = path1 + '/1metadata.csv'
        pd.DataFrame([row_of_file_names]).to_csv(paths,  mode='a', header=False, index=False)
        row_of_file_names = []
    dilution = Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])
    enzymatic_degradation = Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, R0, 1], 1, [0])
    
    cv =.5
    for mu in mean_range:
        #for cv in cv_range:

        pool2.starmap(gillespie_sim, [(mu, cv,alpha,beta,R0 ,C0,yr,param,par,dilution,enzymatic_degradation) for cv in cv_range])

pool2.close()
pool2.join()

