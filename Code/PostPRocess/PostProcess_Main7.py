#!/usr/bin/env python
#         gamma_r40.0
import multiprocessing as mp
safeProcessors = max(1, mp.cpu_count() - 2)
pool2 = mp.Pool(safeProcessors)
import PostProcessing_Functions7Short as Fun
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os


#directory1 = ["/scratch/Infor/__pycache__/PostProcessing/Simulations/Co20.0/","/scratch/Infor/__pycache__/PostProcessing/Simulations/Co5.0/", "/scratch/Infor/__pycache__/PostProcessing/Simulations/Ro0.5/","/scratch/Infor/__pycache__/PostProcessing/Simulations/Ro2.0/"]
directory1 = ["/scratch/Infor/PostProcessing/Simulations/Ro0.5/","/scratch/Infor/PostProcessing/Simulations/Ro2.0/"]
#directory2 = "/scratch/Infor/__pycache__/"
#freqs = np.linespace(1,3,10)
for directory in directory1:

    file_names = np.array(pd.read_csv(directory+'1metadata.csv', header=None))
    heat_map_axes = [16]
    heat_map_axes0 = [1]
    heat_map_matrices = np.zeros([4, file_names.shape[0], file_names.shape[1]])


    directory2 = "scratch/Infor/"
    file_names = np.array(pd.read_csv(directory+'1metadata.csv', header=None))
     
    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap, [(np.genfromtxt(directory3, delimiter=','), mean_axis , cv_axis ,  100, 10, 0) for cv_axis in range(file_names.shape[1])])
    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'Once_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'Once_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'Once_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'Once_amplitude_CV.csv', mode = "a")

    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap,[(np.genfromtxt(directory3, delimiter=','), mean_axis , cv_axis , 100,10, 1) for cv_axis in range(file_names.shape[1])])

    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'UniPass_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'UniPass_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'UniPass_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'UniPass_amplitude_CV.csv', mode = "a")

    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap, [(np.genfromtxt(directory3, delimiter=','), mean_axis , cv_axis , 100,10, 2) for cv_axis in range(file_names.shape[1])])


    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'Triang_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'Triang_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'Triang_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'Triang_amplitude_CV.csv', mode = "a")
    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap, [(np.genfromtxt(directory3, delimiter=','), mean_axis , cv_axis , 100,10, 3) for cv_axis in range(file_names.shape[1])])

    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'Binom_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'Binom_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'Binom_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'Binom_amplitude_CV.csv', mode = "a")
pool2.close()
pool2.join()
'''mean_period_heat_map = Fun.generate_heat_map(mean_period, 'Mean Period', ['CV of Delay', 'Delay Mean (min)'],
                                             np.mean(list(heat_map_matrices[0, :, :])))
mean_amplitude_heat_map = Fun.generate_heat_map(mean_amplitude, 'Mean Amplitude', ['CV of Delay', 'Delay Mean (min)'],
                                                np.mean(list(heat_map_matrices[1, :, :])))
period_CV_heat_map = Fun.generate_heat_map(period_CV, 'CV of Period', ['CV of Delay', 'Delay Mean (min)'],
                                           np.mean(list(heat_map_matrices[2, :, :])))
amplitude_CV_heat_map = Fun.generate_heat_map(amplitude_CV, 'CV of Amplitude', ['CV of Delay', 'Delay Mean (min)'],
                                              np.mean(list(heat_map_matrices[3, :, :])))

normalized_heat_map_matricies = np.zeros(heat_map_matrices.shape)'''

'''
directory = "/scratch/Infor/__pycache__/PostProcessing/Simulations/gamma_r160.0/"
directory2 = "/scratch/Infor/__pycache__/"
file_names = np.array(pd.read_csv(directory+'1metadata.csv', header=None))
heat_map_axes = [16]
heat_map_axes0 = [1]
heat_map_matrices = np.zeros([4, file_names.shape[0], file_names.shape[1]])


def cleanStatsHeatMap(directory2,file_names, mean_axis,cv_axis, heat_map_matrices):
    with open(directory2 + file_names[mean_axis, cv_axis]) as file:
        length = len(list(csv.reader(file)))
    t1 = time.time()
    stats = Fun.all_together_now(np.genfromtxt(file_names[mean_axis, cv_axis], delimiter=','),
                                 int(length*.02), 100, [1/256, 1/32, 7/64, 7/32, 35/128, 7/32, 7/64, 1/32, 1/256])
    heat_map_matrices[:, mean_axis, cv_axis] = stats

    pass

for mean_axis in range(file_names.shape[0]):
    pool2.starmap(cleanStatsHeatMap, [(directory2,file_names, mean_axis,cv_axis, heat_map_matrices) for cv_axis in range(file_names.shape[1])])
        

        

pool2.close()
pool2.join()
mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

mean_period.to_csv(directory + 'mean_period.csv', mode = "a")
mean_amplitude.to_csv(directory + 'mean_amplitude.csv', mode = "a")
period_CV.to_csv(directory + 'period_CV.csv', mode = "a")
amplitude_CV.to_csv(directory + 'amplitude_CV.csv', mode = "a")

mean_period_heat_map = Fun.generate_heat_map(mean_period, 'Mean Period', ['CV of Delay', 'Delay Mean (min)'],
                                             np.mean(list(heat_map_matrices[0, :, :])))
mean_amplitude_heat_map = Fun.generate_heat_map(mean_amplitude, 'Mean Amplitude', ['CV of Delay', 'Delay Mean (min)'],
                                                np.mean(list(heat_map_matrices[1, :, :])))
period_CV_heat_map = Fun.generate_heat_map(period_CV, 'CV of Period', ['CV of Delay', 'Delay Mean (min)'],
                                           np.mean(list(heat_map_matrices[2, :, :])))
amplitude_CV_heat_map = Fun.generate_heat_map(amplitude_CV, 'CV of Amplitude', ['CV of Delay', 'Delay Mean (min)'],
                                              np.mean(list(heat_map_matrices[3, :, :])))

normalized_heat_map_matricies = np.zeros(heat_map_matrices.shape)

import multiprocessing as mp
safeProcessors = max(1, mp.cpu_count() - 2)
pool2 = mp.Pool(safeProcessors)
import PostProcessing_Functions5o as Fun
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import time
import os

directory = "/scratch/Infor/__pycache__/PostProcessing/Simulations/gamma_r40.0/"
directory2 = "/scratch/Infor/__pycache__/"
file_names = np.array(pd.read_csv(directory+'1metadata.csv', header=None))
heat_map_axes = [16]
heat_map_axes0 = [1]
heat_map_matrices = np.zeros([4, file_names.shape[0], file_names.shape[1]])


def cleanStatsHeatMap(directory2,file_names, mean_axis,cv_axis):
    with open(directory2 + file_names[mean_axis, cv_axis]) as file:
        length = len(list(csv.reader(file)))
    t1 = time.time()
    stats = Fun.all_together_now(np.genfromtxt(file_names[mean_axis, cv_axis], delimiter=','),
                                 int(length*.02), 100, [1/256, 1/32, 7/64, 7/32, 35/128, 7/32, 7/64, 1/32, 1/256])
    heat_map_matrices[:, mean_axis, cv_axis] = stats

    pass

for mean_axis in range(file_names.shape[0]):
    pool2.starmap(cleanStatsHeatMap, [(directory2,file_names, mean_axis,cv_axis) for cv_axis in range(file_names.shape[1])])
        

        

pool2.close()
pool2.join()
mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

mean_period.to_csv(directory + 'mean_period.csv')
mean_amplitude.to_csv(directory + 'mean_amplitude.csv')
period_CV.to_csv(directory + 'period_CV.csv')
amplitude_CV.to_csv(directory + 'amplitude_CV.csv')

mean_period_heat_map = Fun.generate_heat_map(mean_period, 'Mean Period', ['CV of Delay', 'Delay Mean (min)'],
                                             np.mean(list(heat_map_matrices[0, :, :])))
mean_amplitude_heat_map = Fun.generate_heat_map(mean_amplitude, 'Mean Amplitude', ['CV of Delay', 'Delay Mean (min)'],
                                                np.mean(list(heat_map_matrices[1, :, :])))
period_CV_heat_map = Fun.generate_heat_map(period_CV, 'CV of Period', ['CV of Delay', 'Delay Mean (min)'],
                                           np.mean(list(heat_map_matrices[2, :, :])))
amplitude_CV_heat_map = Fun.generate_heat_map(amplitude_CV, 'CV of Amplitude', ['CV of Delay', 'Delay Mean (min)'],
                                              np.mean(list(heat_map_matrices[3, :, :])))

normalized_heat_map_matricies = np.zeros(heat_map_matrices.shape)
'''
#
