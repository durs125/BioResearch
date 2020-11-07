#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import math
import scipy
from scipy import signal
def burn_in_time_series(signal, burn_in_time):
    temp_signal = signal
    temp_signal[:, 0] = signal[:, 0] - burn_in_time
    new_start_time = np.where(temp_signal[:, 0] > 0)[0][0]
    new_start_state = np.zeros(temp_signal.shape[1])
    new_start_state[1:] = temp_signal[new_start_time - 1, 1:]
    temp_signal[new_start_time - 1, :] = new_start_state
    burned_in_signal = temp_signal[(new_start_time - 1):, :]
    return burned_in_signal



def uniformly_sample(signal, rate = 0, number_of_samples = 0):
    if rate > 0:
        n = min(np.shape(signal))
        end = signal[np.shape(signal)[0] - 1, 0]
        number_of_samples = int( end / rate )
    elif number_of_samples < 1:
        raise ValueError(("No samples specified or no sampling rate specified") ) 
    
    uniform_sampling = np.zeros([number_of_samples, n], dtype="float32")
    uniform_timestamps = np.linspace(0, end, number_of_samples)
    uniform_sampling[:, 0] = uniform_timestamps
    counter = int(0)
    for index in range(number_of_samples):
        while counter < max(np.shape(signal)):
            if signal[counter, 0] > uniform_timestamps[index]:
                uniform_sampling[index, 1:n] = signal[counter - 1, 1:n]
                break
            counter += 1
    return uniform_sampling


def low_pass_filter(signal, weights):
    if len(weights) % 2 == 1:
        chop = int((len(weights) - 1) / 2)
        filtered_signal = np.zeros([max(signal.shape) - chop * 2, 2])
        weights = np.array(weights)
        for index in range(max(filtered_signal.shape)):
            filtered_signal[index, 0] = signal[index + chop, 0]
            filtered_signal[index, 1] = np.dot(weights, signal[index:(index+chop*2+1), 1])
        return filtered_signal
    else:
        print("Don't make me do this.")
        return "error"
	#[1/256, 1/32, 7/64, 7/32, 35/128, 7/32, 7/64, 1/32, 1/256]

def is_max_in_window(signal, length_of_signal, window_size):
    def window_checker(index):
        if index <= window_size or index >= length_of_signal - window_size:
            return 1 + np.random.uniform()
        elif signal[index, 1] < signal[index - window_size, 1] \
                or signal[index, 1] < signal[index + window_size, 1]:
            return 1 + np.random.uniform()
        else:
            return 0
    return window_checker




def compute_optimal_time_window_WhileLoop(signal, splits_matrix_into = 1):
    '''Split the matrix to evaluate in memory saving chunks.
    Can usually split into 2048 pieces with no issue.'''
    if splits_matrix_into == 1 :
        splits_matrix_into = 2 + (np.shape(signal)[0] >> 11 )
    iterator = splits_matrix_into
    nsample = max(np.shape(signal))
    end_row = int(np.ceil(nsample / (2)) - 1)  
    start_row = int(np.ceil((iterator - 1) * nsample / (2 * splits_matrix_into)))    
    lms = np.zeros([end_row - start_row, nsample], dtype="float16")  #LMS = Local Maxima Scalogram
    #more efficient to initialize the first chunk outside the loop
    
    #make first chunk
    for index in range(start_row, end_row):
        lms[index - start_row, :] = np.array(list(map(is_max_in_window(signal, (end_row-start_row),
								       index + 1), range(nsample))))
    #end make first chunk

    row_sum = np.sum(lms, axis=1)
    gamma = np.where(row_sum == np.amin(row_sum))[0]
    previous_block = lms[0:(gamma[0] + 1), :]
    
    iterator -= 1
    while iterator >= 1:
        end_row = int(np.ceil(iterator * nsample / (2 * splits_matrix_into))) - 1
        start_row = int(np.ceil((iterator - 1) * nsample / (2 * splits_matrix_into)))
        size_combind = end_row - start_row + np.shape(previous_block)[0]
        lms = np.zeros([end_row - start_row + np.shape(previous_block)[0], nsample], dtype="float16")
        
	#make subsequent block
	for index in range(end_row - start_row):
	    lms[index, :] = np.array(list(map(is_max_in_window(signal, (end_row-start_row), 
							       index + 1), range(nsample))))
        #end make subsequent block
                
        lms[range(end_row - start_row, size_combind), :] = previous_block  #attach previous block
        row_sum = np.sum(lms, axis=1)
        gamma = np.where(row_sum == np.amin(row_sum))[0]
        previous_block = lms[0:(gamma[0] + 1), :]  # cut the lms at the optimal search window
        iterator -= 1
    return previous_block ##return value of the row sum to avoid duplication




def compute_optimal_time_window(signal, splits_matrix_into=8):
    """ first half of the peak detection algorithm """
    nsamples = max(np.shape(signal))  
    rows = int(np.ceil(nsamples / 2) - 1)  
    lms = np.zeros((rows, nsamples), dtype="float16")  #LMS = Local Maxima Scalogram

    for x in range(0, rows):
        lms[x, :] = np.array(list(map(is_max_in_window(signal, n, x + 1), range(n))))
    row_sum = np.sum(lms, axis=1)
    gamma = np.where(row_sum == np.amin(row_sum))
    rescaled_lms = np.vsplit(lms, gamma[0] + 1)[0]  # cut the lms at the optimal search window
    return rescaled_lms


def detect_peaks(signal):
    column_sd = np.std(compute_optimal_time_window(signal), axis=0)  # where columns sum to zero is where local maxima occur
    peaks_index = np.where(column_sd == 0)  
    peaks = signal[peaks_index, :]
    peaks = peaks[0, :, :]  # I don't know why this is necessary but it is
    return peaks


def run_statistics(peaks):
    mean_amplitude = np.mean(peaks[:, 1])
    mean_period = np.mean(np.diff(peaks[:, 0]))
    amplitude_coefficient_of_variation = np.std(peaks[:, 1]) / mean_amplitude
    period_coefficient_of_variation = np.std(np.diff(peaks[:, 0])) / mean_period
    return [mean_period, mean_amplitude, period_coefficient_of_variation, amplitude_coefficient_of_variation]

'''
def cbind2(az,bz): #combinds two columns
    azshape =  np.shape(az)
    
    bzshape =  len(bz)  
    fullMon = np.zeros((max(max(azshape),bzshape),2))
    fullMon[0:azshape[0],0] = az
    fullMon[0:bzshape,1] = bz
    return(fullMon)
def cbind(az,bz):#adds a column to an array
    azshape =  np.shape(az)
    bzshape =  length(bz)  
    fulCol = np.zero((max(azshape[0],bzshape),(azshape[1]+bzshape[1])))
    fulCol[0:azshape[0],0:azshape[1]] #= azshape
    fulCol[bzshape[0]:,bzshape[1]:] = bzshape
    return(fulCol)'''


def binomialCoeffSum( n): 
        
        C = [[0]*(n+2) for i in range(0,n+2)] 
        
       
        # Calculate value of Binomial  
        # Coefficient in bottom up manner 
        for i in range(0,n+1): 
            for j in range(0, min(i, n)+1): 
              
                # Base Cases 
                if (j == 0 or j == i): 
                    C[i][j] = 1
       
                # Calculate value using previously 
                # stored values 
                else: 
                    C[i][j] = C[i - 1][j - 1] + C[i - 1][j] 
       
        # Calculating the sum. 
        sums = 0
        for i in range(0,n+1): 
            sums += C[n][i] 
       
        return sums

def makeBinomial(freq1 = 1, freq2 = 1):
    points_ahead = int(freq1/freq2)
    binomVector = np.zeros(2*points_ahead+1)
    bsum = binomialCoeffSum(2*points_ahead)
    for ixz in range(2*points_ahead+1):
        binomVector[ixz] =scipy.special.binom(2*points_ahead,ixz)/bsum
    return(binomVector)
#makeBinomial(2)


def uniform_Coeff( freq1,freq2):# how much overlap?
    n = 18*int(freq1/freq2)+1 

    return np.repeat(1/n,n)




def triang_Coeff( freq1,freq2):# how much overlap?
    n = 4*int(freq1/freq2)+1

    first = np.linspace(1,2*int(freq1/freq2)+1,int(freq1/freq2)+1)
    last = np.linspace(int(2*freq1/freq2),1,int(2*freq1/freq2))
    vecs = np.append(first, last,axis = 0)
    sums = sum(first) +sum(last)
    return vecs/sums

def all_together_now_Binom(signal,  burn_in_time, freq1 = 0, number_of_samples1=0, freq2 = 1,  number_of_samples2 =0):
    weights = makeBinomial(freq1,freq2)
    #print('cleaning signal')
    signal = 1000*np.random.rand(9,3)
    bb = np.where(signal[:,0]<=600)
    signal = signal[bb,:]
    print("Abriged Signal")
    if freq1 < freq2:
        spare = freq1
        freq1 = freq2
        freq2 = spare
    if number_of_samples1 < number_of_samples2:
        spare = number_of_samples1
        number_of_samples1 = number_of_samples2
        number_of_samples2 = spare

    uniformized = uniformly_sample(burn_in_time_series(signal, burn_in_time), freq= freq1, number_of_samples=number_of_samples1)


    clean_signal = low_pass_filter(uniformized, weights)
    
    peaks = detect_peaks(clean_signal)

   
    return peaks

def all_together_now_Triang(signal,  burn_in_time, freq1 = 0, number_of_samples1=0,freq2 = 1,  number_of_samples2=0):
    weights = triang_Coeff(freq1,freq2)
    #print('cleaning signal')
    signal = 1000*np.random.rand(9,3)
    bb = np.where(signal[:,0]<=600)
    signal = signal[bb,:]

    if freq1 < freq2:
        spare = freq1
        freq1 = freq2
        freq2 = spare
    if number_of_samples1 < number_of_samples2:
        spare = number_of_samples1
        number_of_samples1 = number_of_samples2
        number_of_samples2 = spare

    uniformized = uniformly_sample( burn_in_time_series(signal, burn_in_time), freq= freq1, number_of_samples=number_of_samples1)


    clean_signal =  low_pass_filter( uniformized, weights)   
    
    peaks = detect_peaks(clean_signal)

   
    return peaks
def all_together_nowUni(signal,  burn_in_time, freq1 = 0, number_of_samples1=0,freq2 = 1,  number_of_samples2=0): # short Run
    weights = uniform_Coeff(freq1,freq2)
    #print('cleaning signal')
    signal = 1000*np.random.rand(9,3)
    bb = np.where(signal[:,0]<=600) # some of the runs are too long. We only want 600 minute
    signal = signal[bb,:]



    #boxcar = scipy.signal.firwin(sizefil,.5, window = "boxcar") does not work
    if freq1 < freq2:
        spare = freq1
        freq1 = freq2
        freq2 = spare
    if number_of_samples1 < number_of_samples2:
        spare = number_of_samples1
        number_of_samples1 = number_of_samples2
        number_of_samples2 = spare
    #print(signal)
    uniformized = uniformly_sample(burn_in_time_series(signal, burn_in_time), freq= freq1, number_of_samples=number_of_samples1)
    #print("signal cleaned")
    clean_signal =  low_pass_filter( uniformized, weights)   


    #print(clean_signal)
    peaks = detect_peaks(clean_signal)
    #print("done")
    #clean_signal = uniformly_sample(merged, freq = freq2, number_of_samples=number_of_samples2)
#    return low_pass_filter( uniformized,makeBinomial(freq1,freq2))
    #return makeBinomial(freq1,freq2)
   # return (clean_signal)
    return peaks

def all_together_now_Once(signal,  burn_in_time, freq1 = 0, number_of_samples1 = 0,freq2 = 1,  number_of_samples2 = 0):
    #weights = makeBinomial(freq1,freq2)
    #print('cleaning signal')

   
    #print('cleaning signal')
    signal = 1000*np.random.rand(9,3)
    bb = np.where(signal[:,0]<=600) # some of the runs are too long. We only want 600 minute
    signal = signal[bb,:]


    #boxcar = scipy.signal.firwin(sizefil,.5, window = "boxcar") does not work
    if freq1 < freq2:
        spare = freq1
        freq1 = freq2
        freq2 = spare
    if number_of_samples1 < number_of_samples2:
        spare = number_of_samples1
        number_of_samples1 = number_of_samples2
        number_of_samples2 = spare
   # #print(signal)
    uniformized = uniformly_sample(burn_in_time_series(signal, burn_in_time), freq= freq1, number_of_samples=number_of_samples1)


    peaks = detect_peaks(uniformized)

    return peaks

def generate_heat_map(data, title, axis_labels, heat_center):
    f = plt.figure()
    heat_map = sb.heatmap(data, center=heat_center)
    heat_map.set_yticklabels(heat_map.get_yticklabels(), rotation=0)
    heat_map.set_xticklabels(heat_map.get_xticklabels(), rotation=30)
    heat_map.invert_yaxis()
    plt.title(title)
    plt.ylabel(axis_labels[1])
    plt.xlabel(axis_labels[0])
#   file_name = title + '.png'
#   heat_map.get_figure().savefig(file_name)
    #plt.show()
    f.savefig(str(title) +"_heatmap.pdf", bbox_inches='tight')
    return heat_map



def cleanStatsHeatMap(directory2,file_names, mean_axis,cv_axis,burn_in_time, freq1, weights = 0):
    if weights ==0;
        stats = Fun.all_together_now_Once(signal=np.genfromtxt(file_names[mean_axis, cv_axis], delimiter=','),burn_in_time = burn_in_time,  freq = freq1, )
    if weights ==1;
        stats = Fun.all_together_nowUni(signal=np.genfromtxt(file_names[mean_axis, cv_axis], delimiter=','),burn_in_time = burn_in_time,  freq = freq1, )
    if weights ==2;
        stats = Fun.all_together_now_Triang(signal=np.genfromtxt(file_names[mean_axis, cv_axis], delimiter=','),burn_in_time = burn_in_time, freq = freq1, )
    if weights ==3;
        stats = Fun.all_together_now_Binom(signal=np.genfromtxt(file_names[mean_axis, cv_axis], delimiter=','),burn_in_time = burn_in_time,  freq = freq1, )





    heat_map_matrices[:, mean_axis, cv_axis] = stats

    pass



def plot_time_series_with_peaks(signal, burn_in_time, weights):
    burned_in_signal = burn_in_time_series(signal, burn_in_time)
    uniform_signal = uniformly_sample(burned_in_signal, round(max(burned_in_signal.shape)/10))
    clean_signal = low_pass_filter(uniform_signal, weights)
    peaks = detect_peaks(clean_signal)
    curve = plt.plot(clean_signal[0, :], clean_signal[1, :])
    maxima = plt.scatter(peaks[0, :], peaks[1, :], cmap='r')
    plt.show()
    return [curve, maxima]
