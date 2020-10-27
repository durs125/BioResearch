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



def uniformly_sample(signal, freq=0, number_of_samples = 0 ):
    if freq > 0:
        n = min(np.shape(signal))
        end = signal[np.shape(signal)[0]-1,0]
        number_of_samples = int( end*freq )
    elif freq < 0 or number_of_samples < 1:
        raise ValueError(("No samples specified or no sampling rate specified") ) 
    else:
        freq = signal[-1, 0]/number_of_samples

    uniform_sampling = np.zeros([number_of_samples, n])
    uniform_timestamps = np.linspace(0, end, number_of_samples)
    uniform_sampling[:, 0] = uniform_timestamps
    counter = 0
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
            return np.float32(1 + np.random.uniform())
        elif signal[index, 1] < signal[index - window_size, 1] \
                or signal[index, 1] < signal[index + window_size, 1]:
            return np.float32(1 + np.random.uniform())
        else:
            return np.float32(0)
    return window_checker


def compute_optimal_time_window(signal, splits_matrix_into=1):
    """ first half of the peak detection algorithm """
    n = max(np.shape(signal))  # number of samples
    rows = int(np.ceil(n / 2) - 1)  # length of lms matrix
    lms = np.float32(np.zeros((rows, n)))  # everything is a peak until proven guilty
    for x in range(0, rows):  # I wrote a function that creates a function to be called in another function for this...
        lms[x, :] = np.array(list(map(is_max_in_window(signal, n, x + 1), range(n))))
    row_sum = np.sum(lms, axis=1)
    #plt.plot(row_sum, range(1,np.shape(row_sum)+1) )#new, to debug inconsitent output
    #plt.show()
    gamma = np.where(row_sum == np.amin(row_sum))
    rescaled_lms = np.vsplit(lms, gamma[0] + 1)[0]  # cut the lms at the optimal search
    return rescaled_lms


def detect_peaks(signal):
#    print('detecting peaks')
    column_sd = np.std(compute_optimal_time_window(signal), axis=0) #column wise standard deviation
    # where columns sum to zero is where local maxima occur
    peaks_index = np.where(column_sd == 0)  # peaks occur when column-wise standard deviation is zero
    peaks = signal[peaks_index, :]  # select peaks based on their index
    peaks = peaks[0, :, :]  # I don't know why this is necessary but it is
    directory = "/home/david/BioResearch/python/PostProcessing/Simulations/beta0.2/plots/"
 
    pd.DataFrame(np.transpose(peaks)).to_csv((directory + 'peaks2000s.csv'), mode='a', header=False, index=False)
#pd.DataFrame(peaks).transpose().to_csv((directory + 'peaks2000s.csv'), mode='a', header=False, index=False)
    print("One down")
    return peaks  # pick off the peaks with timestamps and return them in a numpy array


def run_statistics(peaks):
   # print('generating statistics')
    mean_amplitude = np.mean(peaks[:, 1])
    mean_period = np.mean(np.diff(peaks[:, 0]))
    amplitude_coefficient_of_variation = np.std(peaks[:, 1]) / mean_amplitude
    period_coefficient_of_variation = np.std(np.diff(peaks[:, 0])) / mean_period
    return [mean_period, mean_amplitude, period_coefficient_of_variation, amplitude_coefficient_of_variation]

'''
def cbind2(az,bz):
    azshape =  np.shape(az)
    
    bzshape =  len(bz)  
    fullMon = np.zeros((max(max(azshape),bzshape),2))
    fullMon[0:azshape[0],0] = az
    fullMon[0:bzshape,1] = bz
    return(fullMon)
def cbind(az,bz):
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
    n = 2*int(.3333333333*freq1/freq2)+1 

    return np.repeat(1/n,n)




def triang_Coeff( freq1,freq2):# how much overlap?
    n = 2*int(.3333333333*freq1/freq2)+1

    first = np.linspace(1,int(.3333333333*freq1/freq2)+1,int(.3333333333*freq1/freq2)+1)
    last = np.linspace(int(.3333333333*freq1/freq2),1,int(.3333333333*freq1/freq2))
    vecs = np.append(first, last,axis = 0)
    sums = sum(first) +sum(last)
    return vecs/sums

def all_together_now(signal,  burn_in_time, freq1 = 0, number_of_samples1=0,freq2 = 0,  number_of_samples2=0):
    weights = makeBinomial(freq1,freq2)
    print('cleaning signal')



    #boxcar = scipy.signal.firwin(sizefil,.5, window = "boxcar") does not work
    if freq1 < freq2:
        spare = freq1
        freq1 = freq2
        freq2 = spare
    if number_of_samples1 < number_of_samples2:
        spare = number_of_samples1
        number_of_samples1 = number_of_samples2
        number_of_samples2 = spare
    uniformized = uniformly_sample(burn_in_time_series(signal, burn_in_time), freq= freq1, number_of_samples=number_of_samples1)
    clean_signal = burn_in_time_series(low_pass_filter(uniformized, makeBinomial(freq1,freq2),
    #clean_signal = uniformly_sample(scipy.signal.lfilter(boxcar,1, freq= freq1, number_of_samples=number_of_samples1)), freq = freq2, number_of_samples=number_of_samples2)
# I need rto make this accesible to other people 
#I need to make this work on more than one input column: do that with  looking at the maximum size of the array, filtering each, and merging the streams
# see if column bind is still around somewhere
    #merged = cbind(uniformized[,0], scipy.signal.lfilter(boxcar,1,uniformized[:,-1]))

    #clean_signal = uniformly_sample(merged, freq = freq2, number_of_samples=number_of_samples2)
#    return low_pass_filter( uniformized,makeBinomial(freq1,freq2))
    #return makeBinomial(freq1,freq2)
    return (clean_signal)
  #  return run_statistics(detect_peaks(clean_signal))



'''def clean_signal(signal,  freq,  number_of_samples=0,burn_in_time, weights):
    clean_signal = low_pass_filter(uniformly_sample(burn_in_time_series(signal,     burn_in_time),freq=freq number_of_samples = number_of_samples), weights)

    return clean_signal'''


def generate_heat_map(data, title, axis_labels, heat_center):
    heat_map = sb.heatmap(data, center=heat_center)
    heat_map.set_yticklabels(heat_map.get_yticklabels(), rotation=0)
    heat_map.set_xticklabels(heat_map.get_xticklabels(), rotation=30)
    heat_map.invert_yaxis()
    plt.title(title)
    plt.ylabel(axis_labels[1])
    plt.xlabel(axis_labels[0])
#   file_name = title + '.png'
#   heat_map.get_figure().savefig(file_name)
    plt.show()
    return heat_map






def plot_time_series_with_peaks(signal, burn_in_time, weights):
    burned_in_signal = burn_in_time_series(signal, burn_in_time)
    uniform_signal = uniformly_sample(burned_in_signal, round(max(burned_in_signal.shape)/10))
    clean_signal = low_pass_filter(uniform_signal, weights)
    peaks = detect_peaks(clean_signal)
    curve = plt.plot(clean_signal[0, :], clean_signal[1, :])
    maxima = plt.scatter(peaks[0, :], peaks[1, :], cmap='r')
    plt.show()
    return [curve, maxima]
