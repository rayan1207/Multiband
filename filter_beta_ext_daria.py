# -*- coding: utf-8 -*-
"""
Created on Sun May 14 22:41:15 2023

@author: Rayan
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
axis_font = {'fontname':'Arial', 'size':'15'}
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)



def summ(dict_of_results: dict, dict_of_w: dict) -> dict:
    result_sum = defaultdict(lambda: np.array([0 + 1j*0]))
    for key in dict_of_results:
        cur_result = dict_of_results[key]
        cur_w = dict_of_w[key]
        temp = cur_w[0]
        for i in range(len(cur_result)):
            if abs(cur_w[i] - temp) <= 0.0001:
                result_sum[key][-1] += cur_result[i]
            else:
                result_sum[key] = np.append(result_sum[key], cur_result[i])
                temp = cur_w[i]
    return result_sum



def filter(filename,extline,beta):
    dict_of_results = defaultdict(lambda: defaultdict(lambda: np.array([])))
    dict_of_w = defaultdict(lambda: defaultdict(lambda: np.array([])))
    data = np.loadtxt(filename)
    for k in range(len(data)):
        string = data[k]
        dict_of_results[string[0]][string[4], string[5]] = np.append(dict_of_results[string[0]][string[4], string[5]], string[2] + 1j*string[3] )
        dict_of_w[string[0]][string[4], string[5]] = np.append(dict_of_w[string[0]][string[4], string[5]], string[1])
    
  
    
    w = defaultdict(lambda: np.array([]))
    w[data[0, 0]] = np.array([data[0, 1]])
    temp_beta = data[0, 0]
    for i in range(len(data[:, 1])):
        if abs(data[i, 0] - temp_beta) < 0.001:
            if abs(data[i, 1] - w[temp_beta][-1]) > 0.0001:
                w[temp_beta] = np.append(w[temp_beta], data[i, 1])
        else:
            w[data[i, 0]] = np.array([data[i, 1]])
            temp_beta = data[i, 0]
    
    print(w)
    
    result_sum = dict()
    for key, value in dict_of_results.items():
        result_sum[key] = summ(value, dict_of_w[key])
        print( dict_of_w[key])
    return [w[beta],result_sum[beta][extline[0],extline[1]]]

ext=[1,1]
beta = 50
graph0 = filter('sto-10_02_g0_v0.txt', ext,beta)
graph1 = filter('sto-10_02_g0_v1.txt', ext,beta)


Jdata = np.loadtxt('sigma11_james.dat')


result_sigma_real = (graph0[1]+graph1[1]).real
result_sigma_imag = (graph0[1]+graph1[1]).imag
print(graph0[1])
print(graph1[1])

freq = graph0[0]


plt.subplot(1, 2, 1)
plt.plot(freq, result_sigma_real,'.-',label='G_AMI')
plt.plot(Jdata[:, 2], Jdata[:, 3],'x-',label='Sym-det')
plt.xlabel('Mfreq',**axis_font)
plt.ylabel('Real Sigma11',**axis_font)
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(freq, result_sigma_imag,'.-',label='G_AMI')
plt.plot(Jdata[:, 2], Jdata[:, 4],'x-',label='Sym-det')
plt.ylabel('Imag sigma11',**axis_font)
plt.xlabel('Mfreq',**axis_font)
plt.tight_layout()
plt.legend()
plt.show()

ext=[2,2]
beta = 50
graph0 = filter('sto-10_02_g0_v0.txt', ext,beta)
graph1 = filter('sto-10_02_g0_v1.txt', ext,beta)
print(graph0[1])
print(graph1[1])

Jdata = np.loadtxt('sigma22_james.dat')


result_sigma_real = (graph0[1]+graph1[1]).real
result_sigma_imag = (graph0[1]+graph1[1]).imag


freq = graph0[0]


plt.subplot(1, 2, 1)
plt.plot(freq, result_sigma_real,'.-',label='G_AMI')
plt.plot(Jdata[:, 2], Jdata[:, 3],'x-',label='Sym-det')
plt.xlabel('Mfreq',**axis_font)
plt.ylabel('Real Sigma22',**axis_font)
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(freq, result_sigma_imag,'.-',label='G_AMI')
plt.plot(Jdata[:, 2], Jdata[:, 4],'x-',label='Sym-det')
plt.ylabel('Imag sigma22',**axis_font)
plt.xlabel('Mfreq',**axis_font)
plt.tight_layout()
plt.legend()
plt.show()