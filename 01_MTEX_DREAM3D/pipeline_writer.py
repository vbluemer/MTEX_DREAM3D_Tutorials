# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:30:50 2023

@author: x
"""
## ----------------------------------------- DESCRIPTION -------------------------------------- ##
## Automatic generation of Dream3D json files. Data can be put into these files without using
## the GUI of Dream3D. 
## ----------------------------------------- \DESCRIPTION ---------------------------------- ##
import json
import math 
import numpy as np


def write_phase(**kwargs):
    min_ESD         = math.exp(kwargs['mu']-kwargs['min_cutoff']*kwargs['sigma'])
    max_ESD         = math.exp(kwargs['mu']+kwargs['max_cutoff']*kwargs['sigma'])
    #mean_ESD        = math.exp(mu+(sigma**2/2))

    bin_number      = int(np.ceil((max_ESD-min_ESD)/0.5))
    bin_list        = list(min_ESD + np.arange(bin_number)*kwargs['bin_size'])


    d['0']['StatsDataArray'][kwargs['phase_ID']]['Bin Count']                                                         = bin_number
    
    d['0']['StatsDataArray'][kwargs['phase_ID']]['BinNumber']                                                         = bin_list
    

    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Distribution']                                          = {}
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Distribution']['Average']                               = kwargs['mu']
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Distribution']['Standard Deviation']                    = kwargs['sigma']
    
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs B Over A Distributions']                             = {}
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs B Over A Distributions']['Alpha']                    = kwargs['b_a_alpha']*bin_number
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs B Over A Distributions']['Beta']                     = kwargs['b_a_beta']*bin_number
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs B Over A Distributions']['Distribution Type']        = 'Beta Distribution'
    
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs C Over A Distributions']                             = {}
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs C Over A Distributions']['Alpha']                    = kwargs['c_a_alpha']*bin_number
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs C Over A Distributions']['Beta']                     = kwargs['c_a_beta']*bin_number
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs C Over A Distributions']['Distribution Type']        = 'Beta Distribution'
    
    
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs Omega3 Distributions']                               = {}
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs Omega3 Distributions']['Alpha']                      = kwargs['featuresize_omega3_mu']*bin_number
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs Omega3 Distributions']['Beta']                       = kwargs['featuresize_omega3_sigma']*bin_number
    d['0']['StatsDataArray'][kwargs['phase_ID']]['FeatureSize Vs Omega3 Distributions']['Distribution Type']          = 'Beta Distribution'
    
    d['0']['StatsDataArray'][kwargs['phase_ID']]['Feature_Diameter_Info']                                             = [kwargs['bin_size'], max_ESD, min_ESD]



## ------------------------ MAIN PROGRAM ------------------------------------## 
with open('minimal_pipeline.json') as f:
    d = json.load(f)
    
    

## ------------------------ USER INPUTS -------------------------------------## 
write_phase(phase_ID        = '1',
            mu              = 1.2,
            sigma           = 0.6,
            bin_size        = 0.8,
            max_cutoff      = 5,
            min_cutoff      = 5,
            b_a_alpha                    = [15],   
            b_a_beta                     = [1.5],   
            c_a_alpha                    = [3.5],  
            c_a_beta                     = [50],   
            featuresize_omega3_mu        = [10],  
            featuresize_omega3_sigma     = [1.8])
## ----------------- END OF USER INPUTS -------------------------------------## 

with open('generated_pipeline.json', 'w', encoding='utf-8') as f:
    json.dump(d, f, ensure_ascii=False, indent=4)

## -------------------- END MAIN PROGRAM ------------------------------------## 
