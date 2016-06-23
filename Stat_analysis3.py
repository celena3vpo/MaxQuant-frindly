#!/usr/bin/python
import numpy as np
import math

# this function accepts an array of peptides with corresponding MS/MS count and
# H/(H+L) values from 1 protein, filters out outliers and calculates
# mean, median, standard deviation, standard error, finds the number of unique peptides 
# and sums up all the MS/MS counts, also assigns number of unique peptides to weight for fitting
def stat_analysis(peptide_array):
    peptides=[];
    msms_counts=0;
    ratios=[];
    for line in peptide_array:
        line.rstrip('\n');
        values = line.split('\t');
        peptides.append(values[0]);
        msms_counts=msms_counts+int(values[1])
        ratios.append(round(float(values[2]),3));
    num_unique_peptides=len(set(peptides))
    if num_unique_peptides>3:
        stats=find_outlier(ratios);
    elif num_unique_peptides==1:
        stats=(ratios[0], ratios[0], 0, 0);
    else:
        stats=mean_sterror(ratios);
    return num_unique_peptides, msms_counts,  stats, num_unique_peptides
        
          
# this part calculates the first and third quartile for interquartile range calculation
# to figure out which data values are out of range, the function accepts an input of a list
# of H/(H+L) ratios, after the outlier elimination (only 1 peptide measurement per protein)
# the function calls onto mean_sterror function       
        
def find_outlier(ratios):
    #ratios.sort();
    #np_ratios=np.ratios;
    p25=np.percentile(ratios,25);
    p75=np.percentile(ratios,75);
    iqr=p75-p25;
    lower_limit=p25-1.5*iqr;
    upper_limit=p75+1.5*iqr;
    if min(ratios) < lower_limit:
        ratios.remove(min(ratios));
        stats=mean_sterror(ratios);
    elif max(ratios) > upper_limit:
        ratios.remove(max(ratios));
        stats=mean_sterror(ratios);
    else:
        stats=mean_sterror(ratios);
    return stats

# this part calculates the mean, median, standard deviation, standard error of
# the heavy light ratio list    
    
def mean_sterror(ratios):
    mu=round(np.median(ratios), 3);
    mn=round(np.average(ratios), 3);
    std=round(np.std(ratios), 3);
    ste=round(std/math.sqrt(len(ratios)-1),3);
    return mu, mn, std, ste
    
    
    