#!/usr/bin/python
# this script is to prepare the MaxQuant evidence file for parameter fitting in Matlab,
# it also calculates various statistics (average, median, standard dev, error) and pulls out important info such as number of unique peptides 
# identified, total sum of spectral counts, which tells us how reliable the data is.
# For the purpose of Matlab curve fitting I only need the protein identifier (systematic name),
# average H/(H+L) over time, standard error to know the spread of the data and the sum of
# spectral counts to use as weights for fitting.The other info can be used as needed later on.

import Stat_analysis3
import sys

f = open(sys.argv[1], 'r');
# open the evidence.txt file in read only mode

file_name = sys.argv[1];
# create a variable name from the file name

filt_name = file_name + '_filtered_out';
analyz_name = file_name + '_analyzed';
matlab_name = file_name + '_matlab';
# create names for various files were the calculated, analyzed or excluded data will
# be stored, which is based on the original file name

filt = open(filt_name+'.txt', 'w');
analyzed = open(analyz_name+'.txt', 'w');
matlab = open(matlab_name+'.txt', 'w')
# now create writable text files with all of the names
# the final text file contains tables that need to go into matlab before grouping

first_line = f.readline();
# gets the first line, column headers

first_line=first_line.rstrip('\r\n');
# gets rid of the newline character

filt.write(first_line);
filt.write('\n');
analyzed.write(first_line);
analyzed.write('\tRatio H/(H+L)\n')
# the data that cannot be used (e.g. contaminant, reverse or no K peptides)
# will be simply copied to _filtered_out.txt file to keep track, while
# all the data that will go into analysis will be written into _analyzed.txt file
# with an aditional column Ratio H/(H+L) that will be calculated from Ratio H/L

# make all lower letters
first_line=first_line.lower();
headers = first_line.split('\t');
#this is for case insensitive matching

# get the column indices of interest
i0=headers.index("sequence");
i1=headers.index("k count");
i2=headers.index("leading razor protein");
i3=headers.index("raw file");
i4=headers.index("ms/ms count");
i5=headers.index("ratio h/l");
i6=headers.index("reverse");
i7=headers.index("potential contaminant");
#i7=headers.index("contaminant");

# this part of the program writes all the unusable data to a _filtered_out.txt 
#file and writes all the usable data into _analysed.txt file and keeps track of all the different raw files
# also we create a dictionary, where the keys are raw files and the values are arrays containing all the peptide
# characteristics to be analyzed from the corresponding raw file
# the dictionary is printed out based on keys into separate txt files which have 4 column entries - "Sequence", 
#"Leading Razor Protein", "MS/MS Count", "Ratio H/total"
raw_files = {};
for lines in f:
    entries=lines.split('\t');
    if entries[i1] == '0':
        filt.write('\t'.join(entries));
    elif entries[i5] == '' or entries[i5] == "NaN":
        filt.write('\t'.join(entries));
    elif entries[i6] == "+":
        filt.write('\t'.join(entries));
    elif entries[i7] == "+":
        filt.write('\t'.join(entries));
    else:
        line = '\t'.join(entries);
        analyzed.write(line.rstrip('\r\n'));
        analyzed.write('\t');
        ratio = float(entries[i5])/(float(entries[i5])+float(1));
        analyzed.write(str(ratio));
        analyzed.write('\n');
        line = (entries[i0]+'\t'+entries[i2]+'\t'+entries[i4]+'\t'+str(ratio)+'\n');
        if entries[i3] not in raw_files:
            raw_files[entries[i3]] = [line];
        else:
            raw_files[entries[i3]].append(line);
f.close();
filt.close();
analyzed.close();

# this part writes the contents of the raw files into separate text files

raw_file_header = ("Sequence"+'\t'+"Leading Razor Protein"+'\t'+"MS/MS Count"+'\t'+"Ratio H/(H+L)"+'\n')

for files in sorted(raw_files):
    f = open(file_name+files+'.txt', 'w');
    f.write(raw_file_header)
    for line in raw_files[files]:
        f.write(line);
    f.close();
# now that each raw file is separate, I need to read it again and organise entries by proteins
# so I create a dictionary of proteins for each file, the keys are the leading razor proteins,
# the values are arrays with corresponding peptides and characteristics

protein_analysis_headers = ("Protein systematic name"+'\t'+"Number of unique peptides"+'\t'+"Total MS/MS counts"+'\t'+"Average H/(H+L)"+'\t'+"Median H/(H+L)"+'\t'+"Standard deviation"+'\t'+"Standard error"+'\t'+"Weight"+'\n')
for files in sorted(raw_files):
    f = open(file_name+files+'.txt', 'r');
    f1 = open(file_name+files+'_processed.txt', 'w');
    f1.write(protein_analysis_headers);
    first_line = f.readline();
    headers = first_line.split('\t');
    i0=headers.index("Sequence");
    i2=headers.index("Leading Razor Protein");
    i4=headers.index("MS/MS Count");
    i8=headers.index("Ratio H/(H+L)\n");
    proteins={};
    for lines in f:
        entries=lines.split('\t');
        line = (entries[i0]+'\t'+entries[i4]+'\t'+entries[i8]);
        if entries[i2] not in proteins:
            proteins[entries[i2]]=[line];
        else:
            proteins[entries[i2]].append(line);
    for protein in proteins:
        l=Stat_analysis3.stat_analysis(proteins[protein]);
        f1.write(protein+'\t');
        f1.write(str(l[0])+'\t'+str(l[1])+'\t'+str(l[2][0])+'\t'+str(l[2][1])+'\t'+str(l[2][2])+'\t'+str(l[2][3])+'\t'+str(l[3])+'\n')
    f1.close()
    
# open a file to write the final analysis in the same directory
f2=open(file_name+'_final_analysis.txt', 'w');
# create an empty list of protein names to be populated with protein names from all files
protein_names=[];    

# create an empty list of dictioneries, that are going to contain protein=key, 
# the rest of descriptors=valus pairs
file_names=[];

for files in sorted(raw_files):
    f3 = open(file_name+files+'_processed.txt', 'r');
    f3.readline(); # get rid of the header row
    # make an empty dictionary named "filename" and add it to the file_name list
    file_name1={};
    # get rid of the new line character, split the line on tabs, add the protein name 
    #to the list of proteins, create a key-value pair from the protein name and the rest 
    #of the characteristics
    for lines in f3:
        lines=lines.rstrip('\n');
        entries=lines.split('\t');
        protein_names.append(entries[0])
        file_name1[entries[0]]=('\t'.join(entries[1:]));
    f3.close();
    file_names.append(file_name1);
# go through the list of protein names and write the protein names and values from each file
# next to each other
protein_analysis_header = ("Protein systematic name"+('\t'+"Number of unique peptides"+'\t'+"Total MS/MS counts"+'\t'+"Average H/(H+L)"+'\t'+"Median H/(H+L)"+'\t'+"Standard deviation"+'\t'+"Standard error"+'\t'+"Weight")*len(raw_files)+'\n')
f2.write(protein_analysis_header);
matlab_header = ("Protein systematic name"+('\t'+"Average H/(H+L)"+'\t'+"Standard error"+'\t'+"Weight")*len(raw_files)+'\n')
matlab.write(matlab_header);

# this part of the code pulls out the protein measuremnt characteristics from each file
# and prints them next to each other in the _final_analysis.txt file for side by side comparison.
# the order of values is based on string sorted file_name, so need to make sure that the sequential raw 
#files going to MaxQuant are sequentially (alphanumerically) named.
unique_protein_names=set(protein_names);
seven_tabs=('\t'*7);
three_tabs=('\t'*3);
for proteins in unique_protein_names:
    f2.write(proteins);
    matlab.write(proteins);
    sorted(file_names, key=str);
    for file_name1 in file_names:
        if proteins not in file_name1:
            f2.write(seven_tabs);
            matlab.write(three_tabs);
            continue;
        values=file_name1[proteins];
        val_entries = values.split('\t');
        f2.write('\t'+values); 
        matlab.write('\t'+val_entries[2]+'\t'+val_entries[5]+'\t'+val_entries[6]);
    f2.write('\n');
    matlab.write('\n');
f2.close();
matlab.close();