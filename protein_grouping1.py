#!/usr/bin/python
import re
import sys
#import os
# note this scrip is in the Python_analysis direcotry, so to provide text input 
#that is in another folder you need to provide the path, e.g "name of folder"/"name of file"

f = open(sys.argv[1], 'r');
file_name = sys.argv[1];
f1=open(file_name+"_nups_grouped.txt", 'w');
f2=open(file_name+"_kaps_grouped.txt", 'w');
f3=open(file_name+"_other_grouped.txt", 'w');

#open the file given as input to the script in read only form
# create an empty dictionary of protein systematic and common names as well as descriptions
protein_desc={};
prot=open("orf_trans_all_2015_Nup145.txt", 'r');
for line in prot:
    if re.match('^>', line):
        line=line.rstrip();
        #print line;
        line=line[1:];
        #print line;
        entries=line.split(' ');
        #print entries[0];
        #print entries[1];
        protein_description=(' '.join(entries[2:]));
        protein_desc[entries[0]]=(entries[1]+'\t'+protein_description);
    else:
        next;
        
header=f.readline();
names= header.split('\t');
new_header=(names[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(names[1:])));
f1.write(new_header);
f2.write(new_header);
f3.write(new_header);
for line in f:
    entries=line.split('\t');
    protein_stats=('\t'.join(entries[1:]));
    #if re.search('nuclear pore complex', protein_desc[entries[0]], flags=re.IGNORECASE):
    if re.search('nuclear pore complex', protein_desc[entries[0]], flags=re.IGNORECASE):
        if entries[0]=='YLR347C' or entries[0]=='YMR308C' or entries[0]=='YHR170W' or entries[0]=='YML107C':
            f2.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
        else:
            f1.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
    elif re.search('karyopherin', protein_desc[entries[0]], flags=re.IGNORECASE):
        if entries[0]=='YDL075W' or entries[0]=='YLR406C':
            f3.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
        else:
            f2.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
    elif re.search('mRNA export', protein_desc[entries[0]], flags=re.IGNORECASE):
        f2.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
    elif re.search('nuclear import', protein_desc[entries[0]], flags=re.IGNORECASE):
        f2.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats); 
    elif re.search('nucleocytoplasmic transport', protein_desc[entries[0]], flags=re.IGNORECASE):
        f2.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
    elif entries[0]=='YDR424C':
        f1.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);     
    else:
        f3.write(entries[0]+'\t'+protein_desc[entries[0]]+'\t'+protein_stats);
f.close();
f1.close();
f2.close();
f3.close();
prot.close();

        
        
