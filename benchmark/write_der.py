#!/usr/bin/python

import sys
import re


# delete blank lines
def del_blank(lines1):
    lines = []
    for line in lines1:
        if line.strip():
            lines.append(line)
    return lines


def search_der(file,flag):
    fp = open(file,'r')
    lines_orig = fp.readlines()
    fp.close()
   
    lines = del_blank(lines_orig);
    
    # find the line number
    LAYER_ID = []
    RAY_ID = [] 
    for i, line in enumerate(lines):
        if 'LAYER' in line:
            LAYER_ID = i
        if flag in line:
            RAY_ID.append(i)
    
    # find the number of layers
    NL = RAY_ID[0]-LAYER_ID-1
    print 'Number of layers: ', NL

    # derivatives for different frequencies 
    for i in range(len(RAY_ID)):
        start_ID = RAY_ID[i]     # write out the header
        end_ID = start_ID + NL -1 + 4
        
        # find the mode
        tmp = lines[start_ID]
        match=re.search('(?<=#)\s*(\d+)',tmp) 
        mode=match.group(0)

        # find the period
        tmp = lines[start_ID+1]
        match=re.search('(T\s*=\s*)(.*)C',tmp)
        period=match.group(2)
        
        # write out
        file1=str(int(float(period))).strip()+'.'+mode.strip()+'.'+flag+'.der.txt'
        
        fp = open(file1,'w+')
        for i in range(start_ID,end_ID+1):
            fp.write('%s' % (lines[i]))
        fp.close()

    return

def main():
    files=['SRDER.TXT','SLDER.TXT']
    flags=['RAYLEIGH','LOVE']
    for i, file in enumerate(files):
        search_der(file,flags[i])

if __name__=='__main__':
    main()
