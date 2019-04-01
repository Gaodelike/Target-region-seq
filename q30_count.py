# -*- coding:utf-8 -*-
from __future__ import division
from pysam import FastqFile
import os
import sys
import re

read1 = sys.argv[1]
read2 = sys.argv[2]
if os.path.splitext(read1.strip())[-1] == '.gz' and os.path.splitext(read2.strip())[-1] == '.gz':
    gz_handle1 = os.popen( 'gunzip -cd %s' % read1)
    gz_handle2 = os.popen( 'gunzip -cd %s' % read2)
else:
    gz_handle1 = read1
    gz_handle2 = read2
basename1 = os.path.basename(read1)
basename1 = re.match('(\S+)_R1\.fastq(\.gz)?',basename1).group(1)
basename2 = os.path.basename(read2)
basename2 = re.match('(\S+)_R2\.fastq(\.gz)?',basename2).group(1)
if basename1 != basename2:
    raise 'Two Read are not mapped!'
cwd = os.getcwd()
out_handle = open('%s/%s.stat'%(cwd,basename1),'w')
out_handle.write('AllReadsNum\tRead1_PE150_ReadsNum\tRead2_PE150_ReadsNum\tUseful_ReadsNum(Read1>=150 and Read2>=150)\tAll_Bases\tRead1_Q30_Bases(PE150)\tRead2_Q30_Bases(PE150)\tQ30_PE_Reads(Q30>50%)\tUseful_Bases(All)\tUseful_Ratio\n')

AllReadsNum = 0
AllBases = 0
Read1_PE150_ReadsNum = 0
Read2_PE150_ReadsNum = 0
Useful_ReadsNum = 0
Read1_Q30_Bases = 0
Read2_Q30_Bases = 0
Q30_PE_Reads = 0
Useful_Bases = 0

def PE150(seq):
    if len(seq) >= 150:
        return True
    else:
        return False

def Q30(qual_list):
    num = 0
    for qual in qual_list:
        if qual >= 30:
            num += 1
    return num

reads2 = FastqFile(gz_handle2)
for read1 in FastqFile(gz_handle1):
    read2 = reads2.next()
    seq1 = read1.sequence
    qual1 = read1.get_quality_array()
    seq2 = read2.sequence
    qual2 = read2.get_quality_array()
    AllReadsNum += 1
    AllBases += len(seq1)
    AllBases += len(seq2)
    R1_150 = PE150(seq1)
    R2_150 = PE150(seq2)
    if R1_150 and R2_150:
        Useful_ReadsNum +=1
        R1_Q30 = Q30(qual1)
        R2_Q30 = Q30(qual2)
        Read1_Q30_Bases += R1_Q30
        Read2_Q30_Bases += R2_Q30
        if ( R1_Q30 / len(seq1) >= 0.5 ) and ( R2_Q30 / len(seq2) >= 0.5 ):
            Q30_PE_Reads += 1
            Useful_Bases += R1_Q30
            Useful_Bases += R2_Q30
    elif R1_150:
        Read1_PE150_ReadsNum += 1
    elif R2_150:
        Read2_PE150_ReadsNum += 1

Useful_Ratio = Useful_Bases / AllBases
out_handle.write('%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n'%(AllReadsNum,Read1_PE150_ReadsNum,Read2_PE150_ReadsNum,Useful_ReadsNum,AllBases,Read1_Q30_Bases,Read2_Q30_Bases,Q30_PE_Reads,Useful_Bases,Useful_Ratio))
out_handle.close()
