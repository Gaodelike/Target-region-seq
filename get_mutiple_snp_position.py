#-*- coding:UTF-8 -*-
import sys
import re
f = open(sys.argv[1])
L = [line for line in f.readlines()]
L_name = L[0::2]#按行读取reference.fasta文件，生成序列名（rsID）的 list
L_sequence = L[1::2]#按行读取reference.fasta文件，生成序列本身的 list
double_snp_A=re.compile('[A]/[A-z]/[A-z]')
double_snp_C=re.compile('[C]/[A-z]/[A-z]')
double_snp_G=re.compile('[G]/[A-z]/[A-z]')
double_snp_T=re.compile('[T]/[A-z]/[A-z]')
single_ref_A=re.compile('[A]/[ACGT]')
single_ref_C=re.compile('[C]/[ACGT]')
single_ref_G=re.compile('[G]/[ACGT]')
single_ref_T=re.compile('[T]/[ACGT]')

def change_lower_base(x):
    L_sequence[x]=L_sequence[x].upper()
    filter_double_A=double_snp_A.sub('a',L_sequence[x])
    filter_double_C=double_snp_C.sub('c',filter_double_A)
    filter_double_G=double_snp_G.sub('g',filter_double_C)
    filter_double_T=double_snp_T.sub('t',filter_double_G)
    filter_A=single_ref_A.sub('a',filter_double_T)
    filter_C=single_ref_C.sub('c',filter_A)
    filter_G=single_ref_G.sub('g',filter_C)
    filter_T=single_ref_T.sub('t',filter_G)
    return filter_T

def get_position(x):
    position_in_amplicon=[]
    for i in xrange(len(change_lower_base(x))):
        a=re.match('[acgt]',change_lower_base(x)[i])
        if bool(a) == True:
            position_in_amplicon.append(i+1)
    return position_in_amplicon

def get_reference(x):
    ref_sequence=change_lower_base(x).upper()
    return ref_sequence

snp_file=open('snp-position','w')
for x in xrange(len(L_name)):
    id=L_name[x].strip()[1:]
    split_id=id.split('_')
    if len(split_id) == 1:
        t='%s\t%s\n'%(id,get_position(x)[0])
    elif len(split_id) >= 2:
        i=0
        s=''
        while i < len(get_position(x)):
            s=s+','+str(get_position(x)[i])
            i+=1
        t='%s\t'%(id)+s[1:]+'\n'
    snp_file.write(t)
snp_file.close()

snp_file2=open('snp-position')
all_position=[line for line in snp_file2.readlines()]
double_position=[]
triple_position=[]
rest_position=[]
for x in range(len(all_position)):
    if len(all_position[x].split('_')) == 2:
        t=all_position[x]
        second_pos=t.split('\t')[0]+'\t'+t.split(',')[1]
        double_position.append(second_pos)
    elif len(all_position[x].split('_')) == 3:
        s=all_position[x]
        second_site=s.split('\t')[0]+'\t'+s.split(',')[1]+'\n'
        third_pos=s.split('\t')[0]+'\t'+s.split(',')[2]
        double_position.append(second_site)
        triple_position.append(third_pos)
    elif len(all_position[x].split('_')) > 3:
        m=all_position[x]
        rest_pos=m.split('\t')[0]+'\t'+''.join(m.split(',')[3:])+'\n'
        rest_position.append(rest_pos)

snp_file2=open('snp-position2','w')
sec_pos=''.join(double_position)
snp_file2.write(sec_pos)
snp_file2.close()

snp_file3=open('snp-position3','w')
thr_pos=''.join(triple_position)
snp_file3.write(thr_pos)
snp_file3.close()

snp_file4=open('snp-position_rest','w')
rt_pos=''.join(rest_position)
snp_file4.write(rt_pos)
snp_file4.close()

w = open(sys.argv[1])
#生成一个去掉斜杠的reference.fasta文件
W = [l for l in w.readlines()]
W[-1]=W[-1].strip('/n')

W_name=W[0::2]
W_seq=W[1::2]
for x in range(len(W_seq)):
    W_seq[x] = W_seq[x].upper()
    W_seq[x] = re.sub('/G','',W_seq[x])
    W_seq[x] = re.sub('/A','',W_seq[x])
    W_seq[x] = re.sub('/T','',W_seq[x])
    W_seq[x] = re.sub('/C','',W_seq[x])
out_put = open(sys.argv[2],'w')#打开一个新文件，用来生成用于BWA 比对的reference.fasta文件（即去掉斜杠）
for i in range(len(W_name)):
    out_put.write(W_name[i])
    out_put.write(W_seq[i])
print 'The information of snp position is in snp-position'
print 'The sequences for alignment is in %s' % (sys.argv[2])
out_put.close()