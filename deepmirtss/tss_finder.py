# -*- coding:utf-8 -*-

import os
import subprocess
import re

import numpy as np

import get_cnn_score as cnn_score

help_info = """
usage: DeepMiRTss_analysis [-h] [--help]
                           [-p] [--polymerase_filename] The RNA polymerase signal file is a necessity for the calculation process.
                           [-k] [--h3k4me3_filename] The h3k4m3 signal file can help to improve the accuracy of prediction .
                           [-e] [--expressed_filename] miRNA expression file, we suggest providing no more than 20 miRNAs at a time which reduce the waiting time.
                           [-n] [--number_alternative_tss] The predicted number of alternative TSSs,and the default value is 3.

"""

file_dir = os.getcwd()
pre_load_dir =file_dir + '/pre_load'
set_inter_mirna = set()
with open(pre_load_dir + '/2014_miriad_inter.bed') as inter_mirna_file:
    for inter_line in inter_mirna_file.readlines():
        inter_name=inter_line.split('\t')[3]
        set_inter_mirna.add(inter_name)
    


class File_Detection(object):
    def __init__(self, pol2_file = None,
    h3k4me3_file = None, express_file = None):
        self.express_file = express_file
        self.pol2_file = pol2_file
        self.h3k4me3_file = h3k4me3_file
        if not self.pol2_file:
            print 'Error,there must be pol2_signal file'
            print help_info
            exit()
        if self.express_file:       
            # try to use the express_file
            with open(self.express_file) as e_file:
                e_line = e_file.readline()
                col_num = len(e_line.split('\t'))
                if col_num != 1:
                    print 'Error:please ensure your expressed miRNA file \
has only one column in one row'
                    exit()
                else:
                    with open(self.express_file) as e_file:
                        e_list = e_file.readlines()
                        self.set_express_mirna = {e_mirna.replace('\n','')
                    .replace('\t','').replace('\r','') for e_mirna in e_list}
        else:
            # use miriad_inter-genic miRNA
            self.set_express_mirna = set_inter_mirna
  
                       
                    
    def return_set_express(self):
        set_inter_express = set_inter_mirna & self.set_express_mirna
        alter_tss = len(set_inter_express)
        num_intr_mirna = len(set_inter_mirna)
        if alter_tss == num_intr_mirna:
            return set_inter_express
        else:
            if alter_tss:
                print 'there is', alter_tss, 'intergenic miRNA in \
your expression file'
                return set_inter_express
            else:
                print 'Error:please ensure the human intergenic \
    miRNA with the correct format in mirbase'
                exit()
     

    def judge_each_row(self,one_row_line):
    # judge the format of bed file with six columns
        one_row_line = one_row_line.strip('\n').strip('\r').strip('\t')
        one_row=one_row_line.split('\t')
        if len(one_row) == 6:
            '''
            if one_row[-1].strip('\n')  in ['+','-','.']:
                pass
            else:
                print 'Error, the last column of your signal file should be strand like + , - or .'
                return False
            '''    
            if re.match('chr[\d|X|Y]',one_row[0]):
                pass
            else:
                print 'Error, the first column of your signal file should be chromosome number like chr1 or chrX.'
                return False
            if int(one_row[2]) > int(one_row[1]):
                pass
            else:
                print 'Error, the third column should be bigger than the second column in your signal file.'
                return False      
            try:
                float(one_row[4])
            except:
                print 'Error, the fifth column should be the float signal value in your signal file.'
                return False           
            return True
        
        else:
            print 'Error: your signal file should be 6 columns'
            return False
                 
    def signal_file_test(self):
    # this help test the signal file to satisfied the program input 
        if self.h3k4me3_file:
            self.h3k4me3_num = 0
            self.pol2_num=0
            for test_h3k4me3_line in open(self.h3k4me3_file):
                self.h3k4me3_num += 1
                if self.judge_each_row(test_h3k4me3_line):        
                    pass
                else:
                    print 'your h3k4me3 signal file has something wrong in line ', self.h3k4me3_num
                    exit()
            for test_pol2_line in open(self.pol2_file):
                self.pol2_num += 1
                if self.judge_each_row(test_pol2_line):        
                    pass
                else:
                    print 'your pol2 signal file has something wrong in line', self.pol2_num      
                    exit()
        else:
            self.pol2_num=0        
            for test_pol2_line in open(self.pol2_file):
                self.pol2_num += 1
                if self.judge_each_row(test_pol2_line):        
                    pass
                else:
                    print 'your pol2 signal file has something wrong in line', self.pol2_num   
                    exit()
        return True
                
class Tss_Finder(object):
    '''
    compute tss by different kinds of inputs combination
    '''
    def __init__(self, set_express, dic_dir):
        self.set_express = set_express
        self.dic_dir = dic_dir
        
    def get_pol2_inter_region(self):
        # find out the inter region by expressed mirna and intersect with pol2
        inter_region_file=open('inter_region.txt', 'w')
        for inter_region_line in open(
            pre_load_dir + '/inter_mirna_region.txt'
        ):
            if inter_region_line.split('\t')[3] in self.set_express:
                inter_region_file.write(inter_region_line)
        inter_region_file.close()
        os.system('bedtools intersect -a %s -b %s/inter_region.txt \
-wa -wb>inter_pol2.txt'%(self.dic_dir['pol2'], file_dir))
        os.remove('inter_region.txt')
        
    def call_tss_with_one_signal(self):
        self.get_pol2_inter_region()
        dict_mi_score = dict()
        for pol2_line in open('inter_pol2.txt'):
            pol2_new_line, mirna_name = self.inter_pol2_region_new_line(pol2_line)
            dict_mi_score.setdefault(mirna_name,[]).append(pol2_new_line.strip('\n')+'\t#\t#\t#\t#\t1\t#\n')
        mirna_tss_three=open('mirna_tss_three.txt','w')                               
        for i in dict_mi_score.keys():
            if len(dict_mi_score[i]) < self.dic_dir['alter_tss']:
                for mi in dict_mi_score[i]:
                    mirna_tss_three.write(mi[0])
            else:
                list_mi_tss=list()
                for mi in dict_mi_score[i]:
                    list_mi_tss.append(float(mi.split('\t')[4]))
                sort_list_mi_tss=sorted(list_mi_tss)
                for mi1 in dict_mi_score[i]:
                    list_sort_mi_tss=[]
                    for num_exp in range(self.dic_dir['alter_tss']):
                        list_sort_mi_tss.append(sort_list_mi_tss[-1 - num_exp])
                    if float(mi1.split('\t')[4]) in list_sort_mi_tss:
                        mirna_tss_three.write(mi1)                         
        mirna_tss_three.close()    
        os.system('sortBed -i mirna_tss_three.txt >sorted_mirna_tss_region.txt')
        os.remove('inter_pol2.txt')
        os.remove('mirna_tss_three.txt')
        print 'top', self.dic_dir['alter_tss'], 'regions \
finished\nprogram is searching tss ...'
        self.search_tss()
        
        print 'Finish, you can find results in \
miRNA_alternative_tss.bed in %s'%file_dir
        os.remove('sorted_mirna_tss_region.txt')        
           
               
    def call_tss_with_two_signal(self):
        self.get_pol2_inter_region()
        dict_mi_score = dict()       
        for pol2_line in open('inter_pol2.txt'):
            inter_pol2_region = open('inter_pol2_region.txt','w')
            pol2_new_line, mirna_name = self.inter_pol2_region_new_line(pol2_line)
            inter_pol2_region.write(pol2_new_line)
            inter_pol2_region.close()
            s=subprocess.Popen(
'bedtools intersect -a inter_pol2_region.txt \
-b %s -wa -wb'%self.dic_dir['h3k4me3'], shell=True, stdout=subprocess.PIPE
        )
            output_value = s.stdout.readlines()
            if output_value <> []:
                dict_mi_score.setdefault(mirna_name,[]).append(output_value)
        mirna_tss_three=open('mirna_tss_three.txt','w')
        for i in dict_mi_score.keys():           
            if len(dict_mi_score[i]) <= self.dic_dir['alter_tss']:
                for mi in dict_mi_score[i]:
                    mirna_tss_three.write(mi[0])
            else:
                list_mi_tss=list()
                for mi in dict_mi_score[i]:
                    list_mi_tss.append(float(mi[0].split('\t')[4]))
                sort_list_mi_tss=sorted(list_mi_tss)
                for mi1 in dict_mi_score[i]:
                    list_sort_mi_tss=[]
                    for num_exp in range(self.dic_dir['alter_tss']):
                        list_sort_mi_tss.append(sort_list_mi_tss[-1 - num_exp])
                    if float(mi1[0].split('\t')[4]) in list_sort_mi_tss:
                        mirna_tss_three.write(mi1[0])                         
        mirna_tss_three.close()    
        os.system('sortBed -i mirna_tss_three.txt >sorted_mirna_tss_region.txt')
        os.remove('inter_pol2.txt')
        os.remove('inter_pol2_region.txt')
        os.remove('mirna_tss_three.txt')

        
        print 'top', self.dic_dir['alter_tss'], 'regions \
finished\nprogram is searching tss ...'
        self.search_tss()
        
        print 'Finish, you can find results in \
miRNA_alternative_tss.bed in %s'%file_dir
        os.remove('sorted_mirna_tss_region.txt')



    def seq_translate_list(self,templet,seq):  
        final_list=[]
        templet=templet.strip()
        seq=seq.strip()
        seq=seq.lower()
        for i in seq:
            for j in templet:
                if i in j:
                    final_list.append(1.0)
                else:
                    final_list.append(0.0)
        return final_list

    def search_tss(self):
        # search the tss according the middle of signal region
        final_table=open('miRNA_alternative_tss.bed','w')
        for i in open('sorted_mirna_tss_region.txt'):
            bcd=open('bcd_best1.bed','w')
            split_line=i.split('\t')
            chrom=split_line[0]
            mirna_name=split_line[3]
            strand=split_line[5]
            pol2_score=float(split_line[4])
            h3k4me3_score=float(split_line[10])
            #tss_site=(int(split_line[1])+int(split_line[2]))/2
            middle=(int(split_line[1])+int(split_line[2]))/2
            split_num=11
            split_length=1000/(split_num-1)
            for i in range(split_num):
                line=chrom+'\t'+str(middle-500+split_length*i-100)+'\t'+str(middle-500+split_length*i+100)+'\t'+'#'+'\t'+'0'+'\t'+strand+'\n'
                bcd.write(line)
            bcd.close()
            os.system('bedtools getfasta -fi %s/hg19.fa \
-bed bcd_best1.bed -s -fo peak_with_fasta_best1.fasta'%pre_load_dir)
            list_a=[]
            for fasta_line in open('peak_with_fasta_best1.fasta'): 
                if fasta_line[0]=='>':
                    pass
                else:      
                    list_a.append(self.seq_translate_list('atcg',fasta_line.rstrip('\n')))
            array_data=np.array(list_a,dtype='float32').reshape(split_num,len(fasta_line.rstrip('\n')),4)
            score,conv_w=cnn_score.get_score(array_data)
            list_score=score[:,1].tolist()
            tss=str(middle-500+split_length*list_score.index(max(list_score)))
            final_line=chrom+'\t'+tss+'\t'+str(int(tss)+1)+'\t'+mirna_name+'\t'+str(pol2_score*h3k4me3_score*max(list_score))+'\t'+strand+'\n'
            #print final_line
            final_table.write(final_line)
        final_table.close()
        os.remove('bcd_best1.bed')
        os.remove('peak_with_fasta_best1.fasta')
                    
    def inter_pol2_region_new_line(self, line):
        # organize the line more comfortable
        split_line = line.split('\t')
        middle = int(int(split_line[1]) + int(split_line[2])) / 2  
        start_site = middle - 500
        end_site = middle + 500
        chr_name = split_line[0]
        mirna_name = split_line[9]
        mi_name = split_line[10]
        strand = split_line[-1].rstrip('\n')
        score = split_line[4]
        mirna_line = chr_name + '\t' + str(start_site) \
+ '\t' + str(end_site) + '\t' + mirna_name + '\t' + score + '\t' + strand + '\n'
        return mirna_line, mirna_name
           
''' 
b = File_Detection(pol2_file = 'a1.bed', h3k4me3_file = None, express_file = None)
set_b = b.return_set_express()
if b.signal_file_test():
    print 'file has been validated, the program is running \nplease wait ... '

#a=Tss_Finder(set_b,{'h3k4me3': 'a1.bed', 'express': None, 'pol2': 'a1.bed', 'alter_tss': 2})
#a.call_tss_with_two_signal()     
a1=Tss_Finder(set_b,{'h3k4me3': None, 'express': None, 'pol2': 'a1.bed', 'alter_tss': 2})
a1.call_tss_with_one_signal()            
'''
