# -* - coding:utf-8 -*-

import os
import sys
import getopt

import tensorflow as tf
import numpy as np
import h5py
import pandas as pd


def seq_translate_list(templet,seq):  
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



def conv_forward_naive(x, w, b, conv_param):
  out = None
  N,C,H,W = x.shape
  F,_,WW,HH = w.shape
  S = conv_param['stride']
  P = conv_param['pad']
  Ho = 1 + (H + 2 * P - HH) / S
  Wo = 1 + (W + 2 * P - WW) / S
  x_pad = np.zeros((N,C,H+2*P,W+2*P))
  x_pad[:,:,P:P+H,P:P+W]=x
  out = np.zeros((N,F,Ho,Wo))
  for f in xrange(F):
    for i in xrange(Ho):
      for j in xrange(Wo):
        # N*C*HH*WW, C*HH*WW = N*C*HH*WW, sum -> N*1
        out[:,f,i,j] = np.sum(x_pad[:, :, i*S : i*S+HH, j*S : j*S+WW] * w[f, :, :, :], axis=(1, 2, 3)) 
    out[:, f, :, :]+=b[f]
  cache = (x, w, b, conv_param)
  return out, cache



class LSTMcell(object):
    def __init__(self, input, D_input, D_cell,param_rnn, h_act=tf.tanh, init_h=None, init_c=None):
        self.input = input
        self.D_input = D_input
        self.D_cell = D_cell
        self.h_act = h_act
        self.type = 'lstm' 
        self.param_rnn=param_rnn
        if init_h is None and init_c is None:
          # If init_h and init_c are not provided, initialize them
          # the shape of init_h and init_c is [n_samples, D_cell]
            self.init_h = tf.matmul(self.input[0,:,:], tf.zeros([self.D_input, self.D_cell]))
            self.init_c = self.init_h
            self.previous = tf.stack([self.init_h, self.init_c])

    def sigmoid_hard(self,x):
        """Hard sigmoid."""
        return tf.minimum(1.0, tf.maximum(0.0, 0.25 * x + 0.5))

    def Step(self, previous_h_c_tuple, current_x):
        prev_h, prev_c = tf.unstack(previous_h_c_tuple)
        i = self.sigmoid_hard(
          tf.matmul(current_x, self.param_rnn[0]) + 
          tf.matmul(prev_h, self.param_rnn[1]) + 
          self.param_rnn[2])

        f = self.sigmoid_hard(
          tf.matmul(current_x, self.param_rnn[6]) + 
          tf.matmul(prev_h, self.param_rnn[7]) + 
          self.param_rnn[8])

        o = self.sigmoid_hard(
          tf.matmul(current_x, self.param_rnn[9]) + 
          tf.matmul(prev_h, self.param_rnn[10]) + 
          self.param_rnn[11])

        c = tf.tanh(
          tf.matmul(current_x, self.param_rnn[3]) + 
          tf.matmul(prev_h, self.param_rnn[4]) + 
          self.param_rnn[5])

        current_c = f*prev_c + i*c

        current_h = o*tf.tanh(current_c)

        return tf.stack([current_h, current_c])      


def RNN(cell, cell_b=None, merge='concat'):
    hstates = tf.scan(fn = cell.Step,
                  elems = cell.input,
                  initializer = cell.previous,
                  name = 'hstates')

    if cell.type == 'lstm':
        hstates = hstates[:,0,:,:]

    if cell_b is not None:
        input_b = tf.reverse(cell.input, axis=[0])
        b_hstates_rev = tf.scan(fn = cell_b.Step,
                  elems = input_b,
                  initializer = cell_b.previous, 
                  name = 'b_hstates')
        if cell_b.type == 'lstm':
            b_hstates_rev = b_hstates_rev[:,0,:,:]
        b_hstates = tf.reverse(b_hstates_rev, axis=[0])
        if merge == 'sum':
            hstates = hstates + b_hstates
        else:
            hstates = tf.concat(values=[hstates, b_hstates], axis=2)
        return hstates

help_info = """
usage: DeepMiRTss_analysis [-h] [--help]
                           [-u] [--upstream] The distance value between tss and analysis region ,and the default value is 0.
"""
def main():
    opts,args = getopt.getopt(
        sys.argv[1:],'-h-u:',
        ['help','upstream=']
    )


    u = 0
    for opt_name, opt_value in opts:
        if opt_name in ('-h', '--help'): 
            print help_info     #这个地方help信息写了一个简单过程，在github说明书上应该写详细，提供-u 可以是任意整数，但生物学意义上应该是给负数即上游，整数即下游，所以最好给负数。不给就是0。
            exit()
        if opt_name in ('-u','--upstream'): 
            u = opt_value

    file_dir = os.getcwd()
    pre_load_dir = file_dir + '/pre_load'

    file_exist_result = os.path.exists('./miRNA_alternative_tss.bed')
    if not file_exist_result:
        print 'miRNA_alternative_tss.bed file must be in the current directory'
        exit()

    set_mirna=set()
    mirna_alternative_tss_file = open('./miRNA_alternative_tss.bed')
    mirna_alternative_tss_list = mirna_alternative_tss_file.readlines()
    for mirna_line in mirna_alternative_tss_list:
        mirna_name = mirna_line.split('\t')[3]
        set_mirna.add(mirna_name)
    dict_mirna={}
    for mirna in set_mirna:
        for mirna_tss_line in mirna_alternative_tss_list:
            if mirna_tss_line.split('\t')[3] == mirna:
                dict_mirna.setdefault(mirna,[]).append(mirna_tss_line)
    mirna_alternative_tss_file.close()
    sort_mirna = open('sort_mirna_withtss.txt','w')
    for mirna_key in dict_mirna.keys():
        if len(dict_mirna[mirna_key]) == 1:
            split_mirna = dict_mirna[mirna_key][0].split('\t')
            new_line = split_mirna[0] + '\t' + split_mirna[1] + '\t' + split_mirna[2] + '\t' + split_mirna[3]+'_1' + '\t' + split_mirna[4] + '\t' + split_mirna[5]
            sort_mirna.write(new_line)
        else:
            list_mi_tss_score=list()
            for line in dict_mirna[mirna_key]:
                list_mi_tss_score.append(float(line.split('\t')[4]))        
            sort_list_mi_tss_score = sorted(list_mi_tss_score)
            for sort_mi_tss_index, sort_mi_tss_score in enumerate(sort_list_mi_tss_score):
                for line in dict_mirna[mirna_key]:
                    if float(line.split('\t')[4]) == sort_mi_tss_score:
                        line_split = line.split('\t')
                        new_line = line_split[0] + '\t' + line_split[1] + '\t' + line_split[2] + '\t' + line_split[3]+'_'+str(sort_mi_tss_index+1) + '\t' + line_split[4] + '\t' + line_split[5]

                        sort_mirna.write(new_line)
    sort_mirna.close()              

    promoter_region = open('promoter_region.bed','w')
    mirna_list = []
    sort_mirna_withtss_file = open('sort_mirna_withtss.txt')
    sort_mirna_withtss_list = sort_mirna_withtss_file.readlines()
    for line in sort_mirna_withtss_list: 
        u = int(u)
        split_line = line.split('\t')
        chr_name = split_line[0]
        strand = split_line[-1].strip('\n').strip('\r')
        mirna_name = split_line[3]

        mirna_score = split_line[4]
        tss_site = int(split_line[1])
        if strand == '+':
            promoter_site = tss_site + u
            new_line = chr_name +'\t' + str(promoter_site-500) +'\t' + str(promoter_site+500) + '\t' + mirna_name + '\t' +mirna_score + '\t' + strand +'\n'
        else:
            promoter_site = tss_site - u
            new_line = chr_name +'\t' + str(promoter_site-500) +'\t' + str(promoter_site+500) + '\t' + mirna_name + '\t' +mirna_score + '\t' + strand +'\n'
        promoter_region.write(new_line)
        mirna_list.append(mirna_name)    
    promoter_region.close()
    os.remove('sort_mirna_withtss.txt')

    genome_file_exist_result = os.path.exists('%s/hg19.fa'%pre_load_dir)
    if not genome_file_exist_result:
        print 'genome file hg19.fa must be in the %s'%pre_load_dir
        exit()

    os.system('bedtools getfasta -fi %s/hg19.fa \
    -bed promoter_region.bed -s -fo promoter_region.bed.fasta'%pre_load_dir)
    os.remove('promoter_region.bed')

    print 'The program is running. It may take a long time.\nPlease be patient...'

    num_seq = 0
    all_seq_list =[]
    seq_file = open('promoter_region.bed.fasta')
    seq_list = seq_file.readlines()
    for seq_line in seq_list:
        if seq_line[0] == '>':
            pass
        else:
            num_seq += 1
            seq_list = seq_translate_list('atcg',seq_line)
            all_seq_list.append(seq_list)
    all_seq_array=np.array(all_seq_list,'float32').reshape(num_seq,4,1,1000)
    seq_file.close()
    os.remove('promoter_region.bed.fasta')

    x_shape = (num_seq, 4, 1, 1000) #n,c,w,h
    w_shape = (320, 4, 1, 26) #f,c,ww,hh
    x=all_seq_array
    danq_model=h5py.File('%s/DanQ_bestmodel.hdf5'%pre_load_dir)
    con_w=danq_model['/layer_0/param_0'][...]
    con_b=danq_model['/layer_0/param_1'][...]
    conv_param = {'stride': 1, 'pad': 0}
    out, _ = conv_forward_naive(x, con_w, con_b, conv_param)
    out_tensor = tf.convert_to_tensor(out,dtype=tf.float32)
    sess=tf.Session()
    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())
    relu_1=tf.nn.relu(out_tensor)
    pool_1 = tf.nn.max_pool(relu_1, ksize=[1, 1, 1, 13], strides=[1, 1, 1, 13], padding='VALID')
    pool_1 =tf.reshape(pool_1,[num_seq,320,75])
    inputs_T= tf.transpose(pool_1, perm=[2,0,1])

    W_i=danq_model['/layer_3/param_0'][...]
    U_i=danq_model['/layer_3/param_1'][...]
    b_i=danq_model['/layer_3/param_2'][...]

    W_c=danq_model['/layer_3/param_3'][...]
    U_c=danq_model['/layer_3/param_4'][...]
    b_c=danq_model['/layer_3/param_5'][...]

    W_f=danq_model['/layer_3/param_6'][...]
    U_f=danq_model['/layer_3/param_7'][...]
    b_f=danq_model['/layer_3/param_8'][...]

    W_o=danq_model['/layer_3/param_9'][...]
    U_o=danq_model['/layer_3/param_10'][...]
    b_o=danq_model['/layer_3/param_11'][...]

    back_W_i=danq_model['/layer_3/param_12'][...]
    back_U_i=danq_model['/layer_3/param_13'][...]
    back_b_i=danq_model['/layer_3/param_14'][...]

    back_W_c=danq_model['/layer_3/param_15'][...]
    back_U_c=danq_model['/layer_3/param_16'][...]
    back_b_c=danq_model['/layer_3/param_17'][...]

    back_W_f=danq_model['/layer_3/param_18'][...]
    back_U_f=danq_model['/layer_3/param_19'][...]
    back_b_f=danq_model['/layer_3/param_20'][...]

    back_W_o=danq_model['/layer_3/param_21'][...]
    back_U_o=danq_model['/layer_3/param_22'][...]
    back_b_o=danq_model['/layer_3/param_23'][...]

    param_rnn_f=[W_i,U_i,b_i,W_c,U_c,b_c,W_f,U_f,b_f,W_o,U_o,b_o]
    param_rnn_b=[back_W_i,back_U_i,back_b_i,back_W_c,back_U_c,back_b_c,back_W_f,back_U_f,back_b_f,back_W_o,back_U_o,back_b_o]

    num_units=320        
    rnn_fcell = LSTMcell(input=inputs_T, D_input=num_units, D_cell=num_units,param_rnn=param_rnn_f)       
    rnn_bcell = LSTMcell(input=inputs_T, D_input=num_units, D_cell=num_units,param_rnn=param_rnn_b) 
    rnn0 = RNN(cell=rnn_fcell, cell_b=rnn_bcell)
    rnn1 = tf.reshape(rnn0, [num_seq,75*num_units*2])
    w_dense=danq_model['/layer_6/param_0'][...]
    b_dense=danq_model['/layer_6/param_1'][...]
    relu_2=tf.nn.relu(tf.matmul(rnn1,w_dense)+b_dense)
    w_dense_1=danq_model['/layer_8/param_0'][...]
    b_dense_1=danq_model['/layer_8/param_1'][...]
    sigmoid_1=tf.nn.sigmoid(tf.matmul(relu_2,w_dense_1)+b_dense_1)
    y = sess.run(sigmoid_1)

    feature_list = []
    for feature_line in open('%s/features.txt'%pre_load_dir):
        feature_split = feature_line.split('\t')
        feature = feature_split[1] + ':' + feature_split[2]
        feature_list.append(feature)
    y_array = np.array(y,'float32')
    y_df = pd.DataFrame(y_array, index = mirna_list, columns = feature_list) 
    y_df.to_csv('y.csv')
    print 'Finish, you can find results in \
y.csv in %s'%file_dir

if __name__ == '__main__':
    main()
