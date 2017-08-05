# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------


#         This file is used to generate DNA sequence scores.
#         CNN parameters are stored in the tmp folder.


#------------------------------------------------------------------------
import time
import commands,os,sys

import tensorflow as tf 
import numpy as np

path = sys.path[0]

number=0





def get_score(array_data):
    x_channels = tf.reshape(array_data, [-1, 200,1,4])
    sess = tf.Session()
    new_saver = tf.train.import_meta_graph('pre_load/model_parameters/my-model.meta.meta')
    new_saver.restore(sess, tf.train.latest_checkpoint('pre_load/model_parameters/'))
    all_vars = tf.trainable_variables()
    conv=tf.nn.conv2d(x_channels,all_vars[0],strides=[1,1,1,1],padding='VALID')
    act=tf.nn.relu(conv+all_vars[1])
    pool_1=tf.nn.max_pool(act,ksize=[1,8,1,1],strides=[1,8,1,1],padding='VALID')
    flatten=tf.reshape(pool_1,[-1,16*((200-16+1)*2/16)])
    fc_1=tf.nn.relu(tf.matmul(flatten,all_vars[2])+all_vars[3])
    fc_2=tf.matmul(fc_1,all_vars[4])+all_vars[5]
    y_=tf.nn.softmax(tf.matmul(fc_2, all_vars[6] )+all_vars[7])
    return sess.run(y_),sess.run(all_vars[0])





def main():
    #seq='tataattata'*20
    num2seq=['A','T','C','G']
    list_a=[]
    list_a.append(seq_translate_list('atcg',seq.rstrip('\n')))
    array_data=np.array(list_a,dtype='float32').reshape(1,len(seq.rstrip('\n')),4)
    seq_model=get_seq_score()
    score,conv_w=seq_model.get_score(array_data)
    return score[:,0]    
    #for i in range(conv_w.shape[0]):
        #print ''.join([num2seq[np.argmax(conv_w[i,0,:,j])] for j in range(conv_w.shape[3])])

    


if __name__=='__main__':
    main()

	
