# -* - coding:utf-8 -*-


import numpy as np
import h5py
import tensorflow as tf
'''from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput'''


class Sequence_Translate(object):
    def __init__(self, input_data, base_template):
        self.input_data = input_data
        self.base_template = base_template

    def test_data(self):
    # test if input length is 1000
        seq_length = len(self.input_data)
        if seq_length != 1000:
            print 'seq_length should be 1000 bps,please check again'
            pass
        else:
            self.seq_translate_list()

    def seq_translate_array(self):
    # translate the seq to one hot according template
        sample_list = []
        self.base_templatet = self.base_template.strip()
        self.input_data = self.input_data.strip()
        self.input_data = self.input_data.lower()
        for i in self.input_data:
            for j in self.base_template:
                if i in j:
                    sample_list.append(1.0)
                else:
                    sample_list.append(0.0)
        # print len(sample_list)
        sample_array = np.array(sample_list, 'float32').reshape(1, 4, 1, 1000)
        return sample_array


class Conv_Layer_Feedforward(object):
    def __init__(self, input_x, param, conv_param):
        self.input_x = input_x
        self.filter_w = param[0]
        self.filter_b = param[1]
        self.conv_param = conv_param

    def conv_forward(self):
    # conv_feed_forward
        self.N, self.C, self.H, self.W = self.input_x.shape
        self.F, self._, self.WW, self.HH = self.filter_w.shape
        self.S = self.conv_param['stride']
        # A few space on both sides are not swept, but usually set to 0
        self.P = self.conv_param['pad']
        self.Ho = 1 + (self.H + 2 * self.P - self.HH) / self.S
        self.Wo = 1 + (self.W + 2 * self.P - self.WW) / self.S
        self.x_pad = np.zeros(
            (self.N, self.C, self.H + 2 * self.P, self.W + 2 * self.P) # 是可以算多个序列的。
        )
        # It doesn't work here. If pad is not 0, it will work.
        self.x_pad[
            :, :, self.P:self.P + self.H, self.P:self.P + self.W
        ] = self.input_x
        # x_pad = np.pad(x, ((0,), (0,), (P,), (P,)), 'constant')
        self.out = np.zeros((self.N, self.F, self.Ho, self.Wo))
        for f in xrange(self.F):
            for i in xrange(self.Ho):
                for j in xrange(self.Wo):
                    # N*C*HH*WW, C*HH*WW = N*C*HH*WW, sum -> N*1
                    self.out[:, f, i, j] = np.sum(
                        self.x_pad[
                            :, :, i * self.S:i * self.S + self.HH,
                            j * self.S:j * self.S + self.WW
                        ]
                        * self.filter_w[f, :, :, :], axis=(1, 2, 3) # 这里的*并不是矩阵的dot操作，而是按位置逐个乘。
                    )
                self.out[:, f, :, :] += self.filter_b[f]
        return self.out
# x_shape = (1, 4, 1, 1000) # n, c, h, c
# w_shape = (320, 4, 1, 26) # f, c, hh, ww


class Lstm_Layer_Feedforward(object):
    def __init__(
        self, input, D_input, D_cell, param_rnn, h_act = tf.tanh, init_h = None, init_c = None
    ):
        self.input = input
        self.D_input = D_input
        self.D_cell = D_cell
        self.h_act = h_act # 这里可以选择LSTM的hidden state的激活函数
        self.type = 'lstm' # 区分gru
        self.param_rnn = param_rnn
        
        # 如果没有提供最初的hidden state和memory cell，会全部初始为0
        if init_h is None and init_c is None:
          # If init_h and init_c are not provided, initialize them
          # the shape of init_h and init_c is [n_samples, D_cell]
            self.init_h = tf.matmul(self.input[0, :, :], tf.zeros([self.D_input, self.D_cell]))
            self.init_c = self.init_h
            self.previous = tf.stack([self.init_h, self.init_c])

    def sigmoid_hard(self, x):
        # Hard sigmoid
        return tf.minimum(1.0, tf.maximum(0.0, 0.25 * x + 0.5))

    def Step(self, previous_h_c_tuple, current_x):
        # 分离上一时刻的hidden state和memory cell的组合
        prev_h, prev_c = tf.unstack(previous_h_c_tuple)
        i = self.sigmoid_hard(
          tf.matmul(current_x, self.param_rnn[0]) + 
          tf.matmul(prev_h, self.param_rnn[1]) + 
          self.param_rnn[2]
        )

        f = self.sigmoid_hard(
          tf.matmul(current_x, self.param_rnn[6]) + 
          tf.matmul(prev_h, self.param_rnn[7]) + 
          self.param_rnn[8]
        )

        o = self.sigmoid_hard(
          tf.matmul(current_x, self.param_rnn[9]) + 
          tf.matmul(prev_h, self.param_rnn[10]) + 
          self.param_rnn[11]
        )

        c = tf.tanh(
          tf.matmul(current_x, self.param_rnn[3]) + 
          tf.matmul(prev_h, self.param_rnn[4]) + 
          self.param_rnn[5]
        )

        current_c = f * prev_c + i * c

        current_h = o * tf.tanh(current_c)

        return tf.stack([current_h, current_c])      


def RNN(cell, cell_b=None, merge='concat'):
    """
    该函数接受的数据需要是[n_steps, n_sample, D_output],
    函数的输出也是[n_steps, n_sample, D_output].
    如果输入数据不是[n_steps, n_sample, D_input],
    使用'inputs_T = tf.transpose(inputs, perm=[1,0,2])'.
    """
    # 正向rnn的计算
    hstates = tf.scan(fn = cell.Step,
                  elems = cell.input,
                  initializer = cell.previous,
                  name = 'hstates')
    # lstm的step经过scan计算后会返回4维tensor，
    # 其中[:,0,:,:]表示hidden state，
    # [:,1,:,:]表示memory cell，这里只需要hidden state
    if cell.type == 'lstm':
        hstates = hstates[:,0,:,:]
    # 如果提供了第二个cell，将进行反向rnn的计算
    if cell_b is not None:
        # 将输入数据变为反向
        input_b = tf.reverse(cell.input, axis=[0])
        # scan计算反向rnn
        b_hstates_rev = tf.scan(fn = cell_b.Step,
                  elems = input_b,
                  initializer = cell_b.previous, # 每个cell自带的初始值
                  name = 'b_hstates')
        if cell_b.type == 'lstm':
            b_hstates_rev = b_hstates_rev[:,0,:,:]
        # 用scan计算好的反向rnn需要再反向回来与正向rnn所计算的数据进行合并
        b_hstates = tf.reverse(b_hstates_rev, axis=[0])
        # 合并方式可以选择直接相加，也可以选择concat
        if merge == 'sum':
            hstates = hstates + b_hstates
        else:
            hstates = tf.concat(values=[hstates, b_hstates], axis=2)
        return hstates




def load_model(sample_sequence):
    sample_sequence=sample_sequence
    template = "atcg"
    seq_translate = Sequence_Translate(sample_sequence, template)
    seq_array = seq_translate.seq_translate_array()
    x = seq_array  
    danq_model = h5py.File('DanQ_bestmodel.hdf5')
    con_w = danq_model['/layer_0/param_0'][...]
    con_b = danq_model['/layer_0/param_1'][...]
    conv_w_b = [con_w, con_b]
    conv_param = {'stride': 1, 'pad': 0}
    conv_layer_feedforward = Conv_Layer_Feedforward(x, conv_w_b, conv_param)
    conv_out = conv_layer_feedforward.conv_forward()
    out_tensor = tf.convert_to_tensor(conv_out, dtype = tf.float32)
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())
    relu_1 = tf.nn.relu(out_tensor)
    pool_1 = tf.nn.max_pool(relu_1, ksize = [1, 1, 1, 13], 
        strides=[1, 1, 1, 13], padding = 'VALID'
    )
    pool_1 = tf.reshape(pool_1,[1,320,75])
    inputs_T = tf.transpose(pool_1, perm = [2, 0, 1]) # [75, 1, 320] 75是step（看danq论文的图）
    W_i = danq_model['/layer_3/param_0'][...]
    U_i = danq_model['/layer_3/param_1'][...]
    b_i = danq_model['/layer_3/param_2'][...]

    W_c = danq_model['/layer_3/param_3'][...]
    U_c = danq_model['/layer_3/param_4'][...]
    b_c = danq_model['/layer_3/param_5'][...]

    W_f = danq_model['/layer_3/param_6'][...]
    U_f = danq_model['/layer_3/param_7'][...]
    b_f =danq_model['/layer_3/param_8'][...]

    W_o = danq_model['/layer_3/param_9'][...]
    U_o = danq_model['/layer_3/param_10'][...]
    b_o = danq_model['/layer_3/param_11'][...]

    back_W_i = danq_model['/layer_3/param_12'][...]
    back_U_i = danq_model['/layer_3/param_13'][...]
    back_b_i = danq_model['/layer_3/param_14'][...]

    back_W_c = danq_model['/layer_3/param_15'][...]
    back_U_c = danq_model['/layer_3/param_16'][...]
    back_b_c = danq_model['/layer_3/param_17'][...]

    back_W_f = danq_model['/layer_3/param_18'][...]
    back_U_f = danq_model['/layer_3/param_19'][...]
    back_b_f = danq_model['/layer_3/param_20'][...]

    back_W_o = danq_model['/layer_3/param_21'][...]
    back_U_o = danq_model['/layer_3/param_22'][...]
    back_b_o = danq_model['/layer_3/param_23'][...]
    param_rnn_f = [
        W_i, U_i, b_i, W_c, U_c, b_c, W_f, U_f, b_f, W_o, U_o, b_o
    ]
    param_rnn_b = [
        back_W_i, back_U_i, back_b_i, back_W_c, back_U_c, back_b_c,
        back_W_f, back_U_f, back_b_f, back_W_o, back_U_o, back_b_o
    ]    
    
    
    
    num_units = 320
    # 对参数的疑问必须考虑，论文给的权重w和u都是（320,320），而w应该是(D_input, self.D_cell)    
    rnn_fcell = Lstm_Layer_Feedforward(input=inputs_T, D_input=num_units, D_cell=num_units,param_rnn=param_rnn_f)   
    rnn_bcell = Lstm_Layer_Feedforward(input=inputs_T, D_input=num_units, D_cell=num_units,param_rnn=param_rnn_b)
    #print rnn_bcell.shape
    rnn0 = RNN(cell=rnn_fcell, cell_b=rnn_bcell)
    rnn0 = tf.reshape(rnn0, [1,75 * num_units * 2])
    
    w_dense=danq_model['/layer_6/param_0'][...]
    b_dense=danq_model['/layer_6/param_1'][...]
    relu_2=tf.nn.relu(tf.matmul(rnn0,w_dense)+b_dense)

    w_dense_1=danq_model['/layer_8/param_0'][...]
    b_dense_1=danq_model['/layer_8/param_1'][...]
    sigmoid_1=tf.nn.sigmoid(tf.matmul(relu_2,w_dense_1)+b_dense_1)
    abc=sess.run(sigmoid_1)
    for s in abc[0]:
        if s>0.1:
            print s
    return abc
    
def main():
    sample_sequence = "tataatcaattataatcaattataatcaattataatcaattataat\
    caattataatcaattataatcaattataatcaattataatcaattataatcaattataatc\
    aattataatcaattataatcaattataatcaattataatcaattataatcaattataatca\
    attataatcaattataatcaattataatcaattataatcaattataatcaattataatcaa\
    ttataatcaattataatcaattataatcaattataatcaattataatcaattataatcaat\
    tataatcaattataatcaattataatcaattataatcaattataatcaattataatcaatt\
    ataatcaattataatcaattataatcaattataatcaattataatcaattataatcaatta\
    taatcaattataatcaattataatcaattataatcaattataatcaattataatcaattat\
    aatcaattataatcaattataatcaattataatcaattataatcaattataatcaattata\
    atcaattataatcaattataatcaattataatcaattataatcaattataatcaattataa\
    tcaattataatcaattataatcaattataatcaattataatcaattataatcaattataat\
    caattataatcaattataatcaattataatcaattataatcaattataatcaattataatc\
    aattataatcaattataatcaattataatcaattataatcaattataatcaattataatca\
    attataatcaattataatcaattataatcaattataatcaattataatcaattataatcaa\
    ttataatcaattataatcaattataatcaattataatcaattataatcaattataatcaat\
    tataatcaattataatcaattataatcaattataatcaattataatcaattataatcaatt\
    ataatcaattataatcaattataatcaattataatcaat"
    abc=load_model(sample_sequence)
    
    
if __name__ == '__main__':
    main()

'''











rnn_fcell = LSTMcell(input=inputs_T, D_input=num_units, D_cell=num_units,param_rnn=param_rnn_f)       
rnn_bcell = LSTMcell(input=inputs_T, D_input=num_units, D_cell=num_units,param_rnn=param_rnn_b) 
#print rnn_bcell.shape
rnn0 = RNN(cell=rnn_fcell, cell_b=rnn_bcell)
#abc=sess.run(rnn0)
#print abc,abc.shape
rnn1 = tf.reshape(rnn0, [1,75*num_units*2])

w_dense=danq_model['/layer_6/param_0'][...]
b_dense=danq_model['/layer_6/param_1'][...]
relu_2=tf.nn.relu(tf.matmul(rnn1,w_dense)+b_dense)

w_dense_1=danq_model['/layer_8/param_0'][...]
b_dense_1=danq_model['/layer_8/param_1'][...]
sigmoid_1=tf.nn.sigmoid(tf.matmul(relu_2,w_dense_1)+b_dense_1)
abc=sess.run(sigmoid_1)
print abc




#print rnn0.shape
'''