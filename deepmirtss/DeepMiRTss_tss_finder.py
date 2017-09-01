# -*- coding:utf-8 -*-

import sys
import getopt

import tss_finder as tsf

help_info = """
usage: DeepMiRTss_analysis [-h] [--help]
                           [-p] [--polymerase_filename] The RNA polymerase signal file is a necessity for the calculation process.
                           [-k] [--h3k4me3_filename] The h3k4m3 signal file can help to improve the accuracy of prediction .
                           [-e] [--expressed_filename] miRNA expression file, we suggest providing no more than 20 miRNAs at a time which reduce the waiting time.
                           [-n] [--number_alternative_tss] The predicted number of alternative TSSs,and the default value is 3.

"""

def main():
    opts,args = getopt.getopt(
        sys.argv[1:],'-h-p:-k:-e:-n:',
        ['help','polymerase_filename=','h3k4me3_filename=',
        'expressed_filename=','number_alternative_tss']
    )
    if not opts:
        print '...\nError, you should use at least a parmeter'
        print help_info
        exit()
    else:
        dic_dir = {}
        pol2_dir = None
        h3k4me3_dir = None
        expressed_dir = None
        alternative_num = 3
        for opt_name, opt_value in opts:
            if opt_name in ('-h', '--help'): 
                print help_info
                exit()
            if opt_name in ('-p','--polymerase_filename'): 
                pol2_dir = opt_value 
            dic_dir['pol2'] = pol2_dir
            if opt_name in ('-k','--h3k4me3_filename'): 
                h3k4me3_dir = opt_value
            dic_dir['h3k4me3'] = h3k4me3_dir
            if opt_name in ('-e','--expressed_filename'):
                expressed_dir = opt_value
            dic_dir['express'] = expressed_dir
            if opt_name in ('-n','--number_alternative_tss'):
                alternative_num = int(opt_value)
            dic_dir['alter_tss'] = alternative_num 

    if dic_dir['alter_tss'] <= 0:
        print '-n parameter must be greater than 0'
        exit()

    if dic_dir['pol2']:
        if dic_dir['h3k4me3']:
            file_detect = tsf.File_Detection(
            pol2_file = dic_dir['pol2'],
            h3k4me3_file = dic_dir['h3k4me3'],
            express_file = dic_dir['express']
        )           
            set_express = file_detect.return_set_express()
            if file_detect.signal_file_test():
                print 'file has been validated,\
    the program is running \nplease wait ... '
                tss_finder=tsf.Tss_Finder(set_express, dic_dir)
                tss_finder.call_tss_with_two_signal()         
        else:
            file_detect = tsf.File_Detection(
            pol2_file = dic_dir['pol2'],
            h3k4me3_file = dic_dir['h3k4me3'],
            express_file = dic_dir['express']
        )           
            set_express = file_detect.return_set_express()
            if file_detect.signal_file_test():
                print 'file has been validated,\
    the program is running \nplease wait ... '
                tss_finder=tsf.Tss_Finder(set_express, dic_dir)
                tss_finder.call_tss_with_one_signal()
    else:
        print 'Error,you should at least input pol2 signal value.'
        exit()

if __name__=='__main__':
    main()
