## INTRODUCTION
This python package is designed to compute, analyze, and visualize the results of alternative TSSs based on the ChIP-Seq data.
## SUPPORTED PLATFORMS
Linux ,MacOS(到时候找一台mac试一试)
## PREREQUISITES
python 2.7.x  
Python packages dependencies:  
--numpy  
--h5py  
--pandas  
--tensorflow  
The pip install -r requirement.txt will try to install these packages automatically. However, please install them manually if, by any reason, the automatic installation fails.  
External packages dependencies:  
--bedtools  
Install instructions:[http://bedtools.readthedocs.io/en/latest/](http://bedtools.readthedocs.io/en/latest/)  
## INSTALLATION
First, create a python virtualenv

```
sudo virtualenv DeepMiRTss_env
```
Second, active the virtualenv and cd into it.

```
sudo source ./DeepMiRTss_env/bin/activate
cd ./DeepMiRTss_env
```
Third, download the software package.  

```
git clone https://github.com/jwanglabwangzhi/DeepMiRTss

```
Fouth, cd to the package directory,run requirements and setup.py script to install.

```
cd ./DeepMiRTss
pip install -r requirements.txt
python install setup.py install
```
Fifth, download the genomic data and DanQ model data and put it into the pre_load folder.

```
hg19:#把两个文件上传到我们的服务器，再把链接放上来
DanQ model:
```
## USAGE
The following 3 commands were provided by the DeepMiRTss package:  
1 DeepMiRTss.tssfinder, 2 DeepMiRTss.analysis, 3 DeepMiRTss.visulization.
### 1. DeepMiRTss.tssfinder  
This command is used to compute the alternative TSSs according to pol2 and histone ChIP-seq data.  
#### Parameters:  
<-h --help>  
<-p --polymerase_filename> [The RNA polymerase signal file ]  
RNA pol2 file is a necessity for the whole calculation process. You can use different kinds of  peak calling tools to obtain the signal file.  
<-k --h3k4me3_filename> [The h3k4me3 signal file]  
The h3k4me3 signal file can help to improve the accuracy of prediction , but it not a necessity.  
<-e --expressed_filename> [miRNA expression file]  
We suggest providing no more than 20 miRNAs at a time which reduce the waiting time.  
<-n --number_alternative_tss>[number of alternative TSSs]    
The predicted number of alternative TSSs, and the default value is 3, you can set any positive number to obtain multiple alternative TSSs.
#### Examples：
```
# If you only have pol2 signal file, you can run the programme to find TSSs.
DeepMiRTss.tssfinder -p ./example/pol2.bed 

# When you have both pol2 and histone signals, you can use the code below. Histone file can help improve the accuracy.
DeepMiRTss.tssfinder -p ./example/pol2.bed -k ./example/h3k4me3.bed

# When you have both pol2 and histone signal file and express file, code like this help you find TSSs based on expression file.
DeepMiRTss.tssfinder -p ./example/pol2.bed -k ./example/h3k4me3.bed -e ./example/express.bed

# You can identify only one TSS for a miRNA, or identify multiple alternative TSSs, as long as you change the -n parameter.
DeepMiRTss.tssfinder -p ./example/pol2.bed -n 1
DeepMiRTss.tssfinder -p ./example/pol2.bed -n 5
```
### 2. DeepMiRTss.analysis  
Make full use of DanQ model to calculate and analyze sequence features of region nearby TSSs.The DanQ model can refer to [``DanQ: a hybrid convolutional and recurrent neural network for predicting the function of DNA sequences''](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw226)
#### Parameters:  
<-h --help>  
<-u --upstream> [The distance value]  
The distance value between TSS and analysis region,and the default value is 0. You can choose the TSS region (TSS-500, TSS+500) as the analysis region which is the default or choose other region as you like.  
#### Examples：
```
# If you analysis TSSs region (TSS-500,  TSS+500), no parameter needs to be set.
DeepMiRTss.analysis
# If the upstream TSS 1000bp is you need, set -u as 1000.
DeepMiRTss.analysis -u 1000
```
### 3. DeepMiRTss.visulization  
We use JavaScript and Ajax to dynamically visualize the features of miRNA alternative TSSs. So you can find the feature value and compare different alternative TSSs at the same time.  
#### No parameter here
#### Examples：
```
DeepMiRTss.visuliation
```

## FAQ  
1 Signal file including RNA pol2 and histone must be in the format '[chr] [start] [end] [name] [score] [strand]'.  
2 If you don't provide -e parameter as a express file in 'DeepMiRTss.tssfinder' command, the program will compute all miRNA which might be last for a long time. 
So we suggest providing no more than 20 miRNAs at a time which reduce the waiting time.  
3 The command 'DeepMiRTss.analysis' must base on the result of command 'DeepMiRTss.tssfinder' which generates a 'miRNA_alternative_tss.bed' file.  
4 The command 'DeepMiRTss.visulization' must base on the result of command 'DeepMiRTss.analysis' which generates a 'y.csv' file.


## CREDITS
This software was developed by mcube bioinformatics group @ NanJing University  
Implemented by Wang Zhi.
## LICENSE
This software is under MIT license.  
See the LICENSE.txt file for details.  
## CONTACT
szxszx@foxmail.com

















