#!/usr/bin/env python
# -*- coding:utf-8 -*-
import cgi,cgitb
import urllib2
import os
import pandas as pd

form=cgi.FieldStorage()
name=form.getvalue('mirna_name')
show_num=10
y_df = pd.read_csv('./y.csv', index_col = 0)
mirna_sort_list = y_df.index.tolist()
feature_list = y_df.columns.tolist()
mir_sort_list = []
local_path=os.getcwd()
pre_load_path=local_path+'/pre_load/'
print 'Content-Type: text/html\n\n'
y_df = pd.read_csv('./y.csv', index_col = 0)
mirna_list = y_df.index.tolist()
feature_list = y_df.columns.tolist()

mir_sort_list = []
mirna_tf_dict={}
for mirna_sort in mirna_list:
    if mirna_sort[:-2] == name:
        mir_sort_list.append(mirna_sort)
mir_sort_list=sorted(mir_sort_list)
if len(mir_sort_list)>3:
    mir_sort_list=mir_sort_list[:3]


if len(feature_list)==919:
    for mirna_index,mirna_name in enumerate(mir_sort_list):
        tf_dic_tmp = {}
        num = 0
        for tf_name_score in zip(y_df.ix[mirna_name].index,y_df.ix[mirna_name]):
            num +=1
            tf_name=tf_name_score[0].split(':')[1].split('.')[0]
            tf_score=tf_name_score[1]
            if tf_name in tf_dic_tmp.keys():
                tf_dic_tmp[tf_name].append(tf_score)
            else:
                tf_dic_tmp[tf_name]=[]
                tf_dic_tmp[tf_name].append(tf_score)
            if num >=690:
                break
        tf_dic_tmp_sort ={} 
        for tf_key in tf_dic_tmp.keys():
            tf_dic_tmp_sort[tf_key] = float(sum(tf_dic_tmp[tf_key])) / len(tf_dic_tmp[tf_key])
        mirna_tf_dict[mirna_name]=tf_dic_tmp_sort

    all_0=[0 for x in range(170)] 
    num_data_array=all_0
    for mirna_name in mir_sort_list:
        new_add_array=mirna_tf_dict[mirna_name].values()
        num_data_array=[x+y for x,y in zip(num_data_array,new_add_array)]

    dict_value_index={}
    value_sort_list=[]
    value_index=0
    for value_data in num_data_array:    
        dict_value_index[value_data]=value_index
        value_index+=1
    for value_key in sorted(dict_value_index.keys(),reverse=True):
        value_sort_list.append(dict_value_index[value_key]) 

    tf_name_show_list=[mirna_tf_dict[name+'_1'].keys()[x] for x in value_sort_list][:show_num]
else:
    for mirna_index,mirna_name in enumerate(mir_sort_list):
        tf_dic_tmp = {}
        num = 0
        for tf_name_score in zip(y_df.ix[mirna_name].index,y_df.ix[mirna_name]):
            num+=1
            tf_name=tf_name_score[0]
            tf_score=tf_name_score[1]
            if tf_name in tf_dic_tmp.keys():
                tf_dic_tmp[tf_name].append(tf_score)
            else:
                tf_dic_tmp[tf_name]=[]
                tf_dic_tmp[tf_name].append(tf_score)
        tf_dic_tmp_sort ={} 
        for tf_key in tf_dic_tmp.keys():
            tf_dic_tmp_sort[tf_key] = float(sum(tf_dic_tmp[tf_key])) / len(tf_dic_tmp[tf_key])
        mirna_tf_dict[mirna_name]=tf_dic_tmp_sort

    all_0=[0 for x in range(len(mirna_tf_dict[name+'_1'].keys()))]
    num_data_array=all_0
    for mirna_name in mir_sort_list:
        new_add_array=mirna_tf_dict[mirna_name].values()
        num_data_array=[x+y for x,y in zip(num_data_array,new_add_array)]

    dict_value_index={}
    value_sort_list=[]
    value_index=0
    for value_data in num_data_array:    
        dict_value_index[value_data]=value_index
        value_index+=1
    for value_key in sorted(dict_value_index.keys(),reverse=True):
        value_sort_list.append(dict_value_index[value_key]) 

    tf_name_show_list=[mirna_tf_dict[name+'_1'].keys()[x] for x in value_sort_list][:show_num]

hjkl1 = '''
<!-- 为ECharts准备一个具备大小（宽高）的Dom --> 
<div id='main1' style='width: 90%;height:400px;'></div> 
<script type='text/javascript'> 
        // 基于准备好的dom，初始化echarts实例 
    var myChart = echarts.init(document.getElementById('main1')); 
    // 指定图表的配置项和数据
    var option = {
        title: {
            text: 'occurrent frequency of TFs in miRNA promoters'
        },
        tooltip: {
            trigger: 'axis'
        },
        legend: {
            right:10,
            data:''' 

hjkl2 = '''            
        },
        toolbox: {
            show : false,
            feature : {
                dataView : {show: true, readOnly: false},
                magicType : {show: true, type: ['line', 'bar']},
                restore : {show: true},
                saveAsImage : {show: true}
            }
        },   
        calculable : true,        
          
        xAxis: [
            {
                type : 'category',            
                data:''' 

hjkl3 = '''                
            }
        ],
        yAxis : [
            {
                type : 'value'
            }
        ],
               
        series: ['''
        
num_mir_sort = len(mir_sort_list)
hjkl4 =""
for mir_sort in mir_sort_list:
    tf_value_show_list=[mirna_tf_dict[mir_sort].values()[x] for x in value_sort_list][:show_num]
    num_mir_sort -= 1
    if num_mir_sort == 0:
        hjkl4 +="   {name:'%s', type: 'bar', data: %s,} "%(mir_sort, str(tf_value_show_list))
    else:
        hjkl4 +="   {name:'%s', type: 'bar', data: %s,}, "%(mir_sort, str(tf_value_show_list))
              
hjkl5 = """        
        ]
    };
    myChart.setOption(option);
</script>
"""

hjkl = hjkl1 + str(mir_sort_list) + hjkl2 + str(tf_name_show_list) + hjkl3 + hjkl4 + hjkl5


wasd1 = """
    <div id='main2' style='width: 100%;height:1200px;'></div>
    <script type='text/javascript'>    
        // 基于准备好的dom，初始化echarts实例
        var myChart2 = echarts.init(document.getElementById('main2'));
        // 指定图表的配置项和数据
        option = {
            backgroundColor: 'white',            
            title:{
                show: true,
                text: '\\n\\n TF-miRNA regulatory network' ,          
                fontSize:20
            },                   

            label: {
                normal: {
                    show: false,
                    textStyle: {
                        fontSize: 12
                    },
                }
            },			
			
			legend:{
				show:true,
				right:30,
				data:
			""" 
wasd2 = """                
            },
            animation: false,
            series: [{                
                type: 'graph',
                layout: 'force',
                symbol: 'circle',
                symbolSize: 50,
                roam: true,
                edgeSymbol: ['circle','arrow'],
                edgeSymbolSize: [0, 8],
                // 连接线上的文字
                focusNodeAdjacency: true, //划过只显示对应关系
                edgeLabel: {
                    normal: {
                        show: false,
                        textStyle: {
                            fontSize: 20
                        },
                        formatter: '{c}'
                    }
                },
                categories: [
				"""

wasd3=				"""
				],
                lineStyle: {
                    normal: {
                        opacity: 1,
                        width: 2,
                        curveness: 0
                    }
                },
                // 圆圈内的文字
                label: {
                    normal: {
                        show: true
                    }
                },
                force: {
                    repulsion: 8000,   //张力，同时也是线长
                    layoutAnimation:false
                },
                data: ["""
                
wasd4 = ""
for tf_name in tf_name_show_list:
        wasd4 += "{name:'%s',value:1,itemStyle: {normal: {color: '#7388AC', label: {position: 'bottom',textStyle: {color: '#7388AC'}}}},},"%tf_name
                         
wasd5 = ""
tf_mir_dict = {}

color_list=['#EA7E52','#74A274','#EFDD79']
for index,mir_name in enumerate(mir_sort_list):
    feature_color = color_list[index] 
    wasd5 += "    {name:'%s',category:1,itemStyle:{normal: {color: '%s',},emphasis:{color: '#000'}}},"%(mir_name, feature_color)

wasd5 = wasd5.strip(',')

wasd6 =""" ],     
    links: ["""

tf_ppi_con_list=[]
for ppi_con in open(pre_load_path+'/tf_170_ppi_connection.txt'):
    split_ppi_con=ppi_con.strip('\n').split('\t')
    if split_ppi_con[0] in tf_name_show_list and split_ppi_con[1] in tf_name_show_list:
        tf_ppi_con_list.append([split_ppi_con[0],split_ppi_con[1]])


wasd7 =""
for tf_name in tf_name_show_list:
    for mir_index,mir_name in enumerate(mir_sort_list):
        line_color=color_list[mir_index]
        wasd7 += "   {source: '%s',target: '%s',value:'',lineStyle:{normal:{color:'%s'}}  },"%(tf_name, mir_name,line_color)

for tf_ppi_con in tf_ppi_con_list:
    wasd7 += "   {source: '%s',target: '%s',value:'',lineStyle:{normal:{color:'#7388AC'}}  },"%(tf_ppi_con[0],tf_ppi_con[1])

wasd8 ="""                                                            ]
            }]
            };
    // 使用刚指定的配置项和数据显示图表。
    myChart2.setOption(option);       
    </script>"""

mir_sort_list.append('TFs')
mir_tf_color=[]
for mir_index,mir_name in enumerate(mir_sort_list):
	if mir_index+1 != len(mir_sort_list):
		cate_dict="{name: '%s',itemStyle:{color:'%s'}},"%(mir_name,color_list[mir_index])
	else:
		cate_dict="{name: '%s',itemStyle:{color:'#7388AC'}}"%mir_name
	wasd2+=cate_dict

if wasd4:
    wasd = wasd1+str(mir_sort_list) + wasd2 +wasd3 + wasd4 + wasd5 + wasd6 + wasd7 +wasd8
else:
    wasd =''

print

if name != None and len(name) > 0:   
    print hjkl
    print wasd

else:
    print "<h3 id=name>You should choose a miRNA !!!</h3>"  
