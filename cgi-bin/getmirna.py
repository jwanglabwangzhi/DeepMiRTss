#!/usr/bin/env python
# -*- coding:utf-8 -*-
import cgi,cgitb
import urllib2

import pandas as pd

form=cgi.FieldStorage()
name=form.getvalue('mirna_name')
thre_value=0.7
thre_add_value = (1.0 - thre_value)/3.0
thre_low = thre_value
thre_middle = thre_low + thre_add_value
thre_high = thre_middle + thre_add_value

print 'Content-Type: text/html\n\n'
y_df = pd.read_csv('./y.csv', index_col = 0)
mirna_sort_list = y_df.index.tolist()
feature_list = y_df.columns.tolist()
mir_sort_list = []
feature_set = set()
for mirna_sort in mirna_sort_list:
    if mirna_sort[:-2] == name:
        mir_sort_list.append(mirna_sort)
        for mirna_sort_value_index in zip(y_df.ix[mirna_sort], y_df.ix[mirna_sort].index):
            if mirna_sort_value_index[0] > thre_low:                     
                feature_set.add(mirna_sort_value_index[1])        
feature_sort_list = sorted([x for x in feature_set])        

hjkl1 = '''
<!-- 为ECharts准备一个具备大小（宽高）的Dom --> 
<div id='main1' style='width: 90%;height:400px;'></div> 
<script type='text/javascript'> 
        // 基于准备好的dom，初始化echarts实例 
    var myChart = echarts.init(document.getElementById('main1')); 
    // 指定图表的配置项和数据
    var option = {
        title: {
            text: 'miRNA with different features \\nin different tss regions'
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
    num_mir_sort -= 1
    if num_mir_sort == 0:
        hjkl4 +="   {name:'%s', type: 'bar', data: %s,} "%(mir_sort, str(y_df.ix[mir_sort,feature_sort_list].tolist()))
    else:
        hjkl4 +="   {name:'%s', type: 'bar', data: %s,}, "%(mir_sort, str(y_df.ix[mir_sort,feature_sort_list].tolist()))
              
hjkl5 = """        
        ]
    };
    myChart.setOption(option);
</script>
"""

hjkl = hjkl1 + str(mir_sort_list) + hjkl2 + str(feature_sort_list) + hjkl3 + hjkl4 + hjkl5


wasd1 = """
    <div id='main2' style='width: 100%;height:1200px;'></div>
    <script type='text/javascript'>    
        // 基于准备好的dom，初始化echarts实例
        var myChart2 = echarts.init(document.getElementById('main2'));
        // 指定图表的配置项和数据
        option = {
            backgroundColor: 'rgb(128, 128, 128)',            
            title:{
                show: true,
                text: '\\n\\n The picture shows features with different intensities' ,          
                fontSize:20
            },                   
            tooltip: {
                show: false
            },
            legend: {
                x: 'left',
                data:'feature'""" 
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
                categories: [{
                    name: ' ',
                    itemStyle: {
                        normal: {
                            color: '#009',
                        }
                    }
                }, {}],
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
                    repulsion: 900,   //张力，同时也是线长
                    layoutAnimation:false
                },
                data: ["""
                
wasd3 = ""
for mir_sort in mir_sort_list:
        wasd3 += "{name:'%s',value:1,itemStyle: {normal: {color: '#f90', label: {position: 'bottom',textStyle: {color: '#f90'}}}},},"%mir_sort
                         
wasd4 = ""
mir_feature_dict = {}
for mir_sort in mir_sort_list:  
    for feature_sort in feature_sort_list:
        new_feature_sort = feature_sort + mir_sort[-2:]
        feature_color = ''
        if thre_low < y_df.ix[mir_sort, feature_sort] <= thre_middle:     
            feature_color ='#999'
            mir_feature_dict.setdefault(mir_sort,[]).append(new_feature_sort)
        elif thre_middle < y_df.ix[mir_sort, feature_sort] <= thre_high:      
            feature_color ='#099'
            mir_feature_dict.setdefault(mir_sort,[]).append(new_feature_sort)
        elif thre_high < y_df.ix[mir_sort, feature_sort]:          
            feature_color ='#009'
            mir_feature_dict.setdefault(mir_sort,[]).append(new_feature_sort)
        if feature_color:
            wasd4 += "    {name:'%s',category:1,itemStyle:{normal: {color: '%s',},emphasis:{color: '#000'}}},"%(new_feature_sort, feature_color)
        else:
            pass   
wasd4 = wasd4.strip(',')

wasd5 =""" ],     
    links: ["""

wasd6 =""
for mir_sort in mir_feature_dict.keys():
    for feature_sort in mir_feature_dict[mir_sort]:
        wasd6 += "   {source: '%s',target: '%s',value:'feature'  },"%(mir_sort, feature_sort)

wasd7 ="""                                                            ]
            }]
            };
    // 使用刚指定的配置项和数据显示图表。
    myChart2.setOption(option);       
    </script>"""

if wasd4:
    wasd = wasd1 + wasd2 + wasd3 + wasd4 + wasd5 + wasd6 + wasd7
else:
    wasd =''

print

if name != None and len(name) > 0:   
    #print str(feature_sort_list)
    if feature_sort_list:
        print hjkl
        print wasd
    else:
        print "<h3 id=name>This miRNA's promoter regulatory region does not meet the requirements of the DanQ model</h3>"
        print "<img src='/image/smile1.PNG'/>"
else:
    print "<h3 id=name>You should choose a miRNA !!!</h3>"   
    print "<img src='/image/smile0.gif'/>"
