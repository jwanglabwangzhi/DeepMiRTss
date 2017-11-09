# -*- coding:utf-8 -*-

import os, sys
import getopt
import re
from BaseHTTPServer import HTTPServer
from CGIHTTPServer import CGIHTTPRequestHandler

import pandas as pd


help_info = """
usage: DeepMiRTss.visulization [-h] [--help]
                           [-t] [--threhold] The threhold value you set to filer\
the features for visulization. Please set it between (0,1). The default value is 0.5.
                           
"""

# 必须判断是否存在'y.csv'文件
def main():
    opts,args = getopt.getopt(
        sys.argv[1:],'-h-t:',
        ['help', 'threshold=']
    )
    thre_value = 0.5
    for opt_name, opt_value in opts:
        if opt_name in ('-h', '--help'): 
            print help_info
            exit()
        if opt_name in ('-t','--threshold'): 
            thre_value = float(opt_value)
    if thre_value> 0 and thre_value < 1:
        pass
    else:
        print 'You should set the threhold between (0,1).'
        exit()


    file_dir = os.getcwd()
    cgi_bin_dir = file_dir + '/cgi-bin'
    with open('%s/getmirna.py'%cgi_bin_dir,'r') as r:
        file_read=r.read()
    with open('%s/getmirna.py'%cgi_bin_dir,'w') as w:
        line_9 = 'thre_value='+str(thre_value)
        pattern = 'thre_value=0(\.\d*)'
        file_read_change=re.sub(pattern, line_9, file_read)
        w.write(file_read_change)


    file_exist_result = os.path.exists('./y.csv')
    if file_exist_result:
        pass
    else:
        print "Error ,visulization system must be based on the result \
    of DeepMiRTss_analysis.py .\nIn other words, you must have 'y.csv' \
    under the current folder"
        exit()

    y_df = pd.read_csv('./y.csv')
    mirna_sort_list = y_df.iloc[:,0].tolist()
    set_mirna = set()
    for mirna_sort in mirna_sort_list:
        set_mirna.add(mirna_sort[:-2])
    mirna_list = [x for x in set_mirna]
    sort_mirna_list = sorted(mirna_list)

    visu_html = open('visulization.html','w')
    visu_html.write('''
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8"><title>miRNA visulization system for TSS regulatory region</title>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" type="text/css" href="css/home.css">

    <script src="js/echarts.min.js"></script>
    </head>

    <script language="javascript">
    function executeScript(html)
        {
            var reg = /<script[^>]*>([^\x00]+)$/i;
              //对整段HTML片段按<\/script>拆分
                var htmlBlock = html.split("<\/script>");
                  for (var i in htmlBlock) 
                    {
                        var blocks;//匹配正则表达式的内容数组，blocks[1]就是真正的一段脚本内容，因为前面reg定义我们用了括号进行了捕获分组
                          if (blocks = htmlBlock[i].match(reg)) 
                            {
                                //清除可能存在的注释标记，对于注释结尾-->可以忽略处理，eval一样能正常工作
                                  var code = blocks[1].replace(/<!--/, '');
                                    try 
                                      {
                                          eval(code) //执行脚本
                                        } 
                                          catch (e) 
                                            {
                                            }
                                        }
                                    }
                                }

    function postToPic()
    {   
        //获取接受返回信息层   
        var pic = document.getElementById("pic");  
        //获取表单对象和用户信息值
        var f = document.getElementById("form-mir");     
        var mirnaName = document.getElementById("select-mirna");
        //接收表单的URL地址
        var url = "/cgi-bin/getmirna.py";
        //需要POST的值，把每个变量都通过&来联接    
        var postStr   = "mirna_name="+mirnaName.value;
        //实例化Ajax
        //var ajax = InitAjax();
        var ajax = false;
        //开始初始化XMLHttpRequest对象
        if(window.XMLHttpRequest) 
        {     //Mozilla 浏览器
            ajax = new XMLHttpRequest();
            if (ajax.overrideMimeType) 
            {    //设置MiME类别
                ajax.overrideMimeType("text/xml");
            }
        }
        else if (window.ActiveXObject) 
        {     // IE浏览器
            try 
            {
                ajax = new ActiveXObject("Msxml2.XMLHTTP");
            } 
            catch (e) 
            {
                try 
                {
                    ajax = new ActiveXObject("Microsoft.XMLHTTP");
                } 
                catch (e) {}
             }
        }
        if (!ajax) 
        {     // 异常，创建对象实例失败
            window.alert("不能创建XMLHttpRequest对象实例.");
            return false;
        }

        //通过Post方式打开连接
        ajax.open("POST", url, true);
        //定义传输的文件HTTP头信息
        ajax.setRequestHeader("Content-Type","application/x-www-form-urlencoded");
        //发送POST数据
        ajax.send(postStr);
        //获取执行状态
        ajax.onreadystatechange = function() 
        { 
               //如果执行状态成功，那么就把返回信息写到指定的层里
               if (ajax.readyState == 4 && ajax.status == 200) 
            { 
                //改一下
                pic.innerHTML = ajax.responseText;
                //pic.innerHTML = ajax.return; 
                executeScript(ajax.responseText);
               } 
        } 
    }
    </script>

    <body>
        <div class="section-01">
            <h3>Welcome to visulization system<br>You can choose one miRNA in the left form and visualize the characteristics of its TSS regulatory region</h3>
        </div>
        <div>
            <div class="row-2">
                <div class="w-col-3">
                    <div class="form-block">
                        <form class="form-mirna" name="form-mir" action="">
                            <div><label for="label-mirna"><h3>Choose a miRNA:</h3></label></div>
                            <div><select class="select-field-mirna" id="select-mirna" multiple="multiple" name="select-mir" required="required">
                                <option value= >Select one...</option>
    '''         
    )                   
    for mirna_name in sort_mirna_list:
        visu_html.write('<option value="%s">%s</option>'%(mirna_name,mirna_name))
    visu_html.write(                                                        
    '''                            
                            </select></div>         
                            <div><input type="button" class="w-button" value="Submit" onClick="postToPic()"></div>
                        </form>
                    </div>
                </div>
                <div class="w-col-9">
                    <div id="pic">
                    <img src="/image/smile0.gif"/>                
                    </div>
                </div>
                            
            </div>
        </div>
    </body>

    <script type="text/javascript" src="js/jquery.min.js"></script>
    <!--[if lte IE 9]><script src="https://cdnjs.cloudflare.com/ajax/libs/placeholders/3.0.2/placeholders.min.js"></script><![endif]-->
    </body>
    </html>
    '''
    )

    visu_html.close()
    port=8080
    httpd=HTTPServer(('',port),CGIHTTPRequestHandler)
    print 'The visual system has been started. You can access \033[0;31m%s\033[0m on your browser'%'127.0.0.1:8080/visulization.html'
    httpd.serve_forever()

if '__name__' == '__main__':
    main()
