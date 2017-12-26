import os
import sys
from setuptools import setup,find_packages

def main():
    # judge the version of python 
    if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:        
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)
    setup(name="DeepMiRTss",
          version="0.1",
          description="miRNA TSS Analysis By Deep-learning method",
          author='Wang Zhi',
          author_email='szxszx@foxmail.com',
          url='https://github.com/jwanglabwangzhi/DeepMiRTss',
          packages=['deepmirtss'],
          include_package_data=True,
          entry_points = {
            'console_scripts':[           
                'DeepMiRTss.tssfinder=deepmirtss.DeepMiRTss_tss_finder:main',
                'DeepMiRTss.analysis=deepmirtss.DeepMiRTss_analysis:main',
                'DeepMiRTss.visualization=deepmirtss.DeepMiRTss_visulization:main'
],          
}
)

if __name__ == '__main__':
    main()
