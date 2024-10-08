from distutils.core import setup

setup(
  name = 'bmdrc',        
  packages = ['bmdrc'],   
  version = '0.0.1.1',      
  description = 'A package for fitting benchmark dose curves to dichotomous data and light photomotor response data',   # Give a short description about your library
  author = ['David Degnan', 'Sara Gosline', 'Katrina Waters'],               
  author_email = ['David.Degnan@pnnl.gov', 'Sara.Gosline@pnnl.gov', 'Katrina.Waters@pnnl.gov'],     
  url = 'https://github.com/PNNL-CompBio/bmdrc', 
  download_url = 'https://github.com/PNNL-CompBio/bmdrc/archive/refs/tags/v_0.0.1.tar.gz',
  install_requires=[           
    "pandas>=2.1.4",
    "numpy>=1.26.2",
    "scipy>=1.11.4",
    "astropy>=6.0.0",
    "argparse",
    "statsmodels>=0.14.1",
    "matplotlib>=3.8.2"  
  ],
)