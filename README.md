This works best with anaconda python

Dependencies
------------
conda  
numpy  
scipy  
matplotlib  
[bs4](https://www.crummy.com/software/BeautifulSoup/)  
svg.path
vtk  
nibabel   
pandas  
[shapeio](https://github.com/bmapdev/shapeio)  
[curvematch](https://github.com/bmapdev/curvematch)  



Installation
------------
To install the dependencies, do

pip install numpy scipy matplotlib bs4 svg.path pandas nibabel  
conda install vtk -y  
pip install git+https://github.com/bmapdev/shapeio git+https://github.com/bmapdev/curvematch  
pip install git+https://github.com/bmapdev/ccshape.git

Execution
----------
For offline plotting on macos, you will have to add a backend for matplotlib
On the terminal do  
`echo "backend: TkAgg" >> ~/.matplotlib/matplotlibrc`   

To get help, type:  
`corpus_callosum_analyze.py -h`

The corpus callosum shape analysis can be run in two different ways.
1. Calculate thickness  
Assuming `toplist.txt` and `botlist.txt` contain the paths to the top and bottom curves, and `subjid.txt` contains the list of unique subject identifiers, the below command will match the top and bottom curves, compute callosal thickness and save them in the directory `ccout`.  
`corpus_callosum_analyze.py subjid.txt toplist.txt botlist.txt -odir ccout`

2. Registration to a template  
Assuming `toplist.txt` and `botlist.txt` contain the paths to the top and bottom curves, and `subjid.txt` contains the list of unique subject identifiers, the below command will register the top and bottom curves to the top and bottom curves of a template, and save the warped outputs in the directory `ccout_register`.  
`corpus_callosum_analyze.py subjid.txt toplist.txt botlist.txt -template curve_template -odir ~/Desktop/ccout_register`



