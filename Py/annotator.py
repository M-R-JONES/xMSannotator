import pandas as pd
import numpy as np
import scipy as sp
import matplotlib
import networkx
import csv

class variables:
    
    #define default processing parameters for clustering
    def __init__(self):
        self.max_mz_diff = 10
        self.max_rt_diff = 10
        self.cormethod = 'pearson'
        self.clustermethod = 'WGCNA'

        
vars = variables()


#run script
data = pd.read_csv('C:/Users/jonesmar/Documents/Git_xMSannotator/xMSannotator/modeldata.csv')



