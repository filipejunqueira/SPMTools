import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
import csv
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")
sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})

# Loading the data from the csv file and then it is calculating the force for each spectra.
folder_path = "/media/captainbroccoli/DATA/2022-07-17/"
prefix = "20220717-163110_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = ".Df(V)_mtrx"
file_number="4_1"
path = f"{folder_path}{prefix}{file_number}{sufix}"
curve_list=[]

curve_trace, curve_retrace = import_spectra(path)
