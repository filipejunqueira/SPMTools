import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
from scipy.signal import savgol_filter

import csv
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")
sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})


folder_path = "/media/captainbroccoli/DATA/2022-07-29/"
prefix = "20220729-095353_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
type = "Aux2(V)"
sufix = "_mtrx"
file_number=["9_1"]

x = np.zeros(512)
y = np.zeros(512)
for number in file_number:

    path = f"{folder_path}{prefix}{number}.{type}{sufix}"

    mtrx_data = access2thematrix.MtrxData()
    data_file = f'{path}'
    traces, message = mtrx_data.open(data_file)
    curve_trace, message = mtrx_data.select_curve(traces[0])
    x = Spec_curve(curve_trace).X
    y = y + Spec_curve(curve_trace).Y

y = y/4


filter_order = 3
filter_window = 5
#y = savgol_filter(y, filter_window, filter_order)

figure, axis = plt.subplots(1,1,figsize=(14, 7), sharex=True)
sns.lineplot(ax=axis, x=x, y=y, color="#3386FF", alpha=1)

if type == "Aux2(V)":
    title = f"dI/dV"
    name_x = "V"
    name_y = "dI/dV"
    unit_x = "[V]"
    unit_y = "[1/Ohm]"

elif type =="Df(V)":
    title = f"Df(V)"
    name_x = "Bias"
    name_y = "Df(V)"
    unit_x = "[V]"
    unit_y = "[Hz]"

axis.set_title(f"{title}")
axis.set_xlabel(f"{name_x}{unit_x}")
axis.set_ylabel(f"{name_y}{unit_y}")
#axis.legend(bbox_to_anchor = (1.01, 1), loc = 'upper left')

pwd = os.getcwd()
name= f"{type}_{file_number}"
dir = "curves"
if os. path.isdir(f"{pwd}/{dir}") == False:
    os.mkdir(f"{pwd}/{dir}")
else:
    pass
plt.savefig(fname=f"{pwd}/{dir}/{name}",formatstr='.eps',facecolor='auto', edgecolor='auto')

print("I owe you nothing")