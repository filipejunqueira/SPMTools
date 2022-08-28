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
#sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})


folder_path = "/media/filipejunqueira/DATA/2022-08-12/"
prefix = "20220812-124820_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
type = "Aux2(V)"
sufix = "_mtrx"
file_number=["37_1"]

x = np.zeros(512)
y = np.zeros(512)

for number in file_number:

    path = f"{folder_path}{prefix}{number}.{type}{sufix}"
    print(path)

    mtrx_data = access2thematrix.MtrxData()
    data_file = f'{path}'
    traces, message = mtrx_data.open(data_file)
    curve_trace, message = mtrx_data.select_curve(traces[0])
    x = Spec_curve(curve_trace).X
    y = y + Spec_curve(curve_trace).Y

y = y/len(file_number)

filter_order = 3
filter_window = 5
#y = savgol_filter(y, filter_window, filter_order)

figure, axis = plt.subplots(1,1,figsize=(14, 7), sharex=True)

if type == "Aux2(V)":
    sns.lineplot(ax=axis, x=x, y=y, color="#3386FF", alpha=1)
    title = f"dI/dV"
    name_x = "V"
    name_y = "dI/dV"
    unit_x = "[V]"
    unit_y = "[1/Ohm]"

elif type =="Df(V)":
    sns.lineplot(ax=axis, x=x, y=y, color="#3386FF", alpha=1)
    title = f"Df(V)"
    name_x = "Bias"
    name_y = "Df(V)"
    unit_x = "[V]"
    unit_y = "[Hz]"

    axis.set_title(f"{title}")
    axis.set_xlabel(f"{name_x}{unit_x}")
    axis.set_ylabel(f"{name_y}{unit_y}")
    #axis.legend(bbox_to_anchor = (1.01, 1), loc = 'upper left')

if type=="Df(Z)":

    path = f"{folder_path}{prefix}{number}.{type}{sufix}"
    print(path)

    mtrx_data = access2thematrix.MtrxData()
    data_file = f'{path}'
    traces, message = mtrx_data.open(data_file)
    curve_trace, message = mtrx_data.select_curve(traces[0])
    curve_retrace, message = mtrx_data.select_curve(traces[1])

    df_ON_trace = Spec_curve(curve_trace).data_framedf
    df_ON_retrace = Spec_curve(curve_retrace).data_framedf
    z = Spec_curve(curve_trace).X

    sns.lineplot(ax=axis, x=z * 10 ** 9, y=df_ON_trace['deltaF'], color=color_map["green"], label="df ON trace")
    sns.lineplot(ax=axis, x=z * 10 ** 9, y=df_ON_retrace['deltaF'], color=color_map["yellow"], label="df ON retrace")
    axis.set_title("df(Z)")
    axis.set_xlabel("Z[nm]")
    axis.set_ylabel("df[Hz]")
    axis.legend(loc=0)

pwd = os.getcwd()
name= f"{type}_{file_number}"
dir = "curves"
if os. path.isdir(f"{pwd}/{dir}") == False:
    os.mkdir(f"{pwd}/{dir}")
else:
    pass
plt.savefig(fname=f"{pwd}/{dir}/{name}",formatstr='.eps',facecolor='auto', edgecolor='auto')

print("I owe you nothing")