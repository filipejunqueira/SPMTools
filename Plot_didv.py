import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
from scipy.signal import savgol_filter

sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")

folder_path = "/media/filipejunqueira/DATA/2022-07-29/"
prefix = "20220729-095353_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
type = "Aux2(V)"
sufix = "_mtrx"
file_number=["3_1"]
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
x = x[:500]
y = y[:500]
figure, axis = plt.subplots(1,1,figsize=(12, 12), sharex=True)

sns.lineplot(ax=axis, x=x, y=y, color="#3386FF", alpha=1)
title = f"dI/dV"
name_x = "V"
name_y = "dI/dV"
unit_x = "[V]"

axis.set_title(f"{title}",fontsize=30)
axis.set_xlabel(f"{name_x}{unit_x}",fontsize=28)
plt.xticks(fontsize=28)
plt.yticks([])

#axis.set_ylabel(f"{name_y}{unit_y}")
pwd = os.getcwd()
name= f"{type}_{file_number}"
dir = "curves"
if os. path.isdir(f"{pwd}/{dir}") == False:
    os.mkdir(f"{pwd}/{dir}")
else:
    pass
plt.savefig(fname=f"{pwd}/{dir}/{name}",formatstr='.png',facecolor='auto', edgecolor='auto')