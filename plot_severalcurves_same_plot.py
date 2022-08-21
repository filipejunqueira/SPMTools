import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
from matplotlib import cm
import csv
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")
sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})

#folder_path = "/media/captainbroccoli/DATA/2022-07-17/"
folder_path = "/media/captainbroccoli/DATA/2022-08-07/"

#prefix = "20220717-163110_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
prefix = "20220807-174532_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"

sufix = ".Df(Z)_mtrx"
csv_list = load_spec_list_from_cvs()
csv_list_temp = []

for file,parameter in csv_list:
    csv_list_temp.append((file,float(parameter)))

#[csv_list_temp.append((item[i][0],int(item[i][1]))) for i,item in enumerate(csv_list)]

sorted_list = sorted(csv_list_temp, key=lambda t: t[1])
palette = sns.color_palette(['blue'], len(np.unique(np.asarray(sorted_list))))

figure, axis = plt.subplots(1,1,figsize=(14, 7), sharex=True)

#color_trace = 3554368 # black/grey to blue
#color_retrace = 16758272 # orange to pink


for idx,item in enumerate(sorted_list):

    path = f"{folder_path}{prefix}{item[0]}{sufix}"
    print(path)
    curve_trace, curve_retrace = import_spectra(path)

    x_trace = curve_trace.X
    y_trace = curve_trace.Y

    x_retrace = curve_retrace.X
    y_retrace = curve_retrace.Y

    #color_trace_int = hex(color_trace)[2:]
    #color_retrace_int = hex(color_retrace)[2:]


    #sns.lineplot(ax=axis, x=x_trace*10**9, y=y_trace,color=f"#{color_trace_int}",alpha=1)
    #sns.lineplot(ax=axis, x=x_retrace*10**9, y=y_retrace,color=f"#{color_retrace_int}",alpha=1)

    sns.lineplot(ax=axis, x=x_trace * 10 ** 9, y=y_trace, color="green", alpha=1)
    sns.lineplot(ax=axis, x=x_retrace * 10 ** 9, y=y_retrace, color=f"yellow", alpha=1)

    #color_trace = color_trace + 12
    #color_retrace = color_retrace +12

    #sns.lineplot(ax=axis, x=x_trace*10**9, y=y_trace)
    #sns.lineplot(ax=axis, x=x_retrace*10**9, y=y_retrace)


title = f"df(Z) - from {sorted_list[0][1]}mV to {sorted_list[-1][1]}mV"
name_x = "Z"
name_y = "df(Z)"
unit_x = "[nm]"
unit_y = "[Hz]"
axis.set_title(f"{title}")
axis.set_xlabel(f"{name_x}{unit_x}")
axis.set_ylabel(f"{name_y}{unit_y}")
#axis.legend(bbox_to_anchor = (1.01, 1), loc = 'upper left')

pwd = os.getcwd()

name=f"{csv_list[0][0]}to{csv_list[-1][0]}"
dir = "curves"
if os. path.isdir(f"{pwd}/{dir}") == False:
    os.mkdir(f"{pwd}/{dir}")
else:
    pass
plt.savefig(fname=f"{pwd}/{dir}/{name}",formatstr='.eps',facecolor='auto', edgecolor='auto')

print("I owe you nothing")

