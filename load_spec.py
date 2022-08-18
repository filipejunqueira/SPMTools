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

#plt.style.use("dark_background")



# Loading the data from the csv file and then it is calculating the force for each spectra.
folder_path = "/media/captainbroccoli/DATA/2022-07-17/"
prefix = "20220717-163110_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = ".Df(Z)_mtrx"
csv_list = load_spec_list_from_cvs()

path_on = []
path_off = []
results = pd.DataFrame()

for idx,item in enumerate(csv_list):

    path_on = f"{folder_path}{prefix}{item[0]}{sufix}"
    path_off = f"{folder_path}{prefix}{item[1]}{sufix}"
    print(path_on)
    print(path_off)


    curve_ON_trace, curve_ON_retrace = import_spectra(path_on)
    curve_OFF_trace, curve_OFF_retrace = import_spectra(path_off)

    # THESE ARE THE 2 ON data frames in which the force will be calculated.
    dfON_trace = curve_ON_trace.data_frame
    dfON_retrace = curve_ON_retrace.data_frame

    # OFF trace and retrace should be the same. Therefore, we will take the average between them.

    dfOFF_trace = curve_OFF_trace.data_frame
    dfOFF_retrace = curve_OFF_retrace.data_frame
    dfOFF_temp = dfOFF_trace["deltaF"].add(dfOFF_retrace["deltaF"]) / 2
    dfOFF = pd.DataFrame({'deltaF': dfOFF_temp, 'Z': dfOFF_trace['Z']})

    # Data frames of the difference of frequency shifts

    df_diff_trace_temp = dfON_trace["deltaF"].subtract(dfOFF["deltaF"])
    df_diff_retrace_temp = dfON_retrace["deltaF"].subtract(dfOFF["deltaF"])

    df_diff_trace = pd.DataFrame({'deltaF': df_diff_trace_temp, 'Z': dfON_trace['Z']})
    df_diff_retrace = pd.DataFrame({'deltaF': df_diff_retrace_temp, 'Z': dfON_retrace['Z']})

    # Parameters for initializing the sjarvis devoconvolution

    A=0.01E-9
    f0 = 25000
    k = 1800

    # force calculation -> ON (2x) OFF(1 averaged) + difference (2x) = 5 forces in total.

    force_ON_trace, z_on, _, _, _ = sjarvis_deconvolution(dfON_trace, A=A, f0=f0, k=k)
    force_ON_retrace, _, _, _, _ = sjarvis_deconvolution(dfON_retrace, A=A, f0=f0, k=k)
    force_OFF, z_off, _, _, _ = sjarvis_deconvolution(dfOFF, A=A, f0=f0, k=k)

    force_diff_trace, _, _,_,_= sjarvis_deconvolution(df_diff_trace, A=A, f0=f0, k=k)
    force_diff_retrace, _, _,_,_= sjarvis_deconvolution(df_diff_retrace, A=A, f0=f0, k=k)

    # Result data frame receiving data.
    results.insert(idx, f"{item[0]}Force_ON_trace", force_ON_trace, True)
    results.insert(idx+1, f"{item[0]}Force_ON_retrace", force_ON_retrace, True)
    results.insert(idx+2, f"{item[0]}Z_ON", z_on, True)
    results.insert(idx+3, f"{item[1]}Force_OFF", force_OFF, True)
    results.insert(idx+4, f"{item[1]}Z_OFF", z_off, True)
    results.insert(idx+5, f"{item[0]}-{item[1]}Force_diff_trace", force_diff_trace, True)
    results.insert(idx+6, f"{item[0]}-{item[1]}Force_diff_retrace", force_diff_retrace, True)

    z_force = z_on
    z_df = df_diff_trace["Z"]
    plot_df(dfON_trace.deltaF,dfON_retrace.deltaF,dfOFF.deltaF, dfON_retrace.Z, name=f"dfVsZ{item[0]}_ON_{item[1]}_OFF", save=True)
    plot_forces_direct(force_ON_trace,force_ON_retrace,force_OFF,z_on, name=f"Force_Full_VsZ{item[0]}_ON_{item[1]}_OFF", save=True)
    plot_forces_short_range(force_diff_trace,force_diff_retrace,z_on, name=f"Force_short_range_VsZ{item[0]}_ON_{item[1]}_OFF", save=True)
    plot_forces_and_df(force_diff_trace,force_diff_retrace,dfON_trace["deltaF"],dfON_retrace["deltaF"], dfOFF["deltaF"],z_force,z_df,name=f"Force_and_dfVsZ{item[0]}_ON_{item[1]}_OFF", save=True)

print("Done, I owe you nothing!")