import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
import csv
sns.set()  # Setting seaborn as default style even if use only matplotlib

# Loading the data from the csv file and then it is calculating the force for each spectra.
folder_path = "/home/captainbroccoli/Documents/2022-06-30/"
prefix = "20220630-152823_Cu(111)--AFM_NonContact_QPlus_AtomManipulation--"
sufix = ".Df(Z)_mtrx"
csv_list = load_spec_list_from_cvs()

path_on = []
path_off = []
results = pd.DataFrame()

for item in csv_list:

    path_on = f"{folder_path}{prefix}{item[0]}{sufix}"
    path_off = f"{folder_path}{prefix}{item[1]}{sufix}"

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

    force_ON_trace, z_on, _, _, _ = sjarvis_deconvolution(dfON_trace, A=0.01E-9, f0=-25000, k=1800)
    force_ON_retrace, z_on, _, _, _ = sjarvis_deconvolution(dfON_retrace, A=0.01E-9, f0=-25000, k=1800)
    force_OFF, z_off, _, _, _ = sjarvis_deconvolution(dfOFF, A=0.01E-9, f0=-25000, k=1800)

    results.insert(0, f"{item[0]}Force_ON_trace", force_ON_trace, True)
    results.insert(1, f"{item[0]}Force_ON_retrace", force_ON_retrace, True)
    results.insert(2, f"{item[1]}Force_OFF", force_OFF, True)
    results.insert(3, f"{item[0]}Z", z_on, True)

    plot_forces_short_range(force_ON_trace,force_ON_retrace,force_OFF,z_on, name=f"ForceVsZ{item[0]}_ON_{item[1]}_OFF", save=True)



# for line in dflist:
#
#     path_on = line[0]
#     path_off = line[1]
#
#     cu_ON = import_spectra(path_on)
#     cd_OFF = import_spectra(path_off)
#
#
#     force_ON_trace = sjarvis_deconvolution(cu_ON.trace)
#     force_ON_retrace = sjarvis_deconvolution(cu_ON.retrace)
#
#     force_OFF_trace = sjarvis_deconvolution(cd_OFF.trace)
#     force_OFF_retrace = sjarvis_deconvolution(cd_OFF.retrace)
#
#     store_dataframe(results,force_ON_trace,force_ON_retrace,force_OFF_trace,force_OFF_retrace)
# ,file_prefix="20220630-152823_Cu(111)--AFM_NonContact_QPlus_AtomManipulation--",file_extension=".Df(Z)_mtrx"
