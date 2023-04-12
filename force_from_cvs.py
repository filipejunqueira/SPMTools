import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
import csv
from scipy.signal import savgol_filter

#sns.set()  # Setting seaborn as default style even if use only matplotlib
#sns.set(style="ticks", context="talk")
#sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})


# Loading the data from the csv file and then it is calculating the force for each spectra.
folder_path = "/media/captainbroccoli/DATA/2022-08-12/"
prefix = "20220812-124820_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_1.Df(Z)_mtrx"
csv_list = load_on_off_spec_list_from_cvs(cvs_name="Cu3spec_on_off_list")

path_on = []
path_off = []
results = pd.DataFrame()


number_of_curves =len(csv_list)
print(f"Number of curves: {number_of_curves}")

counter_off_not_match = 0

for idx,item in enumerate(csv_list):

    path_on = f"{folder_path}{prefix}{item[0]}{sufix}"
    path_off = f"{folder_path}{prefix}{item[1]}{sufix}"


    curve_ON_trace, curve_ON_retrace = import_spectra(path_on)
    curve_OFF_trace, curve_OFF_retrace = import_spectra(path_off)

    print(f"Parsed: {prefix}{item[0]}{sufix}")
    print(f"Parsed: {prefix}{item[1]}{sufix}")

    # THESE ARE THE 2 ON data frames in which the force will be calculated
    dfON_trace = curve_ON_trace.data_framedf
    dfON_retrace = curve_ON_retrace.data_framedf

    # OFF trace and retrace should be the same. Therefore, we will take the average between them.

    dfOFF_trace = curve_OFF_trace.data_framedf
    dfOFF_retrace = curve_OFF_retrace.data_framedf
    dfOFF_temp = dfOFF_trace["deltaF"].add(dfOFF_retrace["deltaF"]) / 2
    dfOFF = pd.DataFrame({'deltaF': dfOFF_temp, 'Z': dfOFF_trace['Z']})

    z_off_check, z_on_check = dfOFF_trace['Z'].to_numpy(), dfON_trace['Z'].to_numpy()
    check = (np.round(z_off_check,10) == np.round(z_on_check,10)).all()

    if check == True:
        pass

    elif check == False:
        counter_off_not_match = counter_off_not_match+1
        #warning
        print(f"{item[0]} and {item[1]} did not match, I found these ranges: {np.round(np.min(z_on_check)*10E9,4)}nm to {np.round(np.max(z_on_check)*10E9,4)}nm vs {np.round(np.min(z_off_check)*10E9,4)}nm to {np.round(np.max(z_off_check)*10E9,4)}nm")
        print(f"Will have to interpolate the OFF curve to match the ON curve")

        #sends the data frame df_OFF which has been averaged containg Z and deltaF to the fit_lennard_jones function.
        #Returns a lambda function
        lj_function = fit_lennard_jones(dfOFF, simple=True)



        new_Z_off = np.zeros(len())
        for z in df_ON_trace['Z']:
            df

        dfOFF = pd.DataFrame({'deltaF': new_df_off, 'Z': new_Z_off})



    # Data frames of the difference of frequency shifts
    df_diff_trace_temp = dfON_trace["deltaF"].subtract(dfOFF["deltaF"])
    df_diff_retrace_temp = dfON_retrace["deltaF"].subtract(dfOFF["deltaF"])


    df_diff_trace = pd.DataFrame({'deltaF': df_diff_trace_temp, 'Z': dfON_trace['Z']})
    df_diff_retrace = pd.DataFrame({'deltaF': df_diff_retrace_temp, 'Z': dfON_retrace['Z']})


    # smooth out curves using savitzky_golay filter.

    '''
    It uses least squares to regress a small window of your data onto a polynomial, 
    then uses the polynomial to estimate the point in the center of the window. 
    Finally the window is shifted forward by one data point and the process repeats. 
    This continues until every point has been optimally adjusted relative to its neighbors. 
    It works great even with noisy samples from non-periodic and non-linear sources.
    '''
    # window parameters are the curve, window size 5, polynomial order 3
    filter_order=3
    filter_window =5
    df_diff_trace_filtered_temp = savgol_filter(df_diff_trace_temp, filter_window, filter_order)
    df_diff_retrace_filtered = savgol_filter(df_diff_retrace_temp, filter_window, filter_order)

    df_diff_trace_filtered = pd.DataFrame({'deltaF': df_diff_trace_filtered_temp, 'Z': dfON_trace['Z']})
    df_diff_retrace_filtered = pd.DataFrame({'deltaF': df_diff_retrace_filtered, 'Z': dfON_retrace['Z']})

    # Parameters for initializing the sjarvis devoconvolution

    A = 0.52E-9
    f0 = 24234
    k = 1800

    # force calculation -> ON (2x) OFF(1 averaged) + difference (2x) = 5 forces in total.

    force_ON_trace, z_on, _, _, _ = sjarvis_deconvolution(dfON_trace, A=A, f0=f0, k=k)
    force_ON_retrace, _, _, _, _ = sjarvis_deconvolution(dfON_retrace, A=A, f0=f0, k=k)
    force_OFF, z_off, _, _, _ = sjarvis_deconvolution(dfOFF, A=A, f0=f0, k=k)

    force_diff_trace, _, _,_,_= sjarvis_deconvolution(df_diff_trace, A=A, f0=f0, k=k)
    force_diff_retrace, _, _,_,_= sjarvis_deconvolution(df_diff_retrace, A=A, f0=f0, k=k)

    force_diff_trace_filtered, _, _,_,_= sjarvis_deconvolution(df_diff_trace_filtered, A=A, f0=f0, k=k)
    force_diff_retrace_filtered, _, _,_,_= sjarvis_deconvolution(df_diff_retrace_filtered, A=A, f0=f0, k=k)


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

    fontsize = 24

    plot_df(dfON_trace.deltaF,dfON_retrace.deltaF,dfOFF.deltaF, dfON_retrace.Z, name=f"dfVsZ{item[0]}_ON_{item[1]}_OFF", save=True, retrace=True, fontsize=fontsize)
    #plot_forces_direct(force_ON_trace,force_ON_retrace,force_OFF,z_on, name=f"Force_Full_VsZ{item[0]}_ON_{item[1]}_OFF", save=True, retrace=False, fontsize=fontsize)
    plot_forces_short_range(force_diff_trace,force_diff_retrace,z_on, name=f"Force_short_range_VsZ{item[0]}_ON_{item[1]}_OFF", save=True, retrace=False, fontsize=fontsize)
    plot_forces_short_range(force_diff_trace_filtered,force_diff_retrace_filtered,z_on, name=f"Force_short_range_VsZ{item[0]}_ON_{item[1]}_OFF_filtered_order{filter_order}", save=True, retrace=False, fontsize=fontsize)
    plot_forces_and_df(force_diff_trace,force_diff_retrace,dfON_trace["deltaF"],dfON_retrace["deltaF"], dfOFF["deltaF"],z_force,z_df,name=f"Force_and_dfVsZ{item[0]}_ON_{item[1]}_OFF", save=True, retrace=False, fontsize=fontsize)
    plot_forces_and_df(force_diff_trace_filtered,force_diff_retrace_filtered,dfON_trace["deltaF"],dfON_retrace["deltaF"], dfOFF["deltaF"],z_force,z_df,name=f"Force_and_dfVsZ_{item[0]}_ON_{item[1]}_OFF_filtered_order{filter_order}", save=True, retrace=False, fontsize=fontsize)



print(f"Out of {number_of_curves} pair of curves, {counter_off_not_match} did not match")
print("Done, I owe you nothing!")
