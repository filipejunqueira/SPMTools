from spmUtils import *
#Layout of the graph:###################################################################################################
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")


# Hex color for graphs
n_component_colors = 3
color_components_initial = Color("red")
color_components_final = Color('green')
color_map_comp = [color.hex_l for color in list(color_components_initial.range_to(color_components_final,n_component_colors))]
color_map = {"data": "#7b40c9","data_retrace": "#90EE90", "best_fit": "#171717", "best_fit_retrace": "#B60005", "color_map_comp": color_map_comp}
color_map2 = {"data": "#90EE90","data_retrace": "#90EE90", "best_fit": "#171717", "best_fit_retrace": "#B60005", "color_map_comp": color_map_comp}

#Folder Structure#######################################################################################################
root_path = "/media/filipejunqueira/DATA/"
project_folder_name = "2022-08-07"
prefix = "20220807-174532_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_mtrx"
project_folder_path = os.path.join(root_path,project_folder_name)
prefix_full_path = os.path.join(project_folder_path,prefix)

#INPUT##################################################################################################################
#Graph
fontsize = 36
marker_size = 260
curve_type = "Aux2(V)"
figsize = (20,14)

# File
file_id =27
n_files = 4
plot_retrace_flag = False
slice_start = 0
slice_end = 500

# Filtering
filter_flag = False
filter_order = 3
filter_window = 3  # This needs to be an odd number


# Fitting
n_curves= 1
fit_type = "lorentzian"
algo = "leastsq"


#fitting background
fit_linear = True
fit_exponential = True
linear_background_fit_cutoff = 1/3
plot_bestfit = True
plot_components = False

# Saving
save_dir_name = "single_curve_graphs"


#####################################################################


###### Script starts here ############################################

file_number_list = create_file_number_list(file_id, n_files)
print(f"Curves {file_number_list} of type: {curve_type} have been averaged ")

if plot_retrace_flag is False:
    x, y = average_curves(prefix_full_path,file_number_list, curve_type, direction=0)
    print(f"Attention: I'm slicing the graph from {slice_start} to {slice_end} points. Make sure this is what you want.")
    x = x[slice_start:slice_end]
    y = y[slice_start:slice_end]
    y_retrace = None
else:
    pass #THIS NEEDS TO BE CHANGED> RIGHT NOW THE AVERAGING IS INSIDE THE PLOT THIS IS BADDLY WRITTEN.

# Slicing the data to be plotted/filtered/fitted.

# Filtering the data
if filter_flag is True:
    y = savgol_filter(y, filter_window, filter_order)
    filtered_str = f"_filter_w{filter_window}o{filter_order}"
    if plot_retrace_flag is True:
        y_retrace = savgol_filter(y_retrace, filter_window, filter_order)
else:
    filtered_str = ""

# Fitting the data to a specified curve.

if plot_bestfit is True:
    results, background, components, init = quick_fit(x, y, linear = fit_linear, exponential = fit_exponential, fit_type=fit_type, lin_cutoff=linear_background_fit_cutoff, n_curves=n_curves, algo=algo)
    y_bestfit = results.best_fit
    #print(results.fit_report())

else:
    y_bestfit = None
    results = None
    pass
#neeed to initiate figure and axis
figure, axis = plt.subplots(1, 1, figsize=figsize, sharex=True)

if plot_components is True:
    plot_single_curve(x, y, axis=axis, y_retrace=y_retrace, curve_type=curve_type,
                      fontsize=fontsize, marker_size=marker_size, y_fit=y_bestfit, results_object=results,
                      color_map=color_map)

else:
    plot_single_curve(x, y, axis=axis, y_retrace=y_retrace, curve_type=curve_type,
                      fontsize=fontsize, marker_size=marker_size, y_fit=y_bestfit,
                      color_map=color_map)

# Building path to save the image:
path_dic = {'curve_type': curve_type, 'special_string': project_folder_name, 'file_id': file_id, 'n_files': n_files ,'filtered_str': filtered_str }
save_plot(path_dic,save_dir_name)
plt.show()

print("I owe you nothing")