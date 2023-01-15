from spmUtils import *
#Layout of the graph:###################################################################################################
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")


# Hex color for graphs
n_component_colors = 3
color_components_initial = Color("#FF9671")
color_components_final = Color('#2C73D2')
color_map_comp = [color.hex_l for color in list(color_components_initial.range_to(color_components_final,n_component_colors))]

color_map_type1 = {"data": "#00C9A7","data_retrace": "#FFC75F", "best_fit": "#171717", "best_fit_retrace": "#D65DB1", "color_map_comp": color_map_comp}
color_map_type2 = {"data": "#D13A28","data_retrace": "#F9F871", "best_fit": "#171717", "best_fit_retrace": "#FF6F91", "color_map_comp": color_map_comp}
color_map_type3 = {"data": "#00C9A7","data_retrace": "#F9F871", "best_fit": "#171717", "best_fit_retrace": "#FF6F91", "color_map_comp": color_map_comp}
color_map_type4 = {"data": "#845EC2","data_retrace": "#FFC75F", "best_fit": "#171717", "best_fit_retrace": "#D65DB1", "color_map_comp": color_map_comp}


generic_gradient = ["#845EC2","#D65DB1","#FF6F91","#FF9671","#FFC75F","#F9F871" ]
matching_gradient = ["#845EC2", "#2C73D2", "#0081CF","#0089BA", "#008E9B", "#008F7A"]
spot_pallet = ["#845EC2","#B39CD0","#FBEAFF","#00C9A7"]
twisted_spot_pallet = ["#845EC2","#00C9A7", "#C4FCEF", "#4D8076"]
collective_pallet = ["#845EC2","#BD38B2","#D13A28"]


#Folder Structure#######################################################################################################
root_path = "/media/captainbroccoli/DATA/"
project_folder_name = "2022-07-17"
prefix = "20220717-163110_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_mtrx"
project_folder_path = os.path.join(root_path,project_folder_name)
prefix_full_path = os.path.join(project_folder_path,prefix)
#INPUT##################################################################################################################
#Graph

marker_size = 260
figsize = (30,16)

# File
curve_type = "Df(V)"
file_id =5
n_files = 1
plot_retrace_flag = False
slice_start = 0
slice_end = 512

# Filtering
filter_flag = False
filter_order = 3
filter_window = 7  # This needs to be an odd number


# Fitting
n_curves= 1
fit_type = "voigt"
algo = "leastsq"
intial_guess = 3.2

#Fitting background
fit_exponential = False
fit_quadratic = False
plot_background = False
linear_background_fit_cutoff = 1
plot_bestfit = False
plot_components = False

#Fix Linear background fitting
fix_background_linear = False

# Saving
save_dir_name = "single_curve_graphs"


#####################################################################
###### Script starts here ############################################

#Color map to use:

match curve_type:
    case "Df(Z)":
        color_map = color_map_type1
        fontsize = 60

    case "Aux2(V)":
        color_map = color_map_type1
        fontsize = 60

    case "Df(V)":
        color_map = color_map_type2
        fontsize = 60

    case "I(Z)":
        color_map = color_map_type1
        fontsize = 60

    case "A(f)":
        color_map = color_map_type1
        fontsize = 60

    case "Phase(f)":
        color_map = color_map_type1
        fontsize = 60

##########################################################################

file_number_list = create_file_number_list(file_id, n_files)
print(f"Curves {file_number_list} of type: {curve_type} have been averaged ")

if plot_retrace_flag is False:
    x, y = average_curves(prefix_full_path,file_number_list, curve_type, direction=0)
    if curve_type == "Aux2(V)":
        y =  y - np.nanmin(y) # re-scale to 0

    print(f"Attention: I'm slicing the graph from {slice_start} to {slice_end} points. Make sure this is what you want.")
    x = x[slice_start:slice_end]
    y = y[slice_start:slice_end]
    y_retrace = None

else:
    x, y = average_curves(prefix_full_path, file_number_list, curve_type, direction=0)
    print(
        f"Attention: I'm slicing the graph from {slice_start} to {slice_end} points. Make sure this is what you want.")
    x = x[slice_start:slice_end]
    y = y[slice_start:slice_end]
    _, y_retrace = average_curves(prefix_full_path, file_number_list, curve_type, direction=1)

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
    results, background, components, init = quick_fit(x, y, linear = fix_background_linear, exponential = fit_exponential, quadratic=fit_quadratic, fit_type=fit_type, lin_cutoff=linear_background_fit_cutoff, n_curves=n_curves, algo=algo, guess= intial_guess)
    y_bestfit = results.best_fit
    print(results.fit_report())

else:
    y_bestfit = None
    results = None
    pass
#neeed to initiate figure and axis
figure, axis = plt.subplots(1, 1, figsize=figsize, sharex=True)

if plot_components is True:
    plot_single_curve(x, y, axis=axis, y_retrace=y_retrace, curve_type=curve_type,
                      fontsize=fontsize, marker_size=marker_size, y_fit=y_bestfit, results_object=results,
                      color_map=color_map, plot_background=plot_background)

else:
    plot_single_curve(x, y, axis=axis, y_retrace=y_retrace, curve_type=curve_type,
                      fontsize=fontsize, marker_size=marker_size, y_fit=y_bestfit,
                      color_map=color_map)

# Building path to save the image:
path_dic = {'curve_type': curve_type, 'special_string': project_folder_name, 'file_id': file_id, 'n_files': n_files ,'filtered_str': filtered_str }
save_plot(path_dic,save_dir_name)
plt.show()

print("I owe you nothing")