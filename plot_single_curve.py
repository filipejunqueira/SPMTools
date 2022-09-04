from spmUtils import *




#Layout of the graph:###################################################################################################
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")


# Hex color for graphs
#color_map_filipe = {"green": "#90EE90", "dark_green": "#95BA61" , "orange": "#FFAB00", "dark_yellow": "#8B8000", "yellow": "#FFDB58", "red": "#B60005", "blue": "#3c59ff", "white": "#FFFFFF.","purple": "#7b40c9", "pink": "#FFB490", "black": "#171717"}
#color_map_comp2 = ["green","blue","#FFDB58","#C500B2"]


n_component_colors = 2
color_components_initial = Color("red")
color_components_final = Color('green')
color_map_comp = [color.hex_l for color in list(color_components_initial.range_to(color_components_final,n_component_colors))]


color_map = {"data": "#7b40c9","data_retrace": "#90EE90", "best_fit": "#171717", "best_fit_retrace": "#B60005", "color_map_comp": color_map_comp}


#Folder Structure#######################################################################################################
root_path = "/media/filipejunqueira/DATA/"
project_folder_name = "2022-08-07"
prefix = "20220807-174532_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_mtrx"
project_folder_path = os.path.join(root_path,project_folder_name)

#INPUT##################################################################################################################
fontsize = 36
marker_size = 260
curve_type = "Aux2(V)"
figsize = (20,14)
file_id =27
n_files = 4
plot_retrace = False
slice_start = 0
slice_end = 500

# Filtering
filter_flag = False
filter_order = 3
filter_window = 3  # This needs to be an odd number

#fitting background
fit_type = "lorentzian"
fit_linear = True
fit_exponential = True
linear_background_fit_cutoff = 1/3
plot_bestfit = True
plot_components = True

#####################################################################


###### Script starts here ############################################

file_number_list = create_file_number_list(file_id, n_files)
print(f"Curves {file_number_list} of type: {curve_type} have been averaged ")

x, y = average_curves(file_number_list, curve_type, direction=0)


# Slicing the data to be plotted/filtered/fitted.
print(f"Attention: I'm slicing the graph from {slice_start} to {slice_end} points. Make sure this is what you want.")
x = x[slice_start:slice_end]
y = y[slice_start:slice_end]
print(f"X size{len(x)}")
print(f"Y size{len(y)}")

# Filtering the data
if filter_flag is True:
    y = savgol_filter(y, filter_window, filter_order)
    filtered_str = f"_filter_w{filter_window}o{filter_order}"
else:
    filtered_str = ""

# Fitting the data to a specified curve.

if plot_bestfit is True:
    results, background, components, init = quick_fit(x, y, linear = fit_linear, exponential = fit_exponential, fit_type=fit_type, n_curves= 1, lin_cutoff=linear_background_fit_cutoff)
    y_bestfit = results.best_fit
    print(f"Y_best size{len(y_bestfit)}")
    print(results.fit_report())

else:
    y_bestfit = None
    results = None
    pass

if plot_components is True:
    plot_single_curve(x, y, plot_retrace=plot_retrace, curve_type=curve_type,
                      fontsize=fontsize, marker_size=marker_size,
                      figsize=figsize, y_fit=y_bestfit, results_object=results,
                      color_map=color_map)
else:
    plot_single_curve(x, y, plot_retrace=plot_retrace, curve_type=curve_type,
                      fontsize=fontsize, marker_size=marker_size,
                      figsize=figsize, y_fit=y_bestfit,
                      color_map=color_map)

# Building path to save the image:
pwd = os.getcwd()

name = f"{curve_type}_{project_folder_name}--{file_id}_1-{n_files}{filtered_str}"
save_dir_name = "single_curve_graphs"
save_dir = os.path.join(pwd, save_dir_name)
# Checks if the folder exists
if os.path.isdir(save_dir) == False:
    os.mkdir(save_dir)
else:
    pass

image_name = os.path.join(pwd, save_dir_name, name)
# Finally saves the image
plt.savefig(fname=f"{image_name}.svg", facecolor='auto', edgecolor='auto', transparent=True)
plt.show()

print("I owe you nothing")