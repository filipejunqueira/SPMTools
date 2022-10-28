from spmUtils import *
#Layout of the graph:###################################################################################################
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")


# Hex color for graphs
#color_map_filipe = {"green": "#90EE90", "dark_green": "#95BA61" , "orange": "#FFAB00", "dark_yellow": "#8B8000", "yellow": "#FFDB58", "red": "#B60005", "blue": "#3c59ff", "white": "#FFFFFF.","purple": "#7b40c9", "pink": "#FFB490", "black": "#171717"}
#color_map_comp2 = ["green","blue","#FFDB58","#C500B2"]

n_component_colors = 10
color_components_initial = Color("red")
color_components_final = Color('yellow')

color_map_comp = [color.hex_l for color in list(color_components_initial.range_to(color_components_final,n_component_colors))]

color_map = {"data": "#90EE90","data_retrace": "#FFFF00", "best_fit": "#171717", "best_fit_retrace": "#B60005", "color_map_comp": color_map_comp}
color_map2 = {"data": "#90EE90","data_retrace": "#90EE90", "best_fit": "#171717", "best_fit_retrace": "#B60005", "color_map_comp": color_map_comp}



#Folder Structure#######################################################################################################
root_path = "/media/captainbroccoli/DATA/"
project_folder_name = "2022-07-17"
prefix = "20220717-163110_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_mtrx"
project_folder_path = os.path.join(root_path,project_folder_name)
prefix_full_path = os.path.join(project_folder_path,prefix)

#INPUT##################################################################################################################
#Graph
fontsize = 36
marker_size = 260
curve_type = "Df(Z)"
figsize = (20,14)

# File
plot_retrace_flag = True
slice_start = 0
slice_end = 500

# Filtering
filter_flag = False
filter_order = 3
filter_window = 3  # This needs to be an odd number

# Fitting
n_curves= 2
fit_type = "lorentzian"
algo = "leastsq"

#fitting background
fit_linear = True
fit_exponential = True
linear_background_fit_cutoff = 1/3
plot_bestfit = True
plot_components = True

# Saving
save_dir_name = "several_curve_graphs"
#####################################################################
###### Script starts here ############################################

csv_list = load_spec_list_from_cvs()
csv_list_temp = []

[csv_list_temp.append((file,int(number),float(parameter))) for file, number, parameter in csv_list]
sorted_list = sorted(csv_list_temp, key=lambda t: t[2])
print(sorted_list)

figure, axis = plt.subplots(1, 1, figsize=figsize, sharex=True)

#COLOR SCHEME for all the curves in the same graph::


color_trace_initial= Color("#013220")
color_trace_final= Color('#008080')

color_retrace_initial= Color("#ff7600")
color_retrace_final= Color("#eeff00")

color_map_trace = [color.hex_l for color in list(color_trace_initial.range_to(color_trace_final,len(sorted_list)))]
color_map_retrace = [color.hex_l for color in list(color_retrace_initial.range_to(color_retrace_final,len(sorted_list)))]

for i, item in enumerate(sorted_list):

    file_id, n_files, parameter = item

    # this modifies the color_map dictionary such that the colors of each curve are different. The gradient is automatic and determined by the size of the list.
    #only availble
    color_map = {"data": color_map_trace[i], "data_retrace": color_map_retrace[i], "best_fit": "#171717", "best_fit_retrace": "#B60005",
                 "color_map_comp": color_map_comp}

    file_number_list = create_file_number_list(file_id, n_files)
    x, y = average_curves(prefix_full_path, file_number_list, curve_type, direction=0)
    print(f"Curves {file_id} of type: {curve_type} have been averaged ")

    if plot_retrace_flag is False:
        plot_single_curve(x, y, axis=axis, curve_type=curve_type, fontsize=fontsize, marker_size=marker_size, color_map=color_map,
                          label_flag=False)
    else:
        _, y_retrace = average_curves(prefix_full_path,file_number_list, curve_type, direction=1)
        plot_single_curve(x, y, axis=axis, curve_type=curve_type, fontsize=fontsize, marker_size=marker_size,
                          y_retrace=y_retrace, color_map=color_map,label_flag=False)


filtered_str =""
# Building path to save the image:
# Building path to save the image:
path_dic = {'curve_type': curve_type, 'special_string': project_folder_name, 'file_id': f"{sorted_list[0][0]}-{sorted_list[-1][0]}", 'n_files': "" ,'filtered_str': filtered_str }
save_plot(path_dic,save_dir_name)
plt.show()

print("I owe you nothing")



