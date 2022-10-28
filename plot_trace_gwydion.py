from spmUtils import *

#Layout of the graph:###################################################################################################
sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")


# First step is to locate file. One method is to manually select the spectroscopy txt .file.
# To do this one must use the following code:
# Depending how it is exported from gwideon;

#file_path = get_path_gui()
file_path = "/media/captainbroccoli/DATA/thesis_figs/traces/583_height_difference.txt"
root_path = file_path.split("/traces/")[0]

trace = np.genfromtxt(file_path, skip_header=3)
units = get_one_line(file_path, 3).decode("utf-8").strip().split("    ")
legend = get_one_line(file_path,1).decode("utf-8").strip().split("             ")
number_traces = int(len(trace[0, :]) / 2)
fontsize = 40
figsize = (24,14)
marker_size=200
figure, axis = plt.subplots(1, 1, figsize=figsize, sharex=True)

plot_title = file_path.split("/traces/")[1].strip(".txt")
colours_palette = ["#00C9A7","#D13A28", "#0057D9", "#5B6C5D","#7b40c9", "#FFC145",]
#plt.title(plot_title)

index = 0

while index < number_traces:

    x = np.array(trace[:, 2 * index])*1000_000_000_000
    y = np.array(trace[:, 2 * index + 1])*1000_000_000_000
    y = y - np.nanmin(y)
    sns.scatterplot(ax=axis, x=x, y=y, alpha=1, edgecolor=colours_palette[index%5], facecolor="None", s=marker_size, label=legend[index])
    sns.lineplot(ax=axis, x=x, y=y, alpha=1, color = colours_palette[index%5])
    #plt.plot(trace[:, 2 * index], trace[:, 2 * index + 1], colours_palette[index%5], label=legend[index], )
    print(index%5)
    index += 1

name_x = "Trace Profile"
name_y = "Z height"


#axis.set_title(f"{plot_title}", fontsize=fontsize)
#axis.set_xlabel(f"{name_x} {units[0]}", fontsize=(fontsize - 2))
#axis.set_ylabel(f"{name_y} {units[1]}", fontsize=(fontsize - 2))
axis.set_xlabel(f"{name_x} [pm]", fontsize=(fontsize - 2))
axis.set_ylabel(f"{name_y} [pm]", fontsize=(fontsize - 2))

plt.xticks(fontsize=(fontsize - 2))
plt.yticks(fontsize=(fontsize - 2), )
plt.grid(color = 'grey', linestyle = '--', linewidth = 0.5, alpha = 0.4)

plt.legend(fontsize=(fontsize - 10))
plt.savefig(f"{root_path}/traces/{plot_title}.png", bbox_inches="tight")

plt.show()
