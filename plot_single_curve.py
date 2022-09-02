import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
from scipy.signal import savgol_filter

import csv

# Layout of the graph:

sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")

fontsize = 30
marker_size = 160
#sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})


type = "Aux2(V)"
figsize = (20,14)
file_id = 37
n_files = 1
plot_retrace = False
end = 512
# If there is no retrace it will give out an error!!!

path = os.path.join(project_folder_path, f"{prefix}{file_id}_{n_files}.{type}{sufix}")
spec = import_spectra(path)
x = spec.X
y = spec.Y
x = x[:end]
y = y[:end]


result, background, components, init = quick_fit(x, y, linear = True, exponential = True, fit_type="voigt", n_curves= 1)
y_bestfit = result.best_fit


plot_single_curve(file_id, n_files, plot_retrace=plot_retrace, type =type, filter = False, fontsize=fontsize,marker_size=marker_size, figsize=figsize, slice_end=end, bestfit=True, y_fit=y_bestfit, results_object=result)

#print(result.fit_report())
result.params.pretty_print()




print("I owe you nothing")