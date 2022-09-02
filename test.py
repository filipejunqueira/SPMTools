import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import access2thematrix
import seaborn as sns
from spmUtils import *
from scipy.signal import savgol_filter
from lmfit import Model
from lmfit.models import LorentzianModel, LinearModel, VoigtModel, ExponentialModel
import csv

root_path = "/media/captainbroccoli/DATA/"
project_folder_name = "2022-08-12"
prefix = "20220812-124820_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_mtrx"


sns.set()  # Setting seaborn as default style even if use only matplotlib
sns.set(style="ticks", context="talk")

fontsize = 30
marker_size = 150
#sns.set_style("darkgrid", {"grid.color": "1", "grid.linestyle": ":"})


type = "Aux2(V)"
figsize = (20,14)
file_id = 37
n_files = 1
plot_retrace = False
end = 512
# If there is no retrace it will give out an error!!!

#plot_single_curve(file_id, n_files, plot_retrace=plot_retrace, type =type, filter = False, fontsize=fontsize,marker_size=marker_size, figsize=figsize, slice_end=end)
path = os.path.join(project_folder_path, f"{prefix}{file_id}_{n_files}.{type}{sufix}")
print(path)
spec = import_spectra(path)
x = spec.X
y = spec.Y
x = x[:end]
y = y[:end]

linear_mod =LinearModel()
parameters = linear_mod.make_params(intercept=y.min(), slope=0)

exp_mod = ExponentialModel(prefix='exp_')
parameters.update(exp_mod.guess(y,x=x))

lorentzian_mod =LorentzianModel()
parameters.update(lorentzian_mod.make_params())

voigt_mod =VoigtModel()
parameters.update(voigt_mod.guess(y,x=x))

model =  lorentzian_mod + exp_mod

init = model.eval(parameters, x=x)
result = model.fit(y, x=x, params=parameters)
print(result.fit_report())

y_bestfit = result.best_fit

plt.plot(x, y, 'o')
plt.plot(x, y_bestfit, '-', label='best fit', )

plt.legend()
plt.show()

exponential = list(map(''.join, itertools.product(*zip("exponential".upper(), "exponential".lower()))))
exponential.append(list(map(''.join, itertools.product(*zip("exp".upper(), "exp".lower())))))


print("I owe you nothing")