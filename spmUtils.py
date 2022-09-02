import numpy as np
import pandas as pd
import access2thematrix
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import hyp2f1 as Fy
from scipy.optimize import curve_fit
import os
import csv
from lmfit.models import (StepModel,ExpressionModel, PolynomialModel, QuadraticModel, LorentzianModel, LinearModel, VoigtModel, ExponentialModel, GaussianModel, SkewedVoigtModel, SkewedGaussianModel, SplitLorentzianModel, PseudoVoigtModel, SkewedGaussianModel,LognormalModel,ExponentialGaussianModel)
from lmfit import Model
import itertools

## DISCLAIMER: Documentation was mostly created using AI! called  Mintilify DocWriter.


# Hex color for graphs
color_map = {"green": "#90EE90", "dark_green": "#95BA61" , "orange": "#FFAB00", "dark_yellow": "#8B8000", "yellow": "#FFDB58", "red": "#FF5733", "blue": "#3c59ff", "white": "#FFFFFF.",
             "purple": "#7b40c9", "pink": "#FFB490", "black": "#171717"}


# THIS NEEDS TO BE CHANGE IF THE PROJECT FOLDER CHANGES!

root_path = "/media/captainbroccoli/DATA/"
project_folder_name = "2022-08-12"
prefix = "20220812-124820_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"
sufix = "_mtrx"

project_folder_path = os.path.join(root_path,project_folder_name)

# USEFUL FOR SINGLE PLOT   ############################################################################################

def create_file_number_list(file_id,n_files):
    """
    This function takes a file id and the number of files in the directory and returns a list of file numbers

    :param file_id: the file number you want to start with
    :param n_files: the number of files you want to create
    Examples: 37 and 5 will generate ['37_1', '37_2', '37_3', '37_4', '37_5']
    """
    file_number_list=[]
    for i in range(n_files):
        file_number_list.append(f"{file_id}_{i+1}")
    return file_number_list


def average_curves(file_number_list,type,direction=0):
    """
    It takes a list of file numbers, a type of file, and a direction (0 or 1) and returns the average of the curves in the
    files

    :param file_number_list: a list of the file numbers you want to average
    :param type: is the type of file you want to average
    :param direction: 0 = forward, 1 = backward, defaults to 0 (optional)
    :return: The x and y values of the averaged curve.
    """
    path = os.path.join(project_folder_path,f"{prefix}{file_number_list[0]}.{type}{sufix}")

    #I need to know how
    #path = f"{project_folder_path}{prefix}{file_number_list[0]}.{type}{sufix}"
    mtrx_data = access2thematrix.MtrxData()
    data_file = f'{path}'
    traces, message = mtrx_data.open(data_file)
    curve_trace, message = mtrx_data.select_curve(traces[direction])
    array_lengh = len(Spec_curve(curve_trace).X)

    x = np.zeros(array_lengh)
    y = np.zeros(array_lengh)


    for number in file_number_list:
        # this path needs to be global.
        path = os.path.join(project_folder_path, f"{prefix}{number}.{type}{sufix}")
        print(f"Loading: {path}")
        mtrx_data = access2thematrix.MtrxData()
        data_file = f'{path}'
        traces, message = mtrx_data.open(data_file)
        curve_trace, message = mtrx_data.select_curve(traces[direction])
        x = Spec_curve(curve_trace).X
        y = y + Spec_curve(curve_trace).Y

    y = y / len(file_number_list)

    return x, y

def quick_fit(x, y, linear = True, exponential = True, fit_type="lorentzian", n_curves=1, algo=""):
    """
    > This function takes in a set of x and y data, and a set of booleans that determine which models to fit to the data. It
    then returns the fit parameters and the initial guess for the fit

    :param x: the x-axis data
    :param y: the data to fit
    :param linear: True/False, defaults to False (optional)
    :param exponential: y = a*exp(b*x), defaults to False (optional)
    :param gaussian: Gaussian/normal distribution, defaults to False (optional)
    :param lorentzian: , defaults to False (optional)
    :param voigt: VoigtModel(), defaults to False (optional)
    :return: The result of the fit and the initial guess.
    """



    def empty(x):
        return 0

    initial_model = Model(empty)
    parameters = initial_model.make_params()

    # generating every possible combination of exponential and linear string inputs
    #exponential = list(map(''.join, itertools.product(*zip("exponential".upper(), "exponential".lower()))))
    #exponential.append(list(map(''.join, itertools.product(*zip("exp".upper(), "exp".lower())))))


    #background

    if linear == True:

        linear_mod = LinearModel(prefix='lin_')
        parameters.update(linear_mod.guess(y, x=x))
        initial_model = initial_model + linear_mod

    if exponential == True:

        exp_mod = ExponentialModel(prefix='exp_')
        parameters.update(exp_mod.guess(y, x=x))
        initial_model = initial_model + exp_mod

    background = initial_model
    final_model = background

    for i in range(n_curves):

        if fit_type == "lorentzian":

            fit_model = LorentzianModel(prefix=f'lor_{i}')
            parameters.update(fit_model.guess(y, x=x))
            final_model = background + fit_model

        elif fit_type == "voigt":

            fit_model = VoigtModel(prefix=f'voi_{i}')
            parameters.update(fit_model.guess(y, x=x))
            final_model = background + fit_model

        elif fit_type == "gaussian":

            fit_model = GaussianModel(prefix=f'gauss_{i}')
            parameters.update(fit_model.guess(y, x=x))
            final_model = background + fit_model


    init = final_model.eval(parameters, x=x)
    result = final_model.fit(y, x=x, params=parameters, method="leastsq")

    return result, background, init



def plot_single_curve(file_id, n_files=1, type ="Aux2(V)", plot_retrace=False, slice_start= 0, slice_end=512, fontsize=36, filter = False, filter_order = 3, filter_window = 5, root_path="", figsize=(14, 14),marker_size=100, bestfit=True, y_fit= None, y_fit_retrace=None):

    file_number_list = create_file_number_list(file_id, n_files)
    filtered_str = ""

    x, y = average_curves(file_number_list, type, direction=0)
    x = x[slice_start:slice_end]
    y = y[slice_start:slice_end]

    print(f"Curves {file_number_list} of type: {type} have been averaged ")
    print(f"Attention: I'm slicing the graph from {slice_start} to {slice_end} points. Make sure this is what you want.")

    if filter == True:
        y = savgol_filter(y, filter_window, filter_order)
        filtered_str = f"filter_w{filter_window}o{filter_order}"

    if plot_retrace==True:
        x_retrace, y_retrace = average_curves(file_number_list, type, direction=1)
        x_retrace = x_retrace[slice_start:slice_end]
        y_retrace = y_retrace[slice_start:slice_end]
        print(f"Attention - Retraces are also being computed")
        if filter == True:
            y = savgol_filter(y, filter_window, filter_order)
            filtered_str = f"filter_w{filter_window}o{filter_order}"

    figure, axis = plt.subplots(1, 1, figsize=figsize, sharex=True)

    if type == "Aux2(V)":

        if plot_retrace == False:
            if bestfit == True:
                sns.lineplot(ax=axis, x=x, y=y_fit, color=f"{color_map['black']}", alpha=1, label="best fit")
            sns.scatterplot(ax=axis, x=x, y=y, alpha=1, edgecolor=f"{color_map['purple']}", facecolor="None", s=marker_size)

        if plot_retrace == True:
            sns.lineplot(ax=axis, x=x, y=y_fit, color=f"{color_map['black']}", alpha=1, label="best fit")
            sns.scatterplot(ax=axis, x=x, y=y, alpha=1, edgecolor=f"{color_map['purple']}", facecolor="None",
                            s=marker_size,label ="trace")
            sns.lineplot(ax=axis, x=x_retrace, y=y_fit_retrace, color=f"{color_map['black']}", alpha=1, label="best fit")
            sns.scatterplot(ax=axis, x=x_retrace, y=y_retrace, alpha=1, edgecolor=f"{color_map['blue']}", facecolor="None",label ="retrace",
                            s=marker_size)

        title = f"dI/dV"
        name_x = "Bias"
        name_y = "dI/dV"
        unit_x = "[V]"
        unit_y = "[arb. units]"

        axis.set_title(f"{title}", fontsize=fontsize)
        axis.set_xlabel(f"{name_x} {unit_x}", fontsize=(fontsize - 2))
        axis.set_ylabel(f"{name_y} {unit_y}", fontsize=(fontsize - 2))
        plt.xticks(fontsize=(fontsize - 2))
        plt.yticks([])

    elif type == "Df(V)":

        if plot_retrace == False:
            sns.lineplot(ax=axis, x=x, y=y, color=f"{color_map['black']}", alpha=0.1)
            sns.scatterplot(ax=axis, x=x, y=y, alpha=1, edgecolor=f"{color_map['black']}", facecolor="None", s=marker_size)

        if plot_retrace == True:
            sns.lineplot(ax=axis, x=x, y=y, color=f"{color_map['black']}", alpha=0.1)
            sns.scatterplot(ax=axis, x=x, y=y, alpha=1, edgecolor=f"{color_map['black']}", facecolor="None",
                            s=marker_size, label ="trace")
            sns.lineplot(ax=axis, x=x_retrace, y=y_retrace, color=f"{color_map['black']}", alpha=0.1)
            sns.scatterplot(ax=axis, x=x_retrace, y=y_retrace, alpha=1, edgecolor=f"{color_map['pink']}", facecolor="None",
                            s=marker_size, label ="retrace")

        title = f"Df(V)"
        name_x = "Bias"
        name_y = "Df(V)"
        unit_x = "[V]"
        unit_y = "[Hz]"

        axis.set_title(f"{title}", fontsize=fontsize)
        axis.set_xlabel(f"{name_x} {unit_x}", fontsize=(fontsize - 2))
        axis.set_ylabel(f"{name_y} {unit_y}", fontsize=(fontsize - 2))
        plt.xticks(fontsize=(fontsize - 2))
        plt.yticks(fontsize=(fontsize - 2))

    elif type == "Df(Z)":

        x, y = average_curves(file_number_list, type, direction=0)

        if plot_retrace == False:
            sns.lineplot(ax=axis, x=x * 10 ** 9, y=y, color=color_map["green"], label="df trace",alpha=1)
            sns.scatterplot(ax=axis, x=x* 10 ** 9, y=y, alpha=0.1, edgecolor=f"{color_map['dark_green']}", facecolor="None", s=marker_size)


        if plot_retrace == True:
            sns.lineplot(ax=axis, x=x * 10 ** 9, y=y, color=color_map["green"], label="df trace", alpha=1)
            sns.scatterplot(ax=axis, x=x * 10 ** 9, y=y, alpha=0.1, edgecolor=f"{color_map['dark_green']}",
                            facecolor="None", s=marker_size)
            x_retrace, y_retrace = average_curves(file_number_list, type, direction=1)
            sns.lineplot(ax=axis, x=x_retrace * 10 ** 9, y=y_retrace, color=color_map["yellow"], label="df retrace", alpha=1)
            sns.scatterplot(ax=axis, x=x_retrace * 10 ** 9, y=y_retrace, alpha=0.1, edgecolor=f"{color_map['dark_yellow']}", facecolor="None",
                            s=marker_size)

        title = f"Df(Z)"
        name_x = "Z"
        name_y = "Df(Z)"
        unit_x = "[nm]"
        unit_y = "[Hz]"

        axis.set_title(f"{title}", fontsize=fontsize)
        axis.set_xlabel(f"{name_x} {unit_x}", fontsize=(fontsize - 2))
        axis.set_ylabel(f"{name_y} {unit_y}", fontsize=(fontsize - 2))
        plt.xticks(fontsize=(fontsize - 2))
        plt.yticks(fontsize=(fontsize - 2))
        axis.legend(loc=0)



    pwd = os.getcwd()
    name = f"{type}_{project_folder_name}--{file_id}_{n_files}{filtered_str}"
    save_dir_name = "single_curve_graphs"
    save_dir= os.path.join(pwd,save_dir_name)
    if os.path.isdir(save_dir) == False:
        os.mkdir(save_dir)
    else:
        pass

    image_name = os.path.join(pwd,save_dir_name,name)
    plt.savefig(fname=f"{image_name}.svg", facecolor='auto', edgecolor='auto', transparent=True)
    plt.show()

        # axis.legend(bbox_to_anchor = (1.01, 1), loc = 'upper left')
    return True



# Object spec curve   ############################################################################################


# The class Spec_curve takes a curve object from the AFM_data class and creates a dataframe with the Z and deltaF values
class Spec_curve:
    def __init__(self, curve):
        self.meta = curve.referenced_by
        self.X = curve.data[0] #Parameter X have the x vector as numpy array
        self.Y = curve.data[1] #Parameter Y have the y vector as numpy array

        # this section makes sure that the X which is the first vector of (X) is always crescent.
        # If it is not it will flip the self.trace or self.retrace in order to make it crescent in X.
        # THis could be in theory anything.

        if self.X[0] < self.X[-1]:
            pass
        elif self.X[0] > self.X[-1]:
            self.X = np.flipud(self.X)
            self.Y = np.flipud(self.Y)

        # THis needs to be changed. But in principle it assumes this is deltaF and Z collums in the dataframe.
        self.data_framedf = pd.DataFrame({'Z': self.X, 'deltaF': self.Y}, columns=['Z', 'deltaF'])


def plot_curve_list(curve_list, X, save=True, title="Title", name_x="X", name_y="Y", unit_x="[X]", unit_y="[Y]",
                    folder_name="curves"):
    figure, axis = plt.subplots(1, 1, figsize=(10, 6), sharex=True)

    for curve in curve_list:
        sns.lineplot(ax=axis, x=X, y=curve)

    axis.set_title(f"{title}")
    axis.set_xlabel(f"{name_x}{unit_x}")
    axis.set_ylabel(f"{name_y}{unit_y}")

    if save == True:
        pwd = os.getcwd()
        if os.path.isdir(f"{pwd}/{folder_name}") == False:
            os.mkdir(f"{pwd}/{folder_name}")
        else:
            pass
        plt.savefig(fname=f"{pwd}/{folder_name}/{name_y}Vs{name_x}", formatstr='.eps', facecolor='auto',
                    edgecolor='auto')
    else:
        pass

    plt.show()


def get_matrix_image(image_path):
    mtrx_data = access2thematrix.MtrxData()
    traces, message = mtrx_data.open(image_path)
    print(message)

    im, message_im = mtrx_data.select_image(traces)

    return im, message, message_im


def plot_df(df_ON_trace, df_ON_retrace, df_OFF, z, save=False, name="dfvsZ", retrace=False, off=True):
    figure, axis = plt.subplots(1, 1, figsize=(10, 7), sharex=True)

    sns.lineplot(ax=axis, x=z * 10 ** 9, y=df_ON_trace, color=color_map["green"], label="df ON trace")
    if retrace == True:
        sns.lineplot(ax=axis, x=z * 10 ** 9, y=df_ON_retrace, color=color_map["yellow"], label="df ON retrace")
    elif retrace == False:
        pass
    if off == True:
        sns.lineplot(ax=axis, x=z * 10 ** 9, y=df_OFF, color=color_map["red"], label="df OFF")
    elif off == False:
        pass

    axis.set_title("df (ON and OFF) vs Z")
    axis.set_xlabel("Z[nm]")
    axis.set_ylabel("df[Hz]")
    axis.legend(loc=0)

    if save == True:
        pwd = os.getcwd()
        directory_path = "force_graphs"
        if os.path.isdir(f"{pwd}/{directory_path}") == False:
            os.mkdir(f"{pwd}/{directory_path}")
        else:
            pass
        plt.savefig(fname=f"{pwd}/{directory_path}/{name}", formatstr='.eps', facecolor='auto', edgecolor='auto')
    else:
        pass

    plt.show()


def plot_forces_direct(Force_ON_trace, Force_ON_retrace, Force_OFF, z, save=False, name="ForceVsZ", retrace=False):
    figure, axis = plt.subplots(1, 1, figsize=(10, 7), sharex=True)

    sns.lineplot(ax=axis, x=z * 10 ** 9, y=Force_ON_trace * 10 ** 9, color=color_map["green"], label="force trace")
    if retrace == True:
        sns.lineplot(ax=axis, x=z * 10 ** 9, y=Force_ON_retrace * 10 ** 9, color=color_map["yellow"],
                 label="force retrace")
        axis.text(0.5, 0.4, f"min F (trace): {np.round(np.min(Force_ON_retrace) * 10 ** 9, 1)}nN",
                  transform=axis.transAxes)
    elif retrace == False:
        pass

    sns.lineplot(ax=axis, x=z * 10 ** 9, y=Force_OFF * 10 ** 9, color=color_map["red"], label="force off")

    axis.set_title("Force(ON) and Force(OFF) Vs Z")
    axis.set_xlabel("Z[nm]")
    axis.set_ylabel("Force [nN]")
    axis.text(0.5, 0.5, f"min F (trace): {np.round(np.min(Force_ON_trace) * 10 ** 9, 1)}nN", transform=axis.transAxes)


    axis.legend(loc=0)

    if save == True:
        pwd = os.getcwd()
        directory_path = "force_graphs"
        if os.path.isdir(f"{pwd}/{directory_path}") == False:
            os.mkdir(f"{pwd}/{directory_path}")
        else:
            pass
        plt.savefig(fname=f"{pwd}/{directory_path}/{name}", formatstr='.eps', facecolor='auto', edgecolor='auto')
    else:
        pass

    plt.show()


def plot_forces_short_range(force_diff_trace, force_diff_retrace, z_on, save=False, name="ForceVsZ", retrace=False):
    figure, axis = plt.subplots(1, 1, figsize=(10, 7), sharex=True)

    sns.lineplot(ax=axis, x=z_on * 10 ** 9, y=force_diff_trace * 10 ** 9, color=color_map["green"],
                 label="force diff trace")
    axis.text(0.5, 0.5, f"min F (trace): {np.round(np.min(force_diff_trace) * 10 ** 9, 1)}nN",
              transform=axis.transAxes)
    if retrace == True:
        sns.lineplot(ax=axis, x=z_on * 10 ** 9, y=force_diff_retrace * 10 ** 9, color=color_map["yellow"],
                 label="force diff retrace")

        axis.text(0.5, 0.4, f"min F (trace): {np.round(np.min(force_diff_retrace) * 10 ** 12, 1)}nN",
                  transform=axis.transAxes)
    elif retrace == False:
        pass

    axis.set_title("Force (ON - OFF) vs Z")
    axis.set_xlabel("Z[nm]")
    axis.set_ylabel("Force [nN]")



    axis.legend(loc=0)
    if save == True:
        pwd = os.getcwd()
        directory_path = "force_graphs"
        if os.path.isdir(f"{pwd}/{directory_path}") == False:
            os.mkdir(f"{pwd}/{directory_path}")
        else:
            pass
        plt.savefig(fname=f"{pwd}/{directory_path}/{name}", formatstr='.eps', facecolor='auto', edgecolor='auto')
    else:
        pass

    plt.show()


def plot_forces_and_df(force_trace, force_retrace, df_trace, df_retrace, df_off, z_force, z_df, save=False,
                       name="ForceVsZ", retrace=False):
    figure, axis = plt.subplots(1, 2, figsize=(20, 7), sharex=True)

    sns.lineplot(ax=axis[0], x=z_force * 10 ** 9, y=force_trace * 10 ** 9, color=color_map["green"],
                 label="force trace")
    axis[0].text(0.5, 0.5,
                 f"min F (trace): {np.round(np.min(force_trace) * 10 ** 9, 1)}nN", transform=axis[0].transAxes)
    if retrace == True:
        sns.lineplot(ax=axis[0], x=z_force * 10 ** 9, y=force_retrace * 10 ** 9, color=color_map["yellow"],
                 label="force retrace")

        axis[0].text(0.5, 0.4,f"min F (retrace): {np.round(np.min(force_retrace) * 10 ** 9, 1)}nN", transform=axis[0].transAxes)
    elif retrace == False:
        pass

    axis[0].set_title("Force (ON - OFF) vs Z")
    axis[0].set_xlabel("Z[nm]")
    axis[0].set_ylabel("Force[nN]")
    axis[0].legend(loc=0)

    sns.lineplot(ax=axis[1], x=z_df * 10 ** 9, y=df_trace, color=color_map["green"], label="df trace")
    sns.lineplot(ax=axis[1], x=z_df * 10 ** 9, y=df_off, color=color_map["red"], label="df off")
    sns.lineplot(ax=axis[1], x=z_df * 10 ** 9, y=df_retrace, color=color_map["yellow"], label="df retrace")

    axis[1].set_title("Df (Ttrace,Retrace,OFF) vs Z")
    axis[1].set_xlabel("Z[nm]")
    axis[1].set_ylabel("df[Hz]")
    axis[1].legend(loc=0)

    if save == True:
        pwd = os.getcwd()
        force_dir_graphs = "force_graphs"
        if os.path.isdir(f"{pwd}/{force_dir_graphs}") == False:
            os.mkdir(f"{pwd}/{force_dir_graphs}")
        else:
            pass
        plt.savefig(fname=f"{pwd}/{force_dir_graphs}/{name}", formatstr='.eps', facecolor='auto', edgecolor='auto')
    else:
        pass

    plt.show()


def import_spectra(path):
    """
    > The function `import_spectra` takes a path to a .mtrx file and returns a `Spec_curve` object for the trace and retrace

    :param path: the path to the file you want to import
    :return: two objects of the class Spec_curve.
    """
    mtrx_data = access2thematrix.MtrxData()
    data_file = f'{path}'
    traces, message = mtrx_data.open(data_file)
    curve_trace, message = mtrx_data.select_curve(traces[0])
    try:
        curve_retrace, message = mtrx_data.select_curve(traces[1])
    except:
        return Spec_curve(curve_trace)
        print("Spectra imported - Only trace exists")
    else:
        return Spec_curve(curve_trace), Spec_curve(curve_retrace)
        print("Spectra imported - Both trace and retrace exists")

    # TODO: Return this as a dictionary instead of an object {trace: curve_trace, retrace: curve_retrace, Z: }


def load_on_off_spec_list_from_cvs(folder_base_path=f"{os.getcwd()}", cvs_name="spec_on_off_list"):
    """
    It opens the spec_on_off_list.csv file, reads the data, and returns a list of tuples

    :param folder_base_path: The path to the folder where the spec_on_off_list.csv file is located, defaults to f"{os.getcwd()}"
    (optional)
    :return: A list of tuples. Those are your ON and OFF number files.
    """
    data_list = []
    with open(f'{folder_base_path}/{cvs_name}.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, skipinitialspace=True)
        for on, off in reader:
            data_list.append((on, off))
    return data_list


def load_spec_list_from_cvs(folder_base_path=f"{os.getcwd()}", cvs_name="spec_list"):
    """
    It opens the spec_on_off_list.csv file, reads the data, and returns a list of tuples

    :param folder_base_path: The path to the folder where the spec_on_off_list.csv file is located, defaults to f"{os.getcwd()}"
    (optional)
    :return: A list of tuples. Those are your ON and OFF number files.
    """
    data_list = []
    with open(f'{folder_base_path}/{cvs_name}.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, skipinitialspace=True)
        for item1, item2 in reader:
            data_list.append((item1, item2))
    return data_list


def fit_lennard_jones(z_on, z_off, df_off, A=0.2E-9, k=1800, f0=25000, simple=False):
    """

    # NOT WORKING AT THE MOMENT!!!
    It fits the Lennard-Jones potential to the OFF trace and returns the fitted OFF trace

    :param z_on: the z-values for the ON trace
    :param z_off: the z-values for the OFF trace
    :param df_off: the force OFF trace
    :param A: tip Amplitude
    :param k: spring constant, defaults to 1800 (optional)
    :param f0: the resonance frequency of the cantilever, defaults to 25000 (optional)
    :param simple: if True, then the fit is done with a simple Lennard-Jones potential. If False, then the fit is done with
    a more complex Lennard-Jones potential, defaults to False (optional)
    :return: The fitted_df_off is the fitted force OFF trace.
    """
    print(f"range z_on [{np.min(z_on)} to {np.max(z_on)}]")
    print(f"range z_off [{np.min(z_off)} to {np.max(z_off)}]")

    #translating z_on first

    u_on = z_on-min(z_on)+10**-9
    u_off =z_off-min(z_off)+10**-9

    print(f"range u_on [{np.min(u_on)} to {np.max(u_on)}]")
    print(f"range u_off [{np.min(u_off)} to {np.max(u_off)}]")
    if simple == False:

        def off_objective(z, a, b, c):
            return a*(1/(z+b)**(-13)) + c

        E = 0.371E-18
        sigma = 0.235E-9
        prefactor = -12 * f0 * E / (sigma * k * A)

        # fitting, here is where the magic happens:
        parameters_opt, parameters_pcov = curve_fit(off_objective, u_off, df_off, p0=(prefactor, sigma, sigma))
        #parameters_opt, parameters_pcov = curve_fit(off_objective, u_off, df_off)

        print(parameters_opt)
        a, b, c = parameters_opt
        fitted_df_off = off_objective(u_on, a, b, c)

        ## plot for checking
        figure, axis = plt.subplots(1, 1, figsize=(10, 7), sharex=True)

        sns.lineplot(ax=axis, x=u_off, y=df_off, color=color_map["green"],
                     label="original")
        sns.lineplot(ax=axis, x=u_on, y=fitted_df_off, color=color_map["red"],
                     label="fitted")

        axis.set_title("DF (OFF and fitted OFF) vs Z")
        axis.set_xlabel("Z[nm]")
        axis.set_ylabel("df [Hz]")
        plt.show()

        return fitted_df_off, parameters_opt, parameters_pcov

    elif simple == True:
        pass


def sjarvis_benchmark(zmin=0.23E-9, zmax=5.000E-9, points=5000, sigma=0.235E-9, E=0.371E-18, f0=32768, k=1800,
                      A=11.8E-12, simple=False, plot=False):
    """
        It takes in a bunch of parameters and returns the minimum force as well as a benchmark for Lennard Jones

        :param zmin: minimum z value to calculate used to generate the z array [in m]
        :param zmax: The maximum distance between the tip and the surface, used to generate the z array [in m]
        :param points: number of points to calculate the mapping function, defaults to 5000 (optional)
        :param sigma: The distance at which the potential energy is zero
        :param E: Young's Modulus
        :param f0: Resonant frequency of the qplus cantilever, defaults to 32768 (optional)
        :param k: cantilever stiffness, defaults to 1800
        :param A: The amplitude of the oscillation
        :param simple: If you want to use the simple approximation of the Lennard Jones potential, set this to True, defaults to
        False (optional)
        :param plot: Boolean, if True, plots the force vs Z and deltaF vs Z, defaults to False (optional)
        :return: The minimum force.
    """

    if simple == False:

        # BENCHMARK Lennard Jones using Hyper-geometric function
        z = np.linspace(zmin, zmax, num=points)
        # Using Lenard Jones Potential - Fy is that fancy function - check equation 10 of the paper [2].
        pre_factor = -12 * f0 * E / (sigma * k * A)

        df_on = pre_factor * (
                ((sigma / z) ** (7)) * (Fy(7, 0.5, 1, (-2 * A / z)) - Fy(7, 1.5, 2, (-2 * A / z)))) - pre_factor * (
                        ((sigma / z) ** (13)) * (Fy(13, 0.5, 1, (-2 * A / z)) - Fy(13, 1.5, 2, (-2 * A / z))))
        ON = pd.DataFrame({'deltaF': df_on, 'Z': z})

    elif simple == True:
        # BENCHMARK Lennard Jones using simple approximation
        z = np.linspace(zmin, zmax, num=points)
        pre_factor = -12 * f0 * E / (sigma * k * A)
        df = pre_factor * ((sigma / z) ** (8)) - pre_factor * ((sigma / z) ** (14))
        ON = pd.DataFrame({'deltaF': df_on, 'Z': z})

    force, delta_f, omega, dz_omega, z_force = sjarvis_deconvolution(df_Z=ON, A=A, f0=f0, k=k)

    # Plotting stuff.
    cmap = sns.light_palette("#77DD77", as_cmap=True)

    # create a line plot for the mapping function
    # Initialise the subplot function using number of rows and columns

    if plot == True:
        figure, axis = plt.subplots(1, 2, figsize=(15, 5), sharex=True)
        sns.lineplot(ax=axis[0], x=z_force, y=force * 10 ** 9, color='#ff6969')
        axis[0].set_title("Force vs Z")
        axis[0].set_xlabel("Z[m]")
        axis[0].set_ylabel("Force [nN]")
        axis[0].text(x=0, y=0, s=f"Min Force {np.round(np.min(force) * 10 ** 9, 2)} nN")

        sns.lineplot(ax=axis[1], x=z_force, y=delta_f[:-1])
        axis[1].set_title("deltaF vs Z")
        axis[1].set_xlabel("Z[m]")
        axis[1].set_ylabel("Df [Hz]")

        plt.show()

    else:
        pass

    min_force = np.min(force)
    return min_force


# Deconvoluting the data.
def sjarvis_deconvolution(df_Z, A=0.01E-9, f0=-25000, k=1800):
    ############### specification: ##################
    '''
    It receives two panda data frames with the following columns:
    It calculates the force from the frequency shift and the distance to the surface

    INPUT: (1 data frame + 3 float numbers).

    df_Z["deltaF", "Z"] - it is assumed that z(i) is closer to the surface than z(i+1) in meters | z(i) < z(i+1)
    If this condition is not met we get errors (sqrt or negative numbers).

    a: is the amplitude in m
    f0: is resonance frequency in Hz
    k: is the stiffness of the cantilever in N/m, defaults to 1800 (optional)

    OUTPUT: [5 numpy arrays]

    force: force in N
    delta_f: that's result of df_ON - df_OFF
    z: distance in m (same size as force).
    omega = delta_f/f0
    dz_omega = derivative of omega (calculated by element wise ratio of diff omega and diff z)

    '''

    # 1) Transform into numpy arrays as it's easier to work with them like this inside the function.
    delta_f = df_Z["deltaF"].to_numpy()
    z = df_Z["Z"].to_numpy()

    # 2) Calculates reduced frequency shift omega = delta_f/f_0.

    omega = delta_f / f0

    # 3) Derivative of the reduced frequency shift dz_omega

    dz_omega = np.ediff1d(omega) / np.ediff1d(z)

    # 4) Because we are umbral/discrete "world" we must adjust the size of Z, omega and delta_f to be the same as dz_omega

    z = z[:-1]
    delta_f = delta_f[:-1]
    omega = omega[:-1]

    # 5) Main loop.

    # Initiate force array
    force = np.zeros(len(z) - 1)

    for j in range(len(z) - 2):
        # 5a) Adjust length to the size of t.
        t = z[j + 1:]  # this is due to a pole at t = z
        omega_temp = omega[j + 1:]
        dz_omega_temp = dz_omega[j + 1:]

        # 5b) calculate integral Eq.(9) in [1] using the trapezoidal method.
        integral = np.trapz(
            y=(1 + np.sqrt(A) / (8 * np.sqrt(np.pi * (t - z[j])))) * omega_temp - ((A ** (3 / 2)) / (np.sqrt(
                2 * (t - z[j])))) * dz_omega_temp, x=t)

        # 5c) corrections for t=z C(j) #All corrections are already divided automatically by f0 since Im using the corrected omega (delta_f/f0)
        correction1 = omega[j] * (z[j + 1] - z[j])
        correction2 = 2 * (np.sqrt(A) / (8 * np.sqrt(np.pi))) * (omega[j] * np.sqrt(z[j + 1] - z[j]))
        correction3 = -(np.sqrt(2)) * (A ** (3 / 2)) * (dz_omega[j] * np.sqrt(z[j + 1] - z[j]))
        # 6) F(i) = 2*k*(correction1 + correction2 + correction3 + integral)
        force[j] = 2 * k * (correction1 + correction2 + correction3 + integral)

    # 7) adjust z to be the same size as F

    z = z[:len(force)]

    # else:
    #     break
    #     #ValueError('the Z values of the off curve are not the same as the on, fix this first!.')
    return force, z, delta_f, omega, dz_omega


############### References: ##################

"""

[1] J. E. Sader and S. P. Jarvis
       "Accurate formulas for interaction force and energy 
       in frequency modulation force spectroscopy"
       Applied Physics Letters, 84, 1801-1803 (2004)  
       https://aip.scitation.org/doi/10.1063/1.1667267

[2] Joachim Welker, Esther Illek and Franz J. Giessibl:
       "Analysis of force-deconvolution methods in frequency-modulation atomic force microscopy"
        https://www.beilstein-journals.org/bjnano/articles/3/27#E7
        In particular equations 10, 12 and 13. 

Other important papers to read about this include:


[3] Adam Sweetman and Andrew Stannard
    Beilstein J Nanotechnol. 2014
    "Uncertainties in forces extracted from non-contact atomic force microscopy measurements by fitting of long-range background forces"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3999863/

[4] Discriminating short-range from van der Waals forces using total force data
    in noncontact atomic force microscopy 
    Phys. Rev. B 89, 235417 – Published 13 June 2014
    https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.235417

[5] Andrew Stannard & Adam M. Sweetman 
    "A Considered Approach to Force Extraction from Dynamic Force Microscopy Measurements"
    First Online: 29 April 2015 - Part of the Advances in Atom and Single Molecule Machines book series (AASMM)
    https://link.springer.com/chapter/10.1007/978-3-319-17401-3_4

[6] John E. Sader, Barry D. Hughes, Ferdinand Huber & Franz J. Giessibl 
    "Interatomic force laws that evade dynamic measurement"
    Nature Nanotechnology volume 13 (2018)
    https://www.nature.com/articles/s41565-018-0277-x

[7] Daniel Heile, Reinhard Olbrich, Michael Reichling, and Philipp Rahe
    Phys. Rev. B 103, 075409 – Published 5 February 2021
    "Alignment method for the accurate and precise quantification of tip-surface forces"
    https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.075409


OTHER: 

[8] Quantitative Measurement of Short-Range Chemical Bonding Forces
    M. A. LantzH. J. HugR. HoffmannP. J. A. van SchendelP. KappenbergerS. MartinA. Baratoffand H.-J. Güntherodt
    DOI: 10.1126/science.1057824 - Science 30 Mar 2001
    https://www.science.org/doi/10.1126/science.1057824

[9] Manipulating Si(100) at 5 K using qPlus frequency modulated atomic force microscopy: Role of defects and dynamics in 
    the mechanical switching of atoms
    A. Sweetman, S. Jarvis, R. Danza, J. Bamidele, L. Kantorovich, and P. Moriarty
    Phys. Rev. B 84, 085426 – Published 25 August 2011
    https://doi.org/10.1103/PhysRevB.84.085426

"""

if __name__ == '__main__':
    main()
