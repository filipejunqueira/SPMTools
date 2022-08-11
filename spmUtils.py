import numpy as np
import pandas as pd
import access2thematrix
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import hyp2f1 as Fy
from scipy.optimize import curve_fit
import os
import csv
sns.set()  # Setting seaborn as default style even if use only matplotlib


def plot_forces_short_range(Force_ON_trace,Force_ON_retrace,Force_OFF,z):

    Force_diference_trace = Force_ON_trace - Force_OFF
    Force_diference_retrace = Force_ON_retrace - Force_OFF

    figure, axis = plt.subplots(1, 2, figsize=(20, 5), sharex=True)

    sns.lineplot(ax=axis[0], x=z, y=Force_diference_trace * 10 ** 9, color='red')

    axis[0].set_title("Force (ON - OFF) vs Z - Trace")
    axis[0].set_xlabel("Z[m]")
    axis[0].set_ylabel("Force [nN]")
    #axis[0].text(x=0, y=0, s=f"Min Force {np.round(np.min(Force_ON_trace) * 10 ** 9, 2)} nN")
    #axis[0].text(x=0, y=0, s=f"Min Force {np.round(np.min(Force_ON_trace) * 10 ** 9, 2)} nN")

    sns.lineplot(ax=axis[1], x=z, y=Force_diference_retrace * 10 ** 9, color='red')
    axis[1].set_title("Force (ON - OFF) vs Z - Retrace")
    axis[1].set_xlabel("Z[m]")
    axis[1].set_ylabel("Force [nN]")

    plt.show()




class Spec_curve:
    def __init__(self, curve):
        self.meta = curve.referenced_by
        self.Z = curve.data[0]
        self.df_Z = curve.data[1]

    # this section makes sure that the z which is the first vector of (Z) is always crescent.
    # If it is not it will flip the self.trace or self.retrace in order to make it crescent in Z.
        if self.Z[0] < self.Z[-1]:
            pass
        elif self.Z[0] > self.Z[-1]:
            self.Z = np.flipud(self.Z)
            self.df_Z = np.flipud(self.df_Z)

        self.data_frame = pd.DataFrame({'Z':self.Z, 'deltaF': self.df_Z}, columns=['Z','deltaF'])

def import_spectra(path):

    mtrx_data = access2thematrix.MtrxData()
    data_file = f'{path}'
    traces, message = mtrx_data.open(data_file)
    curve_trace, message = mtrx_data.select_curve(traces[0])
    curve_retrace, message = mtrx_data.select_curve(traces[1])
    return Spec_curve(curve_trace), Spec_curve(curve_retrace)


def load_spec_list_from_cvs(folder_base_path=f"{os.getcwd()}"):
    data_list = []
    with open(f'{folder_base_path}/spec_list.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, skipinitialspace=True)
        for on, off in reader:
            data_list.append((on, off))
    return data_list



def fit_lennard_jones():
    pass


def sjarvis_benchmark(zmin=0.23E-9, zmax=5.000E-9, points=5000, sigma=0.235E-9, E=0.371E-18, f0=32768, k=1800,
                      A=11.8E-12, simple=False, plot=False):
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


############### References: ##################

'''

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

'''


def sjarvis_deconvolution(df_Z, A=0.01E-9, f0=-25000, k=1800):
    ############### specification: ##################
    '''
    It receives two panda data frames with the following columns:

    INPUT: (1 data frame + 3 float numbers).

    df_Z["deltaF", "Z"] - it is assumed that z(i) is closer to the surface than z(i+1) in meters | z(i) < z(i+1)
    If this condition is not met we get errors (sqrt or negative numbers).

    a: is the amplitude in m
    f0: is resonance frequency in Hz
    k: is the stiffness of the cantilever in N/m

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
            y=(1 + np.sqrt(A) / (8 * np.sqrt(np.pi * (t - z[j])))) * omega_temp - A ** (3 / 2) / np.sqrt(
                2 * (t - z[j])) * dz_omega_temp, x=t)

        # 5c) corrections for t=z C(j) #All corrections are already divided automatically by f0 since Im using the corrected omega (delta_f/f0)
        correction1 = omega[j] * (z[j + 1] - z[j])
        correction2 = 2 * (np.sqrt(A) / (8 * np.sqrt(np.pi))) * (omega[j] * np.sqrt(z[j + 1] - z[j]))
        correction3 = (-np.sqrt(2)) * (A ** (3 / 2)) * (dz_omega[j] * np.sqrt(z[j + 1] - z[j]))
        # 6) F(i) = 2*k*(correction1 + correction2 + correction3 + integral)
        force[j] = 2 * k * (correction1 + correction2 + correction3 + integral)

    # 7) adjust z to be the same size as F

    z = z[:len(force)]
    print(len(z))

    # else:
    #     break
    #     #ValueError('the Z values of the off curve are not the same as the on, fix this first!.')
    return force, z, delta_f, omega, dz_omega


