"""
    File name: calculate_kappa.py
    Author: Matthias Blaschke
    Python Version: 3.9
"""
import codecs
import configparser
import sys

import numpy as np
import tmoutproc as top
import matplotlib
matplotlib.use('Agg') #for cluster usage!
import matplotlib.pyplot as plt
from numpy import inf


#h_bar in Js
h_bar = 1.0545718*10**(-34)
eV2hartree = 0.0367493
ang2bohr = 1.88973
har2J = 4.35974E-18
bohr2m = 5.29177E-11
u2kg = 1.66054E-27
har2pJ = 4.35974e-6
har2meV = 27211.396641308

def calculate_kappa(tau_ph, E, T):
    """Calculates the thermal conductance according to "Tuning the thermal conductance of molecular junctions with interference effects" (https://doi.org/10.1103/PhysRevB.96.245419)

    Args:
     tau_ph (array): Phonon transmission
     E (array): Phonon transmission
     T (float): Temperature in Kelvin

    Returns:
     kappa (np.array): Thermal conductance in Hartree/sK
    """
    #Planck sontant in hartree s
    h = 1.51983003855935e-16
    #Boltzmann constant in hartree/K
    k_B = 3.167*10**(-6)

    prefactor = 1./(h*k_B*T**2)

    beta = 1/(k_B*T)
    exp_ = np.exp(E*beta)
    #can be replaced because contribution is zuro in integrand (limits of exp)
    exp_[exp_ == inf] = 0
    integrand = E ** 2 * tau_ph * exp_ / ((exp_ - 1) ** 2)

    integral = np.trapz(integrand, E)

    return integral * prefactor

if __name__ == '__main__':
    config_path = sys.argv[1]
    cfg = configparser.ConfigParser()
    cfg.read_file(codecs.open(config_path, "r", "utf8"))

    try:
        data_path = str(cfg.get('Data Input', 'data_path'))
        transp_name = str(cfg.get('Data Input', 'transp_name'))
        transp_units = str(cfg.get('Data Input', 'transp_units', fallback="har"))


        # for thermal conducatance
        kappa_int_lower_E = float(cfg.get('Calculation', 'kappa_int_lower_E', fallback=-1))
        kappa_int_upper_E = float(cfg.get('Calculation', 'kappa_int_upper_E', fallback=-1))
        T_min = float(cfg.get('Calculation', 'T_min'))
        T_max = float(cfg.get('Calculation', 'T_max'))
        kappa_grid_points = int(cfg.get('Calculation', 'kappa_grid_points'))

    except configparser.NoOptionError:
        print("Missing option in config file. Check config file!")
        exit(-1)
    except ValueError:
        print("Wrong value in config file. Check config file!")
        exit(-1)

    transport = top.read_plot_data(data_path + "/" + transp_name)

    #Temperature range in Kelvin with given Resolution res
    res = kappa_grid_points
    T = np.linspace(T_min, T_max, res, dtype=np.float64)

    #Energy must be in Hartree! -> Convert otherwise
    if(transp_units=="har"):
        conv_factor = 1
    elif(transp_units=="har/(bohr**2*u)"):
        # convert from sqrt(har/(bohr**2*u)) to 1/s
        factor = np.sqrt(har2J/(bohr2m**2*u2kg))
        # convert to J
        factor = factor * h_bar
        # convert to hartree
        factor = factor / har2J
        conv_factor = factor

    E = np.asarray(transport[0, :][1:len(transport[0, :])]*conv_factor, dtype=np.float64)
    kappa = list()

    if(kappa_int_lower_E != -1 and kappa_int_upper_E !=-1):
        E_lower = kappa_int_lower_E / har2meV
        E_min_index = np.abs(E - E_lower).argmin()
        E_upper = kappa_int_upper_E / har2meV
        E_max_index = np.abs(E - E_upper).argmin()
        for i in range(E_min_index, E_max_index):
            E_i = E[E_min_index:i]
            thau_ph = np.asarray(transport[1, :][E_min_index:i], dtype=np.float64)
            T = 300
            kappa.append(calculate_kappa(thau_ph, E_i, T))
        kappa = np.asarray(kappa)
    else:
        thau_ph = np.asarray(transport[1, :][1:len(transport[1, :])], dtype=np.float64)

        for i in range(0, len(T)):
            kappa.append(calculate_kappa(thau_ph, E, T[i]))
        kappa = np.asarray(kappa)

    har2pJ = 4.35974e-6
    kappa = kappa*har2pJ

    #now plot everything
    plt.tick_params(axis="x", labelsize=15)
    plt.tick_params(axis="y", labelsize=15)
    #plt.legend(fontsize=13)
    plt.rc('xtick', labelsize=15)

    if(kappa_int_lower_E == -1):
        plt.ylabel(r'$\kappa_{\mathrm{ph}}$ ($\mathrm{pw/K}$)', fontsize=17)
        plt.plot(T, kappa)
        plt.xlabel('Temperature ($K$)', fontsize=17)
        # save kappa data
        top.write_plot_data(data_path + "/kappa.dat", (T, kappa), "T [K], kappa [pW/K]")
        plt.savefig(data_path + "/kappa.pdf", bbox_inches='tight')
    else:
        plt.ylabel(r'$\kappa^{\mathrm{c}}_{\mathrm{ph}}$ ($\mathrm{pw/K}$)', fontsize=17)

        plt.plot(E[E_min_index:E_max_index]*har2meV, kappa)
        plt.xlabel('Energy ($\mathrm{meV}$)',fontsize=17)
        # save kappa data
        top.write_plot_data(data_path + "/kappa_c.dat", (E[E_min_index:E_max_index], kappa), "E [har], kappa [pW/K]" + str(kappa_int_lower_E) + "->" + str(kappa_int_upper_E))
        plt.savefig(data_path + "/kappa_c.pdf", bbox_inches='tight')
    plt.show()
