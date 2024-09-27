import codecs
import configparser
import copy
import json
import sys
from functools import partial
import scipy.signal
import tmoutproc as top
import matplotlib

matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker

h_bar = 1.0545718 * 10 ** (-34)
eV2hartree = 0.0367493
ang2bohr = 1.88973
har2J = 4.35974E-18
bohr2m = 5.29177E-11
u2kg = 1.66054E-27
har2pJ = 4.35974e-6


def calculate_g0(w, w_D):
    """Calculates surface greens function according to Markussen, T. (2013). Phonon interference effects in molecular junctions. The Journal of chemical physics, 139(24), 244101 (https://doi.org/10.1063/1.4849178).

    Args:
    w (array_like): Frequency where g0 is calculated
	w_D (float): Debeye frequency

    Returns:
    g0	(array_like) Surface greens function g0
	"""

    def im_g(w):
        if (w <= w_D):
            Im_g = -np.pi * 3.0 * w / (2 * w_D ** 3)
        else:
            Im_g = 0
        return Im_g

    Im_g = map(im_g, w)
    Im_g = np.asarray(list(Im_g))
    Re_g = -np.asarray(np.imag(scipy.signal.hilbert(Im_g)))
    g0 = np.asarray((Re_g + 1.j * Im_g), complex)
    return g0


def calculate_Sigma(w, g0, gamma, M_L, M_C):
    """Calculates self energy according to Markussen, T. (2013). Phonon interference effects in molecular junctions. The Journal of chemical physics, 139(24), 244101  (https://doi.org/10.1063/1.4849178).

    Args:
    w (np.array): frequency
    g0 (np.array): g0
	gamma (float): gamma
	M_L (str): M_L atom type in reservoir
	M_C (str): M_C atom type coupled to reservoir

    Returns:
    sigma_nu (array_like) self energy term
	"""

    # convert to hartree/Bohr**2
    gamma_hb = gamma * (eV2hartree / ang2bohr ** 2)

    M_L = top.atom_weight(M_L, u2kg=False)
    M_C = top.atom_weight(M_C, u2kg=False)

    gamma_prime = gamma_hb / np.sqrt(M_C * M_L)

    g = g0 / (1 + gamma_prime * g0)
    sigma_nu = gamma_prime ** 2 * g

    return sigma_nu





def eval_propagator(e, eigenvalues, eigenvectors, r, s, direction="x", M=12, gamma=4):
    r"""
    Calculates zeroth order propagator  through terminal atoms:
     D_{lr}(E)=\sum_j \frac{C_{1j}C^\dagger_{rj}}{E^2-E_j^2+i\eta}
    Args:
        e: Energy
        eigenvalues: Eigenvalues of dynamical matrix
        eigenvectors: Eigenvectors of dynamical matrix
        r: Index of left terminal atom
        s: Index of right terminal atom

    Returns:
        |D|**2
    """
    G0_sr = 0
    delta = 1E-4
    r_ = r * 3
    s_ = s * 3
    if (direction == "x"):
        r_ = r_
        s_ = s_
    elif (direction == "y"):
        r_ = r_ + 1
        s_ = s_ + 1
    elif (direction == "z"):
        r_ = r_ + 2
        s_ = s_ + 2
    elif (direction == "xy"):
        r_ = r_
        s_ = s_ + 1
    elif (direction == "xz"):
        r_ = r_
        s_ = s_ + 2
    elif (direction == "yz"):
        r_ = r_ + 1
        s_ = s_ + 2
    r, s = r_, s_

    eV2J = 1.60218E-19
    Ang2m = 10E-10
    u2kg = 1.66054E-27
    # convert to hartree/Bohr**2
    gamma_hb = gamma * eV2hartree / ang2bohr ** 2
    # gamma = gamma * eV2J / Ang2m ** 2
    # M = M * u2kg
    h_bar = 1.0545718 * 10 ** (-34)
    J2meV = 6.24150934190e+21
    for i in range(0, len(eigenvalues)):
        numerator = np.dot(eigenvectors[r, i], eigenvectors[s, i])
        denominator = (e + 1.0j * delta) ** 2 - (eigenvalues[i]) ** 2
        tmp = numerator / denominator
        G0_sr += tmp

    D = np.abs(G0_sr)
    return D


def self_energy_correction(eigenvalues, eigenvectors, r, s, M=12, gamma=4):
    r"""
    Calculates zeroth order :
     D_{lr}(E)=\sum_j \frac{C_{1j}C^\dagger_{rj}}{E^2-E_j^2+i\eta}
    Args:
        e: Energy
        eigenvalues: Eigenvalues of dynamical matrix
        eigenvectors: Eigenvectors of dynamical matrix
        r: Index of left terminal atom
        s: Index of right terminal atom

    Returns:
        |D|**2
    """
    r = r * 3
    s = s * 3

    # convert to hartree/Bohr**2
    gamma_hb = gamma * eV2hartree / ang2bohr ** 2
    shifts = [
        np.sqrt(np.abs(gamma_hb) * (np.sum(eigenvectors[r:r + 3, i] ** 2) + np.sum(eigenvectors[s:s + 3, i] ** 2)) / M)
        for i in range(0, len(eigenvalues))]

    return shifts


def eval_propagator_inv(i, para):
    """Calculates Greens Function with given parameters at given frequency w.

	Args:
		i: (int): frequency index
		para: (tuple): frequency w (array), self energy sigma (complex), filename_hessian (str), filename_coord (str), left atom for transport calculation n_l (int), right atom for transport calculation n_r (int), coupling constant Gamma (complex), in_plane (boolean)

	Returns:
		P (array_like): phonon transmission
    """

    w = para[0]
    sigma = para[1]
    r = para[2]
    s = para[3]
    n_l = para[4]
    n_r = para[5]
    gamma = para[6]
    self_energy = para[7]
    self_energy_mode = para[8]
    D = para[9]
    D = copy.copy(D)

    n_atoms = int(D.shape[0] / 3)

    if (self_energy == True):
        # set up self energies
        sigma_L = np.zeros((n_atoms * 3, n_atoms * 3), complex)
        sigma_R = np.zeros((n_atoms * 3, n_atoms * 3), complex)
        # full
        if (self_energy_mode == 0):
            lower = 0
            upper = 3
        # only in plane
        elif (self_energy_mode == 1):
            lower = 1
            upper = 3
        # only out of plane
        elif (self_energy_mode == 2):
            lower = 0
            upper = 1
        # only x
        elif (self_energy_mode == 3):
            lower = 1
            upper = 2
        # only y
        elif (self_energy_mode == 4):
            lower = 2
            upper = 3
        else:
            exit(-1)

        for n_l_ in n_l:
            for u in range(lower, upper):
                sigma_L[n_l_ * 3 + u, n_l_ * 3 + u] = sigma[i]
        for n_r_ in n_r:
            for u in range(lower, upper):
                sigma_R[n_r_ * 3 + u, n_r_ * 3 + u] = sigma[i]
        sigma_i = sigma[i]

        # correct momentum conservation
        # convert to hartree/Bohr**2
        gamma_hb = gamma * eV2hartree / ang2bohr ** 2

        D_save = copy.deepcopy(D)
        for u in range(lower, 3):
            for n_l_ in n_l:
                # remove mass weighting
                K_ = D[n_l_ * 3 + u][n_l_ * 3 + u] * top.atom_weight(M_C)
                # correct momentum
                K_ = K_ - gamma_hb
                # add mass weighting again
                D_ = K_ / top.atom_weight(M_C)
                D[n_l_ * 3 + u][n_l_ * 3 + u] = D_

            for n_r_ in n_r:
                # remove mass weighting
                K_ = D[n_r_ * 3 + u][n_r_ * 3 + u] * top.atom_weight(M_C)
                # correct momentum
                K_ = K_ - gamma_hb
                # add mass weighting again
                D_ = K_ / top.atom_weight(M_C)
                D[n_r_ * 3 + u][n_r_ * 3 + u] = D_
        H = w[i] ** 2 * np.identity(3 * n_atoms) - D - sigma_L - sigma_R
        G = np.linalg.inv(H)
    else:
        G = np.linalg.inv(w[i] ** 2 * np.identity(3 * n_atoms) - D)

    return np.abs(G[r, s])


def eval_mode_type_all(dyn_eigenval, dyn_eigenvec):
    """

    Args:
        dyn_eigenval:
        dyn_eigenvec:

    Returns:

    """
    eigenvalues, eigenvectors = dyn_eigenval, dyn_eigenvec

    x_part = [eigenvectors[0::3, i] for i in range(0, len(eigenvalues))]
    x_part = np.array([np.linalg.norm(x_part[i]) ** 2 for i in range(0, len(eigenvalues))])

    y_part = [eigenvectors[1::3, i] for i in range(0, len(eigenvalues))]
    y_part = np.array([np.linalg.norm(y_part[i]) ** 2 for i in range(0, len(eigenvalues))])

    z_part = [eigenvectors[2::3, i] for i in range(0, len(eigenvalues))]
    z_part = np.array([np.linalg.norm(z_part[i]) ** 2 for i in range(0, len(eigenvalues))])

    in_plane = z_part + y_part
    out_of_plane = x_part
    return x_part, y_part, z_part


def eval_mode_type_atom(dyn_eigenval, dyn_eigenvec, atom):
    """

    Args:
        dyn_eigenval:
        dyn_eigenvec:

    Returns:

    """
    eigenvalues, eigenvectors = dyn_eigenval, dyn_eigenvec

    x_part = [eigenvectors[atom * 3 + 0, i] for i in range(0, len(eigenvalues))]
    x_part = np.array([np.linalg.norm(x_part[i]) ** 2 for i in range(0, len(eigenvalues))])

    y_part = [eigenvectors[atom * 3 + 1, i] for i in range(0, len(eigenvalues))]
    y_part = np.array([np.linalg.norm(y_part[i]) ** 2 for i in range(0, len(eigenvalues))])

    z_part = [eigenvectors[atom * 3 + 2, i] for i in range(0, len(eigenvalues))]
    z_part = np.array([np.linalg.norm(z_part[i]) ** 2 for i in range(0, len(eigenvalues))])

    in_plane = z_part + y_part
    out_of_plane = x_part

    return x_part / (x_part + y_part + z_part), y_part / (x_part + y_part + z_part), z_part / (x_part + y_part + z_part)


if __name__ == '__main__':
    config_path = sys.argv[1]
    calc_path = sys.argv[2]

    cfg = configparser.ConfigParser()
    # cfg.read(config_path)
    cfg.read_file(codecs.open(config_path, "r", "utf8"))

    try:
        data_path = str(cfg.get('Data Input', 'data_path'))
        hessian_name = str(cfg.get('Data Input', 'hessian_name'))
        coord_name = str(cfg.get('Data Input', 'coord_name'))
        filename_hessian = data_path + "/" + hessian_name
        filename_coord = data_path + "/" + coord_name

        # atoms which are coupled to the leads -> self energy
        n_l = np.asarray(str(cfg.get('Calculation', 'n_l')).split(','), dtype=int)
        n_r = np.asarray(str(cfg.get('Calculation', 'n_r')).split(','), dtype=int)

        # atom type in resevoir M_L and molecule M_C
        M_L = str(cfg.get('Calculation', 'M_L'))
        M_C = str(cfg.get('Calculation', 'M_C'))
        # coupling force constant resevoir in eV/Ang**2
        gamma = float(cfg.get('Calculation', 'gamma'))

        # Debeye energy in meV
        E_D = float(cfg.get('Calculation', 'E_D'))
        # Number of grid points
        N = int(cfg.get('Calculation', 'N'))
        # only in plane motion (-> set x and y coupling to zero)
        in_plane = json.loads(str(cfg.get('Calculation', 'in_plane')).lower())
        # mode = int(cfg.get('Calculation', 'mode'))

        # for thermal conducatance
        T_min = float(cfg.get('Calculation', 'T_min'))
        T_max = float(cfg.get('Calculation', 'T_max'))
        kappa_grid_points = int(cfg.get('Calculation', 'kappa_grid_points'))

        # check if g0 should be plotted
        plot_g0 = json.loads(str(cfg.get('Data Output', 'plot_g')).lower())

    except configparser.NoOptionError as e:
        print("Missing option in config file. Check config file!")
        print(e)
        exit(-1)
    except ValueError:
        print("Wrong value in config file. Check config file!")
        exit(-1)

    """
    K = top.create_dynamical_matrix(filename_hessian, filename_coord, t2SI=True)
    K_eigenval,K_eigenvec = np.linalg.eigh(K)
    idx = np.argsort(np.real(K_eigenval))
    eigenValues = K_eigenval[idx]
    K_eigenvec = K_eigenvec[:, idx]
    K_eigenval[K_eigenval < 0] = 0
    K_eigenval = np.array(np.sqrt(K_eigenval))
    K_eigenval = K_eigenval*hbar*6.241506363094e+21
    print(K_eigenval)
    E = np.linspace(1.0,E_D*1.25,600)
    """
    K = top.create_dynamical_matrix(filename_hessian, filename_coord, t2SI=False)
    K_eigenval, K_eigenvec = np.linalg.eigh(K)
    idx = np.argsort(np.real(K_eigenval))
    K_eigenval = K_eigenval[idx]
    K_eigenvec = K_eigenvec[:, idx]
    K_eigenval[K_eigenval < 0] = 0
    K_eigenval = np.array(np.sqrt(K_eigenval))
    # if(self_energy == true):
    #    K_eigenval, K_eigenvec = self_energy_correction()
    # convert to J
    E_D = 1.0 * E_D
    E_D = E_D * 1.60217656535E-22
    # convert to 1/s
    w_D = E_D / h_bar
    # convert to har/(bohr**2*u)
    w_D = w_D / np.sqrt(9.375821464623672e+29)
    N = 1000
    w = np.linspace(w_D * 1E-2, w_D * 1.1, N)
    i = np.linspace(0, N, N, False, dtype=int)
    E = w * np.sqrt(9.375821464623672e+29) * h_bar / (1.60217656535E-22)
    directions = "x"


    r = n_l[0]
    s = n_r[0]


    g0 = calculate_g0(w, w_D)
    Sigma = calculate_Sigma(w, g0, gamma, M_L, M_C)
    self_energies = [False, True]

    for self_energy in self_energies:

        figure, ax = plt.subplots(figsize=(6, 3))

        params = w, Sigma, r * 3, s * 3, [r], [s], gamma, self_energy, 0, K
        D_xx = np.array(list(map(partial(eval_propagator_inv, para=params), i)))

        params = w, Sigma, r * 3 + 1, s * 3 + 1, [r], [s], gamma, self_energy, 0, K
        D_yy = np.array(list(map(partial(eval_propagator_inv, para=params), i)))

        params = w, Sigma, r * 3 + 2, s * 3 + 2, [r], [s], gamma, self_energy, 0, K
        D_zz = np.array(list(map(partial(eval_propagator_inv, para=params), i)))

        params = w, Sigma, r * 3, s * 3 + 1, [r], [s], gamma, self_energy, 0, K
        D_xy = np.array(list(map(partial(eval_propagator_inv, para=params), i)))

        params = w, Sigma, r * 3, s * 3 + 2, [r], [s], gamma, self_energy, 0, K
        D_xz = np.array(list(map(partial(eval_propagator_inv, para=params), i)))

        params = w, Sigma, r * 3 + 1, s * 3 + 2, [r], [s], gamma, self_energy, 0, K
        D_yz = np.array(list(map(partial(eval_propagator_inv, para=params), i)))

        top.write_plot_data(calc_path + "/propagator_Sigma_{}.dat".format(self_energy),
                            (w, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz), "(w, D_xx, D_yy, D_zz, D_xy,D_xz, D_yz)")


        ax.grid(which='minor', alpha=0.7)
        ax.grid(which='major', alpha=0.85)

        ax.plot(E, D_xx, label="xx", alpha=1.0, color=(1.0, 0, 0))
        ax.plot(E, D_yy, label="yy", alpha=1.0, color=(0, 1.0, 0))
        ax.plot(E, D_zz, label="zz", alpha=1.0, color=(0, 0.0, 1.0))
        ax.plot(E, D_xy, label="xy", alpha=1.0, color=(0.6, 0.6, 0))
        ax.plot(E, D_xz, label="xz", alpha=1.0, color=(0.6, 0, 0.6))
        ax.plot(E, D_yz, label="yz", alpha=1.0, color=(0, 0.6, 0.6))

        min = np.min([D_xx, D_yy, D_zz, D_xy, D_xz, D_yz])
        max = np.max([D_xx, D_yy, D_zz, D_xy, D_xz, D_yz])

        x_part_atom, y_part_atom, z_part_atom = eval_mode_type_atom(K_eigenval, K_eigenvec, 0)
        x_part, y_part, z_part = eval_mode_type_all(K_eigenval, K_eigenvec)
        E_LIM = 20

        if (self_energy == True):
            shifts = self_energy_correction(K_eigenval, K_eigenvec, r, s, 12, gamma=gamma)

            K_eigenval = K_eigenval - shifts
        else:
            for j in range(6, len(K_eigenval)):
                if (K_eigenval[j] * np.sqrt(9.375821464623672e+29) * h_bar / (1.60217656535E-22) < np.max(E_LIM)):
                    plt.axvline(K_eigenval[j] * np.sqrt(9.375821464623672e+29) * h_bar / (1.60217656535E-22),
                                color=(x_part[j], y_part[j], z_part[j]), alpha=0.5)
                    if (j % 2 == 0):
                        plt.text(K_eigenval[j] * np.sqrt(9.375821464623672e+29) * h_bar / (1.60217656535E-22) + 0.1,
                                 max * 0.01, str(j), rotation=90,
                                 color=(x_part_atom[j], y_part_atom[j], z_part_atom[j]))
                        pass
                    else:
                        plt.text(K_eigenval[j] * np.sqrt(9.375821464623672e+29) * h_bar / (1.60217656535E-22) + 0.1,
                                 max * 0.1, str(j), rotation=90, color=(x_part_atom[j], y_part_atom[j], z_part_atom[j]))
                        pass
                else:
                    break
        ax.set_yscale("log")

        ax.set_xlim(0, E_LIM)
        index_E_LIM = np.argmin(np.abs(E - E_LIM))

        min = np.min(
            [D_xx[0:index_E_LIM], D_yy[0:index_E_LIM], D_zz[0:index_E_LIM], D_xy[0:index_E_LIM], D_xz[0:index_E_LIM],
             D_yz[0:index_E_LIM]])
        max = np.max(
            [D_xx[0:index_E_LIM], D_yy[0:index_E_LIM], D_zz[0:index_E_LIM], D_xy[0:index_E_LIM], D_xz[0:index_E_LIM],
             D_yz[0:index_E_LIM]])

        ax.set_ylim(min * 0.95, max * 1.05)
        ax.set_ylim(1E-3, 1E5)

        ax.legend(fontsize=14, ncol=3, loc="lower left")
        ax.tick_params(axis='y', labelsize=20)
        ax.tick_params(axis='x', labelsize=20)
        ax.set_xlabel('Phonon Energy ($\mathrm{meV}$)', fontsize=20)
        ax.set_ylabel(r'$|\mathrm{G_{lr}}|$ ($1/\mathrm{har}^2$)', fontsize=20)
        plt.rc('ytick', labelsize=15)
        plt.rc('xtick', labelsize=15)

        ## set y ticks
        y_major = matplotlib.ticker.LogLocator(base=10.0, numticks=5)
        ax.yaxis.set_major_locator(y_major)
        y_minor = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(1.0, 10.0) * 0.1, numticks=10)
        ax.yaxis.set_minor_locator(y_minor)
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        plt.savefig(calc_path + "/propagator_Sigma_{}.pdf".format(self_energy), bbox_inches='tight', transparent=True)
        plt.clf()
