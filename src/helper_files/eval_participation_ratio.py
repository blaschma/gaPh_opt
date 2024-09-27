import configparser
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tmoutproc as top
from utils import find_atoms

hbar = 1.054571817E-34


def eval_participation_ratio(dyn_eigenval, dyn_eigenvec, start_index, end_index):
    """
    Evaluation of participation ratio for degrees of freedom from start_index to end_index
    Args:
        dyn_eigenval: Eigenvalues of dynamical matrix
        dyn_eigenvec: Eigenvectors of dymnamical matrix
        start_index:
        end_index:

    Returns:
        PR
    """
    eigenvalues, eigenvectors = dyn_eigenval, dyn_eigenvec
    # extract degrees of freedom
    PR = [eigenvectors[start_index:end_index, i] for i in range(0, len(eigenvalues))]
    PR = [np.linalg.norm(PR[i]) ** 2 for i in range(0, len(eigenvalues))]
    return PR


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
    return x_part, y_part, z_part


def plot_mode_type_all(freq, x_part, y_part, z_part, filename):
    """

    Args:
        freq:
        in_plane:
        out_of_plane:

    Returns:

    """
    xlim = 100
    argmin = (np.argmin(np.abs(np.array(freq) - xlim)))
    plt.bar(freq[0:argmin], x_part[0:argmin] + y_part[0:argmin] + z_part[0:argmin], label="z part", align='center',
            width=1.0)
    plt.bar(freq[0:argmin], x_part[0:argmin] + z_part[0:argmin], label="y part", align='center', width=1.0)
    plt.bar(freq[0:argmin], x_part[0:argmin], label="x part", align='center', width=1.0)

    plt.legend()
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()
    plt.close()


def plot_mode_type_two_side(freq, E_D, x_part_1, y_part_1, z_part_1, x_part_2, y_part_2, z_part_2, filename):
    """

    Args:
        freq:
        in_plane:
        out_of_plane:

    Returns:

    """
    xlim = E_D
    argmin = (np.argmin(np.abs(np.array(freq) - xlim))) + 1
    width = 0.2
    plt.bar(freq[0:argmin], x_part_1[0:argmin] + y_part_1[0:argmin] + z_part_1[0:argmin], label="z 1", align='center',
            width=width, color="r")
    plt.bar(freq[0:argmin], x_part_1[0:argmin] + z_part_1[0:argmin], label="y 1", align='center', width=width,
            color="g")
    plt.bar(freq[0:argmin], x_part_1[0:argmin], label="x 1", align='center', width=width, color="b")

    plt.bar(freq[0:argmin], -x_part_2[0:argmin] - y_part_2[0:argmin] - z_part_2[0:argmin], label="z 2", align='center',
            width=width, color="r")
    plt.bar(freq[0:argmin], -x_part_2[0:argmin] - z_part_2[0:argmin], label="y 2", align='center', width=width,
            color="g")
    plt.bar(freq[0:argmin], -x_part_2[0:argmin], label="x 2", align='center', width=width, color="b")

    plt.axhline(0, color="black")
    plt.legend()
    plt.xlim(0, E_D * 1.2)
    plt.yscale("symlog")
    plt.xlabel('Phonon Energy ($\mathrm{meV}$)', fontsize=19)
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()
    plt.close()


def plot_participation_ratio(eigenvalues, PR, E_D, filename):
    xlim = E_D
    fig, ax = plt.subplots()
    ax.plot(eigenvalues, np.array(PR), lw=0.3, alpha=0.5, color="b")
    ax.plot(eigenvalues, np.array(PR), marker="x", color="red", lw=0)

    # write_out
    top.write_plot_data(filename.replace(".pdf", ".dat"), (eigenvalues, PR), "EV, PR")

    ax.set_xlabel('Phonon Energy ($\mathrm{meV}$)', fontsize=19)
    ax.set_ylabel(r'PR', fontsize=19)
    ax.set_xlim(0, xlim)
    argmin = (np.argmin(np.abs(np.array(eigenvalues) - xlim)))
    ylim = np.max(PR[0:argmin])
    ax.set_ylim(0, ylim)
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()
    plt.close()


if __name__ == '__main__':
    config_path = sys.argv[1]
    calc_path = sys.argv[2]
    filename_hessian = calc_path + "/hessian"
    filename_coord = calc_path + "/coord"
    coord = top.read_coord_file(filename_coord)

    cfg = configparser.ConfigParser()
    cfg.read(config_path)
    E_D = float(cfg.get('Phonon Calculation', 'E_D'))

    K = top.create_dynamical_matrix(filename_hessian, filename_coord, t2SI=False)
    K_eigenval, K_eigenvec = np.linalg.eig(K)
    K_eigenval[K_eigenval < 0] = 0
    K_eigenval = np.array(np.sqrt(sorted(K_eigenval)))
    K_eigenval = K_eigenval * np.sqrt(9.375821464623672e+29) * hbar / (1.60217656535E-22)

    PR = [eval_participation_ratio(K_eigenval, K_eigenvec, i * 3, i * 3 + 3) for i in
          range(0, int(len(K_eigenval) / 3))]
    for i in range(0, len(PR)):
        atom = coord[3, i]
        plot_participation_ratio(sorted(K_eigenval), PR[i], E_D, calc_path + "/PR_" + str(i) + "_" + atom + ".pdf")

    x_part, y_part, z_part = eval_mode_type_all(K_eigenval, K_eigenvec)
    plot_mode_type_all(K_eigenval, x_part, y_part, z_part, calc_path + "/mode_type.pdf")
    coord_xyz = top.t2x(coord)
    r, s = find_atoms(coord_xyz)
    x_part_1, y_part_1, z_part_1 = eval_mode_type_atom(K_eigenval, K_eigenvec, atom=r)
    x_part_2, y_part_2, z_part_2 = eval_mode_type_atom(K_eigenval, K_eigenvec, atom=s)
    plot_mode_type_two_side(K_eigenval, E_D, x_part_1, y_part_1, z_part_1, x_part_2, y_part_2, z_part_2,
                            calc_path + "/mode_type_terminal_atom.pdf")
