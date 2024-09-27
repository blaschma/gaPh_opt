import sys
from configparser import ConfigParser
import tmoutproc as top
from utils import find_atoms

if __name__ == '__main__':
    xyz_file = str(sys.argv[15])
    xyz_file = top.read_xyz_file(xyz_file)


    config_path = str(sys.argv[1])
    data_path = str(sys.argv[2])
    hessian_name = str(sys.argv[3])
    coord_name = str(sys.argv[4])
    M_L = str(sys.argv[5])
    M_C = str(sys.argv[6])
    gamma = str(sys.argv[7])
    E_D = str(sys.argv[8])
    N = str(sys.argv[9])
    in_plane = str(sys.argv[10])
    T_min = str(sys.argv[11])
    T_max = str(sys.argv[12])
    kappa_grid_points = str(sys.argv[13])
    anchor_mode = int(sys.argv[14])

    if(anchor_mode==0):
        gold_index = find_atoms(xyz_file)
    elif(anchor_mode==1):
        gold_index = find_atoms(xyz_file, "s")
    else:
        raise ValueError('Unknown anchor mode')

    config = ConfigParser()
    config.optionxform = str
    config.add_section('Data Input')
    config.set('Data Input', 'data_path', data_path)
    config.set('Data Input', 'hessian_name', hessian_name)
    config.set('Data Input', 'coord_name', coord_name)

    config.set('Data Input', 'transp_name', "phonon_trans.dat")
    config.set('Data Input', 'transp_units', "har/(bohr**2*u)")

    config.add_section('Calculation')
    config.set('Calculation', 'n_l', str(gold_index[0]))
    config.set('Calculation', 'n_r', str(gold_index[1]))
    config.set('Calculation', "M_L", M_L)
    config.set('Calculation', 'M_C', M_C)
    config.set('Calculation', 'gamma', gamma)
    config.set('Calculation', 'E_D', E_D)
    config.set('Calculation', 'N', N)
    config.set('Calculation', 'in_plane', in_plane)
    config.set('Calculation', 'T_min', T_min)
    config.set('Calculation', 'T_max', T_max)
    config.set('Calculation', 'kappa_grid_points', kappa_grid_points)

    config.set('Calculation', 'kappa_int_lower_E', str(0))
    config.set('Calculation', 'kappa_int_upper_E', str(E_D))



    config.add_section('Data Output')
    config.set('Data Output', 'plot_g', "False")

    # save to a file

    with open(config_path, 'w') as configfile:
        config.write(configfile)