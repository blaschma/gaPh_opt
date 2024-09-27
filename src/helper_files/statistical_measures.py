import configparser
import subprocess
import tmoutproc as top
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import os
import sys
from rdkit import Chem
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
import math



def entropy(string):
    "Calculates the Shannon entropy of a string"

    # get probability of chars in string
    prob = [ float(string.count(c)) / len(string) for c in dict.fromkeys(list(string)) ]

    # calculate the entropy
    entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob ])

    return entropy


def entropy_ideal(length):
    "Calculates the ideal Shannon entropy of a string with given length"

    prob = 1.0 / length

    return -1.0 * length * prob * math.log(prob) / math.log(2.0)


import genome_to_molecule as gtm
from skspatial.objects import Plane, Point, Points


plt.style.use(['science','no-latex'])
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=16)
#plt.rcParams["figure.figsize"] = (3, 3)
plt.rcParams["figure.figsize"] = (3, 4)

"""
This script contains functions to analyze the statistical properties of an evolutionary run.
"""



def count_elements(array):
    """
    Count the number of occurrences of each element in an array
    Args:
        array: np.array

    Returns:
        count_dict: Dictionary with element as key and number of occurrences as value
    """
    count_dict = {}
    for element in array:
        if element in count_dict:
            count_dict[element] += 1
        else:
            count_dict[element] = 1
    return count_dict

def add_dicts(dict1, dict2):
    """
    Add two dictionaries
    Args:
        dict1: Dictionary 1
        dict2: Dictionary 2

    Returns:
        dict1: Dictionary 1 with added values from dictionary 2
    """
    for key, value in dict2.items():
        if key in dict1:
            dict1[key] += value
        else:
            dict1[key] = value
    return dict1



def coord_analysis(path, n_best, blocks, building_blocks=None):
    """
    Analyze properties deducted from coord file for generation in path for n_best individuals
    Args:
        path: Path to generation
        n_best: Upper limit of individuals to process
        building_blocks: Building blocks (full information)

    Returns:
        element_count_dict: Dictionary with element as key and number of occurrences as value
        dihedral_angles: List of dihedral angles
    """
    #load individuals in generation
    individuals = os.listdir(f"{path}")
    individuals = [int(individual) for individual in individuals if os.path.isdir(f"{path}/{individual}")]
    individuals = sorted(individuals)[0:n_best]

    element_count_dict = {}
    dihedral_angles = []
    synthetic_accessibilities = []

    #analyze individuals
    for i, individual in enumerate(individuals):

        if os.path.exists(f"{path}/{individual}/HESSIAN_NOT_CONVERGED") or os.path.exists(f"{path}/{individual}/RELAXATION_NOT_CONVERGED"):
            dihedral_angles.append([])
            synthetic_accessibilities.append(10)
            continue

        try:
            coord = top.read_xyz_file(f"{path}/{individual}/xtbopt.xyz")
        except FileNotFoundError as e:
            print(f"File {path}/{individual}/xtbopt.xyz not found")
            dihedral_angles.append([])
            synthetic_accessibilities.append(10)
            continue
        coord = top.x2t(coord)
        #analyze blocks
        block_informations = [building_blocks[block] for block in blocks[individual]]
        coord_without_anchor = coord[:,2:coord.shape[1]-2]

        #analyze syntetic accessibility
        #todo: proper handling
        if(i == 0):
            sdf_file = f"{path}/{individual}/xtbopt_without_au.sdf"
            if os.path.exists(sdf_file) == False:
                #write coord without anchor to file
                coord_without_gold =  coord.T[coord.T[:, 3] != 'au'].T
                coord_without_gold = top.t2x(coord_without_gold)
                top.write_xyz_file(f"{path}/{individual}/xtbopt_without_au.xyz", coord_without_gold)
                #todo: Do not hard code the path to open babel
                command = f"PATH_TO_OPEN_BABEL/obabel {path}/{individual}/xtbopt_without_au.xyz -O {sdf_file}"
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)

            supplier = Chem.SDMolSupplier(sdf_file)
            mol = supplier[0]
            synthetic_accessibility = sascorer.calculateScore(mol)
            synthetic_accessibilities.append(synthetic_accessibility)


        coord_fragments = []
        lower_index = 0
        dihedral_angle = 0
        dihedral_angle_individual = []
        angle_reasobale_list = []

        #fragmentation of coord file based in genetic information, calculate dihedral angles (absolute values)
        for i, block_information in enumerate(block_informations):

            angle_reasobale = True
            upper_index = block_informations[i].num_atoms
            coord_fragment = coord_without_anchor[:,lower_index:lower_index+upper_index]
            coord_fragments.append(coord_fragment)
            if(block_informations[i].para_pos == block_informations[i].meta_pos):
                angle_reasobale = False
            angle_reasobale_list.append(angle_reasobale)
            if(len(coord_fragments)>1 and angle_reasobale_list[-1]==True and angle_reasobale_list[-2]==True):

                angle_ = dihedral_angle_analsis(coord_fragments[-1], coord_fragments[-2])
                dihedral_angle += np.abs(angle_)
                dihedral_angle_individual.append(angle_)
            lower_index += upper_index

        dihedral_angles.append(dihedral_angle_individual)


        element_count_dict_ = count_elements(coord[3,:])
        element_count_dict = add_dicts(element_count_dict, element_count_dict_)
    return element_count_dict, dihedral_angles, synthetic_accessibilities

def dihedral_angle_analsis(coord1, coord2):
    """
    Calculate dihedral angle between two fragments. The fragments are fitted to a plane and the angle between the two
    planes is calculated. Planes are fitted through the center of mass of the fragments. For the calculation of the
    center of mass, only carbon and hydrogen atoms are considered.
    Args:
        coord1:
        coord2:

    Returns:

    """

    def fit_plane(coord):
        point_list = list()
        mass_list = list()
        for index in range(0,coord.shape[1]):
            type = coord[3, index]
            mass = top.constants.ATOM_DICT_SYM[type][2]
            x = coord[0,index]
            y = coord[1,index]
            z = coord[2,index]
            if (type == "c" or type == "h"):
                point_list.append([x, y, z])
                mass_list.append(mass)
        # force through center of mass
        center_of_mass = np.average(point_list, axis=0, weights=mass_list)
        shifted_points = point_list - center_of_mass


        points = Points(shifted_points)

        try:
            plane = Plane.best_fit(points)
            plane.point = Point(center_of_mass)
            normal = plane.normal
            assert plane.distance_point_signed(Point(center_of_mass)) < 1e-6, "Plane does not go through center of mass"
        except ValueError as e:
            normal = [0,0,0]


        return normal

    normal_1 = fit_plane(coord1)
    normal_1 = normal_1/np.linalg.norm(normal_1)
    normal_2 = fit_plane(coord2)
    normal_2 = normal_2 / np.linalg.norm(normal_2)
    dot_product = np.dot(normal_1, normal_2)
    if(dot_product<0):
        dot_product = dot_product * -1
    angle = np.arccos(dot_product)
    angle = angle/(np.pi*2)*360
    return angle



def plot_values_from_list_of_dicts(dict_list, path, keys_to_exclude={}, file_name="coord_analysis", x_label="Generation", y_label="Value", label_dict={}):
    """
    Plot values from list of dictionaries
    Args:
        dict_list: List of dictionaries
        path: Path where to save the plot
        keys_to_exclude: Keys to exclude from plotting
        file_name: File name of the plot

    Returns:

    """
    keys = set()
    for d in dict_list:
        keys.update(d.keys())

    # Prepare data for plotting
    plot_data = {key: [] for key in keys}
    n_total = []
    for d in dict_list:
        n_total_ = sum(value for key, value in d.items() if key not in keys_to_exclude)
        n_total.append(n_total_)
        for key in keys:
            plot_data[key].append(d.get(key, 0))

    fig, ax1 = plt.subplots()
    n_total = np.array(n_total)

    num_curves = len(keys)
    print(f"num curves {num_curves} for {file_name}")
    colormap = cm.get_cmap('nipy_spectral', num_curves)

    counter = 0
    for key, values in plot_data.items():
        if(key not in keys_to_exclude):
            if(len(label_dict) == 0):
                label = key
            else:
                label = label_dict[key]
            ax1.plot(values/n_total, label=label, color=colormap(counter))
            counter += 1
            #ax1.plot(values, label=key)

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax1.legend()
    #sort legend
    handles, labels = plt.gca().get_legend_handles_labels()
    labels = [int(label) if label.isdigit() else label for label in labels]
    sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: x[1])
    sorted_handles, sorted_labels = zip(*sorted_handles_labels)
    plt.legend(sorted_handles, sorted_labels, frameon=True, fancybox=True)

    plt.savefig( f"{path}/{file_name}.pdf", bbox_inches='tight')
    plt.savefig(f"{path}/{file_name}.svg", bbox_inches='tight')

def analyze_coupling_classes(list_of_lists):
    """
    Analyze coupling classes
    Args:
        list_of_lists: List of coupling lists

    Returns:
        class_dict: Dictionary with the frequency of number of meta couplings for given generation
    """


    # determine maximum number of couplings
    max_len = 0
    for individual in list_of_lists:
        len_ = len(individual)
        if(len_ > max_len):
            max_len = len_

    class_dict_list = []
    # determine class for each dict in every generation
    class_dict = {}
    class_dict[0] = 0
    for individual in list_of_lists:
        count_dict = count_elements(individual)
        if 1 in count_dict:
            if count_dict[1] in class_dict:
                class_dict[count_dict[1]] += 1
            else:
                class_dict[count_dict[1]] = 1
        else:
            class_dict[0] += 1

    #sum all keys in class_dict to check if all individuals were labeled
    sum = 0
    for key, value in class_dict.items():
        sum += value
    assert sum == len(list_of_lists), "Not all individuals were labeled"

    return class_dict

def read_population(path):
    """
    Read the population from the path and return the blocks and couplings
    Args:
        path:

    Returns:
        blocks: List of blocks (list of lists)
        couplings: List of couplings (list of ints)
    """
    blocks = []
    couplings = []
    raw_generatic_data = []
    with open(f"{path}/summary.dat", "r") as file:
        for line in file:
            tmp = line.strip().split(", ")
            tmp = [str(tmp[i]).replace("'", "").replace(" ", "") for i in range(0, len(tmp))]
            raw_generatic_data_ = ''.join(tmp)
            raw_generatic_data.append(raw_generatic_data_)
            blocks_ = tmp[1::2]
            blocks_ = [int(block.split('#', 1)[0]) for block in blocks_]
            blocks_ = np.asarray(blocks_, dtype=int)
            blocks.append(blocks_)

            couplings_ = tmp[0::2]
            couplings_ = [int(coupling) for coupling in couplings_]
            couplings_ = np.asarray(couplings_, dtype=int)
            couplings.append(couplings_)
    return blocks, couplings, raw_generatic_data

def analyze_end_groups(blocks, n_best):
    """
    Analyze the end groups of the blocks. The first n_best individuals are considered
    Args:
        blocks: List of Lists

    Returns:
        starts: List of starts
        ends: List of ends
    """
    starts = []
    ends = []
    for individual_blocks in blocks[0:n_best]:
        start = individual_blocks[0]
        end = individual_blocks[-1]
        starts.append(start)
        ends.append(end)
    return starts, ends

def plot_values_from_list_of_lists_containing_a_list(list_of_lists, path, file_name="angle_analysis", x_label="Generation", y_label="Value"):
    """
    Plot values from list of lists containing a list
    Args:
        list_of_lists: List of lists. Each list contains a list
    """
    fig, ax1 = plt.subplots()
    values_to_plot = [np.sum(value[0]) for value in list_of_lists]
    sums = []
    for list in list_of_lists:
        sums_sublist = [np.sum(sub_list) for i, sub_list in enumerate(list)]
        sums.append(sums_sublist)
    assert np.all([values_to_plot[i] == sums[i][0] for i in range(0, len(sums))]), "Wrong sum"

    std_deviation_values = [np.std(sums[i]) for i in range(0, len(sums))]

    deviation_values = [np.mean(value[1:]) for value in list_of_lists]
    x_values = np.arange(len(deviation_values))


    plt.errorbar(x_values, values_to_plot, yerr=std_deviation_values, ecolor='orangered')
    top.write_plot_data(f"{path}/{file_name}.dat", [x_values, values_to_plot, std_deviation_values])

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(f"{path}/{file_name}.pdf", bbox_inches='tight')
    plt.savefig(f"{path}/{file_name}.svg", bbox_inches='tight')

def plot_values_from_list_of_lists(list_of_lists, path, file_name="angle_analysis", x_label="Generation", y_label="Value"):
    """
    Plot the values from a list of lists. List of lists is a list of lists of values
    Args:
        list_of_lists: List of lists
    """
    fig, ax1 = plt.subplots()
    values_to_plot = [list[0] for list in list_of_lists]
    x_values = np.arange(len(values_to_plot))

    plt.plot(x_values,values_to_plot)

    top.write_plot_data(f"{path}/{file_name}.dat", [x_values, values_to_plot])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(f"{path}/{file_name}.pdf", bbox_inches='tight')
    plt.savefig(f"{path}/{file_name}.svg", bbox_inches='tight')

def analyze_shanon_entropy(genetic_string, n_best):
    """
    Analyze the shanon entropy of the genetic strings in generation
    Args:
        genetic_string: List of strings
        n_best: Number of individuals to consider

    Returns:
        shanon_entropy: List of shanon entropies
    """
    shanon_entropies = []
    for i in range(0, n_best):
        shanon_entropy = entropy(genetic_string[i])
        shanon_entropies.append(shanon_entropy)
    return shanon_entropies



if __name__ == '__main__':
    # sys.argv[1] path to process
    # sys.argv[2] config path

    calculation_path = sys.argv[1]
    config_path = sys.argv[2]

    cfg = configparser.ConfigParser()
    cfg.read(config_path, encoding='utf-8')
    building_block_path = cfg.get('Building Procedure', 'building_block_path')

    building_blocks = gtm.load_building_blocks(building_block_path)

    generations = os.listdir(f"{calculation_path}/generation_data")
    generations = [int(generations) for generations in generations if os.path.isdir(f"{calculation_path}/generation_data/{generations}")]
    generations = sorted(generations)
    element_count_dicts = []
    couplings_stat_concat = []
    coupling_class_dict_list = []
    starts_stat = []
    ends_stat = []
    dihedral_angles = []
    synthetic_accessibility = []
    raw_genetic_data = []
    shanon_entropies = []
    n_best = 20
    for generation in generations:
        print(generation)
        generation_path = f"{calculation_path}/generation_data/{generation}"

        #extract blocks and couplings from individuals in generation
        blocks, couplings, raw_genetic_data_ = read_population(generation_path)
        raw_genetic_data.append(raw_genetic_data_)
        starts, ends = analyze_end_groups(blocks, n_best)
        #start stat
        starts_stat_ = count_elements(starts)
        starts_stat.append(starts_stat_)
        #end stat
        ends_stat_ = count_elements(ends)
        ends_stat.append(ends_stat_)
        #couplings stat concat

        coupling_class_dict = analyze_coupling_classes(couplings[0:n_best])
        coupling_class_dict_list.append(coupling_class_dict)
        couplings_concat = np.concatenate(couplings[0:n_best])
        couplings_stat_ = count_elements(couplings_concat)
        couplings_stat_concat.append(couplings_stat_)
        #entropy stat
        shanon_entropies_ = analyze_shanon_entropy(raw_genetic_data_, n_best)
        shanon_entropies.append(shanon_entropies_)

        #analyze properties from coord file
        element_count_dict, dihedral_angles_, synthetic_accessibility_ = coord_analysis(generation_path, n_best, blocks, building_blocks)
        element_count_dicts.append(element_count_dict)
        dihedral_angles.append(dihedral_angles_)
        synthetic_accessibility.append(synthetic_accessibility_)



    couplings_label_dict = {0: 'para', 1: 'meta'}
    blocks_label_dict = {0: 3, 1: 4, 2: 5, 3: 1, 4: 8, 5: 6, 6: 9, 7: 10, 8: 2, 9: 7}
    subs_label_dicts = {'h': 'H', 'f': 'F', 'c': 'C', 'br': 'Br', 'i': 'I', 'cl': 'Cl'}
    plot_values_from_list_of_dicts(starts_stat, calculation_path,keys_to_exclude = [], file_name="starts_statistics", y_label="Relative Frequency", label_dict=blocks_label_dict)
    plot_values_from_list_of_dicts(ends_stat, calculation_path, keys_to_exclude=[], file_name="ends_statistics", y_label="Relative Frequency", label_dict=blocks_label_dict)
    plot_values_from_list_of_dicts(couplings_stat_concat, calculation_path, keys_to_exclude=[], file_name="couplings_statistics", y_label="Relative Frequency", label_dict=couplings_label_dict)
    plot_values_from_list_of_dicts(coupling_class_dict_list, calculation_path, keys_to_exclude=[], file_name="coupling_class_statistics", y_label="Relative Frequency")

    #todo: check if this is correct
    plot_values_from_list_of_dicts(element_count_dicts, calculation_path,keys_to_exclude = ['s', 'c', 'au'], file_name="subs_statistics", y_label="Relative Frequency", label_dict=subs_label_dicts)

    plot_values_from_list_of_lists_containing_a_list(dihedral_angles, calculation_path, file_name="dihedral_angles", y_label="Cumulative Dihedral Angle (°)", x_label="Generation")
    plot_values_from_list_of_lists_containing_a_list(shanon_entropies, calculation_path, file_name="shanon_entropy", y_label="Shannon Entropy", x_label="Generation")
    plot_values_from_list_of_lists(synthetic_accessibility, calculation_path, file_name="synthetic_accessibility", y_label="Synthetic Accessibility", x_label="Generation")