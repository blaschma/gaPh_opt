import numpy as np


def find_atoms(coord, atom="au"):
    """
    Finds index of gold atoms in coord file stored in coord. index_left<index_right from building procedure. Other atoms can be found as well. The first two occurrences are chosen.

    Args:
        param1 (np.ndarray): coord file loaded with top.read_coord_file
    Returns:
        (index_left, index_right)
    """

    index_1 = -1
    index_2 = -1
    for i in range(0, coord.shape[1]):
        # if line is au atom
        if ((coord[0, i]).lower() == atom):
            if (index_1 == -1):
                index_1 = i
            if (index_1 != -1):
                index_2 = i
    if (index_1 != -1 and index_2 != -1):
        min = np.min((index_1, index_2))
        max = np.max((index_1, index_2))
        return (min, max)
    else:
        print("atoms " + atom + " not found")
        return (-1, -1)


def read_generation(config_path, calculation_path):
    """
    Invokes the next generation. The number of generation is read from calculation_path/generation.dat

    Args:
        param1 (String): Path to config file
        param2 (String): Path to calculation

    Returns:

    """
    try:
        filename = calculation_path + "/generation.dat"
        file = open(filename)
        for line in file:
            generation = int(line)

    except (OSError, ValueError) as e:
        print("generation file cannot be found or generation file is faulty " + str(e))

    return generation

if __name__ == '__main__':
    pass
