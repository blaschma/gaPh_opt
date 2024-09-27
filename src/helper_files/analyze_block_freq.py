import configparser
import sys
from collections import Counter
import numpy as np
import matplotlib
from matplotlib import cm

matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import copy
from more_itertools import sort_together
import genome_to_molecule as gtm
import scienceplots

plt.style.use(['science','no-latex'])
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=16)
plt.rcParams["figure.figsize"] = (3, 4)

def read_population(generation, config_path, calculation_path):
    """
    read generation individuals and their fitness from calculation path and extracts building blocks

    Args:
        param1 (int): number of generation which should be read
        param2 (String): Path to config file
        param3 (String): Path to calculation
    Returns:
        (population, fitness_values)

    """
    try:
        filename_population = calculation_path + "/generation_data/" + str(generation) + "/summary.dat"
        file_population = open(filename_population)
        filename_fitness = calculation_path + "/generation_data/" + str(generation) + "/fitness.dat"
        file_fitness = open(filename_fitness)
    except OSError as e:
        print("Cannot open file " + str(e))
        return -1, -1

    # read population
    population = list()
    for line in file_population:
        tmp = line.strip().split(", ")
        tmp = [str(tmp[i]).replace("'", "").replace(" ", "") for i in range(0, len(tmp))]
        population.extend(tmp)

    #filter blocks
    blocks = list()
    for i in range(0, len(population)):
        # skip coupling indices
        if (str(population[i]).find("#") == -1):
            continue
        delimiter_index = population[i].find("#")
        block = population[i][0:delimiter_index]
        blocks.append(int(block.replace("'", "").replace("\"", "")))

    # read fitness values
    fitness_value = list()
    for line in file_fitness:
        try:
            tmp = float(line)
        except ValueError as e:
            print("Faulty fitness file " + str(e))
            return -1, -1
        fitness_value.append(tmp)

    return blocks, fitness_value


if __name__ == '__main__':
    # param1 (String): Path to config file

    config_path = sys.argv[1]

    cfg = configparser.ConfigParser()
    cfg.read(config_path)
    building_block_path = cfg.get('Building Procedure', 'building_block_path')
    calculation_path = cfg.get('Basics', 'calculation_path')

    building_blocks = gtm.load_building_blocks(building_block_path)
    n_blocks = len(building_blocks)

    n_generations = int(cfg.get('Genetic Algorithm', 'generation_limit'))+1
    hist = np.zeros((n_generations, n_blocks), dtype=int)
    n_handles = 10

    for i in range(0, n_generations):
        population = read_population(i, config_path, calculation_path)[0]
        #generation not yet processed
        if(population == -1):
            hist = hist[0:i,:]
            n_generations = i
            break
        dictio = dict(Counter(population[0:int(2 * len(population) / 3)]))
        for j in range(0, n_blocks):
            try:
                hist[i, j] = int(dictio[j])
            except KeyError:
                pass

    fig, ax1 = plt.subplots()

    freq_last_gen = list()
    block_id = list()
    frequencies = list()
    for j in range(0, n_blocks):

        block_count = np.asarray(copy.deepcopy(hist[:, j]))
        freq = list()

        for i in range(0, n_generations):
            freq.append(float(block_count[i]) / np.sum(hist[i, :]))
        freq_last_gen.append(freq[len(freq) - 1])
        block_id.append(j)
        frequencies.append(freq)


        ax1.plot(freq)
    blocks_label_dict = {0: 3, 1: 4, 2: 5, 3: 1, 4: 8, 5: 6, 6: 9, 7: 10, 8: 2, 9: 7}
    # find most important blocks
    colormap = cm.get_cmap('nipy_spectral', n_blocks)
    freq_last_gen, block_id,frequencies = sort_together([freq_last_gen, block_id, frequencies], reverse=True)
    for i in range(0, n_blocks):
        if(i<= n_handles):
            ax1.plot(frequencies[i], label=blocks_label_dict[block_id[i]], color=colormap(i))
        else:
            ax1.plot(frequencies[i], color=colormap(i))


    ax1.tick_params(axis='y')
    ax1.tick_params(axis='x')
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Relative frequency')
    ax1.grid()
    ax1.legend(frameon=True, fancybox=True, ncol=2)
    #sort legend by label
    handles, labels = ax1.get_legend_handles_labels()
    labels = [int(label) for label in labels]
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax1.legend(handles, labels, frameon=True, fancybox=True, ncol=int(n_handles/2))
    plt.tick_params(axis="x")
    plt.tick_params(axis="y")
    plt.savefig(calculation_path + "/block_frequency.pdf", bbox_inches='tight', transparent=True)
    plt.savefig(calculation_path + "/block_frequency.svg", bbox_inches='tight', transparent=True)
