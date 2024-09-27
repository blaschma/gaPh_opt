import os
import fnmatch
import os.path
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import scienceplots
plt.style.use(['science','no-latex'])
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=16)
plt.rcParams["figure.figsize"] = (6, 4)

def read_fitness_info(generation, config_path, calculation_path):
	"""
    read current generation individuals and their fitness from calculation path

    Args:
    	param1 (int): number of generation which should be read
		param2 (String): Path to config file
		param3 (String): Path to calculation              
    Returns:
    	(population, fitness_values)
            
    """
	try:
		filename_population = calculation_path + "/curr_population.dat"	
		file_population = open(filename_population)
		filename_fitness = calculation_path + "/generation_data/" + str(generation) + "/fitness.dat"	
		file_fitness = open(filename_fitness)
	except OSError as e:
		print("Cannot open file " + str(e))
		return -1,-1, -1
	#read n_unique
	unique_file = -1
	for file in os.listdir(calculation_path + "/generation_data/" + str(generation)):
			if fnmatch.fnmatch(file, '*n_unique.dat'):
				unique_file = file
				#print("unique_file " + str(unique_file))
	
	if(unique_file != -1):
		unique_file = calculation_path + "/generation_data/" + str(generation) + "/" +unique_file
		unique_file = open(unique_file, "r")
		for line in unique_file:
			unique_line = line
		n_unique = int(line)

	#read fitness values
	fitness_value = list()
	for line in file_fitness:
		#print(line)
		try:
			tmp = float(line)
		except ValueError as e:
			print("Faulty fitness file " + str(e))
			return -1,-1
		fitness_value.append(tmp)
		
	if(unique_file != -1):
		return fitness_value, n_unique
	else:
		return fitness_value, -1


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

	except (OSError,ValueError) as e:
		print("generation file cannot be found or generation file is faulty " + str(e))

	return generation





if __name__ == '__main__':
	# sys.argv[1] path to process
	# sys.argv[2] config path
	# sys.argv[3] figure path
	calculation_path = sys.argv[1]
	config_path = sys.argv[2]
	#how many generations have been calculated
	generation = read_generation(config_path, calculation_path)
	if(generation == None):
		exit()
	
	generations_to_check = np.arange(0, generation)

	#load fitness values
	fitness_values = list()
	fitness_means=list()
	std_deviation=list()
	n_unique = list()
	for i in generations_to_check:
		fitness_value = read_fitness_info(i, config_path, calculation_path)
		fitness_values.append(fitness_value[0])
		n_unique.append(fitness_value[1])
		if(fitness_value[1]==-1):
			fitness_means.append(np.mean(fitness_value[1]))
			std_deviation.append(np.std(np.asarray(fitness_value[1])))			
		else:
			print("else")
			fitness_means.append(np.mean(fitness_value[0][0:fitness_value[1]]))
			std_deviation.append(np.std(np.asarray(fitness_value[0][0:fitness_value[1]])))
	print(fitness_means)
	print(std_deviation)
	fig, ax = plt.subplots(1)
	num_individuals = len(fitness_value[0])
	#color & plotting stuff
	phi = np.linspace(0, 2*np.pi, num_individuals)
	rgb_cycle = np.vstack((            # Three sinusoids
    .5*(1.+np.cos(phi          )), # scaled to [0,1]
    .5*(1.+np.cos(phi+2*np.pi/3)), # 120° phase shifted.
    .5*(1.+np.cos(phi-2*np.pi/3)))).T # Shape = (60,3)
	#print(rgb_cycle)
	markersize = 5
	for xe, ye in zip(generations_to_check, fitness_values):
		if(n_unique[xe] == -1):
			for i in range(0, len(ye)):
				if ([xe] != 0):
					ax.scatter([xe], ye[i], c=rgb_cycle[i], s=markersize, marker="x")
				else:
					ax.scatter([xe], ye[i], c=rgb_cycle[i], s=markersize, marker="o")
		else:
			for i in range(0, n_unique[xe]):
				if([xe] != 0):
					ax.scatter([xe], ye[i], c=rgb_cycle[i], s=markersize, marker="x")
				else:
					ax.scatter([xe], ye[i], c=rgb_cycle[i], s=markersize, marker="o")
			for i in range(n_unique[xe], len(ye)):
				ax.scatter([xe], ye[i], c=rgb_cycle[i], s=markersize, marker="o")
		#ax.scatter([xe] * len(ye), ye, c=rgb_cycle, s=num_individuals, marker="x")
		
	#print(generations_to_check)
	#print(fitness_means)
	#print(std_deviation)
	ax.plot(generations_to_check, fitness_means, color="black", lw = 3)
	#ax.plot(generations_to_check, fitness_means-np.asarray(std_deviation), color="blue",linestyle='dashed')
	#ax.plot(generations_to_check, fitness_means+np.asarray(std_deviation), color="blue",linestyle='dashed')
	ax.set_xlabel('Generation')
	ax.set_ylabel('Fitness values')
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	#ax.set_yscale('symlog')
	ax.set_ylim((0.0,None))
	ax.set_xlim((0.0, 100))
	plt.savefig( sys.argv[3] + "/fitness_values.pdf", bbox_inches='tight', transparent=True)
	plt.savefig(sys.argv[3] + "/fitness_values.svg", bbox_inches='tight', transparent=True)
	plt.savefig(sys.argv[3] + "/fitness_values.png", bbox_inches='tight', transparent=True, dpi=300)