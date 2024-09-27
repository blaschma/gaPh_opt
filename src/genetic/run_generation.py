import configparser
import os
import os.path
from os import path
import sys
import genetic_algorithm as ga
from functools import partial


def run_generation(generation : int, config_path, calculation_path):
	"""
    Runs evolution. Inits calculation dirs and invokes evaluation of first generation

    Args:
    		param1 (int) : number of generation
            param2 (String) : Path to config file
            param3 (String) : Path to calculation

    Returns:
            error codes (int)
    """

    #load specification from config file 
	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	population_size = int(cfg.get('Genetic Algorithm', 'population_size'))
	generation_limit = int(cfg.get('Genetic Algorithm', 'generation_limit'))	
	n_blocks_max = int(cfg.get('Genetic Algorithm', 'n_blocks_max'))
	n_blocks_min = int(cfg.get('Genetic Algorithm', 'n_blocks_min'))
	evaluation_methods_path = cfg.get('Genetic Algorithm', 'evaluation_methods_path')
	genetic_algorithm_path = cfg.get('Basics', 'genetic_algorithm_path')
	n_elitism = int(cfg.get('Genetic Algorithm', 'n_elitism'))

	#check if default setting for evaluation methods
	if(evaluation_methods_path == "default"):
		#make sure problem_specification path is found
		sys.path.append(os.path.realpath(genetic_algorithm_path + "/src/"))
		sys.path.append(os.path.realpath('..'))
		from problem_specification import bbEV 
		ev = bbEV.bbEv(generation, 0, config_path, calculation_path) #-> todo individual processing
		pass
		
	#load specific evaluation methods
	else:
		#make sure problem_specification path is found
		sys.path.append(os.path.realpath(genetic_algorithm_path + "/src/"))
		sys.path.append(os.path.realpath('..'))
		from problem_specification import bbEV_tournament
		ev = bbEV_tournament.bbEv_tournament(generation, 0, config_path, calculation_path) #-> todo individual processing
		pass

	#first generation
	if(generation == 0):
		population, family_register = ga.run_generation(
		populate_func=partial(
			ev.generate_population, size=population_size, n_blocks_min=n_blocks_min, n_blocks_max=n_blocks_max
			),
		fitness_func=ev.fitness,
		selection_func=ev.selection_pair,			
		crossover_func=ev.crossover,			
		mutation_func=ev.mutation,			
		population=0,
		population_size=population_size,
		generation=generation,
		fitness_value=0,
		n_elitism=n_elitism)
		print("Generation zero " + str(generation))

	#every other generation
	else:

		#generation-1 because prevois generation should be read
		population, fitness_value = read_population(generation-1,config_path, calculation_path)
		print(fitness_value)

		population, family_register = ga.run_generation(
		populate_func=partial(
			ev.generate_population, size=population_size, n_blocks_min=n_blocks_min, n_blocks_max=n_blocks_max
			),
		fitness_func=ev.fitness,
		selection_func=ev.selection_pair,			
		crossover_func=ev.crossover,			
		mutation_func=ev.mutation,			
		population=population,
		population_size=population_size,
		generation=generation,
		fitness_value=fitness_value,
		n_elitism=n_elitism)

	#todo move problem specific code to problem_specification
	write_generation(population,generation, config_path, calculation_path)
	write_genomes_to_archive(population, generation, calculation_path, config_path)
	write_family_register(family_register, generation, calculation_path)


def write_genomes_to_archive(population, generation, calculation_path, config_path):
	"""
	Writes genomes to archive. First of all it is checked if genome is already in archige

	Args:
		param1 (Polulation): population
		param2 (int): generation
		param3 (String): path to population
		param4 (String): path to config file

	Returns:
		

	"""
	#find existing archive
	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	archive_path = cfg.get('Basics', 'archive_archive_path')

	if(os.path.exists(archive_path)==False):
		print("create file")
		f = open(archive_path, "w")
		f.close()

	#read existing archinve
	with open(archive_path, "r") as archive_file:
		archive_population = list()
		archive_paths = list()
		for line in archive_file:
			if(len(line)<3):
				continue
			line = line.strip().split("	")
			tmp = line[0].replace("[", "").replace("]", "")
			tmp = tmp.split(",")
			tmp = [str(tmp[i]).replace("'", "").replace(" ", "") for i in range(0,len(tmp))]
			archive_population.append(tmp)
			archive_paths.append(line[1])


	#check which file to add and add the file
	with open(archive_path, "a") as archive_file:
		for i in range(len(population)):
			#print(population[i])
			#print(archive_population)
			if (population[i] in archive_population) == False:
				#print("adding to archive " + str(population[i]))
				path_of_individual = calculation_path + "/generation_data/"+ str(generation) + "/" + str(i)
				archive_file.write(str(population[i]) + "	" + path_of_individual + "\n")




def write_generation(population, generation, config_path, calculation_path):
	"""
    write current generation to calculation_path/current_generation.dat and to generation_data/generation/generation_summary.dat.
    First file contains only information about the population. The second file contains the fittness value, too. In addition a file 
    calculation_path/generation.dat which contains the current number of generations

    Args:
		param1 (Population): population to write
		param2 (int): Generation
		param3 (String): Path to config file
		param4 (String): Path to calculation
    Returns:
            
    """

	try:
		first_file_path = calculation_path + "/curr_population.dat"
		second_file_path = calculation_path + "/generation_data/" + str(generation)
		first_file = open(first_file_path, "w")
		generation_file = open(calculation_path + "/generation.dat", "w")
		if(path.exists(second_file_path) == False):
			os.mkdir(second_file_path)	
		second_file = open(second_file_path + "/summary.dat", "w")	
	except OSError as e:
		print("log file cannot be opened " + str(e))
		return 1	
		
	for i in range(0, len(population)):		
		first_file.write(str(population[i]).replace('[', "").replace("]", "")+ "\n")
		second_file.write(str(population[i]).replace('[', "").replace("]", "")+ "\n")

	generation_file.write(str(generation))
	generation_file.close()

	first_file.close()
	second_file.close()

def read_population(generation, config_path, calculation_path):
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
		return -1,-1

	#read population
	population = list()
	for line in file_population:
		#print(line)
		tmp = line.strip().split(", ")
		tmp = [str(tmp[i]).replace("'", "").replace(" ", "") for i in range(0,len(tmp))]
		population.append(tmp)

	#read fitness values
	fitness_value = list()
	for line in file_fitness:
		try:
			tmp = float(line)
		except ValueError as e:
			print("Faulty fitness file " + str(e))
			return -1,-1
		fitness_value.append(tmp)
		

	return population, fitness_value

def write_family_register(family_register, generation, calculation_path):
	"""
    Invokes the next generation. The number of generation is read from calculation_path/generation.dat

    Args:
    	param1 (String): family_register (created during evolution)
		param2 (int): generation    
		param3 (String): Path to calculation         
    Returns:
            
    """
    #find number of unique individuals and write to file if possible
	n_unique = -1
	try:
		index = family_register.find("individuals")
		n_unique = int(family_register[index+len("individuals"): len(family_register)-1])
	except ValueError as e:
		print("Can not find or cast number of unique_individuals " + str(e))
	if(n_unique != -1):
		unique_individuals_path = calculation_path + "/generation_data/" + str(generation) + "/" + str(generation) + "_n_unique.dat"
		try:
			unique_individuals_file = open(unique_individuals_path, "w")
			unique_individuals_file.write(str(n_unique))
			unique_individuals_file.close()
		except OSError as e:
			print("Cannot open file " + str(e))

	family_register_path = calculation_path + "/generation_data/" + str(generation) + "/" + str(generation) + "_family_register.dat"
	
	try:
		family_register_file = open(family_register_path, "w")
		family_register_file.write(family_register)
		family_register_file.close()
	except OSError as e:
		print("Cannot open file " + str(e))



def next_generation(config_path, calculation_path):
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

	#check if generation limit is reached
	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	generation_limit = int(cfg.get('Genetic Algorithm', 'generation_limit'))
	if(generation>=generation_limit):
		print("Generation limit reached")
		return 0


	#check if generation was correctly processed
	file_to_check = calculation_path + "/generation_data/" + str(generation) + "/fitness.dat"
	if(os.path.exists(file_to_check) == False):
		print("Generation " + str(generation) + " was not processed correctly")
		return -1

	#increase number of genrations
	generation += 1
	run_generation(generation, config_path, calculation_path)


if __name__ == '__main__':
	pass


