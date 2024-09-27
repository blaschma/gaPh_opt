import codecs
import json
import numpy as np
import subprocess
import tmoutproc as top
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import configparser
import os
import sys
from rdkit import Chem
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

def determine_syntetic_acessability(gen_dir):
	syntetic_accessibilities = []
	generation = os.listdir(gen_dir)
	individuals = [int(individual) for individual in generation if  os.path.isdir(f"{gen_dir}/{individual}")]
	individuals = sorted(individuals)
	# analyze individuals
	for i, individual in enumerate(individuals):
		try:
			coord = top.read_xyz_file(f"{gen_dir}/{individual}/xtbopt.xyz")
		except FileNotFoundError as e:
			print(f"File {path}/{individual}/xtbopt.xyz not found")

			syntetic_accessibilities.append(10)
			continue
		coord = top.x2t(coord)
		# analyze blocks

		# analyze syntetic accessibility
		# write coord without anchor to file
		coord_without_gold = coord.T[coord.T[:, 3] != 'au'].T
		coord_without_gold = top.t2x(coord_without_gold)
		top.write_xyz_file(f"{path}/{individual}/xtbopt_without_au.xyz", coord_without_gold)
		#todo: Do not hard code the path to open babel
		command = f"PATH_TO_OPEN_BABEL/obabel {path}/{individual}/xtbopt_without_au.xyz -O {path}/{individual}/xtbopt_without_au.sdf"
		result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
		sdf_file = f"{path}/{individual}/xtbopt_without_au.sdf"
		supplier = Chem.SDMolSupplier(sdf_file)
		mol = supplier[0]
		synthetic_accessibility = sascorer.calculateScore(mol)
		syntetic_accessibilities.append(synthetic_accessibility)
	return syntetic_accessibilities

def load_transmission_data(gen_dir):
	"""
	Loads transmission data of all individuals in given dir. 

	Args:
		param1 (String): dir to process

	Returns:
		list()

	"""	

	dirs = os.listdir(gen_dir)
	dirs = [int(i) for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	dirs = [str(i) for i in dirs]
	
	kappas = list()

	#read kappa files and append maximum kappa value to list
	for i in range(0,len(dirs)):
		filename = gen_dir + dirs[i] + "/kappa.dat"
		dat_Content = top.read_plot_data(filename)

		kappa_max = np.nanmax(dat_Content[1,:])

		#if KILL_SIGNAL_SET
		if(os.path.isfile(gen_dir + dirs[i] + "/KILL_SIGNAL_SET")):
			print("KILL_SIGNAL_SET", gen_dir + dirs[i])

		kappas.append(kappa_max)

	return kappas



def eval_fitness_function(input_values, config_path, filter=-1):
	"""
	Evaluates fitness function read from coord file
	Args:
		input_values: (list) input_values for fitness function
		config_path: (String) path to config file

	Returns:
		fitness: (list) list of corresponding fitness values
	"""

	#Read config file
	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	fitness_function = str(cfg.get('Genetic Algorithm', 'fitness_function')).replace("\"", "")
	subs_penalty = json.loads(str(cfg.get('Genetic Algorithm', 'subs_penalty')).lower())

	if(subs_penalty==False):
		kappas = input_values[0]
		sa = input_values[1]
		fitness_values = list()
		for i in range(0, len(kappas)):
			x = kappas[i]
			z = sa[i]
			#handle not converged calculation
			if(x == -1):
				fitness = 0
			else:
				fitness = eval(fitness_function)
			if(filter != -1):
				if(filter[i]==True):
					fitness = 0


			fitness_values.append(fitness)

		return fitness_values

	if(subs_penalty==True):
		kappas, subs_cost, sa = input_values
		fitness_values = list()
		for i in range(0, len(kappas)):
			x = kappas[i]
			y = subs_cost[i]
			z = sa[i]
			# handle not converged calculation
			if (x == -1):
				fitness = 0
			else:
				fitness = eval(fitness_function)
			if (filter != -1):
				if (filter[i] == True):
					fitness = 0

			fitness_values.append(fitness)

		return fitness_values

def eval_filter(path):
	gen_dir = path
	population = list()
	file_population = open(gen_dir + "/summary.dat")
	for line in file_population:
		tmp = line.strip().split(", ")
		tmp = [str(tmp[i]).replace("'", "").replace(" ", "") for i in range(0, len(tmp))]
		population.append(tmp)
	population_filtered = [population[i][1::2] for i in range(0,len(population))]
	filter = [all(c == '2#' for c in population_filtered[i]) for i in range(0,len(population_filtered))]*1

	return filter

def eval_subs_cost(gen_dir, config_file):

	def atom_to_subs(atom):
		if (atom == "Br"):
			return "B"
		elif (atom == "Cl"):
			return "C"
		elif (atom == "H"):
			return "H"
		elif (atom == "F"):
			return "F"
		elif (atom == "I"):
			return "I"
		else:
			return "H"

	def eval_cost(coord, substituent, subs_cost):
		cost = 0
		for i in range(0, coord.shape[1]):
			# if line is au atom
			a = atom_to_subs(coord[0, i])
			if(a in substituent):
				cost += subs_cost[substituent.index(a)]
		return cost



	#load subs cost and subs
	cfg = configparser.ConfigParser()
	cfg.read(config_file)
	subs_cost = list(np.array(str(cfg.get('Genetic Algorithm', 'substituent_cost')).split(','), dtype=float))
	substituent = list(np.array(str(cfg.get('Genetic Algorithm', 'substituent')).split(','), dtype=str))
	assert len(subs_cost) == len(substituent)

	#load cost from config file
	dirs = os.listdir(gen_dir)
	dirs = [int(i) for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	dirs = [str(i) for i in dirs]

	cost = list()
	# read coord files and eval cost
	for i in range(0, len(dirs)):
		filename = gen_dir + dirs[i] + "/coord.xyz"
		coord = top.read_xyz_file(filename)
		cost.append(eval_cost(coord, substituent, subs_cost))

	return cost


def write_fittness(fittness, path):
	"""
	Write fitness values to file

	Args:
		param1 (List): fittness

	Returns:
		

	"""
	with open(path + "/fitness.dat", "w") as file:
		for i in range(len(fittness)):
			file.write(str(fittness[i])+"\n")




if __name__ == '__main__':
	# sys.argv[1] path to process
	# sys.argv[2] config path

	#"""
	#load all data
	path=sys.argv[1]
	config_path = sys.argv[2]

	cfg = configparser.ConfigParser()
	cfg.read_file(codecs.open(config_path, "r", "utf8"))
	subs_penalty = json.loads(str(cfg.get('Genetic Algorithm', 'subs_penalty')).lower())


	if(subs_penalty == True):
		kappas = load_transmission_data(path)
		subs_cost = eval_subs_cost(path, config_path)
		syntetic_accesability = determine_syntetic_acessability(path)
		fitness = eval_fitness_function((kappas,subs_cost, syntetic_accesability), config_path)

	else:
		kappas = load_transmission_data(path)
		syntetic_accesability = determine_syntetic_acessability(path)
		fitness = eval_fitness_function((kappas, syntetic_accesability), config_path)

	write_fittness(fitness,path)

