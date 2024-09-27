import configparser
import json

from problem_specification.evaluation_methods import evaluation_methods
import numpy as np
from random import choices, randint, randrange, random
from collections import namedtuple
from typing import List, Callable, Tuple
from helper_files import genome_to_molecule as gtm
import copy

class bbEv_tournament(evaluation_methods.Evaluation):
	Genome = List[str]
	Population = List[Genome]
	FitnessFunc = Callable[[Genome, int, int], int]
	PopulateFunc = Callable[[], Population]
	SelectionFunc = Callable[[Population, FitnessFunc],Tuple[Genome, Genome]]
	CrossoverFunc = Callable[[Genome,Genome], Tuple[Genome, Genome]]
	MutationFunc = Callable[[Genome], Genome]

	def __init__(self, generation, individual, config_path, calculation_path):
		super().__init__(generation, individual)
		self.Building_Block = namedtuple('Building_Block', ['abbrev', 'num_atoms','origin', 'para_pos','para_angle', 'meta_pos','meta_angle', 'ortho_pos','ortho_angle','subs_index','fixed_left','complexity', 'path'])
		self.Coupling = namedtuple('Coupling', ['abbrev'])
		self.generation = generation
		self.individual = individual
		self.calculation_path = calculation_path
		self.config_path = config_path

		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		building_block_path = cfg.get('Building Procedure', 'building_block_path')
		self.building_blocks=gtm.load_building_blocks(building_block_path)

		self.force_symmetry = json.loads(str(cfg.get('Genetic Algorithm', 'force_symmetry')).lower())

		self.para = self.Coupling(abbrev="p")
		self.meta = self.Coupling(abbrev="m")

		self.couplings = [self.para, self.meta]

		self.substituent = np.asarray(str(cfg.get('Genetic Algorithm', 'substituent')).split(','), dtype=str)
		self.substituent_prob = np.asarray(str(cfg.get('Genetic Algorithm', 'substituent_prob')).split(','), dtype=float)





	def generate_genome(self, min_length:int, max_length: int) -> Genome:
		"""
		Generates genome mit maximum lengt of max_lengt (-> size random). Genome is generated from available couplings and building blocks

		Args:
			param1 (int) : maximum length of 
	                
	                

		Returns:
			Genome
	        """

	    #coupling must contain at least one block and two couplings
		assert min_length < max_length, "Genome length settings seem to be funny"
		assert min_length >= 1, "Min length must be an least 1.."


		num_building_blocks = randrange(min_length,max_length+1)
		indices_building_blocks = np.linspace(0,len(self.building_blocks),len(self.building_blocks),endpoint=False,dtype=int)   
		indices_couplings = np.linspace(0,len(self.couplings),len(self.couplings),endpoint=False,dtype=int)   
		selected_building_blocks = choices(indices_building_blocks, k=num_building_blocks)
		selected_couplings = choices(indices_couplings, k=num_building_blocks+2)

		genome=list()

		#add coupling to anchor
		genome.append(selected_couplings[0])
		for i in range(0, num_building_blocks):
			n_subs = len(self.building_blocks[selected_building_blocks[i]].subs_index)
			selected_subs = choices(self.substituent, weights=self.substituent_prob, k=n_subs)
			subs_string = ''.join(map(str, selected_subs))

			genome_string = str(selected_building_blocks[i]) + "#" + subs_string
			
			genome.append(genome_string)

			#Proper handling of problems if no difference between para and meta
			next_coupling = selected_couplings[i+1]
			genome.append(next_coupling)

		genome[0]=str(0)
		genome = [str(i) for i in genome]

		if(self.force_symmetry == True):
			genome = self.symmetrization(genome)

		genome = self.remove_coupling_problems(genome)

		return genome

		

	def generate_population(self, size: int,n_blocks_min:int, n_blocks_max: int) -> Population:
		print("tournament")
		return [self.generate_genome(n_blocks_min,n_blocks_max) for _ in range(size)]

	def fitness(self, genome: Genome, generation: int, individual: int) -> float:
		genome_copy = copy.deepcopy(genome)
		gtm.process_genome(generation,individual,genome_copy,self.calculation_path)
		return 0




	def selection_pair(self, population: Population, fitness_value) -> Population:	

		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		k = int(cfg.get('Genetic Algorithm', 'n_tournament_selection'))
		if(k<=2 or k >= len(population)):
			print("No suitable choice. Setting to 2")
			k = 2
		zipped_lists = zip(fitness_value, population) 
		sorted_pairs = sorted(zipped_lists)
		tuples = zip(*sorted_pairs)
		fitness_value, population = [ list(tuple) for tuple in  tuples]
		population.reverse()

		return choices(
			population=population[0:k],			
			k=2
		)

	def remove_coupling_problems(self, genome:Genome) -> Genome:
		"""
		Removes problems arising from blocks which have no difference between para and meta. The first coupling is set to 0 because there is also no difference between para and meta
		Args:
			genome:

		Returns:

		"""
		for j in range(len(genome)-1):
			# extract subs
			delimiter_index = genome[j].find("#")
			if(delimiter_index == -1):
				continue
			building_block = int(genome[j][0])
			block = int(building_block)
			if(self.building_blocks[block].para_pos==self.building_blocks[block].meta_pos):
				genome[j + 1] = str(0)
		#set first coupling to zero
		genome[0] = str(0)

		return genome


	def crossover(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		a,b = self.single_point_crossover(a, b)
		if (self.force_symmetry == True):
			a = self.symmetrization(a)
			b = self.symmetrization(b)
		a = self.remove_coupling_problems(a)
		b = self.remove_coupling_problems(b)
		return self.single_point_crossover(a,b)

	def single_point_crossover_old(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		
		length_a = len(a)
		length_b = len(b)
		if(length_a<=1 or length_b<=1):
			return a, b

		#ensure that coupling and blocks alter
		cut_a = randrange(length_a)	
		if(cut_a%2 == 0):
			cut_b = randrange(int(length_b/2))
			cut_b = 2*cut_b
		else:
			cut_b = randrange(int(length_b/2))
			cut_b = 2*cut_b+1

		return a[0:cut_a] + b[cut_b:length_b], b[0:cut_b] + a[cut_a:length_a]

	def single_point_crossover(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		
		length_a = len(a)
		length_b = len(b)
		minimim_length = np.min((length_a, length_b))
		if(length_a<=1 or length_b<=1):
			return a, b
			
		cut = randrange(minimim_length)			

		return a[0:cut] + b[cut:length_b], b[0:cut] + a[cut:length_a]

	def mutation(self, genome: Genome, num: int=2, probability: float = 0.5) -> Genome:		
		method = randrange(5)
		if(method == 0):
			mutated_genome = self.building_block_mutation(genome, probability)
		elif(method == 1):
			mutated_genome = self.coupling_mutation(genome, probability)
		elif(method == 2):
			mutated_genome = self.insert_mutation(genome, probability)
		elif(method == 3):
			mutated_genome = self.truncate_mutation(genome, probability)
		else:
			mutated_genome = self.substituent_mutation(genome, probability)


		if(self.force_symmetry == True):
			mutated_genome = self.symmetrization(mutated_genome)

		mutated_genome = self.remove_coupling_problems(mutated_genome)

		return mutated_genome

	def building_block_mutation(self, genome: Genome, probability: float = 0.5) -> Genome:	
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'block_mutation_prob'))

		mutated_genome = list()
		if(random()<probability and len(genome)>=3):
			new_block = randrange(len(self.building_blocks))
			#add subs
			n_subs = len(self.building_blocks[new_block].subs_index)
			selected_subs = choices(self.substituent, weights=self.substituent_prob, k=n_subs)
			subs_string = ''.join(map(str, selected_subs))
			genome_string = str(new_block) + "#" + subs_string

			#find position in genome
			block_to_mutate = randrange(int((len(genome)-1)/2))
			block_to_mutate = block_to_mutate+block_to_mutate+2
			
			
			mutated_genome.extend(genome[0:block_to_mutate-1])
			mutated_genome.append(genome_string)
			mutated_genome.extend(genome[block_to_mutate:len(genome)])
			print("block mutation!" + str(mutated_genome))
			return mutated_genome
		return genome

	def coupling_mutation(self, genome: Genome, probability: float = 0.5) -> Genome:	
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'coupling_mutation_prob'))

		mutated_genome = list()
		if(random()<probability and len(genome)>=3):
			new_coupling = randrange(len(self.couplings))
			coupling_to_mutate = randrange(0,int((len(genome)+1)/2))
			coupling_to_mutate = coupling_to_mutate+coupling_to_mutate+1
			
			#genome = genome[0:coupling_to_mutate-1] + couplings[new_coupling].abbrev + genome[coupling_to_mutate:len(genome)]
			
			mutated_genome.extend(genome[0:coupling_to_mutate-1])
			mutated_genome.append(str(new_coupling))
			mutated_genome.extend(genome[coupling_to_mutate:len(genome)])
			print("coupling mutation! " + str(mutated_genome))
			return mutated_genome
		return genome

	def insert_mutation(self, genome: Genome, probability: float = 0.5):
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'insert_mutation_prob'))
		n_blocks_max = float(cfg.get('Genetic Algorithm', 'n_blocks_max'))

		mutated_genome = list()
		if(random()<probability and len(genome)/2 < n_blocks_max):

			new_coupling = randrange(len(self.couplings))
			new_block = randrange(len(self.building_blocks))
			# add subs
			n_subs = len(self.building_blocks[new_block].subs_index)
			selected_subs = choices(self.substituent, weights=self.substituent_prob, k=n_subs)
			subs_string = ''.join(map(str, selected_subs))
			genome_string = str(new_block) + "#" + subs_string

			insert_coupling=randrange(0,int((len(genome)+1)/2))
			insert_coupling=insert_coupling+insert_coupling+1
			to_add_at_end = genome[insert_coupling:len(genome)]
			mutated_genome.extend(genome[0:insert_coupling])

			mutated_genome.append(genome_string)
			mutated_genome.append(str(new_coupling))

			mutated_genome.extend(to_add_at_end)
			print("insert mutation! " + str(mutated_genome))
			return mutated_genome

		return genome

	def truncate_mutation(self, genome: Genome, probability: float = 0.5):
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'truncate_mutation_prob'))
		n_blocks_min = float(cfg.get('Genetic Algorithm', 'n_blocks_min'))

		mutated_genome = list()
		n_blocks = int((len(genome)-1)/2)
		if(random()<probability and n_blocks > n_blocks_min):

			block_to_truncate = randrange(int((len(genome)-1)/2))	
			block_to_truncate = block_to_truncate+block_to_truncate+2

			mutated_genome.extend(genome[0:block_to_truncate-1])
			mutated_genome.extend(genome[block_to_truncate+1:len(genome)])
			print("truncate mutation!",str(genome), str(mutated_genome))
			return mutated_genome

		return genome

	def substituent_mutation(self, genome: Genome, probability: float = 0.5):
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'substituent_mutation_prob'))


		mutated_genome = list()
		if(random()<probability):
			index_to_mutate = randrange(int((len(genome) - 1) / 2))
			index_to_mutate = index_to_mutate + index_to_mutate + 1
			part_to_alter = genome[index_to_mutate]

			#extract subs
			delimiter_index = part_to_alter.find("#")
			if(delimiter_index==-1):
				print("something is wrong", part_to_alter)
			building_block = int(part_to_alter[0])
			block = int(building_block)

			#check if there are any subs
			if(len(self.building_blocks[building_block].subs_index)==0):
				return genome
			subsituents = part_to_alter[delimiter_index + 1:len(part_to_alter)]
			sub_to_alter = randrange(len(subsituents))
			selected_sub = choices(self.substituent, weights=self.substituent_prob, k=1)[0]
			subsituents = subsituents[:sub_to_alter] + selected_sub + subsituents[sub_to_alter + 1:]
			genome_string = str(block) + "#" + subsituents


			mutated_genome.extend(genome[0:index_to_mutate])
			mutated_genome.append(genome_string)
			mutated_genome.extend(genome[index_to_mutate+1:len(genome)])
			print("substituent mutation!",str(genome), str(mutated_genome))
			return mutated_genome

		return genome

	def symmetrization(self, genome) -> Genome:
		"""
		Forces symmetry regarding the building blocks. First genome part (left to center) is assumed to be dominant -> building blocks are symmetrized to right part. Couplings are not changed. To reduce the bias the genome is inverted at the beginning randomly (at back at the end) to change the dominant part to the right
		Args:
			param1 () : maximum length of


		Returns:
			Genome
	        """
		genome = copy.deepcopy(genome)
		reversed = False
		#Let the right part be dominant in 50 percent of cases
		if (random() < 0.5):
			genome.reverse()
			reversed = True
			print("Reverse!", genome)

		genome_symmetrized = list()
		n_couplings = len(genome) - int(len(genome) / 2)
		genome_subs = [genome[i][2:len(genome[i])] for i in range(0, len(genome))]
		genome_blocks = [genome[i][0:1] for i in range(0,len(genome)) if i%2==1]



		# symmetry on block in the middle
		if (n_couplings % 2 == 0):
			blocks_to_keep = genome_blocks[:int(len(genome_blocks) / 2) + 1]
			blocks_to_keep_ = genome_blocks[:len(blocks_to_keep) - 1][::-1]
			blocks_symmetrized = blocks_to_keep + blocks_to_keep_

			for i in range(0, len(genome)):
				#add couplings
				if (i % 2 == 0):
					genome_symmetrized.append(genome[i])
				#add blocks
				else:
					block_to_add = blocks_symmetrized.pop(0)
					n_subs_observed = len(genome_subs[i])
					n_subs_expected = len(self.building_blocks[int(block_to_add)].subs_index)
					if(n_subs_expected == n_subs_observed):
						pass
					#fill up
					elif(n_subs_observed<n_subs_expected):
						n_subs_missing = n_subs_expected-n_subs_observed
						additional_subs = choices(self.substituent, weights=self.substituent_prob, k=n_subs_missing)
						additional_subs = ''.join(map(str, additional_subs))
						genome_subs[i] = genome_subs[i] + additional_subs
					#delete too much subs
					elif(n_subs_observed>n_subs_expected):
						genome_subs[i] = genome_subs[i][:n_subs_expected]
					genome_symmetrized.append(block_to_add + "#" + genome_subs[i])


		# symmetry on coupling in the middle
		else:
			blocks_to_keep = genome_blocks[:int(len(genome_blocks) / 2)]
			blocks_to_keep_ = blocks_to_keep[::-1]
			blocks_symmetrized = blocks_to_keep + blocks_to_keep_

			for i in range(0, len(genome)):
				# add couplings
				if (i % 2 == 0):
					genome_symmetrized.append(genome[i])
				# add blocks
				else:
					block_to_add = blocks_symmetrized.pop(0)
					n_subs_observed = len(genome_subs[i])
					n_subs_expected = len(self.building_blocks[int(block_to_add)].subs_index)
					if (n_subs_expected == n_subs_observed):
						pass
					# fill up
					elif (n_subs_observed < n_subs_expected):
						n_subs_missing = n_subs_expected - n_subs_observed
						additional_subs = choices(self.substituent, weights=self.substituent_prob, k=n_subs_missing)
						additional_subs = ''.join(map(str, additional_subs))
						genome_subs[i] = genome_subs[i] + additional_subs
					# delete too much subs
					elif (n_subs_observed > n_subs_expected):
						genome_subs[i] = genome_subs[i][:n_subs_expected]
					genome_symmetrized.append(block_to_add + "#" + genome_subs[i])

		if(reversed==True):
			genome_symmetrized.reverse()
			print("back to normal Reverse!", genome_symmetrized)

		return genome_symmetrized





