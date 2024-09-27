import sys
import os
import os.path
#sys.path.append(os.path.realpath('..'))
from typing import List, Callable, Tuple

class Evaluation:
	Genome = List[str]
	Population = List[Genome]
	FitnessFunc = Callable[[Genome, int, int], int]
	PopulateFunc = Callable[[], Population]
	SelectionFunc = Callable[[Population, FitnessFunc],Tuple[Genome, Genome]]
	CrossoverFunc = Callable[[Genome,Genome], Tuple[Genome, Genome]]
	MutationFunc = Callable[[Genome], Genome]

	def __init__(self, generation, individual):
		pass
		self.generation = generation
		self.individual = individual


	def generate_genome(self, max_length: int) -> Genome:
		"""
		Generates genome mit maximum lengt of max_lengt (-> size random). Genome is generated from available couplings and building blocks

		Args:
			param1 (int) : maximum length of 
	                
	                

		Returns:
			Genome
	        """
		
		return genome

		

	def generate_population(self, size: int,n_blocks_min: int, n_blocks_max: int) -> Population:
		"""
		Generates population of size size withm genomes of length n_blocks_max

		Args:
			param1 (int) : population size
			param2 (int) : n_blocks_min
			param3 (int) : n_blocks_max
	                
	                

		Returns:
			Genome
		"""
		return 0

	def fitness(self, genome: Genome, generation: int, individual: int) -> float:
		"""
		Evaluates finess

		Args:
			param1 (Genome) : genome
			param2 (int) : generation
			param3 (int) : individual
	                
	                

		Returns:
			(float) fittness
		"""
		
		return 0.0




	def selection_pair(self, population: Population, fitness_value) -> Population:
		"""
		Select pair of population

		Args:
			param1 (Population) : Population
			param2 (Callable) : fitness_function
			
		Returns:
			(Population) population
		"""		
		return 0



	def mutation(self, genome: Genome, num: int=1, probability: float = 0.5) -> Genome:
		"""
		Mutation of genome (num of mutations with given probability )

		Args:
			param1 (Genome) : genome
			param2 (int) : number of mutations
			param3 (float) : probability of mutation           

		Returns:
			(Genome) population
		"""		
		return 0

	def crossover(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		"""
		Crossover of genome a and b

		Args:
			param1 (Genome) : genome a
			param2 (Genome) : genome b       

		Returns:
			(Genome,Genome) 
		"""		
		return ([0],[0])


	