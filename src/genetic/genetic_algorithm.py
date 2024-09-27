from typing import List, Callable, Tuple
import copy


Genome = List[int]
Population = List[Genome]
FitnessFunc = Callable[[Genome, int, int], float]
PopulateFunc = Callable[[], Population]
SelectionFunc = Callable[[Population, FitnessFunc],Tuple[Genome, Genome]]
CrossoverFunc = Callable[[Genome,Genome], Tuple[Genome, Genome]]
MutationFunc = Callable[[Genome], Genome]



def run_generation(
		populate_func: PopulateFunc,
		fitness_func: FitnessFunc,		
		selection_func : SelectionFunc,
		crossover_func : CrossoverFunc,
		mutation_func: MutationFunc,
		population : Population,
		population_size:int,
		generation: int,
		fitness_value,
		n_elitism: int
	) :


	

	#initialize first population
	if(generation == 0):
		population = populate_func()
		#invoke fitness evaluation
		population_for_fitness_eval = copy.deepcopy(population)
		for i in range(0, population_size):
			fitness_func(population_for_fitness_eval[i],generation,i)
		return population, ""


	#sort population
	zipped_lists = zip(fitness_value, population)
	sorted_pairs = sorted(zipped_lists, reverse=True)
	tuples = zip(*sorted_pairs)
	fitness_value, population = [ list(tuple) for tuple in  tuples]


	#family register
	family_register=""

	#take best individuals of generation and..
	assert n_elitism>0
	assert n_elitism<population_size
	next_generation = population[0:n_elitism]

	for i in range(0,n_elitism):
		family_register += str(population[i]) + "\n"

	#mutate best n_elistim
	for j in range(0, n_elitism):
		genome_ = population[j]
		print(genome_)
		mutated_elite = mutation_func(population[j])
		print(str(genome_) + " mutated to " + str(mutated_elite))
		family_register += str(genome_) + " mutated to " + str(mutated_elite) + "\n"
		next_generation += [mutated_elite]

	#...fill generation with mutated and cross over children
	for j in range(int(len(population)/2)-int(n_elitism)):
		#select parent according to selection function
		parents = selection_func(population, fitness_value)#->todo too inefficent
		#combine features of parents to generate offspring
		print("parents[0] " + str(parents[0]) + " parents[1] " + str(parents[1]))
		offspring_a, offspring_b = crossover_func(parents[0], parents[1])
		offspring_a_save = str(offspring_a)
		offspring_b_save = str(offspring_b)

		offspring_a = mutation_func(offspring_a)
		offspring_b = mutation_func(offspring_b)

		#handle family register
		offspring_a_mutation=""
		if(str(offspring_a)!=offspring_a_save):
			offspring_a_mutation=str(offspring_a)
		offspring_b_mutation=""
		if(str(offspring_b)!=offspring_b_save):
			offspring_b_mutation=str(offspring_b)
		family_register+=offspring_a_save + " parents " + str(parents[0]) + "&" +  str(parents[1]) + " mutation " + str(offspring_a_mutation) + "\n"
		family_register+=offspring_b_save + " parents " + str(parents[0]) + "&" +  str(parents[1]) + " mutation " + str(offspring_b_mutation) + "\n"

		#add offspring to generation
		next_generation += [offspring_a, offspring_b]
	population = next_generation


	#find unique individuals
	individuals = list()
	for i in range(len(population)):
		if((population[i] in individuals)==False):
			individuals.append(population[i])
	unique_individuals = len(individuals)

	#fill rest of generation randomly
	missing_individuals = population_size-unique_individuals
	individuals_to_add = populate_func()[0:missing_individuals]
	family_register += "unique_individuals " + str(unique_individuals) + "\n"
	population = individuals
	population+=individuals_to_add

	#invoke evaluation of  new population
	population_for_fitness_eval = copy.deepcopy(population)
	for i in range(0, population_size):
		fitness_func(population_for_fitness_eval[i],generation,i)
	
	return population, family_register

