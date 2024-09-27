import tmoutproc as top
import sys

__HYDROGEN_MASS__ = 1.00794080

def modify_mass_to_element(coord, indices, element):
	for index in indices:
		coord[3, index] = element
	return coord

def determine_indices(input_file):

	indices = list()
	with open(input_file, 'r', encoding='utf-8') as file:
		for line in file:
			if(line.startswith("modify mass")):
				start_index = line.index(': ')
				end_index = line.index(',')
				index = line[start_index+1:end_index]
				try:
					index_int = int(index)
				except ValueError as verr:
					print(f"Error casting {index} to int")
					exit(-1)
				#xtb starts counting at 1. Python at 0
				if(index_int<=0):
					print(f"Index {index_int} cannot be correct. Check your input file")
					exit(-1)
				index_int = index_int-1
				indices.append(index_int)
	if(len(indices) == 0):
		print("No indices found! Check your input file!")
	return indices



if __name__ == '__main__':

	coord = top.read_coord_file(sys.argv[1])
	indices = determine_indices(input_file=sys.argv[2])
	modify_mass_to_element(coord, indices, "h")
	top.write_coord_file(sys.argv[3], coord)



	