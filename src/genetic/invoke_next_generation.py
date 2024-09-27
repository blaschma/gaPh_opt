import run_generation as rg 
import sys

if __name__ == '__main__':
	# argv[1] = config path
	# argv[2] = calculation_path
	rg.next_generation(sys.argv[1], sys.argv[2])
