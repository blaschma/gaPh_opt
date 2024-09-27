import sys

import tmoutproc as top
import numpy as np
from utils import find_atoms


def bend_gold(coord):
	"""
	Places gold atoms along sulfur-sulfur axis
	Args:
		coord: coord file loaded with top.load_coord_file

	Returns: coord

	"""
	gold = find_atoms(coord, "au")
	sulfur = find_atoms(coord, "s")

	s_1_x = coord[1, sulfur[0]]
	s_1_y = coord[2, sulfur[0]]
	s_1_z = coord[3, sulfur[0]]

	s_2_x = coord[1, sulfur[1]]
	s_2_y = coord[2, sulfur[1]]
	s_2_z = coord[3, sulfur[1]]

	sulfur_axis = [s_2_x-s_1_x, s_2_y-s_1_y, s_2_z-s_1_z]
	sulfur_axis = sulfur_axis/np.linalg.norm(sulfur_axis)

	#au s bond length in ang
	au_s_length = 2.240000254

	#bend left gold
	coord[1, gold[0]] = s_1_x - sulfur_axis[0] * au_s_length
	coord[2, gold[0]] = s_1_y - sulfur_axis[1] * au_s_length
	coord[3, gold[0]] = s_1_z - sulfur_axis[2] * au_s_length

	#bend right gold
	coord[1, gold[1]] = s_2_x + sulfur_axis[0] * au_s_length
	coord[2, gold[1]] = s_2_y + sulfur_axis[1] * au_s_length
	coord[3, gold[1]] = s_2_z + sulfur_axis[2] * au_s_length

	return coord


if __name__ == '__main__':
	# argv[1] : coord path: where coord and fixed file is located and where result is stored
	# argv[2] : config path
	# argv[3] : output coord file

	# load coord file
	coord = top.read_xyz_file(sys.argv[1])

	coord = bend_gold(coord)

	top.write_xyz_file(sys.argv[3], coord)