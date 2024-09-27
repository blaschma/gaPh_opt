import codecs
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import tmoutproc as top
import sys
import configparser
from utils import find_atoms


if __name__ == '__main__':
	
	#argv[1] : coord path: where coord and fixed file is located and where result is stored
	#argv[2] : config path
	#argv[3] : output coord file
	#load coord file

	cfg = configparser.ConfigParser()
	cfg.read_file(codecs.open(sys.argv[2], "r", "utf8"))
	anchor_mode = int(cfg.get('Building Procedure', 'anchor'))

	coord = top.read_xyz_file(sys.argv[1])
	#find achors
	print("anchor mode",anchor_mode)
	if(anchor_mode==0):
		fixed_index = find_atoms(coord, 'au')
	elif(anchor_mode==1):
		fixed_index = find_atoms(coord, 's')
	else:
		raise ValueError('Unknown anchor mode')

	#rotate and shift
	coord = top.align_molecule(coord, [0,0,1], [fixed_index[0], fixed_index[1]])

	#write coord file
	top.write_xyz_file(sys.argv[3], coord)


