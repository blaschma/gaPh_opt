import tmoutproc as top
import sys

if __name__ == '__main__':

	coord = top.read_xyz_file(sys.argv[1])
	coord = top.x2t(coord)
	top.write_coord_file(sys.argv[2], coord)