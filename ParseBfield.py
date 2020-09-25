

file_directory="./"
bfield_file="bfield.0057883"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *


time,ntblks,nlblks,coord_logR,coord_theta,coord_phi,B=read_bfield_file(file_directory,bfield_file)
