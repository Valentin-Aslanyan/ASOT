

solar_Radius=1.0*7E10
CNT_file="./arms.cnt"
start_directory="./PFLS/0044295/"
end_directory="./PFLS/0059184/"
connection_filename="./PFLS/Connection_18_0044295_0059184.npy"
start_time=4200.0
end_time=5300.0
delta_t=10.0


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *


save_connection_map(connection_filename,CNT_file,start_directory,end_directory,start_time,end_time,delta_t,solar_Radius)




