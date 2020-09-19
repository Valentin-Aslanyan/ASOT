

solar_Radius=1.0*7E10
target_R=1.01
target_theta=95.0		#degrees
target_phi=-10.0
file_directory="./"
flicks_file="flicks.0000000"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *


def is_block_in_limits(coord_logR,coord_theta,coord_phi,target_R,target_theta,target_phi):
	in_limits=False
	R_low=np.exp(min(coord_logR))
	R_high=np.exp(max(coord_logR))
	th_low=min(coord_theta)
	th_high=max(coord_theta)
	ph_low=min(coord_phi)
	ph_high=max(coord_phi)
	if R_low<=target_R and R_high>=target_R and th_low<=target_theta/180.0*np.pi and th_high>=target_theta/180.0*np.pi and ph_low<=target_phi/180.0*np.pi and ph_high>=target_phi/180.0*np.pi:
		in_limits=True
	return in_limits


def spherical_to_cartesian(R,theta,phi):
	x=R*np.sin(theta)*np.cos(phi+np.pi)
	y=R*np.sin(theta)*np.sin(phi+np.pi)
	z=R*np.cos(theta)
	return np.array([x,y,z])


def vertex_distance(v1,v2):
	return np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)


time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)

for idx in range(len(coord_logR)):
	if is_block_in_limits(coord_logR[idx],coord_theta[idx],coord_phi[idx],target_R,target_theta,target_phi):
		print("Delta R [Mm]:    ",abs(np.exp(coord_logR[idx][1])-np.exp(coord_logR[idx][0]))*solar_Radius*1E-8)	
		print("Delta theta [o]: ",abs(coord_theta[idx][1]-coord_theta[idx][0])/np.pi*180.0)
		print("Delta phi [o]:   ",abs(coord_phi[idx][1]-coord_phi[idx][0])/np.pi*180.0)
		print("")
		print("  A---------(1)----------B ")
		print("  |\                    /| ")
		print("  | \                  / | ")
		print("  | (2)             (3)  | ")
		print("  |   \              /   | ")
		print("  |    \            /    | ")
		print("  |     E----(4)---F     | ")
		print("  |     |          |     | ")
		print("  |     |          |     | ")
		print(" (5)   (6) Low R  (7)   (8)")
		print("  |     |          |     | ")
		print("  |     |          |     | ")
		print("  |     G----(9)---H     | ")
		print("  |    /            \    | ")
		print("  |   /              \   | ")
		print("  | (10)            (11) | ")
		print("  | /                  \ | ")
		print("  |/                    \| ")
		print("  C---------(12)---------D ")
		print("")
		print("Full block")
		R_low=np.exp(min(coord_logR[idx]))
		R_high=np.exp(max(coord_logR[idx]))
		th_low=min(coord_theta[idx])
		th_high=max(coord_theta[idx])
		ph_low=min(coord_phi[idx])
		ph_high=max(coord_phi[idx])
		vertex_A=spherical_to_cartesian(R_high,th_high,ph_low)*solar_Radius
		vertex_B=spherical_to_cartesian(R_high,th_high,ph_high)*solar_Radius
		vertex_C=spherical_to_cartesian(R_high,th_low,ph_low)*solar_Radius
		vertex_D=spherical_to_cartesian(R_high,th_low,ph_high)*solar_Radius
		vertex_E=spherical_to_cartesian(R_low,th_high,ph_low)*solar_Radius
		vertex_F=spherical_to_cartesian(R_low,th_high,ph_high)*solar_Radius
		vertex_G=spherical_to_cartesian(R_low,th_low,ph_low)*solar_Radius
		vertex_H=spherical_to_cartesian(R_low,th_low,ph_high)*solar_Radius
		print("(1):  {:.3e}".format(vertex_distance(vertex_A,vertex_B)*1E-8))
		print("(2):  {:.3e}".format(vertex_distance(vertex_A,vertex_E)*1E-8))
		print("(3):  {:.3e}".format(vertex_distance(vertex_B,vertex_F)*1E-8))
		print("(4):  {:.3e}".format(vertex_distance(vertex_E,vertex_F)*1E-8))
		print("(5):  {:.3e}".format(vertex_distance(vertex_A,vertex_C)*1E-8))
		print("(6):  {:.3e}".format(vertex_distance(vertex_E,vertex_G)*1E-8))
		print("(7):  {:.3e}".format(vertex_distance(vertex_F,vertex_H)*1E-8))
		print("(8):  {:.3e}".format(vertex_distance(vertex_B,vertex_D)*1E-8))
		print("(9):  {:.3e}".format(vertex_distance(vertex_G,vertex_H)*1E-8))
		print("(10): {:.3e}".format(vertex_distance(vertex_C,vertex_G)*1E-8))
		print("(11): {:.3e}".format(vertex_distance(vertex_D,vertex_H)*1E-8))
		print("(12): {:.3e}".format(vertex_distance(vertex_C,vertex_D)*1E-8))





