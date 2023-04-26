"""
Basic functions for reading and basic processing of ARMS data

Spherical coordinate definitions:
'polar' - textbook definition: theta in [0,pi) | theta=0 is north pole | phi in [-pi,pi)  
'control' - from arms.cnt: theta in [0,1) | theta=0 is north pole | phi in [-1,1)
'carrington' - Mercator map on the sun: theta in [-90,90) | theta=90 is north pole | phi in [-180,180)
'flicks' - as written to flicks files: theta in [-pi/2,pi/2) | theta=pi/2 is north pole | phi in [-pi,pi)
"""

import numpy as np
import struct
import os


DEG2RAD=np.pi/180.0
RAD2DEG=180.0/np.pi


def change_angular_coords(theta_in,phi_in,from_type=None,to_type=None):
	"""
	Supply input coordinates in float or numpy array and string for input/output type
	"""
	if from_type==None or to_type==None or type(from_type)!=str or type(to_type)!=str:
		print("Error, specify input/output coordinates. Choices: 'polar','control','carrington','flicks'")
		return theta_in,phi_in
	elif from_type.lower()==to_type.lower():
		return theta_in,phi_in
	elif from_type.lower()=='polar':
		if to_type.lower()=='control':
			theta_out=theta_in/np.pi
			phi_out=phi_in/np.pi
		elif to_type.lower()=='carrington':
			theta_out=90.0-theta_in*RAD2DEG
			phi_out=phi_in*RAD2DEG
		elif to_type.lower()=='flicks':
			theta_out=np.pi*0.5-theta_in
			phi_out=phi_in
		else:
			print("Error, output coordinate type not recognized. Choices: 'polar','control','carrington','flicks'")
			theta_out=theta_in
			phi_out=phi_in
		return theta_out,phi_out
	elif from_type.lower()=='control':
		if to_type.lower()=='polar':
			theta_out=np.pi*theta_in
			phi_out=np.pi*phi_in
		elif to_type.lower()=='carrington':
			theta_out=180.0*(0.5-theta_in)
			phi_out=180.0*phi_in
		elif to_type.lower()=='flicks':
			theta_out=np.pi*(0.5-theta_in)
			phi_out=np.pi*phi_in
		else:
			print("Error, output coordinate type not recognized. Choices: 'polar','control','carrington','flicks'")
			theta_out=theta_in
			phi_out=phi_in
		return theta_out,phi_out
	elif from_type.lower()=='carrington':
		if to_type.lower()=='polar':
			theta_out=np.pi*(0.5-theta_in/180.0)
			phi_out=DEG2RAD*phi_in
		elif to_type.lower()=='control':
			theta_out=0.5-theta_in/180.0
			phi_out=phi_in/180.0
		elif to_type.lower()=='flicks':
			theta_out=DEG2RAD*theta_in
			phi_out=DEG2RAD*phi_in
		else:
			print("Error, output coordinate type not recognized. Choices: 'polar','control','carrington','flicks'")
			theta_out=theta_in
			phi_out=phi_in
		return theta_out,phi_out
	elif from_type.lower()=='flicks':
		if to_type.lower()=='polar':
			theta_out=0.5-theta_in
			phi_out=phi_in
		elif to_type.lower()=='control':
			theta_out=0.5-theta_in/np.pi
			phi_out=phi_in/np.pi
		elif to_type.lower()=='carrington':
			theta_out=RAD2DEG*theta_in
			phi_out=RAD2DEG*phi_in
		else:
			print("Error, output coordinate type not recognized. Choices: 'polar','control','carrington','flicks'")
			theta_out=theta_in
			phi_out=phi_in
		return theta_out,phi_out
	else:
		print("Error, input coordinate type not recognized. Choices: 'polar','control','carrington','flicks'")
		return theta_in,phi_in


#Change flicks theta in [-pi/2,pi/2] --> [0,pi]
def read_flicks_file(file_directory,flicks_file):
	"""
	Requires header file in directory and selected flicks file
	nlblks - number of leaf blocks
	Outputs arrays of block coordinates (each one a pair of upper/lower limits) for leaf blocks and data at all stored grid points
	"""
	#Read header file
	hdrfile=open(os.path.join(file_directory,"flicks.hdr"),"r")
	for idx in range(2):
		hdrfile.readline()
	line=hdrfile.readline()
	if 'spherical exponential' in line.lower():
		flicks_type="SphE"
	elif 'cartesian' in line.lower():
		flicks_type="Cart"
	else:
		print('Warning! Coordinate type "'+line+'" not recognized')
		flicks_type="Unknown"
	line=hdrfile.readline()
	n1p=int(line[:2])
	line=hdrfile.readline()
	n2p=int(line[:2])
	line=hdrfile.readline()
	n3p=int(line[:2])
	for idx in range(2):
		hdrfile.readline()

	nvar=0
	for idx in range(12):
		line=hdrfile.readline()
		if len(line)>0:
			nvar+=int(line[:1])
	hdrfile.close()


	#Read main serial flicks file
	flicksfile=open(os.path.join(file_directory,flicks_file),"rb")

	#Use struct.unpack, big endian => use '>' at start of format string
	#Note! Every time Fortran write command is invoked, there are 4 bytes at the start and end with the number of bytes being written
	#Use file.read(4) or file.read(8) to skip this padding
	flicksfile.read(25)
	time=struct.unpack('>f', flicksfile.read(4))[0]
	flicksfile.read(8)
	ntblks=struct.unpack('>i', flicksfile.read(4))[0]
	nlblks=struct.unpack('>i', flicksfile.read(4))[0]
	newgrd=struct.unpack('>i', flicksfile.read(4))[0]
	flicksfile.read(4)

	if flicks_type=="SphE":
		coord_logR=np.zeros((nlblks,2))
		coord_theta=np.zeros((nlblks,2))
		coord_phi=np.zeros((nlblks,2))
		data=np.zeros((nlblks,n3p,n2p,n1p,nvar))

		idx_leaf=0
		for idx_blk in range(ntblks):
			flicksfile.read(4)
			iputwrk=struct.unpack('>'+21*'i', flicksfile.read(21*4))
			flicksfile.read(8)
			rputwrk=struct.unpack('>'+6*'f', flicksfile.read(6*4))
			flicksfile.read(4)
			if iputwrk[2]==1:
				coord_logR[idx_leaf,:]=[rputwrk[0],rputwrk[1]]
				coord_theta[idx_leaf,:]=[rputwrk[2]+np.pi*0.5,rputwrk[3]+np.pi*0.5]
				coord_phi[idx_leaf,:]=[rputwrk[4],rputwrk[5]]
				for idx_p in range(n3p):
					for idx_t in range(n2p):
						for idx_r in range(n1p):
							flicksfile.read(4)
							data[idx_leaf,idx_p,idx_t,idx_r,:]=struct.unpack('>'+nvar*'f', flicksfile.read(nvar*4))
							flicksfile.read(4)
				idx_leaf+=1


		flicksfile.close()
		return time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data

	elif flicks_type=="Cart":
		coord_x=np.zeros((nlblks,2))
		coord_y=np.zeros((nlblks,2))
		coord_z=np.zeros((nlblks,2))
		data=np.zeros((nlblks,n3p,n2p,n1p,nvar))

		idx_leaf=0
		for idx_blk in range(ntblks):
			flicksfile.read(4)
			iputwrk=struct.unpack('>'+21*'i', flicksfile.read(21*4))
			flicksfile.read(8)
			rputwrk=struct.unpack('>'+6*'f', flicksfile.read(6*4))
			flicksfile.read(4)
			if iputwrk[2]==1:
				coord_x[idx_leaf,:]=[rputwrk[0],rputwrk[1]]
				coord_y[idx_leaf,:]=[rputwrk[2],rputwrk[3]]
				coord_z[idx_leaf,:]=[rputwrk[4],rputwrk[5]]
				for idx_z in range(n3p):
					for idx_y in range(n2p):
						for idx_x in range(n1p):
							flicksfile.read(4)
							data[idx_leaf,idx_z,idx_y,idx_x,:]=struct.unpack('>'+nvar*'f', flicksfile.read(nvar*4))
							flicksfile.read(4)
				idx_leaf+=1


		flicksfile.close()
		return time,ntblks,nlblks,newgrd,coord_x,coord_y,coord_z,data
	else:
		return 0.0,1,1,0,np.zeros((1,2)),np.zeros((1,2)),np.zeros((1,2)),np.zeros((1,1,1,1,1))


def get_flicks_coordinate_type(file_directory):
	"""
	Requires header file in directory
	"""
	#Read header file
	hdrfile=open(os.path.join(file_directory,"flicks.hdr"),"r")
	for idx in range(2):
		hdrfile.readline()
	line=hdrfile.readline()
	if 'spherical exponential' in line.lower():
		flicks_type="SphE"
	elif 'cartesian' in line.lower():
		flicks_type="Cart"
	else:
		print('Warning! Coordinate type "'+line+'" not recognized')
		flicks_type="Unknown"
	hdrfile.close()
	return flicks_type


#Dummy function for replicating header with added variable
def duplicate_flicks_header(input_directory,output_directory,new_variable_name):
	infile=open(os.path.join(input_directory,"flicks.hdr"),"r")
	outfile=open(os.path.join(output_directory,"flicks.hdr"),"w")
	for line in infile:
		print(line,end="",file=outfile)
	infile.close()
	print(new_variable_name+"\n",end="",file=outfile)
	outfile.close()


#Dummy function for replicating footer with added variable
def duplicate_flicks_footer(input_directory,output_directory,new_variable_name,new_variable_min,new_variable_max):
	infile=open(os.path.join(input_directory,"flicks.ftr"),"r")
	outfile=open(os.path.join(output_directory,"flicks.ftr"),"w")
	for line in infile:
		print(line,end="",file=outfile)
	infile.close()
	print("{: .5E}Min ".format(new_variable_min)+new_variable_name+"\n",end="",file=outfile)
	print("{: .5E}Max ".format(new_variable_max)+new_variable_name+"\n",end="",file=outfile)
	outfile.close()


#Duplicate flicks file while appending/overwriting data inside; filenames must be full paths; ensure header is correct
def duplicate_flicks_file(input_filename,output_filename,nvar_old,data):
	nlblks,n3p,n2p,n1p,nvar_new=np.shape(data)
	infile=open(input_filename,"rb")
	outfile=open(output_filename,"wb")
	outfile.write(infile.read(37))
	ntblks=struct.unpack('>i', infile.read(4))[0]
	outfile.write(struct.pack('>i',ntblks))
	outfile.write(infile.read(12))

	idx_leaf=0
	for idx_blk in range(ntblks):
		outfile.write(infile.read(4))
		iputwrk_raw=infile.read(21*4)
		outfile.write(iputwrk_raw)
		iputwrk=struct.unpack('>'+21*'i', iputwrk_raw)
		outfile.write(infile.read(36))
		if iputwrk[2]==1:
			for idx_p in range(n3p):
				for idx_t in range(n2p):
					for idx_r in range(n1p):
						infile.read(8+nvar_old*4)
						outfile.write(struct.pack('>i',nvar_new*4))
						outfile.write(struct.pack('>'+nvar_new*'f',*data[idx_leaf,idx_p,idx_t,idx_r,:]))
						outfile.write(struct.pack('>i',nvar_new*4))
			idx_leaf+=1
	infile.close()
	outfile.close()


def get_flicks_time(flicks_file):
	"""
	Simulation time from flicks file
	"""
	flicksfile=open(flicks_file,"rb")
	flicksfile.read(25)
	time=struct.unpack('>f', flicksfile.read(4))[0]
	flicksfile.close()
	return time


def get_flicks_grid_dimensions(data):
	"""
	Dimensions from a data array read in from flicks or bfield
	for bfield nvar==3
	"""
	nlblk=np.shape(data)[0]
	n1p=np.shape(data)[3]
	n2p=np.shape(data)[2]
	n3p=np.shape(data)[1]
	nvar=np.shape(data)[4]
	return nlblk,n1p,n2p,n3p,nvar


def read_bfield_file(file_directory,bfield_file):
	"""
	Requires header file in directory and selected bfield file
	nlblks - number of leaf blocks
	Outputs arrays of block coordinates (each one a pair of upper/lower limits) for leaf blocks and data at all stored grid points
	"""
	#Read header file
	hdrfile=open(os.path.join(file_directory,"bfield.hdr"),"r")
	hdrfile.readline()
	line=hdrfile.readline()
	if 'spherical exponential' not in line.lower():
		print('Warning! Coordinate type "'+line+'" not recognized')
	line=hdrfile.readline()
	n1p=int(line[:2])
	line=hdrfile.readline()
	n2p=int(line[:2])
	line=hdrfile.readline()
	n3p=int(line[:2])
	hdrfile.close()

	#Read main bfield file
	bfieldfile=open(os.path.join(file_directory,bfield_file),"rb")
	bfieldfile.read(4)
	time=struct.unpack('>d', bfieldfile.read(8))[0]
	bfieldfile.read(8)
	ntblks=struct.unpack('>i', bfieldfile.read(4))[0]
	nlblks=struct.unpack('>i', bfieldfile.read(4))[0]
	bfieldfile.read(4)

	coord_logR=np.zeros((nlblks,2))
	coord_theta=np.zeros((nlblks,2))
	coord_phi=np.zeros((nlblks,2))
	data=np.zeros((nlblks,n3p,n2p,n1p,3))

	idx_leaf=0
	for idx in range(ntblks):
		bfieldfile.read(4)
		iputwrk=struct.unpack('>'+35*'i', bfieldfile.read(35*4))
		bfieldfile.read(8)
		rputwrk=struct.unpack('>'+6*'d', bfieldfile.read(6*8))
		bfieldfile.read(4)
		if iputwrk[2]==1:
			coord_logR[idx_leaf,:]=[rputwrk[0],rputwrk[1]]
			coord_theta[idx_leaf,:]=[rputwrk[2]+np.pi*0.5,rputwrk[3]+np.pi*0.5]
			coord_phi[idx_leaf,:]=[rputwrk[4],rputwrk[5]]
			for idx_p in range(n3p):
				for idx_t in range(n2p):
					for idx_r in range(n1p):
						bfieldfile.read(4)
						data[idx_leaf,idx_p,idx_t,idx_r,:]=struct.unpack('>'+3*'d', bfieldfile.read(24))
						bfieldfile.read(4)
			idx_leaf+=1


	bfieldfile.close()
	return time,ntblks,nlblks,coord_logR,coord_theta,coord_phi,data


def is_block_in_limits(coord_logR,coord_theta,coord_phi,target_phi,R_limits):
	"""
	Check if a block straddles a given target_phi
	R_limits can be 2-element array or None 
	"""
	in_limits=False
	if R_limits!=None:
		if coord_phi[0]<=target_phi/180.0*np.pi and coord_phi[1]>=target_phi/180.0*np.pi and np.exp(coord_logR[1])<=R_limits[1] and np.exp(coord_logR[0])>=R_limits[0]:
			in_limits=True
	else:
		if coord_phi[0]<=target_phi/180.0*np.pi and coord_phi[1]>=target_phi/180.0*np.pi:
			in_limits=True
	return in_limits


def interp_pointpair(x,x0,x1,y0,y1):
	"""
	Linear interpolation
	"""
	if x0==x1:
		return y0
	else:
		return y0+(x-x0)*(y1-y0)/(x1-x0)


def read_fcube_file(infile_name,data_indices='all'):
	"""
	fcube - flicks data interpolated to regular grid
	data_indices must be a list or 'all'
	data_indices corresponds to elements of flicks data, e.g. 0 is rho usually
	"""
	infile=open(infile_name,"rb")
	n1p=struct.unpack('>i', infile.read(4))[0]
	n2p=struct.unpack('>i', infile.read(4))[0]
	n3p=struct.unpack('>i', infile.read(4))[0]
	nvar=struct.unpack('>i', infile.read(4))[0]
	if type(data_indices)==str:
		if data_indices.lower()=='all':
			data_indices_actual=list(range(nvar))
		else:
			data_indices_actual=[]
	elif type(data_indices)==list:
		data_indices_actual=[]
		for idx in data_indices:
			if type(idx)==int and idx<nvar:
				data_indices_actual.append(idx)
			elif type(idx)==float and int(idx)<nvar:
				data_indices_actual.append(int(idx))
	else:
		data_indices_actual=[]
	gridx=np.array(struct.unpack('>'+n1p*'f', infile.read(n1p*4)))
	gridy=np.array(struct.unpack('>'+n2p*'f', infile.read(n2p*4)))
	gridz=np.array(struct.unpack('>'+n3p*'f', infile.read(n3p*4)))
	grid_data=np.zeros((n1p,n2p,n3p,len(data_indices_actual)),dtype=np.float32)
	idx_v2=0
	for idx_v in range(nvar):
		for idx3 in range(n3p):
			for idx2 in range(n2p):
				data_raw=infile.read(n1p*4)
				if idx_v in data_indices_actual:
					grid_data[:,idx2,idx3,idx_v2]=struct.unpack('>'+n1p*'f',data_raw)
		if idx_v in data_indices_actual:
			idx_v2+=1
	infile.close()
	return gridx,gridy,gridz,grid_data


def read_bcube_directory(target_directory,timestep):
	"""
	bcube - bfield data interpolated to regular grid
	"""
	coord_R=np.fromfile(os.path.join(target_directory,'xcoords'), dtype='>f8').astype('float64')
	coord_th=np.fromfile(os.path.join(target_directory,'ycoords'), dtype='>f8').astype('float64')
	coord_ph=np.fromfile(os.path.join(target_directory,'zcoords'), dtype='>f8').astype('float64')
	bcube=np.fromfile(os.path.join(target_directory,"bcube."+str(timestep)), dtype='>f8').reshape((len(coord_R), len(coord_th), len(coord_ph), 3), order='F')
	return coord_R,coord_th,coord_ph,bcube


def phi_slice(target_phi,coord_logR,coord_theta,coord_phi,data):
	"""
	Get data in a single plane of given phi
	target_phi must be in degrees, e.g. Carrington coordinates
	Routine will find blocks which straddle target_phi and interpolate linearly
	"""
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	target_phi_actual=target_phi*DEG2RAD
	num_blocks=len(coord_phi[:,0])
	bounding_lblk=[]
	for idx in range(num_blocks):
		if coord_phi[idx,0]<coord_phi[idx,1]:
			if coord_phi[idx,0]<=target_phi_actual and coord_phi[idx,1]>=target_phi_actual:
				bounding_lblk.append(idx)
		else:
			if coord_phi[idx,0]>=target_phi_actual and coord_phi[idx,1]<=target_phi_actual:
				bounding_lblk.append(idx)

	new_logR=np.zeros((len(bounding_lblk)*n2p*n1p))
	new_theta=np.zeros((len(bounding_lblk)*n2p*n1p))
	new_data=np.zeros((len(bounding_lblk)*n2p*n1p,nvar))

	for idx in range(len(bounding_lblk)):
		idx_phi=int(np.floor((target_phi_actual-coord_phi[bounding_lblk[idx],0])/(coord_phi[bounding_lblk[idx],1]-coord_phi[bounding_lblk[idx],0])*(n3p-1)))
		if idx_phi==n3p-1:
			idx_phi=n3p-2
		for idx_t in range(n2p):
			for idx_r in range(n1p):
				new_logR[idx*n2p*n1p+idx_t*n1p+idx_r]=idx_r/(n1p-1)*(coord_logR[bounding_lblk[idx],1]-coord_logR[bounding_lblk[idx],0])+coord_logR[bounding_lblk[idx],0]
				new_theta[idx*n2p*n1p+idx_t*n1p+idx_r]=idx_t/(n2p-1)*(coord_theta[bounding_lblk[idx],1]-coord_theta[bounding_lblk[idx],0])+coord_theta[bounding_lblk[idx],0]
				phi_low=idx_phi/(n3p-1)*(coord_phi[bounding_lblk[idx],1]-coord_phi[bounding_lblk[idx],0])+coord_phi[bounding_lblk[idx],0]
				phi_high=(idx_phi+1)/(n3p-1)*(coord_phi[bounding_lblk[idx],1]-coord_phi[bounding_lblk[idx],0])+coord_phi[bounding_lblk[idx],0]
				for idx_v in range(nvar):
					new_data[idx*n2p*n1p+idx_t*n1p+idx_r,idx_v]=interp_pointpair(target_phi_actual,phi_low,phi_high,data[bounding_lblk[idx],idx_phi,idx_t,idx_r,idx_v], data[bounding_lblk[idx],idx_phi+1,idx_t,idx_r,idx_v])
	return new_logR,new_theta,new_data


def R_slice(target_R,coord_logR,coord_theta,coord_phi,data):
	"""
	Get data in a single plane of given phi
	target_R must be in linear space (as opposed to logR)
	Routine will find blocks which straddle target_R and interpolate linearly
	"""
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	target_logR_actual=np.log(target_R)
	num_blocks=len(coord_logR[:,0])
	bounding_lblk=[]
	for idx in range(num_blocks):
		if coord_logR[idx,0]<coord_logR[idx,1]:
			if coord_logR[idx,0]<=target_logR_actual and coord_logR[idx,1]>=target_logR_actual:
				bounding_lblk.append(idx)
		else:
			if coord_logR[idx,0]>=target_logR_actual and coord_logR[idx,1]<=target_logR_actual:
				bounding_lblk.append(idx)

	new_theta=np.zeros((len(bounding_lblk)*n3p*n2p))
	new_phi=np.zeros((len(bounding_lblk)*n3p*n2p))
	new_data=np.zeros((len(bounding_lblk)*n3p*n2p,nvar))

	for idx in range(len(bounding_lblk)):
		idx_logR=int(np.floor((target_logR_actual-coord_logR[bounding_lblk[idx],0])/(coord_logR[bounding_lblk[idx],1]-coord_logR[bounding_lblk[idx],0])*(n1p-1)))
		if idx_logR==n1p-1:
			idx_logR=n1p-2
		for idx_p in range(n3p):
			for idx_t in range(n2p):
				new_theta[idx*n3p*n2p+idx_p*n2p+idx_t]=idx_t/(n2p-1)*(coord_theta[bounding_lblk[idx],1]-coord_theta[bounding_lblk[idx],0])+coord_theta[bounding_lblk[idx],0]
				new_phi[idx*n3p*n2p+idx_p*n2p+idx_t]=idx_p/(n3p-1)*(coord_phi[bounding_lblk[idx],1]-coord_phi[bounding_lblk[idx],0])+coord_phi[bounding_lblk[idx],0]
				logR_low=idx_logR/(n1p-1)*(coord_logR[bounding_lblk[idx],1]-coord_logR[bounding_lblk[idx],0])+coord_logR[bounding_lblk[idx],0]
				logR_high=(idx_logR+1)/(n1p-1)*(coord_logR[bounding_lblk[idx],1]-coord_logR[bounding_lblk[idx],0])+coord_logR[bounding_lblk[idx],0]
				for idx_v in range(nvar):
					new_data[idx*n3p*n2p+idx_p*n2p+idx_t,idx_v]=interp_pointpair(target_logR_actual,logR_low,logR_high,data[bounding_lblk[idx],idx_p,idx_t,idx_logR,idx_v], data[bounding_lblk[idx],idx_p,idx_t,idx_logR+1,idx_v])
	return new_theta,new_phi,new_data


def z_slice(target_z,coord_x,coord_y,coord_z,data):
	"""
	Get data in a single plane of given z
	Routine will find blocks which straddle target_z and interpolate linearly
	"""
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	num_blocks=len(coord_x[:,0])
	bounding_lblk=[]
	for idx in range(num_blocks):
		if coord_z[idx,0]<coord_z[idx,1]:
			if coord_z[idx,0]<=target_z and coord_z[idx,1]>=target_z:
				bounding_lblk.append(idx)
		else:
			if coord_z[idx,0]>=target_z and coord_z[idx,1]<=target_z:
				bounding_lblk.append(idx)

	new_x=np.zeros((len(bounding_lblk)*n1p*n2p))
	new_y=np.zeros((len(bounding_lblk)*n1p*n2p))
	new_data=np.zeros((len(bounding_lblk)*n1p*n2p,nvar))

	for idx in range(len(bounding_lblk)):
		idx_z=int(np.floor((target_z-coord_z[bounding_lblk[idx],0])/(coord_z[bounding_lblk[idx],1]-coord_z[bounding_lblk[idx],0])*(n3p-1)))
		if idx_z==n3p-1:
			idx_z=n3p-2
		for idx_y in range(n2p):
			for idx_x in range(n1p):
				new_x[idx*n1p*n2p+idx_x*n2p+idx_y]=idx_x/(n1p-1)*(coord_x[bounding_lblk[idx],1]-coord_x[bounding_lblk[idx],0])+coord_x[bounding_lblk[idx],0]
				new_y[idx*n1p*n2p+idx_x*n2p+idx_y]=idx_y/(n2p-1)*(coord_y[bounding_lblk[idx],1]-coord_y[bounding_lblk[idx],0])+coord_y[bounding_lblk[idx],0]
				z_low=idx_z/(n3p-1)*(coord_z[bounding_lblk[idx],1]-coord_z[bounding_lblk[idx],0])+coord_z[bounding_lblk[idx],0]
				z_high=(idx_z+1)/(n3p-1)*(coord_z[bounding_lblk[idx],1]-coord_z[bounding_lblk[idx],0])+coord_z[bounding_lblk[idx],0]
				for idx_v in range(nvar):
					new_data[idx*n1p*n2p+idx_x*n2p+idx_y,idx_v]=interp_pointpair(target_z,z_low,z_high,data[bounding_lblk[idx],idx_z,idx_y,idx_x,idx_v], data[bounding_lblk[idx],idx_z,idx_y,idx_z+1,idx_v])
	return new_x,new_y,new_data


def get_full_unstructured_cart_grid_noedge(coord_logR,coord_theta,coord_phi,data):
	"""
	Returns all flicks data in a large, unstructured Cartesian grid, e.g. X[i],Y[i],Z[i],rho[i] for all grid points
	Do not save the upper edge points (i.e. very highest r, theta, phi points)
	Otherwise, the redundant points at the upper edge of every block, which overlap with other blocks, must not be included
	"""
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	new_logR=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_theta=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_phi=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_X=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_Y=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_Z=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_data=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1),nvar))

	for idx in range(nlblk):
		for idx_p in range(n3pm1):
			for idx_t in range(n2pm1):
				for idx_r in range(n1pm1):
					new_logR[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r]=idx_r/n1pm1*(coord_logR[idx,1]-coord_logR[idx,0])+coord_logR[idx,0]
					new_theta[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r]=idx_t/n2pm1*(coord_theta[idx,1]-coord_theta[idx,0])+coord_theta[idx,0]
					new_phi[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r]=idx_p/n3pm1*(coord_phi[idx,1]-coord_phi[idx,0])+coord_phi[idx,0]
					for idx_v in range(nvar):
						new_data[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r,idx_v]=data[idx,idx_p,idx_t,idx_r,idx_v]
	
	for idx in range(nlblk*(n3p-1)*(n2p-1)*(n1p-1)):
		R=np.exp(new_logR[idx])
		sin_theta=np.sin(new_theta[idx])
		new_X[idx]=R*sin_theta*np.cos(new_phi[idx])
		new_Y[idx]=R*sin_theta*np.sin(new_phi[idx])
		new_Z[idx]=R*np.cos(new_theta[idx])

	return new_logR,new_theta,new_phi,new_X,new_Y,new_Z,new_data


def get_min_sep_sph_grid_noedge(coord_logR,coord_theta,coord_phi,data):
	"""
	Not finished
	Get regular rho grid using Python routines
	"""
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	new_logR=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_theta=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_phi=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_X=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_Y=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_Z=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1)))
	new_data=np.zeros((nlblk*(n3p-1)*(n2p-1)*(n1p-1),nvar))
	min_sep_logR=abs(coord_logR[0,1]-coord_logR[0,0])
	min_sep_theta=abs(coord_theta[0,1]-coord_theta[0,0])
	min_sep_phi=abs(coord_phi[0,1]-coord_phi[0,0])

	for idx in range(nlblk):
		for idx_p in range(n3pm1):
			for idx_t in range(n2pm1):
				for idx_r in range(n1pm1):
					new_logR[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r]=idx_r/n1pm1*(coord_logR[idx,1]-coord_logR[idx,0])+coord_logR[idx,0]
					new_theta[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r]=idx_t/n2pm1*(coord_theta[idx,1]-coord_theta[idx,0])+coord_theta[idx,0]
					new_phi[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r]=idx_p/n3pm1*(coord_phi[idx,1]-coord_phi[idx,0])+coord_phi[idx,0]
					min_sep_logR=min(min_sep_logR,abs(coord_logR[idx,1]-coord_logR[idx,0]))
					min_sep_theta=min(min_sep_theta,abs(coord_theta[idx,1]-coord_theta[idx,0]))
					min_sep_phi=min(min_sep_phi,abs(coord_phi[idx,1]-coord_phi[idx,0]))
					for idx_v in range(nvar):
						new_data[idx*n3pm1*n2pm1*n1pm1+idx_p*n2pm1*n1pm1+idx_t*n1pm1+idx_r,idx_v]=data[idx,idx_p,idx_t,idx_r,idx_v]
	

	min_sep_logR/=n1pm1
	min_sep_theta/=n2pm1
	min_sep_phi/=n3pm1
	final_logR=np.linspace(min(new_logR),max(new_logR),num=(max(new_logR)-min(new_logR))/min_sep_logR)
	final_theta=np.linspace(min(new_theta),max(new_theta),num=(max(new_theta)-min(new_theta))/min_sep_theta)
	final_phi=np.linspace(min(new_phi),max(new_phi),num=(max(new_phi)-min(new_phi))/min_sep_phi)
		
	from scipy.interpolate import griddata
	print(len(final_logR),len(final_theta),len(final_phi))
	grid_logR3,grid_theta3,grid_phi3=np.meshgrid(final_logR,final_theta,final_phi)
	logRho_grid=griddata((new_logR,new_theta,new_phi), np.log(new_data[:,0]), (grid_logR3,grid_theta3,grid_phi3), method='linear')#,fill_value=min(logRho_cart))
	np.save("Regular_grid_logR.np",final_logR)
	np.save("Regular_grid_theta.np",final_theta)
	np.save("Regular_grid_phi.np",final_phi)
	np.save("Regular_grid_rho.np",logRho_grid)
	
	return grid_logR3,grid_theta3,grid_phi3,logRho_grid


def Trilinear_interpolation(x,y,z,x0,x1,y0,y1,z0,z1,f000,f001,f010,f011,f100,f101,f110,f111):
	"""
	Must be a cuboid
	f000=f(x0,y0,z0)
	f100=f(x1,y0,z0) etc
	"""
	x_shift=x-x0
	y_shift=y-y0
	z_shift=z-z0
	delta_x=x1-x0
	delta_y=y1-y0
	delta_z=z1-z0
	a1=(f100-f000)/delta_x
	a2=(f010-f000)/delta_y
	a3=(f001-f000)/delta_z
	a4=(f110-f100-f010+f000)/delta_y/delta_x
	a5=(f110-f100-f001+f000)/delta_z/delta_x
	a6=(f011-f001-f010+f000)/delta_z/delta_y
	a7=(f111-f110-f101-f011+f001+f010+f100-f000)/delta_x/delta_y/delta_z
	return f000+a1*x_shift+a2*y_shift+a3*z_shift+a4*x_shift*y_shift+a5*x_shift*z_shift+a6*y_shift*z_shift+a7*x_shift*y_shift*z_shift


def ref_idx_search(arr,target_value,prev_idx):
	"""
	Find index of value in array
	Assumes array arr is monotonically increasing
	"""
	idx=prev_idx
	if arr[idx]<=target_value:
		while idx<len(arr)-2 and arr[idx+1]<=target_value:
			idx+=1
	else:
		while idx>0 and arr[idx-1]>target_value:
			idx-=1
	return idx


def o1_field_trace(pos_in,X,Y,Z,B,dt):
	"""
	First order field line integration, Cartesian coordinates
	Probably doesn't work
	"""

	delta_X=X[1]-X[0]
	delta_Y=Y[1]-Y[0]
	delta_Z=Z[1]-Z[0]
	
	x_0=pos_in[0]
	y_0=pos_in[1]
	z_0=pos_in[2]
	idx_x=min(max(int(np.floor((x_0-X[0])/delta_X)),0),len(X)-2)
	idx_y=min(max(int(np.floor((y_0-Y[0])/delta_Y)),0),len(Y)-2)
	idx_z=min(max(int(np.floor((z_0-Z[0])/delta_Z)),0),len(Z)-2)
	B_1=Trilinear_interpolation(x_0,y_0,z_0,X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_x,idx_y,idx_z,:], B[idx_x,idx_y,idx_z+1,:], B[idx_x,idx_y+1,idx_z,:], B[idx_x,idx_y+1,idx_z+1,:], B[idx_x+1,idx_y,idx_z,:], B[idx_x+1,idx_y,idx_z+1,:], B[idx_x+1,idx_y+1,idx_z,:], B[idx_x+1,idx_y+1,idx_z+1,:])

	return pos_in+dt*B_1/np.sqrt(B_1[0]*B_1[0]+B_1[1]*B_1[1]+B_1[2]*B_1[2])


def rk4_field_trace_cartesian(pos_in,X,Y,Z,B,dt):
	"""
	RK4 field line integration, Cartesian coordinates
	"""

	delta_X=X[1]-X[0]
	delta_Y=Y[1]-Y[0]
	delta_Z=Z[1]-Z[0]
	
	x_0=pos_in[0]
	y_0=pos_in[1]
	z_0=pos_in[2]
	idx_x=min(max(int(np.floor((x_0-X[0])/delta_X)),0),len(X)-2)
	idx_y=min(max(int(np.floor((y_0-Y[0])/delta_Y)),0),len(Y)-2)
	idx_z=min(max(int(np.floor((z_0-Z[0])/delta_Z)),0),len(Z)-2)
	B_1=Trilinear_interpolation(x_0,y_0,z_0,X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_x,idx_y,idx_z,:], B[idx_x,idx_y,idx_z+1,:], B[idx_x,idx_y+1,idx_z,:], B[idx_x,idx_y+1,idx_z+1,:], B[idx_x+1,idx_y,idx_z,:], B[idx_x+1,idx_y,idx_z+1,:], B[idx_x+1,idx_y+1,idx_z,:], B[idx_x+1,idx_y+1,idx_z+1,:])
	k_1=dt*B_1/np.sqrt(B_1[0]*B_1[0]+B_1[1]*B_1[1]+B_1[2]*B_1[2])

	x_1=pos_in[0]+k_1[0]*0.5
	y_1=pos_in[1]+k_1[1]*0.5
	z_1=pos_in[2]+k_1[2]*0.5
	idx_x=min(max(int(np.floor((x_1-X[0])/delta_X)),0),len(X)-2)
	idx_y=min(max(int(np.floor((y_1-Y[0])/delta_Y)),0),len(Y)-2)
	idx_z=min(max(int(np.floor((z_1-Z[0])/delta_Z)),0),len(Z)-2)
	B_2=Trilinear_interpolation(x_1,y_1,z_1,X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_x,idx_y,idx_z,:], B[idx_x,idx_y,idx_z+1,:], B[idx_x,idx_y+1,idx_z,:], B[idx_x,idx_y+1,idx_z+1,:], B[idx_x+1,idx_y,idx_z,:], B[idx_x+1,idx_y,idx_z+1,:], B[idx_x+1,idx_y+1,idx_z,:], B[idx_x+1,idx_y+1,idx_z+1,:])
	k_2=dt*B_2/np.sqrt(B_2[0]*B_2[0]+B_2[1]*B_2[1]+B_2[2]*B_2[2])
	
	x_2=pos_in[0]+k_2[0]*0.5
	y_2=pos_in[1]+k_2[1]*0.5
	z_2=pos_in[2]+k_2[2]*0.5
	idx_x=min(max(int(np.floor((x_2-X[0])/delta_X)),0),len(X)-2)
	idx_y=min(max(int(np.floor((y_2-Y[0])/delta_Y)),0),len(Y)-2)
	idx_z=min(max(int(np.floor((z_2-Z[0])/delta_Z)),0),len(Z)-2)
	B_3=Trilinear_interpolation(x_2,y_2,z_2,X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_x,idx_y,idx_z,:], B[idx_x,idx_y,idx_z+1,:], B[idx_x,idx_y+1,idx_z,:], B[idx_x,idx_y+1,idx_z+1,:], B[idx_x+1,idx_y,idx_z,:], B[idx_x+1,idx_y,idx_z+1,:], B[idx_x+1,idx_y+1,idx_z,:], B[idx_x+1,idx_y+1,idx_z+1,:])
	k_3=dt*B_3/np.sqrt(B_3[0]*B_3[0]+B_3[1]*B_3[1]+B_3[2]*B_3[2])

	x_3=pos_in[0]+k_3[0]
	y_3=pos_in[1]+k_3[1]
	z_3=pos_in[2]+k_3[2]
	idx_x=min(max(int(np.floor((x_3-X[0])/delta_X)),0),len(X)-2)
	idx_y=min(max(int(np.floor((y_3-Y[0])/delta_Y)),0),len(Y)-2)
	idx_z=min(max(int(np.floor((z_3-Z[0])/delta_Z)),0),len(Z)-2)
	B_4=Trilinear_interpolation(x_3,y_3,z_3,X[idx_x],X[idx_x+1],Y[idx_y],Y[idx_y+1],Z[idx_z],Z[idx_z+1],B[idx_x,idx_y,idx_z,:], B[idx_x,idx_y,idx_z+1,:], B[idx_x,idx_y+1,idx_z,:], B[idx_x,idx_y+1,idx_z+1,:], B[idx_x+1,idx_y,idx_z,:], B[idx_x+1,idx_y,idx_z+1,:], B[idx_x+1,idx_y+1,idx_z,:], B[idx_x+1,idx_y+1,idx_z+1,:])
	k_4=dt*B_4/np.sqrt(B_4[0]*B_4[0]+B_4[1]*B_4[1]+B_4[2]*B_4[2])

	return pos_in+(k_1+2.0*k_2+2.0*k_3+k_4)/6.0


def rk4_field_trace_spherical(pos_in,R,theta,phi,B,dt,idx_r_p,idx_t_p,idx_p_p):
	"""
	RK4 field line integration, spherical coordinates, 'polar' coordinate system
	Currently sub-optimal due to switching back and forth to Cartesian coordinates
	TODO - fully spherical
	"""

	r_0=pos_in[0]
	t_0=pos_in[1]
	p_0=pos_in[2]
	sin_th=np.sin(t_0)
	cos_th=np.cos(t_0)
	sin_ph=np.sin(p_0)
	cos_ph=np.cos(p_0)
	x_0=r_0*sin_th*cos_ph
	y_0=r_0*sin_th*sin_ph
	z_0=r_0*cos_th
	if t_0>np.pi:
		t_0-=np.pi
	elif t_0<0.0:
		t_0+=np.pi
	if p_0>np.pi:
		p_0-=2.0*np.pi
	elif p_0<-np.pi:
		p_0+=2.0*np.pi
	idx_r=ref_idx_search(R,r_0,idx_r_p)
	idx_t=ref_idx_search(theta,t_0,idx_t_p)
	idx_p=ref_idx_search(phi,p_0,idx_p_p)
	B_0=Trilinear_interpolation(r_0,t_0,p_0,R[idx_r],R[idx_r+1],theta[idx_t],theta[idx_t+1],phi[idx_p],phi[idx_p+1],B[idx_r,idx_t,idx_p,:], B[idx_r,idx_t,idx_p+1,:], B[idx_r,idx_t+1,idx_p,:], B[idx_r,idx_t+1,idx_p+1,:], B[idx_r+1,idx_t,idx_p,:], B[idx_r+1,idx_t,idx_p+1,:], B[idx_r+1,idx_t+1,idx_p,:], B[idx_r+1,idx_t+1,idx_p+1,:])
	B_0=B_0/np.sqrt(B_0[0]*B_0[0]+B_0[1]*B_0[1]+B_0[2]*B_0[2])
	k_1_x=dt*(B_0[0]*sin_th*cos_ph+B_0[1]*cos_th*cos_ph-B_0[2]*sin_ph)
	k_1_y=dt*(B_0[0]*sin_th*sin_ph+B_0[1]*cos_th*sin_ph+B_0[2]*cos_ph)
	k_1_z=dt*(B_0[0]*cos_th-B_0[1]*sin_th)

	x_1=x_0+k_1_x*0.5
	y_1=y_0+k_1_y*0.5
	z_1=z_0+k_1_z*0.5
	r_1=np.sqrt(x_1**2+y_1**2+z_1**2)
	t_1=np.arccos(z_1/r_1)
	p_1=np.arctan2(y_1,x_1)
	sin_th=np.sin(t_1)
	cos_th=np.cos(t_1)
	sin_ph=np.sin(p_1)
	cos_ph=np.cos(p_1)
	idx_r=ref_idx_search(R,r_1,idx_r)
	idx_t=ref_idx_search(theta,t_1,idx_t)
	idx_p=ref_idx_search(phi,p_1,idx_p)
	B_1=Trilinear_interpolation(r_1,t_1,p_1,R[idx_r],R[idx_r+1],theta[idx_t],theta[idx_t+1],phi[idx_p],phi[idx_p+1],B[idx_r,idx_t,idx_p,:], B[idx_r,idx_t,idx_p+1,:], B[idx_r,idx_t+1,idx_p,:], B[idx_r,idx_t+1,idx_p+1,:], B[idx_r+1,idx_t,idx_p,:], B[idx_r+1,idx_t,idx_p+1,:], B[idx_r+1,idx_t+1,idx_p,:], B[idx_r+1,idx_t+1,idx_p+1,:])
	B_1=B_1/np.sqrt(B_1[0]*B_1[0]+B_1[1]*B_1[1]+B_1[2]*B_1[2])
	k_2_x=dt*(B_1[0]*sin_th*cos_ph+B_1[1]*cos_th*cos_ph-B_1[2]*sin_ph)
	k_2_y=dt*(B_1[0]*sin_th*sin_ph+B_1[1]*cos_th*sin_ph+B_1[2]*cos_ph)
	k_2_z=dt*(B_1[0]*cos_th-B_1[1]*sin_th)
	
	x_2=x_0+k_2_x*0.5
	y_2=y_0+k_2_y*0.5
	z_2=z_0+k_2_z*0.5
	r_2=np.sqrt(x_2**2+y_2**2+z_2**2)
	t_2=np.arccos(z_2/r_2)
	p_2=np.arctan2(y_2,x_2)
	sin_th=np.sin(t_2)
	cos_th=np.cos(t_2)
	sin_ph=np.sin(p_2)
	cos_ph=np.cos(p_2)
	idx_r=ref_idx_search(R,r_2,idx_r)
	idx_t=ref_idx_search(theta,t_2,idx_t)
	idx_p=ref_idx_search(phi,p_2,idx_p)
	B_2=Trilinear_interpolation(r_2,t_2,p_2,R[idx_r],R[idx_r+1],theta[idx_t],theta[idx_t+1],phi[idx_p],phi[idx_p+1],B[idx_r,idx_t,idx_p,:], B[idx_r,idx_t,idx_p+1,:], B[idx_r,idx_t+1,idx_p,:], B[idx_r,idx_t+1,idx_p+1,:], B[idx_r+1,idx_t,idx_p,:], B[idx_r+1,idx_t,idx_p+1,:], B[idx_r+1,idx_t+1,idx_p,:], B[idx_r+1,idx_t+1,idx_p+1,:])
	B_2=B_2/np.sqrt(B_2[0]*B_2[0]+B_2[1]*B_2[1]+B_2[2]*B_2[2])
	k_3_x=dt*(B_2[0]*sin_th*cos_ph+B_2[1]*cos_th*cos_ph-B_2[2]*sin_ph)
	k_3_y=dt*(B_2[0]*sin_th*sin_ph+B_2[1]*cos_th*sin_ph+B_2[2]*cos_ph)
	k_3_z=dt*(B_2[0]*cos_th-B_2[1]*sin_th)

	x_3=x_0+k_3_x
	y_3=y_0+k_3_y
	z_3=z_0+k_3_z
	r_3=np.sqrt(x_3**2+y_3**2+z_3**2)
	t_3=np.arccos(z_3/r_3)
	p_3=np.arctan2(y_3,x_3)
	sin_th=np.sin(t_3)
	cos_th=np.cos(t_3)
	sin_ph=np.sin(p_3)
	cos_ph=np.cos(p_3)
	idx_r=ref_idx_search(R,r_3,idx_r)
	idx_t=ref_idx_search(theta,t_3,idx_t)
	idx_p=ref_idx_search(phi,p_3,idx_p)
	B_3=Trilinear_interpolation(r_3,t_3,p_3,R[idx_r],R[idx_r+1],theta[idx_t],theta[idx_t+1],phi[idx_p],phi[idx_p+1],B[idx_r,idx_t,idx_p,:], B[idx_r,idx_t,idx_p+1,:], B[idx_r,idx_t+1,idx_p,:], B[idx_r,idx_t+1,idx_p+1,:], B[idx_r+1,idx_t,idx_p,:], B[idx_r+1,idx_t,idx_p+1,:], B[idx_r+1,idx_t+1,idx_p,:], B[idx_r+1,idx_t+1,idx_p+1,:])
	B_3=B_3/np.sqrt(B_3[0]*B_3[0]+B_3[1]*B_3[1]+B_3[2]*B_3[2])
	k_4_x=dt*(B_3[0]*sin_th*cos_ph+B_3[1]*cos_th*cos_ph-B_3[2]*sin_ph)
	k_4_y=dt*(B_3[0]*sin_th*sin_ph+B_3[1]*cos_th*sin_ph+B_3[2]*cos_ph)
	k_4_z=dt*(B_3[0]*cos_th-B_3[1]*sin_th)

	x_4=x_0+(k_1_x+2.0*k_2_x+2.0*k_3_x+k_4_x)/6.0
	y_4=y_0+(k_1_y+2.0*k_2_y+2.0*k_3_y+k_4_y)/6.0
	z_4=z_0+(k_1_z+2.0*k_2_z+2.0*k_3_z+k_4_z)/6.0
	r_4=np.sqrt(x_4**2+y_4**2+z_4**2)
	t_4=np.arccos(z_4/r_4)
	p_4=np.arctan2(y_4,x_4)

	return np.array([r_4,t_4,p_4]),idx_r,idx_t,idx_p


def find_flicks_idx(new_logR,new_theta,new_phi,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks):
	"""
	Find index of block in ARMS grid, use previous index as start
	"""
	def coord_bounded(x,x1,x2):
		if x<=x1 and x>=x2:
			return True
		elif x>=x1 and x<=x2:
			return True
		else:
			return False

	block_found=False
	if coord_bounded(new_logR,coord_logR[idx_flicks,0],coord_logR[idx_flicks,1]) and coord_bounded(new_theta,coord_theta[idx_flicks,0],coord_theta[idx_flicks,1]) and coord_bounded(new_phi,coord_phi[idx_flicks,0],coord_phi[idx_flicks,1]):
		block_found=True
	else:
		for idx in range(nlblks):
			if coord_bounded(new_logR,coord_logR[idx,0],coord_logR[idx,1]) and coord_bounded(new_theta,coord_theta[idx,0],coord_theta[idx,1]) and coord_bounded(new_phi,coord_phi[idx,0],coord_phi[idx,1]):
				idx_flicks=idx
				block_found=True
				break
	return idx_flicks,block_found


def rk4_field_trace_flicks_SphE(pos_in,coord_logR,coord_theta,coord_phi,B_flicks,dt,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks):
	"""
	RK4 field line integration, Spherical Exponential flicks grid
	Currently sub-optimal due to switching back and forth to Cartesian coordinates
	TODO - fully spherical
	"""
	r_0=pos_in[0]
	log_r_0=np.log(r_0)
	t_0=pos_in[1]
	p_0=pos_in[2]
	sin_th=np.sin(t_0)
	cos_th=np.cos(t_0)
	sin_ph=np.sin(p_0)
	cos_ph=np.cos(p_0)
	x_0=r_0*sin_th*cos_ph
	y_0=r_0*sin_th*sin_ph
	z_0=r_0*cos_th
	if t_0>np.pi:
		t_0-=np.pi
	elif t_0<0.0:
		t_0+=np.pi
	if p_0>np.pi:
		p_0-=2.0*np.pi
	elif p_0<-np.pi:
		p_0+=2.0*np.pi

	idx_flicks,block_found=find_flicks_idx(log_r_0,t_0,p_0,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_0-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_0-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_0-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]

		B_0=Trilinear_interpolation(log_r_0,t_0,p_0,logRl,logRr,thetal,thetar,phil,phir,B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		B_0=B_0/np.sqrt(B_0[0]*B_0[0]+B_0[1]*B_0[1]+B_0[2]*B_0[2])
		k_1_x=dt*(B_0[0]*sin_th*cos_ph+B_0[1]*cos_th*cos_ph-B_0[2]*sin_ph)
		k_1_y=dt*(B_0[0]*sin_th*sin_ph+B_0[1]*cos_th*sin_ph+B_0[2]*cos_ph)
		k_1_z=dt*(B_0[0]*cos_th-B_0[1]*sin_th)

		x_1=x_0+k_1_x*0.5
		y_1=y_0+k_1_y*0.5
		z_1=z_0+k_1_z*0.5
		r_1=np.sqrt(x_1**2+y_1**2+z_1**2)
		log_r_1=np.log(r_1)
		t_1=np.arccos(z_1/r_1)
		p_1=np.arctan2(y_1,x_1)
		sin_th=np.sin(t_1)
		cos_th=np.cos(t_1)
		sin_ph=np.sin(p_1)
		cos_ph=np.cos(p_1)

	if block_found:
		idx_flicks,block_found=find_flicks_idx(log_r_1,t_1,p_1,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_1-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_1-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_1-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		B_1=Trilinear_interpolation(log_r_1,t_1,p_1,logRl,logRr,thetal,thetar,phil,phir,B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		B_1=B_1/np.sqrt(B_1[0]*B_1[0]+B_1[1]*B_1[1]+B_1[2]*B_1[2])
		k_2_x=dt*(B_1[0]*sin_th*cos_ph+B_1[1]*cos_th*cos_ph-B_1[2]*sin_ph)
		k_2_y=dt*(B_1[0]*sin_th*sin_ph+B_1[1]*cos_th*sin_ph+B_1[2]*cos_ph)
		k_2_z=dt*(B_1[0]*cos_th-B_1[1]*sin_th)
		
		x_2=x_0+k_2_x*0.5
		y_2=y_0+k_2_y*0.5
		z_2=z_0+k_2_z*0.5
		r_2=np.sqrt(x_2**2+y_2**2+z_2**2)
		log_r_2=np.log(r_2)
		t_2=np.arccos(z_2/r_2)
		p_2=np.arctan2(y_2,x_2)
		sin_th=np.sin(t_2)
		cos_th=np.cos(t_2)
		sin_ph=np.sin(p_2)
		cos_ph=np.cos(p_2)
	if block_found:
		idx_flicks,block_found=find_flicks_idx(log_r_2,t_2,p_2,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_2-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_2-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_2-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		B_2=Trilinear_interpolation(log_r_2,t_2,p_2,logRl,logRr,thetal,thetar,phil,phir,B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		B_2=B_2/np.sqrt(B_2[0]*B_2[0]+B_2[1]*B_2[1]+B_2[2]*B_2[2])
		k_3_x=dt*(B_2[0]*sin_th*cos_ph+B_2[1]*cos_th*cos_ph-B_2[2]*sin_ph)
		k_3_y=dt*(B_2[0]*sin_th*sin_ph+B_2[1]*cos_th*sin_ph+B_2[2]*cos_ph)
		k_3_z=dt*(B_2[0]*cos_th-B_2[1]*sin_th)

		x_3=x_0+k_3_x
		y_3=y_0+k_3_y
		z_3=z_0+k_3_z
		r_3=np.sqrt(x_3**2+y_3**2+z_3**2)
		log_r_3=np.log(r_3)
		t_3=np.arccos(z_3/r_3)
		p_3=np.arctan2(y_3,x_3)
		sin_th=np.sin(t_3)
		cos_th=np.cos(t_3)
		sin_ph=np.sin(p_3)
		cos_ph=np.cos(p_3)
	if block_found:
		idx_flicks,block_found=find_flicks_idx(log_r_3,t_3,p_3,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_3-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_3-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_3-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		B_3=Trilinear_interpolation(log_r_3,t_3,p_3,logRl,logRr,thetal,thetar,phil,phir,B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],B_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],B_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], B_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		B_3=B_3/np.sqrt(B_3[0]*B_3[0]+B_3[1]*B_3[1]+B_3[2]*B_3[2])
		k_4_x=dt*(B_3[0]*sin_th*cos_ph+B_3[1]*cos_th*cos_ph-B_3[2]*sin_ph)
		k_4_y=dt*(B_3[0]*sin_th*sin_ph+B_3[1]*cos_th*sin_ph+B_3[2]*cos_ph)
		k_4_z=dt*(B_3[0]*cos_th-B_3[1]*sin_th)

		x_4=x_0+(k_1_x+2.0*k_2_x+2.0*k_3_x+k_4_x)/6.0
		y_4=y_0+(k_1_y+2.0*k_2_y+2.0*k_3_y+k_4_y)/6.0
		z_4=z_0+(k_1_z+2.0*k_2_z+2.0*k_3_z+k_4_z)/6.0
		r_4=np.sqrt(x_4**2+y_4**2+z_4**2)
		t_4=np.arccos(z_4/r_4)
		p_4=np.arctan2(y_4,x_4)
	else:
		r_4=0.0
		t_4=0.0
		p_4=0.0

	return np.array([r_4,t_4,p_4]),idx_flicks,block_found


def rk4_field_trace_flicks_Cart(pos_in,coord_x,coord_y,coord_z,B_flicks,dt,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks):
	"""
	RK4 field line integration, Cartesian flicks grid
	"""
	x_0=pos_in[0]
	y_0=pos_in[1]
	z_0=pos_in[2]

	idx_flicks,block_found=find_flicks_idx(x_0,y_0,z_0,coord_x,coord_y,coord_z,nlblks,idx_flicks)
	if block_found:
		idx_x=int(np.floor((x_0-coord_x[idx_flicks,0])*n1pm1/(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])))
		if idx_x==n1pm1:
			idx_x=n1pm1-1
		idx_y=int(np.floor((y_0-coord_y[idx_flicks,0])*n2pm1/(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])))
		if idx_y==n2pm1:
			idx_y=n2pm1-1
		idx_z=int(np.floor((z_0-coord_z[idx_flicks,0])*n3pm1/(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])))
		if idx_z==n3pm1:
			idx_z=n3pm1-1

		xl=idx_x/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		xr=(idx_x+1)/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		yl=idx_y/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		yr=(idx_y+1)/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		zl=idx_z/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		zr=(idx_z+1)/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]

		B_0=Trilinear_interpolation(x_0,y_0,z_0,xl,xr,yl,yr,zl,zr,B_flicks[idx_flicks,idx_z,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z+1,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y,idx_x+1,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x+1,:])
		B_0=B_0/np.sqrt(B_0[0]*B_0[0]+B_0[1]*B_0[1]+B_0[2]*B_0[2])
		k_1_x=dt*B_0[0]
		k_1_y=dt*B_0[1]
		k_1_z=dt*B_0[2]

		x_1=x_0+k_1_x*0.5
		y_1=y_0+k_1_y*0.5
		z_1=z_0+k_1_z*0.5
	if block_found:
		idx_flicks,block_found=find_flicks_idx(x_1,y_1,z_1,coord_x,coord_y,coord_z,nlblks,idx_flicks)
	if block_found:
		idx_x=int(np.floor((x_1-coord_x[idx_flicks,0])*n1pm1/(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])))
		if idx_x==n1pm1:
			idx_x=n1pm1-1
		idx_y=int(np.floor((y_1-coord_y[idx_flicks,0])*n2pm1/(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])))
		if idx_y==n2pm1:
			idx_y=n2pm1-1
		idx_z=int(np.floor((z_1-coord_z[idx_flicks,0])*n3pm1/(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])))
		if idx_z==n3pm1:
			idx_z=n3pm1-1

		xl=idx_x/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		xr=(idx_x+1)/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		yl=idx_y/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		yr=(idx_y+1)/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		zl=idx_z/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		zr=(idx_z+1)/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		B_1=Trilinear_interpolation(x_1,y_1,z_1,xl,xr,yl,yr,zl,zr,B_flicks[idx_flicks,idx_z,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z+1,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y,idx_x+1,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x+1,:])
		B_1=B_1/np.sqrt(B_1[0]*B_1[0]+B_1[1]*B_1[1]+B_1[2]*B_1[2])
		k_2_x=dt*B_1[0]
		k_2_y=dt*B_1[1]
		k_2_z=dt*B_1[2]
		
		x_2=x_0+k_2_x*0.5
		y_2=y_0+k_2_y*0.5
		z_2=z_0+k_2_z*0.5
	if block_found:
		idx_flicks,block_found=find_flicks_idx(x_2,y_2,z_2,coord_x,coord_y,coord_z,nlblks,idx_flicks)
	if block_found:
		idx_x=int(np.floor((x_2-coord_x[idx_flicks,0])*n1pm1/(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])))
		if idx_x==n1pm1:
			idx_x=n1pm1-1
		idx_y=int(np.floor((y_2-coord_y[idx_flicks,0])*n2pm1/(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])))
		if idx_y==n2pm1:
			idx_y=n2pm1-1
		idx_z=int(np.floor((z_2-coord_z[idx_flicks,0])*n3pm1/(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])))
		if idx_z==n3pm1:
			idx_z=n3pm1-1

		xl=idx_x/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		xr=(idx_x+1)/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		yl=idx_y/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		yr=(idx_y+1)/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		zl=idx_z/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		zr=(idx_z+1)/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		B_2=Trilinear_interpolation(x_2,y_2,z_2,xl,xr,yl,yr,zl,zr,B_flicks[idx_flicks,idx_z,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z+1,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y,idx_x+1,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x+1,:])
		B_2=B_2/np.sqrt(B_2[0]*B_2[0]+B_2[1]*B_2[1]+B_2[2]*B_2[2])
		k_3_x=dt*B_2[0]
		k_3_y=dt*B_2[1]
		k_3_z=dt*B_2[2]

		x_3=x_0+k_3_x
		y_3=y_0+k_3_y
		z_3=z_0+k_3_z
	if block_found:
		idx_flicks,block_found=find_flicks_idx(x_3,y_3,z_3,coord_x,coord_y,coord_z,nlblks,idx_flicks)
	if block_found:
		idx_x=int(np.floor((x_3-coord_x[idx_flicks,0])*n1pm1/(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])))
		if idx_x==n1pm1:
			idx_x=n1pm1-1
		idx_y=int(np.floor((y_3-coord_y[idx_flicks,0])*n2pm1/(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])))
		if idx_y==n2pm1:
			idx_y=n2pm1-1
		idx_z=int(np.floor((z_3-coord_z[idx_flicks,0])*n3pm1/(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])))
		if idx_z==n3pm1:
			idx_z=n3pm1-1

		xl=idx_x/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		xr=(idx_x+1)/n1pm1*(coord_x[idx_flicks,1]-coord_x[idx_flicks,0])+coord_x[idx_flicks,0]
		yl=idx_y/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		yr=(idx_y+1)/n2pm1*(coord_y[idx_flicks,1]-coord_y[idx_flicks,0])+coord_y[idx_flicks,0]
		zl=idx_z/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		zr=(idx_z+1)/n3pm1*(coord_z[idx_flicks,1]-coord_z[idx_flicks,0])+coord_z[idx_flicks,0]
		B_3=Trilinear_interpolation(x_3,y_3,z_3,xl,xr,yl,yr,zl,zr,B_flicks[idx_flicks,idx_z,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z+1,idx_y,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x,:],B_flicks[idx_flicks,idx_z,idx_y,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y,idx_x+1,:],B_flicks[idx_flicks,idx_z,idx_y+1,idx_x+1,:], B_flicks[idx_flicks,idx_z+1,idx_y+1,idx_x+1,:])
		B_3=B_3/np.sqrt(B_3[0]*B_3[0]+B_3[1]*B_3[1]+B_3[2]*B_3[2])
		k_4_x=dt*B_3[0]
		k_4_y=dt*B_3[1]
		k_4_z=dt*B_3[2]

		x_4=x_0+(k_1_x+2.0*k_2_x+2.0*k_3_x+k_4_x)/6.0
		y_4=y_0+(k_1_y+2.0*k_2_y+2.0*k_3_y+k_4_y)/6.0
		z_4=z_0+(k_1_z+2.0*k_2_z+2.0*k_3_z+k_4_z)/6.0
	else:
		x_4=0.0
		y_4=0.0
		z_4=0.0

	return np.array([x_4,y_4,z_4]),idx_flicks,block_found


def rk4_velocity_trace_flicks(pos_in,coord_logR,coord_theta,coord_phi,v_flicks,dt,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks,solar_Radius):
	"""
	Advect passive particle in flicks velocity grid (streamline)
	Currently sub-optimal due to switching back and forth to Cartesian coordinates
	TODO - fully spherical
	"""
	r_0=pos_in[0]
	log_r_0=np.log(r_0)
	t_0=pos_in[1]
	p_0=pos_in[2]
	sin_th=np.sin(t_0)
	cos_th=np.cos(t_0)
	sin_ph=np.sin(p_0)
	cos_ph=np.cos(p_0)
	x_0=r_0*sin_th*cos_ph*solar_Radius
	y_0=r_0*sin_th*sin_ph*solar_Radius
	z_0=r_0*cos_th*solar_Radius
	if t_0>np.pi:
		t_0-=np.pi
	elif t_0<0.0:
		t_0+=np.pi
	if p_0>np.pi:
		p_0-=2.0*np.pi
	elif p_0<-np.pi:
		p_0+=2.0*np.pi

	idx_flicks,block_found=find_flicks_idx(log_r_0,t_0,p_0,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_0-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_0-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_0-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]

		v_0=Trilinear_interpolation(log_r_0,t_0,p_0,logRl,logRr,thetal,thetar,phil,phir,v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		k_1_x=dt*(v_0[0]*sin_th*cos_ph+v_0[1]*cos_th*cos_ph-v_0[2]*sin_ph)
		k_1_y=dt*(v_0[0]*sin_th*sin_ph+v_0[1]*cos_th*sin_ph+v_0[2]*cos_ph)
		k_1_z=dt*(v_0[0]*cos_th-v_0[1]*sin_th)

		x_1=x_0+k_1_x*0.5
		y_1=y_0+k_1_y*0.5
		z_1=z_0+k_1_z*0.5
		r_1=np.sqrt(x_1**2+y_1**2+z_1**2)
		log_r_1=np.log(r_1/solar_Radius)
		t_1=np.arccos(z_1/r_1)
		p_1=np.arctan2(y_1,x_1)
		sin_th=np.sin(t_1)
		cos_th=np.cos(t_1)
		sin_ph=np.sin(p_1)
		cos_ph=np.cos(p_1)

	if block_found:
		idx_flicks,block_found=find_flicks_idx(log_r_1,t_1,p_1,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_1-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_1-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_1-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		v_1=Trilinear_interpolation(log_r_1,t_1,p_1,logRl,logRr,thetal,thetar,phil,phir,v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		k_2_x=dt*(v_1[0]*sin_th*cos_ph+v_1[1]*cos_th*cos_ph-v_1[2]*sin_ph)
		k_2_y=dt*(v_1[0]*sin_th*sin_ph+v_1[1]*cos_th*sin_ph+v_1[2]*cos_ph)
		k_2_z=dt*(v_1[0]*cos_th-v_1[1]*sin_th)
		
		x_2=x_0+k_2_x*0.5
		y_2=y_0+k_2_y*0.5
		z_2=z_0+k_2_z*0.5
		r_2=np.sqrt(x_2**2+y_2**2+z_2**2)
		log_r_2=np.log(r_2/solar_Radius)
		t_2=np.arccos(z_2/r_2)
		p_2=np.arctan2(y_2,x_2)
		sin_th=np.sin(t_2)
		cos_th=np.cos(t_2)
		sin_ph=np.sin(p_2)
		cos_ph=np.cos(p_2)
	if block_found:
		idx_flicks,block_found=find_flicks_idx(log_r_2,t_2,p_2,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_2-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_2-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_2-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		v_2=Trilinear_interpolation(log_r_2,t_2,p_2,logRl,logRr,thetal,thetar,phil,phir,v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		k_3_x=dt*(v_2[0]*sin_th*cos_ph+v_2[1]*cos_th*cos_ph-v_2[2]*sin_ph)
		k_3_y=dt*(v_2[0]*sin_th*sin_ph+v_2[1]*cos_th*sin_ph+v_2[2]*cos_ph)
		k_3_z=dt*(v_2[0]*cos_th-v_2[1]*sin_th)

		x_3=x_0+k_3_x
		y_3=y_0+k_3_y
		z_3=z_0+k_3_z
		r_3=np.sqrt(x_3**2+y_3**2+z_3**2)
		log_r_3=np.log(r_3/solar_Radius)
		t_3=np.arccos(z_3/r_3)
		p_3=np.arctan2(y_3,x_3)
		sin_th=np.sin(t_3)
		cos_th=np.cos(t_3)
		sin_ph=np.sin(p_3)
		cos_ph=np.cos(p_3)
	if block_found:
		idx_flicks,block_found=find_flicks_idx(log_r_3,t_3,p_3,coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
	if block_found:
		idx_logR=int(np.floor((log_r_3-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((t_3-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((p_3-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		v_3=Trilinear_interpolation(log_r_3,t_3,p_3,logRl,logRr,thetal,thetar,phil,phir,v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],v_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],v_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], v_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		k_4_x=dt*(v_3[0]*sin_th*cos_ph+v_3[1]*cos_th*cos_ph-v_3[2]*sin_ph)
		k_4_y=dt*(v_3[0]*sin_th*sin_ph+v_3[1]*cos_th*sin_ph+v_3[2]*cos_ph)
		k_4_z=dt*(v_3[0]*cos_th-v_3[1]*sin_th)

		x_4=x_0+(k_1_x+2.0*k_2_x+2.0*k_3_x+k_4_x)/6.0
		y_4=y_0+(k_1_y+2.0*k_2_y+2.0*k_3_y+k_4_y)/6.0
		z_4=z_0+(k_1_z+2.0*k_2_z+2.0*k_3_z+k_4_z)/6.0

		r_4=np.sqrt(x_4**2+y_4**2+z_4**2)
		t_4=np.arccos(z_4/r_4)
		p_4=np.arctan2(y_4,x_4)
	else:
		r_4=0.0
		t_4=0.0
		p_4=0.0

	return np.array([r_4/solar_Radius,t_4,p_4]),idx_flicks,block_found


def field_line_cartesian(position0,X,Y,Z,B,R_min,R_max,max_steps=1E6,step_size=1E-2):
	"""
	Obtain back-and-forth field line, regular Cartesian coordinates
	"""
	max_steps=int(max_steps)
	R2_min=R_min*R_min
	R2_max=R_max*R_max
	forward_temp=np.zeros((max_steps,3))
	forward_temp[0,:]=position0
	forward_num=max_steps
	for idx in range(1,max_steps):
		forward_temp[idx,:]=rk4_field_trace(forward_temp[idx-1,:],X,Y,Z,B,step_size)
		current_R2=forward_temp[idx,0]*forward_temp[idx,0]+forward_temp[idx,1]*forward_temp[idx,1]+forward_temp[idx,2]*forward_temp[idx,2]
		if current_R2>R2_max or current_R2<R2_min:
			forward_num=idx+1
			break
	forward_curve=np.zeros((forward_num))
	forward_curve=forward_temp[:forward_num]

	backward_temp=np.zeros((max_steps,3))
	backward_temp[0,:]=position0
	backward_num=max_steps
	for idx in range(1,max_steps):
		backward_temp[idx,:]=rk4_field_trace(backward_temp[idx-1,:],X,Y,Z,-B,step_size)
		current_R2=backward_temp[idx,0]*backward_temp[idx,0]+backward_temp[idx,1]*backward_temp[idx,1]+backward_temp[idx,2]*backward_temp[idx,2]
		if current_R2>R2_max or current_R2<R2_min:
			backward_num=idx+1
			break
	backward_curve=np.zeros((backward_num))
	backward_curve=backward_temp[:backward_num]

	return forward_curve,backward_curve


def field_line_spherical(position0,R,theta,phi,B,R_min,R_max,max_steps=1E6,step_size=1E-2):
	"""
	Obtain back-and-forth field line, regular spherical coordinates, 'polar' coordinate system
	"""
	max_steps=int(max_steps)
	backward_temp=np.zeros((max_steps,3))
	backward_temp[0,:]=position0
	backward_num=max_steps
	idx_r=0
	idx_t=0
	idx_p=0
	for idx in range(1,max_steps):
		backward_temp[idx,:],idx_r,idx_t,idx_p=rk4_field_trace_spherical(backward_temp[idx-1,:],R,theta,phi,-B,step_size,idx_r,idx_t,idx_p)
		if backward_temp[idx,0]>R_max or backward_temp[idx,0]<R_min:
			backward_num=idx+1
			break
	forward_temp=np.zeros((max_steps,3))
	forward_temp[0,:]=position0
	forward_num=max_steps
	for idx in range(1,max_steps):
		forward_temp[idx,:],idx_r,idx_t,idx_p=rk4_field_trace_spherical(forward_temp[idx-1,:],R,theta,phi,B,step_size,idx_r,idx_t,idx_p)
		if forward_temp[idx,0]>R_max or forward_temp[idx,0]<R_min:
			forward_num=idx+1
			break
	final_curve=np.zeros((forward_num+backward_num-3,3))
	final_curve[:backward_num-2,:]=backward_temp[backward_num-2:0:-1,:]
	final_curve[backward_num-2:,:]=forward_temp[0:forward_num-1,:]

	return final_curve


#For Spherical Exponential, interpret the arguments as (position0,coord_logR,coord_theta,coord_phi,B_flicks,R_min,R_max...)
def field_line_flicks(position0,coord_x,coord_y,coord_z,B_flicks,z_min,z_max,nlblks,n1pm1,n2pm1,n3pm1,coord_type="SphE",max_steps=1E6,step_size=1E-2):
	"""
	Obtain back-and-forth field line, flicks grid
	"""
	max_steps=int(max_steps)
	backward_temp=np.zeros((max_steps,3))
	backward_temp[0,:]=position0
	backward_num=0
	if coord_type.lower()=="sphe":
		idx_flicks,block_found=find_flicks_idx(np.log(position0[0]),position0[1],position0[2],coord_x,coord_y,coord_z,nlblks,0)
		idx=1
		while block_found and idx<max_steps:
			backward_temp[idx,:],idx_flicks,block_found=rk4_field_trace_flicks_SphE(backward_temp[idx-1,:],coord_x,coord_y,coord_z,-B_flicks,step_size,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks)
			if backward_temp[idx,0]>z_max or backward_temp[idx,0]<z_min:
				backward_num=idx+1
				idx=max_steps+1
			idx+=1
	elif coord_type.lower()[:4]=="cart":
		idx_flicks,block_found=find_flicks_idx(position0[0],position0[1],position0[2],coord_x,coord_y,coord_z,nlblks,0)
		idx=1
		while block_found and idx<max_steps:
			backward_temp[idx,:],idx_flicks,block_found=rk4_field_trace_flicks_Cart(backward_temp[idx-1,:],coord_x,coord_y,coord_z,-B_flicks,step_size,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks)
			if backward_temp[idx,2]>z_max or backward_temp[idx,2]<z_min:
				break
			idx+=1
		backward_num=idx
	forward_temp=np.zeros((max_steps,3))
	forward_temp[0,:]=position0
	forward_num=0
	if coord_type.lower()=="sphe":
		idx_flicks,block_found=find_flicks_idx(np.log(position0[0]),position0[1],position0[2],coord_x,coord_y,coord_z,nlblks,0)
		idx=1
		while block_found and idx<max_steps:
			forward_temp[idx,:],idx_flicks,block_found=rk4_field_trace_flicks_SphE(forward_temp[idx-1,:],coord_x,coord_y,coord_z,B_flicks,step_size,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks)
			if forward_temp[idx,0]>z_max or forward_temp[idx,0]<z_min:
				forward_num=idx+1
				idx=max_steps+1
			idx+=1
	elif coord_type.lower()[:4]=="cart":
		idx_flicks,block_found=find_flicks_idx(position0[0],position0[1],position0[2],coord_x,coord_y,coord_z,nlblks,0)
		idx=1
		while block_found and idx<max_steps:
			forward_temp[idx,:],idx_flicks,block_found=rk4_field_trace_flicks_Cart(forward_temp[idx-1,:],coord_x,coord_y,coord_z,B_flicks,step_size,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks)
			if forward_temp[idx,2]>z_max or forward_temp[idx,2]<z_min:
				break
			idx+=1
		forward_num=idx
	if forward_num+backward_num>3:
		final_curve=np.zeros((forward_num+backward_num-3,3))
		final_curve[:backward_num-2,:]=backward_temp[backward_num-2:0:-1,:]
		final_curve[backward_num-2:,:]=forward_temp[0:forward_num-1,:]
	else:
		final_curve=np.array([position0])

	return final_curve


#Assume data_flicks is a scalar field with the same format as B_flicks
def data_along_field_line_flicks(position0,coord_logR,coord_theta,coord_phi,B_flicks,data_flicks,R_min,R_max,nlblks,n1pm1,n2pm1,n3pm1,max_steps=1E6,step_size=1E-2):
	"""
	Obtain back-and-forth field line, flicks grid
	Also returns interpolated flicks data
	"""
	max_steps=int(max_steps)
	shape_data=len(data_flicks[0,0,0,0,:])
	backward_temp=np.zeros((max_steps,3+shape_data))
	backward_temp[0,:3]=position0
	backward_num=0
	idx_flicks,block_found=find_flicks_idx(np.log(position0[0]),position0[1],position0[2],coord_logR,coord_theta,coord_phi,nlblks,0)

	if block_found:
		idx_logR=int(np.floor((np.log(position0[0])-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
		if idx_logR==n1pm1:
			idx_logR=n1pm1-1
		idx_theta=int(np.floor((position0[1]-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
		if idx_theta==n2pm1:
			idx_theta=n2pm1-1
		idx_phi=int(np.floor((position0[2]-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
		if idx_phi==n3pm1:
			idx_phi=n3pm1-1

		logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
		thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
		phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
		phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]

		backward_temp[0,3:]=Trilinear_interpolation(np.log(position0[0]),position0[1],position0[2],logRl,logRr,thetal,thetar,phil,phir,data_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],data_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],data_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], data_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],data_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], data_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],data_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], data_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])


	idx=1
	while block_found and idx<max_steps:
		backward_temp[idx,:3],idx_flicks,block_found=rk4_field_trace_flicks_SphE(backward_temp[idx-1,:3],coord_logR,coord_theta,coord_phi,-B_flicks,step_size,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks)
		if backward_temp[idx,0]>0.0:
			idx_flicks,block_found=find_flicks_idx(np.log(backward_temp[idx,0]),backward_temp[idx,1],backward_temp[idx,2],coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
			if block_found:
				idx_logR=int(np.floor((np.log(backward_temp[idx,0])-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
				if idx_logR==n1pm1:
					idx_logR=n1pm1-1
				idx_theta=int(np.floor((backward_temp[idx,1]-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
				if idx_theta==n2pm1:
					idx_theta=n2pm1-1
				idx_phi=int(np.floor((backward_temp[idx,2]-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
				if idx_phi==n3pm1:
					idx_phi=n3pm1-1

				logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
				logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
				thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
				thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
				phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
				phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]

				backward_temp[idx,3:]=Trilinear_interpolation(np.log(backward_temp[idx,0]),backward_temp[idx,1],backward_temp[idx,2],logRl,logRr,thetal,thetar,phil,phir,data_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],data_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],data_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], data_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],data_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], data_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],data_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], data_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		if backward_temp[idx,0]>R_max or backward_temp[idx,0]<R_min:
			backward_num=idx+1
			idx=max_steps+1
		idx+=1
	forward_temp=np.zeros((max_steps,3+shape_data))
	forward_temp[0,:3]=position0
	forward_temp[0,3:]=backward_temp[0,3:]
	forward_num=0
	idx_flicks,block_found=find_flicks_idx(np.log(position0[0]),position0[1],position0[2],coord_logR,coord_theta,coord_phi,nlblks,0)
	idx=1
	while block_found and idx<max_steps:
		forward_temp[idx,:3],idx_flicks,block_found=rk4_field_trace_flicks_SphE(forward_temp[idx-1,:3],coord_logR,coord_theta,coord_phi,B_flicks,step_size,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks)
		if forward_temp[idx,0]>0.0:
			idx_flicks,block_found=find_flicks_idx(np.log(forward_temp[idx,0]),forward_temp[idx,1],forward_temp[idx,2],coord_logR,coord_theta,coord_phi,nlblks,idx_flicks)
			if block_found:
				idx_logR=int(np.floor((np.log(forward_temp[idx,0])-coord_logR[idx_flicks,0])*n1pm1/(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])))
				if idx_logR==n1pm1:
					idx_logR=n1pm1-1
				idx_theta=int(np.floor((forward_temp[idx,1]-coord_theta[idx_flicks,0])*n2pm1/(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])))
				if idx_theta==n2pm1:
					idx_theta=n2pm1-1
				idx_phi=int(np.floor((forward_temp[idx,2]-coord_phi[idx_flicks,0])*n3pm1/(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])))
				if idx_phi==n3pm1:
					idx_phi=n3pm1-1

				logRl=idx_logR/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
				logRr=(idx_logR+1)/n1pm1*(coord_logR[idx_flicks,1]-coord_logR[idx_flicks,0])+coord_logR[idx_flicks,0]
				thetal=idx_theta/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
				thetar=(idx_theta+1)/n2pm1*(coord_theta[idx_flicks,1]-coord_theta[idx_flicks,0])+coord_theta[idx_flicks,0]
				phil=idx_phi/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]
				phir=(idx_phi+1)/n3pm1*(coord_phi[idx_flicks,1]-coord_phi[idx_flicks,0])+coord_phi[idx_flicks,0]

				forward_temp[idx,3:]=Trilinear_interpolation(np.log(forward_temp[idx,0]),forward_temp[idx,1],forward_temp[idx,2],logRl,logRr,thetal,thetar,phil,phir,data_flicks[idx_flicks,idx_phi,idx_theta,idx_logR,:],data_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR,:],data_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR,:], data_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR,:],data_flicks[idx_flicks,idx_phi,idx_theta,idx_logR+1,:], data_flicks[idx_flicks,idx_phi+1,idx_theta,idx_logR+1,:],data_flicks[idx_flicks,idx_phi,idx_theta+1,idx_logR+1,:], data_flicks[idx_flicks,idx_phi+1,idx_theta+1,idx_logR+1,:])
		if forward_temp[idx,0]>R_max or forward_temp[idx,0]<R_min:
			forward_num=idx+1
			idx=max_steps+1
		idx+=1
	if forward_num+backward_num>3:
		final_curve=np.zeros((forward_num+backward_num-3,len(backward_temp[0,:])))
		final_curve[:backward_num-2,:]=backward_temp[backward_num-2:0:-1,:]
		final_curve[backward_num-2:,:]=forward_temp[0:forward_num-1,:]
	else:
		final_curve=np.array([np.concatenate((position0,np.zeros((shape_data))))])

	return final_curve


def particle_advance_flicks(num_points,R_points,theta_points,phi_points,coord_logR,coord_theta,coord_phi,v_flicks,dt,nlblks,n1pm1,n2pm1,n3pm1,solar_Radius):
	"""
	Obtain passive particle trajectory from flicks file
	"""
	R_new=np.zeros((num_points))
	theta_new=np.zeros((num_points))
	phi_new=np.zeros((num_points))
	for idx in range(num_points):
		idx_flicks,block_found=find_flicks_idx(R_points[idx],theta_points[idx],phi_points[idx],coord_logR,coord_theta,coord_phi,nlblks,0)
		temp,idx_flicks,block_found=rk4_velocity_trace_flicks(np.array([R_points[idx],theta_points[idx],phi_points[idx]]),coord_logR,coord_theta,coord_phi,v_flicks,dt,nlblks,n1pm1,n2pm1,n3pm1,idx_flicks,solar_Radius)
		if block_found:
			R_new[idx]=temp[0]
			theta_new[idx]=temp[1]
			phi_new[idx]=temp[2]
		else:
			R_new[idx]=R_points[idx]
			theta_new[idx]=theta_points[idx]
			phi_new[idx]=phi_points[idx]
	return R_new,theta_new,phi_new


def parse_trajectories_file(infile_path):
	"""
	Saved passive particle trajectories
	"""
	filesize=os.path.getsize(infile_path)
	infile=open(infile_path,"rb")
	num_particles=struct.unpack("i",infile.read(4))[0]
	num_times=int(np.floor((filesize-4)/(12*num_particles+4)))
	times=np.zeros((num_times),dtype='float32')
	trajectories=np.zeros((3,num_particles,num_times),dtype='float32')
	for idx_t in range(num_times):
		times[idx_t]=struct.unpack("f",infile.read(4))[0]
		for idx_p in range(num_particles):
			temp=struct.unpack("fff",infile.read(12))
			trajectories[0,idx_p,idx_t]=temp[0]
			trajectories[1,idx_p,idx_t]=temp[1]
			trajectories[2,idx_p,idx_t]=temp[2]
	return times,trajectories


def parse_QSL_folder(folder_path):
	"""
	Load in regular grid and B field from folder from QSLsquasher run
	"""
	dims_file=open(os.path.join(folder_path,"dim.txt"),"r")
	X_dim=int(dims_file.readline())
	Y_dim=int(dims_file.readline())
	Z_dim=int(dims_file.readline())
	dims_file.close()
	X_file=open(os.path.join(folder_path,"xs0.dat"),"rb")
	data=X_file.read(8*X_dim)
	X=np.array(struct.unpack("d"*X_dim,data))
	X_file.close()
	Y_file=open(os.path.join(folder_path,"ys0.dat"),"rb")
	data=Y_file.read(8*Y_dim)
	Y=np.array(struct.unpack("d"*Y_dim,data))
	Y_file.close()
	Z_file=open(os.path.join(folder_path,"zs0.dat"),"rb")
	data=Z_file.read(8*Z_dim)
	Z=np.array(struct.unpack("d"*Z_dim,data))
	Z_file.close()
	B=np.zeros((Z_dim,Y_dim,X_dim,3))
	B_X_file=open(os.path.join(folder_path,"bx0.dat"),"rb")
	B_Y_file=open(os.path.join(folder_path,"by0.dat"),"rb")
	B_Z_file=open(os.path.join(folder_path,"bz0.dat"),"rb")
	for idx_Z in range(Z_dim):
		for idx_Y in range(Y_dim):
			for idx_X in range(X_dim):
				B[idx_Z,idx_Y,idx_X,0]=struct.unpack("d",B_Z_file.read(8))[0]
				B[idx_Z,idx_Y,idx_X,1]=struct.unpack("d",B_Y_file.read(8))[0]
				B[idx_Z,idx_Y,idx_X,2]=struct.unpack("d",B_X_file.read(8))[0]
	B_X_file.close()
	B_Y_file.close()
	B_Z_file.close()
	return X,Y,Z,B


def parse_QSL_Rbinfile(filename,coords_spherical=True):
	"""
	Load in QSL output file
	Note - uses special version of QSLsquasher adapted to output binary data
	"""
	filesize=os.path.getsize(filename)
	num_points=filesize//(28)
	infile=open(filename,"rb")
	r=np.zeros((num_points))
	theta=np.zeros((num_points))
	phi=np.zeros((num_points))
	Q=np.zeros((num_points))
	for idx in range(num_points):
		data=struct.unpack('dddf',infile.read(28))
		r[idx]=data[2]
		theta[idx]=data[1]
		phi[idx]=data[0]
		Q[idx]=data[3]
	infile.close()

	theta_dim=1
	idx=1
	while idx<num_points:
		if phi[idx]==phi[0]:
			theta_dim+=1
			idx+=1
		else:
			idx=num_points
	phi_dim=num_points//theta_dim


	R_actual=r[0]
	if coords_spherical:
		theta=theta.reshape((phi_dim,theta_dim))+np.pi*0.5
		phi=phi.reshape((phi_dim,theta_dim))-np.pi
	else:
		theta=theta.reshape((phi_dim,theta_dim))
		phi=phi.reshape((phi_dim,theta_dim))
	Q=Q.reshape((phi_dim,theta_dim))
	return R_actual,theta,phi,Q


def parse_QSL_RbinfileExpansion(filename):
	"""
	Load in QSL output file with expansion factor calculation
	Note - uses special version of QSLsquasher adapted to output binary data
	"""
	filesize=os.path.getsize(filename)
	num_points=filesize//(44)
	infile=open(filename,"rb")
	r=np.zeros((num_points))
	theta=np.zeros((num_points))
	phi=np.zeros((num_points))
	Q=np.zeros((num_points))
	Bs=np.zeros((num_points))
	Be=np.zeros((num_points))
	for idx in range(num_points):
		data=struct.unpack('=dddfdd',infile.read(44))
		r[idx]=data[2]
		theta[idx]=data[1]
		phi[idx]=data[0]
		Q[idx]=data[3]
		Bs[idx]=data[4]
		Be[idx]=data[5]
	infile.close()

	theta_dim=1
	idx=1
	while idx<num_points:
		if phi[idx]==phi[0]:
			theta_dim+=1
			idx+=1
		else:
			idx=num_points
	phi_dim=num_points//theta_dim


	R_actual=r[0]
	theta=theta.reshape((phi_dim,theta_dim))+np.pi*0.5
	phi=phi.reshape((phi_dim,theta_dim))-np.pi
	Q=Q.reshape((phi_dim,theta_dim))
	Bs=Bs.reshape((phi_dim,theta_dim))
	Be=Be.reshape((phi_dim,theta_dim))
	return R_actual,theta,phi,Q,Bs,Be


def linear_interpolation_2pt(x0,y0,x1,y1,x_t):
	y_t=(y1-y0)/(x1-x0)*(x_t-x0)+y1
	return y_t


def get_Polarity_Inversion_Line(phi,theta,B):
	"""
	Get line where B=0, assuming single solution exists for every phi
	"""
	phi_PIL=phi[:]
	theta_PIL=np.zeros((len(phi_PIL)))
	size_t=len(B[:,0])
	for idx_p in range(len(phi_PIL)):
		for idx_t in range(1,size_t):
			if B[idx_t,idx_p]*B[idx_t-1,idx_p]<=0.0:
				theta_PIL[idx_p]=linear_interpolation_2pt(B[idx_t-1,idx_p],theta[idx_t-1],B[idx_t,idx_p],theta[idx_t],0.0)
	return theta_PIL,phi_PIL


def quadratic_three_points(x,y):
	a=(x[0]*(y[2]-y[1])+x[1]*(y[0]-y[2])+x[2]*(y[1]-y[0]))/(x[0]-x[1])/(x[0]-x[2])/(x[1]-x[2])
	b=(y[1]-y[0])/(x[1]-x[0])-a*(x[0]+x[1])
	c=y[0]-a*x[0]*x[0]-b*x[0]
	return a,b,c


def quadratic_minimum(a,b,c):
	return -0.5*b/a


def v_Alfven(B,rho,mu=1.0): #CGS
	return 2.18E11*B/np.sqrt(rho*mu)
mass_Hydrogen_grams=1.673E-24 #CGS


def v_magnetosonic(B,rho,T_e,angle_rads,bool_fast=True,mu=1.0,Z=1.0): #CGS
	"""
	Fast/slow magnetosonic wave
	angle_rads = angle between B and k in radians
	"""
	gamma=5.0/3.0
	c_s=9.79E5*np.sqrt(gamma*Z*T_e/mu)
	v_A=v_Alfven(B,rho,mu=mu)
	radical=np.sqrt(1.0-(2.0*np.cos(angle_rads)*c_s*v_A)**2/(v_A*v_A+c_s*c_s)**2)
	if bool_fast: #Fast wave
		v_ms2=0.5*(v_A*v_A+c_s*c_s)*(1.0+radical)
	else: #slow wave
		v_ms2=0.5*(v_A*v_A+c_s*c_s)*(1.0-radical)
	return np.sqrt(v_ms2)


def parse_sflow_CNT(cnt_file_path):
	"""
	Read in surface flow parameters from control file
	Not all options implemented yet
	"""
	if os.path.isfile(cnt_file_path):
		Time_profiles=[]	#0 - uniform, 1 - linear, 2 - "cos", 3 - sin
		parameters=[]
		cntfile=open(cnt_file_path,"r")
		num_lines=sum(1 for line in cntfile)
		cntfile.close()
		lines_read=0
		cntfile=open(cnt_file_path,"r")
		line=cntfile.readline(); lines_read+=1
		while lines_read<num_lines:
			line=cntfile.readline(); lines_read+=1
			if "separable boundary flow" in line.lower():
				line=cntfile.readline(); lines_read+=1
				if "r"==line.split()[0].lower()[0]: #Radial
					while len(line.split())>0:
						if "ti" in line.split()[0].lower():
							if "u"==line.split()[1].lower()[0]:
								Time_profiles.append(0)
							elif "l"==line.split()[1].lower()[0]:
								Time_profiles.append(1)
							elif "c"==line.split()[1].lower()[0]:
								Time_profiles.append(2)
							elif "s"==line.split()[1].lower()[0]:
								Time_profiles.append(3)
						line=cntfile.readline(); lines_read+=1
				else:
					print("Unknown type: ",line)
					while len(line.split())>0:
						line=cntfile.readline(); lines_read+=1
			
			if "nsflow" in line.lower():
				line=cntfile.readline(); lines_read+=1
				nsflow=int(line.split()[0])
				for idx in range(nsflow):
					line=cntfile.readline(); lines_read+=1
					line=cntfile.readline(); lines_read+=1
					sfl1l=float(line.split()[0])
					sfl1r=float(line.split()[1])
					sfl1c=float(line.split()[2])
					k1sfl=float(line.split()[3])
					v0sfl1=float(line.split()[4])
					line=cntfile.readline(); lines_read+=1
					line=cntfile.readline(); lines_read+=1
					sfl2l=float(line.split()[0])
					sfl2r=float(line.split()[1])
					sfl2c=float(line.split()[2])
					k2sfl=float(line.split()[3])
					v0sfl2=float(line.split()[4])
					line=cntfile.readline(); lines_read+=1
					line=cntfile.readline(); lines_read+=1
					sfl3l=float(line.split()[0])
					sfl3r=float(line.split()[1])
					sfl3c=float(line.split()[2])
					k3sfl=float(line.split()[3])
					v0sfl3=float(line.split()[4])
					line=cntfile.readline(); lines_read+=1
					line=cntfile.readline(); lines_read+=1
					tlsfl=float(line.split()[0])
					trsfl=float(line.split()[1])
					tcsfl=float(line.split()[2])
					ktsfl=float(line.split()[3])
					parameters.append([sfl1l,sfl1r,sfl1c,k1sfl,v0sfl1,sfl2l,sfl2r,sfl2c,k2sfl,v0sfl2,sfl3l,sfl3r,sfl3c,k3sfl,v0sfl3,tlsfl,trsfl,tcsfl,ktsfl])
		cntfile.close()
		return np.array(Time_profiles,dtype=np.int32), np.array(parameters,dtype=np.float32)
	else:
		print("Warning, .cnt file not found!")
		return np.array([],dtype=np.int32), np.array([[]],dtype=np.float32)


def ARMS_Gauss(x,x_c,c):
	return np.exp(-c*(x-x_c)*(x-x_c))


def ARMS_DGauss(x,x_c,c):
	return -2.0*c*(x-x_c)*np.exp(-c*(x-x_c)*(x-x_c))


def surface_move_rk_4(R,theta,phi,t,delta_t,Time_profiles,sflow_parameters):
	"""
	Advect point on photosphere according to ARMS flow profiles
	ARMS coordinates
	"""
	k_1_th=0.0
	k_1_ph=0.0
	for idx_prof in range(len(sflow_parameters[:,0])):
		if Time_profiles[idx_prof]>1:
			sfl1l=sflow_parameters[idx_prof,0]
			sfl1r=sflow_parameters[idx_prof,1]
			sfl1c=sflow_parameters[idx_prof,2]
			k1sfl=sflow_parameters[idx_prof,3]
			v0sfl1=sflow_parameters[idx_prof,4]
			sfl2l=sflow_parameters[idx_prof,5]
			sfl2r=sflow_parameters[idx_prof,6]
			sfl2c=sflow_parameters[idx_prof,7]
			k2sfl=sflow_parameters[idx_prof,8]
			v0sfl2=sflow_parameters[idx_prof,9]
			sfl3l=sflow_parameters[idx_prof,10]
			sfl3r=sflow_parameters[idx_prof,11]
			sfl3c=sflow_parameters[idx_prof,12]
			k3sfl=sflow_parameters[idx_prof,13]
			v0sfl3=sflow_parameters[idx_prof,14]
			tlsfl=sflow_parameters[idx_prof,15]
			trsfl=sflow_parameters[idx_prof,16]
			tcsfl=sflow_parameters[idx_prof,17]
			ktsfl=sflow_parameters[idx_prof,18]
			if t>=tlsfl and t<=trsfl:
				A_time=0.0
				if Time_profiles[idx_prof]==2:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=0.5*(1.0-np.cos(t_p*(t-tcsfl)))
				if Time_profiles[idx_prof]==3:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=np.sin(t_p*(t-tcsfl))
				v_theta=v0sfl2*ARMS_Gauss(theta,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_DGauss(phi,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.sin(np.pi*theta)/np.pi*A_time
				v_phi=v0sfl3*ARMS_DGauss(theta,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_Gauss(phi,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.pi*A_time
				k_1_th+=v_theta/R/np.pi
				k_1_ph+=v_phi/R/np.pi/np.sin(np.pi*theta)

	k_2_th=0.0
	k_2_ph=0.0
	for idx_prof in range(len(sflow_parameters[:,0])):
		if Time_profiles[idx_prof]>1:
			sfl1l=sflow_parameters[idx_prof,0]
			sfl1r=sflow_parameters[idx_prof,1]
			sfl1c=sflow_parameters[idx_prof,2]
			k1sfl=sflow_parameters[idx_prof,3]
			v0sfl1=sflow_parameters[idx_prof,4]
			sfl2l=sflow_parameters[idx_prof,5]
			sfl2r=sflow_parameters[idx_prof,6]
			sfl2c=sflow_parameters[idx_prof,7]
			k2sfl=sflow_parameters[idx_prof,8]
			v0sfl2=sflow_parameters[idx_prof,9]
			sfl3l=sflow_parameters[idx_prof,10]
			sfl3r=sflow_parameters[idx_prof,11]
			sfl3c=sflow_parameters[idx_prof,12]
			k3sfl=sflow_parameters[idx_prof,13]
			v0sfl3=sflow_parameters[idx_prof,14]
			tlsfl=sflow_parameters[idx_prof,15]
			trsfl=sflow_parameters[idx_prof,16]
			tcsfl=sflow_parameters[idx_prof,17]
			ktsfl=sflow_parameters[idx_prof,18]
			if t+0.5*delta_t>=tlsfl and t+0.5*delta_t<=trsfl:
				A_time=0.0
				if Time_profiles[idx_prof]==2:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=0.5*(1.0-np.cos(t_p*(t+0.5*delta_t-tcsfl)))
				if Time_profiles[idx_prof]==3:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=np.sin(t_p*(t+0.5*delta_t-tcsfl))
				v_theta=v0sfl2*ARMS_Gauss(theta+0.5*delta_t*k_1_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_DGauss(phi+0.5*delta_t*k_1_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.sin(np.pi*(theta+0.5*delta_t*k_1_th))/np.pi*A_time
				v_phi=v0sfl3*ARMS_DGauss(theta+0.5*delta_t*k_1_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_Gauss(phi+0.5*delta_t*k_1_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.pi*A_time
				k_2_th+=v_theta/R/np.pi
				k_2_ph+=v_phi/R/np.pi/np.sin(np.pi*(theta+0.5*delta_t*k_1_th))

	k_3_th=0.0
	k_3_ph=0.0
	for idx_prof in range(len(sflow_parameters[:,0])):
		if Time_profiles[idx_prof]>1:
			sfl1l=sflow_parameters[idx_prof,0]
			sfl1r=sflow_parameters[idx_prof,1]
			sfl1c=sflow_parameters[idx_prof,2]
			k1sfl=sflow_parameters[idx_prof,3]
			v0sfl1=sflow_parameters[idx_prof,4]
			sfl2l=sflow_parameters[idx_prof,5]
			sfl2r=sflow_parameters[idx_prof,6]
			sfl2c=sflow_parameters[idx_prof,7]
			k2sfl=sflow_parameters[idx_prof,8]
			v0sfl2=sflow_parameters[idx_prof,9]
			sfl3l=sflow_parameters[idx_prof,10]
			sfl3r=sflow_parameters[idx_prof,11]
			sfl3c=sflow_parameters[idx_prof,12]
			k3sfl=sflow_parameters[idx_prof,13]
			v0sfl3=sflow_parameters[idx_prof,14]
			tlsfl=sflow_parameters[idx_prof,15]
			trsfl=sflow_parameters[idx_prof,16]
			tcsfl=sflow_parameters[idx_prof,17]
			ktsfl=sflow_parameters[idx_prof,18]
			if t+0.5*delta_t>=tlsfl and t+0.5*delta_t<=trsfl:
				A_time=0.0
				if Time_profiles[idx_prof]==2:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=0.5*(1.0-np.cos(t_p*(t+0.5*delta_t-tcsfl)))
				if Time_profiles[idx_prof]==3:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=np.sin(t_p*(t+0.5*delta_t-tcsfl))
				v_theta=v0sfl2*ARMS_Gauss(theta+0.5*delta_t*k_2_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_DGauss(phi+0.5*delta_t*k_2_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.sin(np.pi*(theta+0.5*delta_t*k_2_th))/np.pi*A_time
				v_phi=v0sfl3*ARMS_DGauss(theta+0.5*delta_t*k_2_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_Gauss(phi+0.5*delta_t*k_2_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.pi*A_time
				k_3_th+=v_theta/R/np.pi
				k_3_ph+=v_phi/R/np.pi/np.sin(np.pi*(theta+0.5*delta_t*k_2_th))

	k_4_th=0.0
	k_4_ph=0.0
	for idx_prof in range(len(sflow_parameters[:,0])):
		if Time_profiles[idx_prof]>1:
			sfl1l=sflow_parameters[idx_prof,0]
			sfl1r=sflow_parameters[idx_prof,1]
			sfl1c=sflow_parameters[idx_prof,2]
			k1sfl=sflow_parameters[idx_prof,3]
			v0sfl1=sflow_parameters[idx_prof,4]
			sfl2l=sflow_parameters[idx_prof,5]
			sfl2r=sflow_parameters[idx_prof,6]
			sfl2c=sflow_parameters[idx_prof,7]
			k2sfl=sflow_parameters[idx_prof,8]
			v0sfl2=sflow_parameters[idx_prof,9]
			sfl3l=sflow_parameters[idx_prof,10]
			sfl3r=sflow_parameters[idx_prof,11]
			sfl3c=sflow_parameters[idx_prof,12]
			k3sfl=sflow_parameters[idx_prof,13]
			v0sfl3=sflow_parameters[idx_prof,14]
			tlsfl=sflow_parameters[idx_prof,15]
			trsfl=sflow_parameters[idx_prof,16]
			tcsfl=sflow_parameters[idx_prof,17]
			ktsfl=sflow_parameters[idx_prof,18]
			if t+delta_t>=tlsfl and t+delta_t<=trsfl:
				A_time=0.0
				if Time_profiles[idx_prof]==2:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=0.5*(1.0-np.cos(t_p*(t+delta_t-tcsfl)))
				if Time_profiles[idx_prof]==3:
					t_p=2.0*np.pi*ktsfl/(trsfl-tlsfl)
					A_time=np.sin(t_p*(t+delta_t-tcsfl))
				v_theta=v0sfl2*ARMS_Gauss(theta+delta_t*k_3_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_DGauss(phi+delta_t*k_3_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.sin(np.pi*(theta+delta_t*k_3_th))/np.pi*A_time
				v_phi=v0sfl3*ARMS_DGauss(theta+delta_t*k_3_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))*ARMS_Gauss(phi+delta_t*k_3_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/np.pi*A_time
				k_4_th+=v_theta/R/np.pi
				k_4_ph+=v_phi/R/np.pi/np.sin(np.pi*(theta+delta_t*k_3_th))

	return np.array([theta+(k_1_th+2.0*k_2_th+2.0*k_3_th+k_4_th)*delta_t/6.0,phi+(k_1_ph+2.0*k_2_ph+2.0*k_3_ph+k_4_ph)*delta_t/6.0])


def find_N_max_idx(arr_in,N):
	arr_in=arr_in.flatten()
	values=np.zeros((N))
	indices=np.zeros((N),dtype="int32")
	for idx in range(len(arr_in)):
		if arr_in[idx]>values[0]:
			values[0]=arr_in[idx]
			indices[0]=idx
			sorted_idx=np.argsort(values)
			values=values[sorted_idx]
			indices=indices[sorted_idx]
	return indices[::-1],values[::-1]


def tidy_max_idx(indices,values,dimension):
	total_N=0
	for idx in range(len(values)):
		if values[idx]>0:
			total_N+=1
	indices_new=np.zeros((total_N,2),dtype="int32")
	for idx in range(total_N):
		indices_new[idx,0]=indices[idx]//dimension
		indices_new[idx,1]=indices[idx] % dimension
	return indices_new


def swap_two_angles(theta1,theta2):
	if theta1>theta2:
		temp1=theta1
		theta1=theta2
		theta2=temp1
	return theta1,theta2


def categorize_reconnected(CNT_file,start_file,end_file,start_time,end_time,delta_t,solar_Radius):
	"""
	Categorize connectivity - open/closed and reconnected/not
	TODO: adjust for R!=1.0
	"""
	from scipy.interpolate import interp1d
	Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)
	sfl2l=[]
	sfl2r=[]
	sfl3l=[]
	sfl3r=[]
	for idx in range(len(sflow_parameters)):
		if Time_profiles[idx]>0:
			sfl2l.append(sflow_parameters[idx][5])
			sfl2r.append(sflow_parameters[idx][6])
			sfl3l.append(sflow_parameters[idx][10])
			sfl3r.append(sflow_parameters[idx][11])
	num_parameters=len(sfl2l)
	sfl2l_m=min(sfl2l)
	sfl2r_m=max(sfl2r)
	sfl3l_m=min(sfl3l)
	sfl3r_m=max(sfl3r)

	#connection: 0 - open, 1 - closed
	R_st,theta_grid_st,phi_grid_st,Q_st=parse_QSL_Rbinfile(start_file)
	connection_start=0.5*(1.0+np.sign(Q_st))
	R_end,theta_grid_end,phi_grid_end,Q_end=parse_QSL_Rbinfile(end_file)
	connection_end=0.5*(1.0+np.sign(Q_end))

	#connection_map=0 +1 if start closed, +2 if end closed, +4 if end connection NOT consistent with flow (if flow advects open/closed field lines)
	connection_map=(connection_start+2.0*connection_end).astype('int32')

	idx_t_L=np.argmin(abs(1.0-theta_grid_st[0,:]/np.pi-sfl2l_m))
	idx_t_U=np.argmin(abs(1.0-theta_grid_st[0,:]/np.pi-sfl2r_m))
	idx_p_L=np.argmin(abs(phi_grid_st[:,0]/np.pi-sfl3l_m))
	idx_p_U=np.argmin(abs(phi_grid_st[:,0]/np.pi-sfl3r_m))
	idx_t_L,idx_t_U=swap_two_angles(idx_t_L,idx_t_U)
	idx_p_L,idx_p_U=swap_two_angles(idx_p_L,idx_p_U)
	theta_func=interp1d(1.0-theta_grid_st[0,:]/np.pi,np.array(range(0,len(theta_grid_st[0,:]))))
	phi_func=interp1d(phi_grid_st[:,0]/np.pi,np.array(range(0,len(phi_grid_st[:,0]))))

	# Extended cage 0 - phi, 1 - theta, 2 - idx_phi, 3 - idx_theta, 4 - start connection, 5 - end connection, 6 - is end connection consistent with flow? 0==yes, 1==no
	cage=[[],[],[],[],[],[],[]]
	for idx_t in range(idx_t_L,idx_t_U+1):
		for idx_p in range(idx_p_L,idx_p_U+1):
			inside_flow=False
			for idx_f in range(num_parameters):
				if 1.0-theta_grid_st[0,idx_t]/np.pi>=sfl2l[idx_f] and 1.0-theta_grid_st[0,idx_t]/np.pi<=sfl2r[idx_f] and phi_grid_st[idx_p,0]/np.pi>=sfl3l[idx_f] and phi_grid_st[idx_p,0]/np.pi<=sfl3r[idx_f]:
					inside_flow=True
			if inside_flow:
				cage[0].append(phi_grid_st[idx_p,0]/np.pi)
				cage[1].append(1.0-theta_grid_st[0,idx_t]/np.pi)
				cage[2].append(idx_p)
				cage[3].append(idx_t)
				cage[4].append(connection_start[idx_p,idx_t])
				cage[5].append(connection_end[idx_p,idx_t])
				cage[6].append(0.0)
	cage=np.array(cage)


	delta_t=abs(delta_t)
	t=end_time
	num_t=int((t-start_time)/delta_t)
	for idx_t in range(1,num_t+1):
		for idx_c in range(len(cage[0,:])):
			new_pos=surface_move_rk_4(solar_Radius,cage[1,idx_c],cage[0,idx_c],t,-delta_t,Time_profiles,sflow_parameters)
			cage[0,idx_c]=new_pos[1]
			cage[1,idx_c]=new_pos[0]
		t-=delta_t

	for idx_c in range(len(cage[0,:])):
		idx_p=int(np.round(phi_func(cage[0,idx_c])))
		idx_t=int(np.round(theta_func(cage[1,idx_c])))
		cage[6,idx_c]=0.5-0.5*np.sign(Q_end[int(cage[2,idx_c]),int(cage[3,idx_c])])*np.sign(Q_st[idx_p,idx_t])
		if not np.isnan(cage[6,idx_c]):
			connection_map[int(cage[2,idx_c]),int(cage[3,idx_c])]+=int(4.0*cage[6,idx_c])
	return R_st,theta_grid_st,phi_grid_st,Q_st,Q_end,connection_map


def save_connection_map(connection_filename,CNT_file,start_file,end_file,start_time,end_time,delta_t,solar_Radius):
	"""
	Calculate and store connectivity map
	"""
	R_st,theta_grid_st,phi_grid_st,Q_st,Q_end,connection_map=categorize_reconnected(CNT_file,start_file,end_file,start_time,end_time,delta_t,solar_Radius)
	connection_map=np.array(connection_map,dtype=np.int32)
	len_p=len(connection_map[:,0])
	len_t=len(connection_map[0,:])
	outfile=open(connection_filename,'wb')
	outfile.write(struct.pack('i',len_p))
	outfile.write(struct.pack('i',len_t))
	for idx_p in range(len_p):
		outfile.write(struct.pack('f',phi_grid_st[idx_p,0]))
	for idx_t in range(len_t):
		outfile.write(struct.pack('f',theta_grid_st[0,idx_t]))
	for idx_p in range(len_p):
		outfile.write(bytes(connection_map[idx_p,:]))
	outfile.close()


def load_connection_map(connection_filename):
	infile=open(connection_filename,"rb")
	len_p=struct.unpack('i',infile.read(4))[0]
	len_t=struct.unpack('i',infile.read(4))[0]
	phi_grid=np.zeros((len_p))
	theta_grid=np.zeros((len_t))
	connection_map=np.zeros((len_p,len_t),dtype=np.int32)
	for idx_p in range(len_p):
		phi_grid[idx_p]=struct.unpack('f',infile.read(4))[0]
	for idx_t in range(len_t):
		theta_grid[idx_t]=struct.unpack('f',infile.read(4))[0]
	theta_grid,phi_grid=np.meshgrid(theta_grid,phi_grid)
	for idx_p in range(len_p):
		connection_map[idx_p,:]=np.array(struct.unpack('i'*len_t,infile.read(4*len_t)))
	infile.close()
	return theta_grid,phi_grid,connection_map


def surfaceline_smooth_curve(line_points):
	"""
	Smooth otherwise blocky open/closed contour using interpolation
	"""
	num_pts=len(line_points[0,:])
	idx1=0
	while idx1<num_pts-1:
		#Horizontal
		if line_points[1,idx1]==line_points[1,idx1+1]:
			idx2=idx1+1
			while idx2<num_pts and line_points[1,idx1]==line_points[1,idx2]:
				idx2+=1
			idx2=min(idx2+1,num_pts-1)
			width=idx2-idx1
			for idx3 in range(1,width):
				line_points[1,idx1+idx3]=(line_points[1,idx2]-line_points[1,idx1])*idx3/width+line_points[1,idx1]
			idx1+=max(width,1)
		#Vertical
		elif line_points[0,idx1]==line_points[0,idx1+1]:
			idx2=idx1+1
			while idx2<num_pts and line_points[0,idx1]==line_points[0,idx2]:
				idx2+=1
			idx2=min(idx2+1,num_pts-1)
			height=idx2-idx1
			for idx3 in range(1,height):
				line_points[0,idx1+idx3]=(line_points[0,idx2]-line_points[0,idx1])*idx3/height+line_points[0,idx1]
			idx1+=max(height,1)
		else:
			idx1+=1
	return line_points


def surfaceline_smooth_Npt(line_points,N):
	"""
	Smooth otherwise blocky open/closed contour by N-point averaging
	"""
	num_pts=len(line_points[0,:])
	line_points_out=np.copy(line_points)
	if N%2==0:
		N+=1
	N_half=(N-1)//2
	for idx in range(N_half,num_pts-N_half):
		line_points_out[0,idx]=sum(line_points_out[0,idx-N_half:idx+N_half+1])/N
	return line_points_out


def surfaceline_inside_region(point_phi,point_theta,num_regions,region_phi,region_theta):
	"""
	Return index of point inside given region
	"""
	inside_region=False
	for idx in range(num_regions):
		if point_phi>=region_phi[idx][0] and point_phi<=region_phi[idx][1] and point_theta>=region_theta[idx][0] and point_theta<=region_theta[idx][1]:
			inside_region=True
			break
	if not inside_region:
		idx=num_regions
	return idx


def surfaceline_print_cage(line_cage,num_cage,label,num_cols=3):
	if num_cage==0:
		print("")
		print(label+'=([[],[]])')
		print("")
	else:
		print("")
		print(label+'=([[',end='')
		for idx in range(num_cage-1):
			print(line_cage[0,idx],end=',')
			if idx %num_cols==0:
				print("")
		print(line_cage[0,num_cage-1],end='],')
		print("")
		print('[',end='')
		for idx in range(num_cage-1):
			print(line_cage[1,idx],end=',')
			if idx %num_cols==0:
				print("")
		print(line_cage[1,num_cage-1],end=']])')
		print("")


def surfaceline_write_binary_cage(outfilename,static_cages,dynamic_cages):
	"""
	Write open/closed boundary points, "cages" to file
	static_cages, dynamic_cages are lists of numpy arrays
	dynamic_cages will be advected by surface flows, static_cages are assumed not to move
	"""
	outfile=open(outfilename,"wb")
	outfile.write(struct.pack('i',len(static_cages)))
	for idx in range(len(static_cages)):
		outfile.write(struct.pack('i',len(static_cages[idx][0,:])))
	outfile.write(struct.pack('i',len(dynamic_cages)))
	for idx in range(len(dynamic_cages)):
		outfile.write(struct.pack('i',len(dynamic_cages[idx][0,:])))
	for idx in range(len(static_cages)):
		arr_temp=np.array(static_cages[idx][0,:],dtype=np.float32)
		outfile.write(bytes(arr_temp))
		arr_temp=np.array(static_cages[idx][1,:],dtype=np.float32)
		outfile.write(bytes(arr_temp))
	for idx in range(len(dynamic_cages)):
		arr_temp=np.array(dynamic_cages[idx][0,:],dtype=np.float32)
		outfile.write(bytes(arr_temp))
		arr_temp=np.array(dynamic_cages[idx][1,:],dtype=np.float32)
		outfile.write(bytes(arr_temp))
	outfile.close()


def surfaceline_read_binary_cage(infilename):
	"""
	Read open/closed boundary points, "cages" to file
	static_cages, dynamic_cages are lists of numpy arrays
	dynamic_cages will be advected by surface flows, static_cages are assumed not to move
	"""
	infile=open(infilename,"rb")
	static_cages=[]
	num_static_cage=struct.unpack('i',infile.read(4))[0]
	len_static_cages=np.zeros((num_static_cage),dtype=np.int32)
	for idx in range(num_static_cage):
		len_static_cages[idx]=struct.unpack('i',infile.read(4))[0]
	dynamic_cages=[]
	num_dynamic_cage=struct.unpack('i',infile.read(4))[0]
	len_dynamic_cages=np.zeros((num_dynamic_cage),dtype=np.int32)
	for idx in range(num_dynamic_cage):
		len_dynamic_cages[idx]=struct.unpack('i',infile.read(4))[0]
	for idx in range(num_static_cage):
		cage_temp=np.zeros((2,len_static_cages[idx]),dtype=np.float32)
		cage_temp[0,:]=np.array(struct.unpack('f'*len_static_cages[idx],infile.read(4*len_static_cages[idx])))
		cage_temp[1,:]=np.array(struct.unpack('f'*len_static_cages[idx],infile.read(4*len_static_cages[idx])))
		static_cages.append(cage_temp)
	for idx in range(num_dynamic_cage):
		cage_temp=np.zeros((2,len_dynamic_cages[idx]),dtype=np.float32)
		cage_temp[0,:]=np.array(struct.unpack('f'*len_dynamic_cages[idx],infile.read(4*len_dynamic_cages[idx])))
		cage_temp[1,:]=np.array(struct.unpack('f'*len_dynamic_cages[idx],infile.read(4*len_dynamic_cages[idx])))
		dynamic_cages.append(cage_temp)
	infile.close()
	return static_cages,dynamic_cages


def surfaceline_advance_cage(cage_in,start_t,plot_t,delta_t,CNT_file,solar_Radius):
	"""
	Advect open/closed boundary points, "cages", by surface flows
	"""
	Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)
	end_t=plot_t[-1]
	cage=np.zeros((np.shape(cage_in)[0],np.shape(cage_in)[1]))
	cage_out=np.zeros((np.shape(cage_in)[0],np.shape(cage_in)[1],len(plot_t)))
	cage[0,:]=cage_in[0,:]/180.0
	cage[1,:]=1.0-cage_in[1,:]/180.0
	t=start_t
	idx_t_out=0
	while t<end_t:
		if t<plot_t[idx_t_out]-delta_t:
			delta_t_current=delta_t
			for idx_c in range(len(cage[0,:])):
				new_pos=surface_move_rk_4(solar_Radius,cage[1,idx_c],cage[0,idx_c],t,delta_t_current,Time_profiles,sflow_parameters)
				cage[0,idx_c]=new_pos[1]
				cage[1,idx_c]=new_pos[0]
		else:
			delta_t_current=plot_t[idx_t_out]-t
			for idx_c in range(len(cage[0,:])):
				new_pos=surface_move_rk_4(solar_Radius,cage[1,idx_c],cage[0,idx_c],t,delta_t_current,Time_profiles,sflow_parameters)
				cage[0,idx_c]=new_pos[1]
				cage[1,idx_c]=new_pos[0]
			cage_out[0,:,idx_t_out]=cage[0,:]*180.0
			cage_out[1,:,idx_t_out]=(1.0-cage[1,:])*180.0
			idx_t_out+=1
		t+=delta_t_current
	return cage_out


def read_fieldline_file(target_file):
	"""
	Read output of ARMS fieldline integrator
	Note: fieldline routine has been modified output binary data
	"""
	infile=open(target_file,"rb")
	is_spherical=struct.unpack('>i',infile.read(4))[0]
	x_periodic,y_periodic,z_periodic=struct.unpack('>iii',infile.read(12))
	fieldlines_saved=struct.unpack('>i',infile.read(4))[0]
	n1,n2,n3=struct.unpack('>iii',infile.read(12))
	rlo=np.array(struct.unpack('>ddd',infile.read(24)))
	rhi=np.array(struct.unpack('>ddd',infile.read(24)))
	nsteps=struct.unpack('>i',infile.read(4))[0]
	infile.read(4)
	step_start,step_max,error_bound=struct.unpack('>ddd',infile.read(24))
	turnmap=struct.unpack('>i',infile.read(4))[0]
	infile.read(4)
	if n1==0:
		nn1=n2
		nn2=n3
		R=rlo[0]
		theta=np.linspace(rlo[1],rhi[1],num=nn1+1)
		phi=np.linspace(rlo[2],rhi[2],num=nn2+1)
		data_available=True
	elif n2==0:
		nn1=n1
		nn2=n3
		R=np.linspace(rlo[0],rhi[0],num=nn1+1)
		theta=rlo[1]
		phi=np.linspace(rlo[2],rhi[2],num=nn2+1)
		data_available=True
	elif n3==0:
		nn1=n1
		nn2=n2
		R=np.linspace(rlo[0],rhi[0],num=nn1+1)
		theta=np.linspace(rlo[1],rhi[1],num=nn2+1)
		phi=rlo[2]
		data_available=True
	else:
		R=np.array([0.0])
		theta=np.array([0.0])
		phi=np.array([0.0])
		print("Error: One output dimension must be zero")
		data_available=False
	if data_available:
		B_ends=np.array(struct.unpack('>'+'d'*(nn1+1)*(nn2+1)*6,infile.read(8*(nn1+1)*(nn2+1)*6))).reshape(nn1+1,nn2+1,2,3)
		if fieldlines_saved:
			flines=np.array(struct.unpack('>'+'d'*(nn1+1)*(nn2+1)*(nsteps*2-1)*3,infile.read(8*(nn1+1)*(nn2+1)*(nsteps*2-1)*3))).reshape(nn1+1,nn2+1,(nsteps*2-1),3)
		else:
			flines=np.zeros((1,1,1,3))
		mapping=np.array(struct.unpack('>'+'d'*(nn1+1)*(nn2+1)*6,infile.read(8*(nn1+1)*(nn2+1)*6))).reshape(nn1+1,nn2+1,2,3)
		if n1==0:
			Q=np.array(struct.unpack('>'+'d'*(nn1+1)*(nn2+1),infile.read(8*(nn1+1)*(nn2+1)))).reshape(nn1+1,nn2+1)
		else:
			Q=np.zeros((1,1))
	else:
		B_ends=np.zeros((1,1,2,3))
		flines=np.zeros((1,1,1,3))
		mapping=np.zeros((1,1,2,3))
		Q=np.zeros((1,1))
	infile.close()
	return R,theta,phi,B_ends,flines,mapping,Q


def fieldlines_to_QSL(in_filename,out_filename):
	"""
	Convert fieldline binary data to QSL binary format
	Note: both programs have been altered
	"""
	R,theta,phi,B_ends,flines,mapping,Q=read_fieldline_file(in_filename)
	outfile=open(out_filename,"wb")
	for idx_p in range(len(phi)):
		for idx_t in range(len(theta)):
			outfile.write(struct.pack("=dddf",(phi[idx_p]+180.0)*DEG2RAD,theta[idx_t]*DEG2RAD,R,Q[idx_t,idx_p]))
	outfile.close()


def read_fieldline_ff_file(target_file):
	"""
	Read output of ARMS fieldline_ff integrator
	ff == from file
	Given that this does not use a regular grid, Q is not computed
	"""
	infile=open(target_file,"rb")
	is_spherical=struct.unpack('>i',infile.read(4))[0]
	x_periodic,y_periodic,z_periodic=struct.unpack('>iii',infile.read(12))
	fieldlines_saved=struct.unpack('>i',infile.read(4))[0]
	nsteps=struct.unpack('>i',infile.read(4))[0]
	step_start,step_max,error_bound=struct.unpack('>ddd',infile.read(24))
	turnmap=struct.unpack('>i',infile.read(4))[0]
	nn_tot=struct.unpack('>i',infile.read(4))[0]
	gridfl=np.array(struct.unpack('>'+'d'*nn_tot*3,infile.read(8*nn_tot*3))).reshape(nn_tot,3)
	B_ends=np.array(struct.unpack('>'+'d'*nn_tot*6,infile.read(8*nn_tot*6))).reshape(nn_tot,2,3)
	if fieldlines_saved:
		flines=np.array(struct.unpack('>'+'d'*nn_tot*(nsteps*2-1)*3,infile.read(8*nn_tot*(nsteps*2-1)*3))).reshape(nn_tot,(nsteps*2-1),3)
	else:
		flines=np.zeros((1,1,3))
	mapping=np.array(struct.unpack('>'+'d'*nn_tot*6,infile.read(8*nn_tot*6))).reshape(nn_tot,2,3)
	infile.close()
	return gridfl,B_ends,flines,mapping


#Return possible processors and steps
def get_probe_files(probe_directory):
	probe_files = [f for f in os.listdir(probe_directory) if os.path.isfile(os.path.join(probe_directory, f)) and os.path.getsize(os.path.join(probe_directory, f))>0 and 'probe.' in f and len(f.split('.'))==3 and f.split('.')[1].isnumeric() and f.split('.')[2].isnumeric()]
	probe_files.sort()
	probe_procs=[]
	probe_steps=[]
	for idx in range(len(probe_files)):
		if os.path.getsize(os.path.join(probe_directory,probe_files[idx]))>12:
			proc=int(probe_files[idx].split('.')[1])
			step=int(probe_files[idx].split('.')[2])
			if proc not in probe_procs:
				probe_procs.append(proc)
			if step not in probe_steps:
				probe_steps.append(step)
	return probe_procs,probe_steps


def get_probe_blocks(probe_directory,proc,step,ndim=3):
	"""
	probe version of ARMS: return data from selected grid points at every simulation timestep
	Return available block numbers and coordinates for given processor/timestep
	"""
	filename=os.path.join(probe_directory,"probe."+"{:04d}".format(proc)+".{:07d}".format(step))
	infile=open(filename,"rb")
	infile.read(4)
	num_blocks=struct.unpack('>i', infile.read(4))[0]
	infile.read(4)
	probe_blocks=np.zeros((num_blocks),dtype=np.int32)
	bnd_box=np.zeros((num_blocks,2,ndim),dtype=np.float64)
	for idx in range(num_blocks):
		infile.read(4)
		probe_blocks[idx]=struct.unpack('>i', infile.read(4))[0]
		infile.read(4)
		for j in range(ndim):
			for i in range(2):
				infile.read(4)
				bnd_box[idx,i,j]=struct.unpack('>d', infile.read(8))[0]
				infile.read(4)
	infile.close()
	return probe_blocks,bnd_box


def get_probe_block_location(probe_directory,step,target_R,target_theta,target_phi,ndim=3):
	"""
	probe version of ARMS: return data from selected grid points at every simulation timestep
	Return processor+block number bounding given point (if available)

	target_R [Rsun]
	target_theta [deg], -90->+90
	target_phi [deg], -180->+180
	"""
	probe_files = [f for f in os.listdir(probe_directory) if os.path.isfile(os.path.join(probe_directory, f)) and os.path.getsize(os.path.join(probe_directory, f))>0 and 'probe.' in f and len(f.split('.'))==3 and f.split('.')[1].isnumeric() and f.split('.')[2].isnumeric()]
	probe_files.sort()
	block_found=False
	target_proc=0
	target_block=0
	target_points=[0,0]
	for idx in range(len(probe_files)):
		if os.path.getsize(os.path.join(probe_directory,probe_files[idx]))>12:
			proc_curr=int(probe_files[idx].split('.')[1])
			step_curr=int(probe_files[idx].split('.')[2])
			if step_curr==step:
				probe_blocks,bnd_box=get_probe_blocks(probe_directory,proc_curr,step_curr,ndim=ndim)
				bnd_box[:,:,0]=np.exp(bnd_box[:,:,0])
				bnd_box[:,:,1]=(0.5-bnd_box[:,:,1])*180.0
				bnd_box[:,:,2]=bnd_box[:,:,2]*180.0
				for idx3 in range(len(bnd_box[:,0,0])):
					if (bnd_box[idx3,0,0]>=target_R and bnd_box[idx3,1,0]<=target_R) or (bnd_box[idx3,0,0]<=target_R and bnd_box[idx3,1,0]>=target_R):
						if (bnd_box[idx3,0,1]>=target_theta and bnd_box[idx3,1,1]<=target_theta) or (bnd_box[idx3,0,1]<=target_theta and bnd_box[idx3,1,1]>=target_theta):
							if (bnd_box[idx3,0,2]>=target_phi and bnd_box[idx3,1,2]<=target_phi) or (bnd_box[idx3,0,2]<=target_phi and bnd_box[idx3,1,2]>=target_phi):
								block_found=True
								target_proc=proc_curr
								target_block=probe_blocks[idx3]
								n2pm1=8
								n3pm1=8
								target_points[0]=int(np.round((target_theta-bnd_box[idx3,0,1])*n2pm1/(bnd_box[idx3,1,1]-bnd_box[idx3,0,1])))
								if target_points[0]==n2pm1:
									target_points[0]=n2pm1-1
								target_points[1]=int(np.round((target_phi-bnd_box[idx3,0,2])*n3pm1/(bnd_box[idx3,1,2]-bnd_box[idx3,0,2])))
								if target_points[1]==n3pm1:
									target_points[1]=n3pm1-1
		if block_found:
			break
	return block_found,target_proc,target_block,target_points
		

def read_probe_file(probe_directory,proc,step,ndim=3):
	"""
	probe version of ARMS: return data from selected grid points at every simulation timestep

	probe_data [(time step),(block index),(see below),(theta index),(phi index)]:
	[:,:,0,:,:] ! Bx
	[:,:,1,:,:] ! By
	[:,:,2,:,:] ! Bz
	[:,:,3,:,:] ! mass density
	[:,:,4,:,:] ! x-momentum
	[:,:,5,:,:] ! y-momentum
	[:,:,6,:,:] ! z-momentum
	[:,:,7,:,:] ! energy density

	NOTE! Python vs Fortran indices!

	Bounds:
	bnd_box[ (Number of Blocks) , (Lower, upper bound) , (Dimension) ] - all blocks
	Lower/upper:
	bnd_box[:,0,:]<leaf_bnd_box[:,1,:] (Always)
	bnd_box[:,:,0] ! r
	bnd_box[:,:,1] ! theta
	bnd_box[:,:,2] ! phi
	NOTE! theta=[0,1] (reversed), phi=[-1,1]
	"""

	filename=os.path.join(probe_directory,"probe."+"{:04d}".format(proc)+".{:07d}".format(step))
	filesize=os.path.getsize(filename)
	#TODO - further checks if file is malformed
	infile=open(filename,"rb")
	infile.read(4)
	num_blocks=struct.unpack('>i', infile.read(4))[0]
	infile.read(4)
	probe_blocks=np.zeros((num_blocks),dtype=np.int32)
	bnd_box=np.zeros((num_blocks,2,ndim),dtype=np.float64)
	for idx in range(num_blocks):
		infile.read(4)
		probe_blocks[idx]=struct.unpack('>i', infile.read(4))[0]
		infile.read(4)
		for j in range(ndim):
			for i in range(2):
				infile.read(4)
				bnd_box[idx,i,j]=struct.unpack('>d', infile.read(8))[0]
				infile.read(4)

	num_steps=(filesize-12-108*num_blocks)//(16+num_blocks*4104)
	probe_times=np.zeros((num_steps),dtype=np.float64)
	probe_data=np.zeros((num_steps,num_blocks,8,8,8),dtype=np.float64)
	for idx2 in range(num_steps):
		infile.read(4)
		probe_times[idx2]=struct.unpack('>d', infile.read(8))[0]
		infile.read(4)
		for idx in range(num_blocks):
			infile.read(4)
			data_temp=np.array(struct.unpack('>'+'d'*512, infile.read(4096)))
			probe_data[idx2,idx,:,:,:]=data_temp.reshape((8,8,8))
			infile.read(4)
	infile.close()
	return probe_blocks,probe_times,probe_data


def read_probe_singleblock(probe_directory,proc,block,ndim=3):
	"""
	probe version of ARMS: return data from selected grid points at every simulation timestep
	Routine extracts and stitches together all times of data for given processor+block number
	"""
	probe_files = [f for f in os.listdir(probe_directory) if os.path.isfile(os.path.join(probe_directory, f)) and os.path.getsize(os.path.join(probe_directory, f))>0 and 'probe.' in f and len(f.split('.'))==3 and f.split('.')[1].isnumeric() and f.split('.')[2].isnumeric()]
	probe_files.sort()
	probe_steps=[]
	for idx in range(len(probe_files)):
		if os.path.getsize(os.path.join(probe_directory,probe_files[idx]))>12:
			proc_curr=int(probe_files[idx].split('.')[1])
			step_curr=int(probe_files[idx].split('.')[2])
			if proc_curr==proc:
				probe_steps.append(step_curr)
	if len(probe_steps)==0:
		return np.zeros((1),dtype=np.float64),np.zeros((1,8,8,8),dtype=np.float64)	#No valid timesteps, return zeros
	elif len(probe_steps)==1:
		probe_blocks_curr,probe_times_curr,probe_data_curr=read_probe_file(probe_directory,proc,probe_steps[0],ndim=ndim)
		if block in probe_blocks_curr:
			idx2=np.where(probe_blocks_curr==block)[0][0]
			return probe_times_curr,probe_data_curr[:,idx2,:,:,:]
		else:
			return np.zeros((1),dtype=np.float64),np.zeros((1,8,8,8),dtype=np.float64)	#Block cannot be found, return zeros
	else:
		probe_steps.sort()
		times_list=[]
		data_list=[]
		for idx in range(len(probe_steps)):
			probe_blocks_curr,probe_times_curr,probe_data_curr=read_probe_file(probe_directory,proc,probe_steps[idx],ndim=ndim)
			if block in probe_blocks_curr:
				idx2=np.where(probe_blocks_curr==block)[0][0]
				times_list.append(probe_times_curr)
				data_list.append(probe_data_curr[:,idx2,:,:,:])
		if len(times_list)==0:
			return np.zeros((1),dtype=np.float64),np.zeros((1,8,8,8),dtype=np.float64)	#Block cannot be found, return zeros
		else:
			timesteps_total=0
			for idx in range(len(probe_steps)-1):
				timesteps_total+=min(len(times_list[idx]),probe_steps[idx+1]-probe_steps[idx])
			timesteps_total+=len(times_list[-1])
			probe_times=np.zeros((timesteps_total),dtype=np.float64)
			probe_data=np.zeros((timesteps_total,8,8,8),dtype=np.float64)
			idx2=0
			for idx in range(len(probe_steps)-1):
				len_curr=min(len(times_list[idx]),probe_steps[idx+1]-probe_steps[idx])
				probe_times[idx2:idx2+len_curr]=times_list[idx][:len_curr]
				probe_data[idx2:idx2+len_curr,:,:,:]=data_list[idx][:len_curr,:,:,:]
				idx2+=len_curr
			probe_times[idx2:]=times_list[-1][:]
			probe_data[idx2:,:,:,:]=data_list[-1][:,:,:,:]
			return probe_times,probe_data


def read_arms_prestart(folder_path,timestep,nfaces=6,nchild=8,ndim=3):
	"""
	Read parallel restart files in a folder; blocks will be padded with zeros
	TODO for future: replace padding zeros with existing or interpolated data
	
	fvar1 (python indices):
	[:,:,:,:,0] ! mass density
	[:,:,:,:,1] ! x-momentum
	[:,:,:,:,2] ! y-momentum
	[:,:,:,:,3] ! z-momentum
	[:,:,:,:,4] ! energy density
	[:,:,:,:,5] ! pressure
	[:,:,:,:,6] ! temperature

	Bounds:
	bnd_box[ (Number of Blocks) , (Lower, upper bound) , (Dimension) ] - all blocks
	leaf_bnd_box[...] - just leaf blocks
	Lower/upper:
	leaf_bnd_box[:,0,:]<leaf_bnd_box[:,1,:] (Always)
	leaf_bnd_box[:,:,0] ! r
	leaf_bnd_box[:,:,1] ! theta
	leaf_bnd_box[:,:,2] ! phi
	NOTE! theta=[0,1] (reversed), phi=[-1,1]
	"""

	#Get list of files
	timestep_str="{:07}".format(timestep)
	rst_files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f)) and 'arms.rst.' in f and timestep_str in f]
	rst_files.sort()

	#Check processors are sequential
	if len(rst_files)==0:
		print("Error: no rst files found")
		processors_consistent=False
	elif len(rst_files)==1:
		if int(rst_files[0].split(".")[2])==0:
			processors_consistent=True
		else:
			print("Error: files missing")
			processors_consistent=False
	else:
		if int(rst_files[0].split(".")[2])==0:
			processors_consistent=True
			for idx in range(1,len(rst_files)):
				if int(rst_files[idx-1].split(".")[2])!=int(rst_files[idx].split(".")[2])-1:
					processors_consistent=False
			if not processors_consistent:
				print("Error: files missing")
		else:
			print("Error: files missing")
			processors_consistent=False
	#################################

	if processors_consistent:
		#First pass: loop over files to check consistency, get total number of leaf blocks
		lunkvar=np.zeros((len(rst_files)),dtype=np.int32)
		lfacevar=np.zeros((len(rst_files)),dtype=np.int32)
		lflags=np.zeros((len(rst_files)),dtype=np.int32)
		time=np.zeros((len(rst_files)),dtype=np.float64)
		deltat=np.zeros((len(rst_files)),dtype=np.float64)
		jstep=np.zeros((len(rst_files)),dtype=np.int32)
		nxb=np.zeros((len(rst_files)),dtype=np.int32)
		nyb=np.zeros((len(rst_files)),dtype=np.int32)
		nzb=np.zeros((len(rst_files)),dtype=np.int32)
		lnblocks=np.zeros((len(rst_files)),dtype=np.int32)
		num_leafblocks=0
		for idx in range(len(rst_files)):
			infile=open(os.path.join(folder_path,rst_files[idx]),"rb")

			sdf_type=infile.read(4)
			if sdf_type!=b'SDF1':
				print("Warning: wrong SDF type, "+sdf_type)
			lunkvar[idx]=struct.unpack('>i', infile.read(4))[0]
			lfacevar[idx]=struct.unpack('>i', infile.read(4))[0]
			lflags[idx]=struct.unpack('>i', infile.read(4))[0]
			time[idx]=struct.unpack('>d', infile.read(8))[0]
			deltat[idx]=struct.unpack('>d', infile.read(8))[0]
			jstep[idx]=struct.unpack('>i', infile.read(4))[0]
			nxb[idx]=struct.unpack('>i', infile.read(4))[0]
			nyb[idx]=struct.unpack('>i', infile.read(4))[0]
			nzb[idx]=struct.unpack('>i', infile.read(4))[0]
			lnblocks[idx]=struct.unpack('>i', infile.read(4))[0]

			for lb in range(lnblocks[idx]):
				lrefine=struct.unpack('>i', infile.read(4))[0]
				nodetype=struct.unpack('>i', infile.read(4))[0]
				which_child=struct.unpack('>i', infile.read(4))[0]
				infile.read(4)
				infile.read(4)
				infile.read(nfaces*2*4)
				infile.read(nchild*2*4)
				infile.read(lflags[idx]*4)
				infile.read(ndim*2*8)
				if nodetype==1:
					num_leafblocks+=1
					infile.read(lunkvar[idx]*nzb[idx]*nyb[idx]*nxb[idx]*8)
					if lfacevar[idx]==1:
						infile.read(nzb[idx]*nyb[idx]*(nxb[idx]+1)*8)
						infile.read(nzb[idx]*(nyb[idx]+1)*nxb[idx]*8)
						infile.read((nzb[idx]+1)*nyb[idx]*nxb[idx]*8)
			infile.close()
		field_output=False
		files_consistent=True
		for idx in range(1,len(rst_files)):
			if lunkvar[idx-1]!=lunkvar[idx] or lfacevar[idx-1]!=lfacevar[idx] or lflags[idx-1]!=lflags[idx] or time[idx-1]!=time[idx] or deltat[idx-1]!=deltat[idx] or jstep[idx-1]!=jstep[idx] or nxb[idx-1]!=nxb[idx] or nyb[idx-1]!=nyb[idx] or nzb[idx-1]!=nzb[idx]:
				files_consistent=False

		if files_consistent:
			if lfacevar[0]==1:
				field_output=True

			#Second pass: allocate arrays and loop over files to get data
			n1=nxb[0]
			n2=nyb[0]
			n3=nzb[0]
			n1p=n1+1
			n2p=n2+1
			n3p=n3+1
			num_totalblocks=sum(lnblocks)
			lblock_idx=np.zeros((len(rst_files),max(lnblocks)),dtype=np.int32)
			lrefine=np.zeros((num_totalblocks),dtype=np.int32)
			nodetype=np.zeros((num_totalblocks),dtype=np.int32)
			which_child=np.zeros((num_totalblocks),dtype=np.int32)
			parent=np.zeros((2,num_totalblocks),dtype=np.int32)
			neigh=np.zeros((2,nfaces,num_totalblocks),dtype=np.int32)
			leaf_neigh=np.zeros((2,nfaces,num_leafblocks),dtype=np.int32)
			child=np.zeros((2,nchild,num_totalblocks),dtype=np.int32)
			bflags=np.zeros((lflags[0],num_totalblocks),dtype=np.int32)
			bnd_box=np.zeros((num_totalblocks,2,ndim),dtype=np.float64)
			leaf_bnd_box=np.zeros((num_leafblocks,2,ndim),dtype=np.float64)
			fvar1=np.zeros((num_leafblocks,n3p,n2p,n1p,lunkvar[0]),dtype=np.float64)
			bvar=np.zeros((num_leafblocks,n3p,n2p,n1p,3),dtype=np.float64)
			idx_block=0	#Normal block
			idx_lblock=0	#Any leafblock
			for idx in range(len(rst_files)):
				infile=open(os.path.join(folder_path,rst_files[idx]),"rb")
				infile.read(9*4+2*8)
				for lb in range(lnblocks[idx]):
					lrefine[idx_block]=struct.unpack('>i', infile.read(4))[0]
					nodetype[idx_block]=struct.unpack('>i', infile.read(4))[0]
					which_child[idx_block]=struct.unpack('>i', infile.read(4))[0]
					parent[0,idx_block]=struct.unpack('>i', infile.read(4))[0]
					parent[1,idx_block]=struct.unpack('>i', infile.read(4))[0]
					for j in range(nfaces):
						for i in range(2):
							neigh[i,j,idx_block]=struct.unpack('>i', infile.read(4))[0]
					for j in range(nchild):
						for i in range(2):
							child[i,j,idx_block]=struct.unpack('>i', infile.read(4))[0]
					for i in range(lflags[0]):
						bflags[i,idx_block]=struct.unpack('>i', infile.read(4))[0]
					for j in range(ndim):
						for i in range(2):
							bnd_box[idx_block,i,j]=struct.unpack('>d', infile.read(8))[0]
					if nodetype[idx_block]==1:
						lblock_idx[idx,lb]=idx_lblock
						leaf_neigh[:,:,idx_lblock]=neigh[:,:,idx_block]
						leaf_bnd_box[idx_lblock,:,:]=bnd_box[idx_block,:,:]
						for l in range(lunkvar[0]):
							for k in range(n3):
								for j in range(n2):
									for i in range(n1):
										fvar1[idx_lblock,k,j,i,l]=struct.unpack('>d', infile.read(8))[0]
						if field_output:
							for k in range(n3):
								for j in range(n2):
									for i in range(n1p):
										bvar[idx_lblock,k,j,i,0]=struct.unpack('>d', infile.read(8))[0]
							for k in range(n3):
								for j in range(n2p):
									for i in range(n1):
										bvar[idx_lblock,k,j,i,1]=struct.unpack('>d', infile.read(8))[0]
							for k in range(n3p):
								for j in range(n2):
									for i in range(n1):
										bvar[idx_lblock,k,j,i,2]=struct.unpack('>d', infile.read(8))[0]
						idx_lblock+=1
					idx_block+=1
				infile.close()
			data_good=True
		else: #not files_consistent
			data_good=False
			num_totalblocks=1
			num_leafblocks=1
			leaf_neigh=np.zeros((2,nfaces,num_leafblocks),dtype=np.int32)
			leaf_bnd_box=leaf_bnd_box=np.zeros((num_leafblocks,2,ndim),dtype=np.float64)
			fvar1=np.zeros((num_leafblocks,1,1,1,1),dtype=np.float64)
			bvar=np.zeros((num_leafblocks,1,1,1,3),dtype=np.float64)
	else: #not processors_consistent
		data_good=False
		num_totalblocks=1
		num_leafblocks=1
		leaf_neigh=np.zeros((2,nfaces,num_leafblocks),dtype=np.int32)
		leaf_bnd_box=leaf_bnd_box=np.zeros((num_leafblocks,2,ndim),dtype=np.float64)
		fvar1=np.zeros((num_leafblocks,1,1,1,1),dtype=np.float64)
		bvar=np.zeros((num_leafblocks,1,1,1,3),dtype=np.float64)
	return data_good,num_totalblocks,num_leafblocks,lblock_idx,leaf_neigh,leaf_bnd_box,fvar1,field_output,bvar




