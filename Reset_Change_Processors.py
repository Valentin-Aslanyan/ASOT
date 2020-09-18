

new_processors=64
folder_path="./PDMP_drestarted"
output_path="./PDMP_new"
timestep=505219


import sys
import os
import numpy as np
import struct

nfaces=6
nchild=8
ndim=3

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
				infile.read(lunkvar[idx]*nzb[idx]*nyb[idx]*nxb[idx]*8)
				if lfacevar[idx]==1:
					infile.read(nzb[idx]*nyb[idx]*(nxb[idx]+1)*8)
					infile.read(nzb[idx]*(nyb[idx]+1)*nxb[idx]*8)
					infile.read((nzb[idx]+1)*nyb[idx]*nxb[idx]*8)
		infile.close()
	n1=nxb[0]
	n2=nyb[0]
	n3=nzb[0]
	n1p=n1+1
	n2p=n2+1
	n3p=n3+1
	num_totalblocks=sum(lnblocks)
	lrefine=np.zeros((num_totalblocks),dtype=np.int32)
	nodetype=np.zeros((num_totalblocks),dtype=np.int32)
	which_child=np.zeros((num_totalblocks),dtype=np.int32)
	parent=np.zeros((2,num_totalblocks),dtype=np.int32)
	neigh=np.zeros((2,nfaces,num_totalblocks),dtype=np.int32)
	child=np.zeros((2,nchild,num_totalblocks),dtype=np.int32)
	bnd_box=np.zeros((num_totalblocks,2,ndim),dtype=np.float64)
	bflags=np.zeros((lflags[0],num_totalblocks),dtype=np.int32)

	processor_idx=np.zeros((num_totalblocks),dtype=np.int32)
	lblock_idx=np.zeros((num_totalblocks),dtype=np.int32)
	new_blockperproc=int(np.ceil(num_totalblocks/new_processors))
	lnblocks_new=np.ones((new_processors),dtype=np.int32)*new_blockperproc
	lnblocks_new[-1]=num_totalblocks-new_blockperproc*(new_processors-1)
	old2glbidx=np.zeros((max(lnblocks)+1,len(rst_files)),dtype=np.int32)
	glbidx2old=np.zeros((num_totalblocks,2),dtype=np.int32)
	glbidx2new=np.zeros((num_totalblocks,2),dtype=np.int32)
	new2glbidx=np.zeros((new_blockperproc,new_processors),dtype=np.int32)
	idx_block=0
	new_proc=0
	new_block=1
	for idx in range(len(rst_files)):
		for lb in range(lnblocks[idx]):
			processor_idx[idx_block]=idx
			lblock_idx[idx_block]=lb+1
			old2glbidx[lb+1,idx]=idx_block
			glbidx2old[idx_block,0]=lb+1
			glbidx2old[idx_block,1]=idx
			glbidx2new[idx_block,0]=new_block
			glbidx2new[idx_block,1]=new_proc
			idx_block+=1
			if new_block<new_blockperproc:
				new_block+=1
			else:
				new_block=1
				new_proc+=1

	field_output=False
	files_consistent=True
	for idx in range(1,len(rst_files)):
		if lunkvar[idx-1]!=lunkvar[idx] or lfacevar[idx-1]!=lfacevar[idx] or lflags[idx-1]!=lflags[idx] or time[idx-1]!=time[idx] or deltat[idx-1]!=deltat[idx] or jstep[idx-1]!=jstep[idx] or nxb[idx-1]!=nxb[idx] or nyb[idx-1]!=nyb[idx] or nzb[idx-1]!=nzb[idx]:
			files_consistent=False

	if files_consistent:
		if lfacevar[0]==1:
			field_output=True
		out_start=sdf_type+struct.pack('>i', lunkvar[0])+struct.pack('>i', lfacevar[0])+struct.pack('>i', lflags[0])+struct.pack('>d', time[0])+struct.pack('>d', deltat[0])+struct.pack('>i', jstep[0])+struct.pack('>i', nxb[0])+struct.pack('>i', nyb[0])+struct.pack('>i', nzb[0])

		#Second pass: allocate arrays and loop over files to get data
		idx_block=0	#Normal block
		new_proc=0
		new_block=1
		outfile=open(os.path.join(output_path,'arms.rst.'+'{:04d}'.format(new_proc)+'.'+'{:07d}'.format(timestep)),"wb")
		outfile.write(out_start)
		outfile.write(struct.pack('>i', lnblocks_new[new_proc]))
		for idx in range(len(rst_files)):
			infile=open(os.path.join(folder_path,rst_files[idx]),"rb")
			infile.read(9*4+2*8)
			for lb in range(lnblocks[idx]):
				lrefine[idx_block]=struct.unpack('>i', infile.read(4))[0]
				outfile.write(struct.pack('>i', lrefine[idx_block]))
				nodetype[idx_block]=struct.unpack('>i', infile.read(4))[0]
				outfile.write(struct.pack('>i', nodetype[idx_block]))
				which_child[idx_block]=struct.unpack('>i', infile.read(4))[0]
				outfile.write(struct.pack('>i', which_child[idx_block]))
				parent[0,idx_block]=struct.unpack('>i', infile.read(4))[0]
				parent[1,idx_block]=struct.unpack('>i', infile.read(4))[0]
				if parent[0,idx_block]>=0 or parent[1,idx_block]>=0:
					parent_new=glbidx2new[old2glbidx[parent[0,idx_block],parent[1,idx_block]],:]
				else:
					parent_new=parent[:,idx_block]
				outfile.write(struct.pack('>i', parent_new[0]))
				outfile.write(struct.pack('>i', parent_new[1]))
				for j in range(nfaces):
					for i in range(2):
						neigh[i,j,idx_block]=struct.unpack('>i', infile.read(4))[0]
					if neigh[0,j,idx_block]>=0 or neigh[1,j,idx_block]>=0:
						neigh_new=glbidx2new[old2glbidx[neigh[0,j,idx_block],neigh[1,j,idx_block]],:]
					else:
						neigh_new=neigh[:,j,idx_block]
					outfile.write(struct.pack('>i', neigh_new[0]))
					outfile.write(struct.pack('>i', neigh_new[1]))
				for j in range(nchild):
					for i in range(2):
						child[i,j,idx_block]=struct.unpack('>i', infile.read(4))[0]
					if child[0,j,idx_block]>=0 or child[1,j,idx_block]>=0:
						child_new=glbidx2new[old2glbidx[child[0,j,idx_block],child[1,j,idx_block]],:]
					else:
						child_new=child[:,j,idx_block]
					outfile.write(struct.pack('>i', child_new[0]))
					outfile.write(struct.pack('>i', child_new[1]))
				for i in range(lflags[0]):
					bflags[i,idx_block]=struct.unpack('>i', infile.read(4))[0]
					outfile.write(struct.pack('>i', bflags[i,idx_block]))
				for j in range(ndim):
					for i in range(2):
						bnd_box[idx_block,i,j]=struct.unpack('>d', infile.read(8))[0]
						outfile.write(struct.pack('>d', bnd_box[idx_block,i,j]))
				if nodetype[idx_block]==1:
					for l in range(lunkvar[0]):
						for k in range(n3):
							for j in range(n2):
								for i in range(n1):
									outfile.write(infile.read(8))
					if field_output:
						for k in range(n3):
							for j in range(n2):
								for i in range(n1p):
									outfile.write(infile.read(8))
						for k in range(n3):
							for j in range(n2p):
								for i in range(n1):
									outfile.write(infile.read(8))
						for k in range(n3p):
							for j in range(n2):
								for i in range(n1):
									outfile.write(infile.read(8))
				idx_block+=1
				if new_block<new_blockperproc:
					new_block+=1
				else:
					new_block=1
					new_proc+=1
					outfile.close()
					outfile=open(os.path.join(output_path,'arms.rst.'+'{:04d}'.format(new_proc)+'.'+'{:07d}'.format(timestep)),"wb")
					outfile.write(out_start)
					outfile.write(struct.pack('>i', lnblocks_new[new_proc]))
			infile.close()
		outfile.close()
	else:
		print("Files inconsistent")

else:
	print("Processors inconsistent")
			





