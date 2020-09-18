

base_directory="./"


import os
import sys

if len(sys.argv)<2:
	print("Specify subfolder as command line argument")
else:
	target_directory=os.path.join(base_directory,sys.argv[1])

	prt_files = [f for f in os.listdir(target_directory) if os.path.isfile(os.path.join(target_directory, f)) and '.prt' in f]

	flicks_steps=[]
	bfield_steps=[]
	for fname in prt_files:
		infile=open(os.path.join(target_directory, fname),"r")
		for line in infile:
			if "flicks" in line.lower():
				flicks_steps.append([int(line.split()[4]),float(line.split()[7])])
			if "bfield" in line.lower():
				bfield_steps.append([int(line.split()[4]),float(line.split()[7])])
		infile.close()

	flicks_steps.sort()
	bfield_steps.sort()

	if len(flicks_steps)>0:
		print(flicks_steps[0][0],flicks_steps[0][1])
		for idx in range(1,len(flicks_steps)):
			if flicks_steps[idx][0]!=flicks_steps[idx-1][0]:
				print(flicks_steps[idx][0],flicks_steps[idx][1])
