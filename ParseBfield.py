

import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *


#Read header file
hdrfile=open("flicks.hdr","r")
for idx in range(3):
	hdrfile.readline()
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
flicksfile=open("bfield.0004486","rb")

print(flicksfile.read(4))

#time_raw=flicksfile.read(4)
#print(time_raw,struct.unpack('>i', time_raw),struct.unpack('<i', time_raw),struct.unpack('>f', time_raw),struct.unpack('<f', time_raw))
time=struct.unpack('>d', flicksfile.read(8))[0]
print(time)
print(flicksfile.read(8))
ntotl=struct.unpack('>i', flicksfile.read(4))[0]
nleaf=struct.unpack('>i', flicksfile.read(4))[0]
print(ntotl,nleaf)
print(flicksfile.read(8))
print(struct.unpack('>i', b'\x00\x00\x00\x8c')[0])

for idx in range(1):#range(ntotl):
	iputwrk=struct.unpack('>'+35*'i', flicksfile.read(35*4))
	rputwrk=struct.unpack('>'+6*'f', flicksfile.read(6*4))
	print(iputwrk)
	print(rputwrk)


flicksfile.close()
