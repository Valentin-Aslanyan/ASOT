
"""
Print Python lists of surface flow parameters in a format for the ARMS control file
"""

rlsfl=7.000e+10
rrsfl=7.020e+10
rcsfl=7.010e+10
krsfl=0.0
vrsfl=0.0

rotation_params=[
#   tlsfl    trsfl    tcsfl       ktsfl       vtsfl      plsfl    prsfl     pcsfl     kpsfl        vpsfl        tilsfl        tirsfl     ticsfl      ktisfl
[   0.319,   0.364,   0.342,   +4.00e+01,   -2.50e+05,   0.636,   0.681,   -0.101,   +4.00e+01,   +2.50e+05,    0.0e+3,       8.0e+3,    0.0e+3,     4.0],
[   0.339,   0.383,   0.361,   +4.00e+01,   +2.50e+05,   0.617,   0.661,   -0.098,   +4.00e+01,   -2.50e+05,    1.0e+3,       9.0e+3,    1.0e+3,     4.0],
[   0.334,   0.378,   0.356,   +4.00e+01,   +2.50e+05,   0.622,   0.666,   -0.065,   +4.00e+01,   -2.50e+05,    0.0e+3,       8.0e+3,    0.0e+3,     4.0],
[   0.337,   0.381,   0.359,   +4.00e+01,   -2.50e+05,   0.619,   0.663,   -0.043,   +4.00e+01,   +2.50e+05,    0.5e+3,       8.5e+3,    0.5e+3,     4.0],
[   0.352,   0.397,   0.375,   +4.00e+01,   -1.50e+05,   0.603,   0.648,   -0.020,   +4.00e+01,   +1.50e+05,    0.8e+3,       8.8e+3,    0.8e+3,     4.0],
[   0.372,   0.417,   0.395,   +4.00e+01,   +1.50e+05,   0.583,   0.628,   -0.022,   +4.00e+01,   -1.50e+05,    1.2e+3,       9.2e+3,    1.2e+3,     4.0],
[   0.383,   0.428,   0.406,   +4.00e+01,   -2.00e+05,   0.572,  0.617,   -0.062,   +4.00e+01,   +2.00e+05,    0.2e+3,       8.2e+3,    0.2e+3,     4.0],
[   0.374,   0.418,   0.396,   +4.00e+01,   +2.00e+05,   0.582,   0.626,   0.030,   +4.00e+01,   -2.00e+05,    0.0e+3,       8.0e+3,    0.0e+3,     4.0],
[   0.386,   0.430,   0.408,   +4.00e+01,   +2.50e+05,   0.570,   0.614,   -0.082,   +4.00e+01,   -2.50e+05,    0.1e+3,       8.1e+3,    0.1e+3,     4.0],
[   0.417,   0.461,   0.439,   +4.00e+01,   -2.50e+05,   0.539,   0.583,   -0.087,   +4.00e+01,   +2.50e+05,    1.3e+3,       9.3e+3,    1.3e+3,     4.0],
[   0.427,   0.472,   0.449,   +4.00e+01,   -2.00e+05,   0.528,   0.573,   -0.072,   +4.00e+01,   +2.00e+05,    0.0e+3,       8.0e+3,    0.0e+3,     4.0],
[   0.432,   0.477,   0.455,   +4.00e+01,   +2.00e+05,   0.523,   0.568,   -0.049,   +4.00e+01,   -2.00e+05,    0.0e+3,       8.0e+3,    0.0e+3,     4.0],
[   0.425,   0.470,   0.448,   +4.00e+01,   -1.00e+05,   0.530,   0.575,   -0.010,   +4.00e+01,   +1.00e+05,    0.7e+3,       8.7e+3,    0.7e+3,     4.0],
[   0.376,   0.420,   0.398,   +4.00e+01,   +1.00e+05,   0.580,   0.624,   0.174,   +4.00e+01,   -1.00e+05,    0.3e+3,       8.3e+3,    0.3e+3,     4.0],
[   0.364,   0.409,   0.387,   +4.00e+01,   -2.50e+05,   0.591,   0.636,   0.123,   +4.00e+01,   +2.50e+05,    0.2e+3,       8.2e+3,    0.2e+3,     4.0],
[   0.395,   0.440,   0.417,   +4.00e+01,   +2.50e+05,   0.560,   0.605,   0.099,   +4.00e+01,   -2.50e+05,    0.5e+3,       8.5e+3,    0.5e+3,     4.0],
[   0.395,   0.439,   0.417,   +4.00e+01,   +2.50e+05,   0.561,  0.605,   0.036,   +4.00e+01,   -2.50e+05,    1.0e+3,       9.0e+3,    1.0e+3,     4.0],
[   0.358,   0.403,   0.381,   +4.00e+01,   -2.50e+05,   0.597,   0.642,   0.035,   +4.00e+01,   +2.50e+05,    0.0e+3,       8.0e+3,    0.0e+3,     4.0]]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *


print("% Flow Options")
print("")
print("")
print("Separable Boundary Flow")
print("Radial Left")
print("No Bn Diffusion")
print("Time Uniform")
print("")
for idx in range(len(rotation_params)):
	print("Separable Boundary Flow")
	print("Radial Left")
	print("Theta Csc Theta")
	print("Theta Gauss  Theta")
	print("Theta DGauss Phi")
	print("Phi   DGauss Theta")
	print("Phi   Gauss  Phi")
	print("No Bn Diffusion")
	print("Time Cos")
	print("")
print("")

print("% Flow Data")
print("")
print("    nsflow")
print(str(len(rotation_params)+1).rjust(10))
print("     rlsfl        rrsfl       rcsfl    krsfl     vrsfl")
print(" 7.000e+10    7.020e+10   7.010e+10      0.0       0.0")
print("     tlsfl        trsfl       tcsfl    ktsfl     vtsfl")
print("       0.0          1.0         0.5      0.0       0.0")
print("     plsfl        prsfl       pcsfl    kpsfl     vpsfl")
print("      -1.0          1.0         1.0      0.0       0.0")
print("    tilsfl       tirsfl      ticsfl   ktisfl")
print("    0.0e+4       1.0e+6      0.0e+4      0.0")
for idx in range(len(rotation_params)):
	print("     rlsfl        rrsfl       rcsfl    krsfl     vrsfl")
	print("{:.3e}".format(rlsfl).rjust(10)+"{:.3e}".format(rrsfl).rjust(13)+"{:.3e}".format(rcsfl).rjust(12)+str(krsfl).rjust(9)+str(vrsfl).rjust(10))
	print("     tlsfl        trsfl       tcsfl    ktsfl     vtsfl")
	print("{:.3f}".format(rotation_params[idx][0]).rjust(10)+"{:.3f}".format(rotation_params[idx][1]).rjust(13)+"{:.3f}".format(rotation_params[idx][2]).rjust(12)+"{:.3f}".format(rotation_params[idx][3]).rjust(9)+"{:.1e}".format(rotation_params[idx][4]).rjust(10))
	print("     plsfl        prsfl       pcsfl    kpsfl     vpsfl")
	print("{:.3f}".format(rotation_params[idx][5]).rjust(10)+"{:.3f}".format(rotation_params[idx][6]).rjust(13)+"{:.3f}".format(rotation_params[idx][7]).rjust(12)+"{:.3f}".format(rotation_params[idx][8]).rjust(9)+"{:.1e}".format(rotation_params[idx][9]).rjust(10))
	print("{:.3e}".format(rotation_params[idx][10]).rjust(10)+"{:.3e}".format(rotation_params[idx][11]).rjust(13)+"{:.3e}".format(rotation_params[idx][12]).rjust(12)+"{:.3f}".format(rotation_params[idx][13]).rjust(9))



