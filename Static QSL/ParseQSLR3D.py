

filename="./PFLS/0078272/qsl.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
from mayavi import mlab


R_actual,theta,phi,Q=parse_QSL_Rbinfile(filename)

Q_grid=np.sign(Q)*np.log(abs(Q))
Q_grid[np.isinf(Q_grid)]=np.nan


X_grid=np.sin(theta)*np.cos(phi)
Y_grid=np.sin(theta)*np.sin(phi)
Z_grid=np.cos(theta)


mlab.figure("Q")
sphere_mesh=mlab.mesh(X_grid[::-1,:],Y_grid[::-1,:],Z_grid[:,:],scalars=Q_grid[:,:],colormap='RdBu',vmin=-10,vmax=10)
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True
mlab.colorbar(orientation="vertical")

mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)
mlab.show()




