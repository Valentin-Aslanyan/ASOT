

import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


flux_integral=np.array([3.1962443672447566e+19, 3.73786660915174e+19, 4.4703635907033235e+19, 4.2929237532945424e+19, 4.124097829476898e+19, 3.6546849892660736e+19, 2.5017664821596987e+19, 2.187182702162757e+19, 8.596540892811663e+18])
time=np.array(range(len(flux_integral)))*1000.0


plt.figure()
plt.plot(time,flux_integral)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\Phi$ [G cm$^2$]",fontsize=20)
plt.xlabel(r"$t$ [s]",fontsize=20)


plt.show()
