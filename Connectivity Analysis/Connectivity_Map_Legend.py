
"""
Colormap for reconnected or not reconnected connectivities; use these colors for reference

O   Open
C   Closed

Connectivity types: partial = full/2 (rounded down)
0  O (final state)   No rec         Blue            #ABD2E5     (171/255,210/255,229/255)
1  C (final state)   No rec         Red             #F4A683     (244/255,166/255,131/255)
2  O (final state)   Reconnected    Grey            #666666     (102/255,102/255,102/255)
3  C (final state)   Reconnected    Brass           #8F9100     (143/255,145/255,0/255)
4  Not connected to solar surface   Orange          #FF8000     (255/255,128/255,0)

Connectivity types: full (not used)
0  O->O No rec		#7FC97F
1  C->O No rec		#BEAED4
2  O->C No rec		#FDC086
3  C->C No rec		#FFFF99
4  O->O Reconnected	#386CB0
5  C->O Reconnected	#F0027F
6  O->C Reconnected	#BF5B17
7  C->C Reconnected	#666666
8  Not connected to solar surface	#FF8000

Outdated, published previously, Accent_new:
Connectivity types: partial = full/2 (rounded down)
0  O (final state)   No rec         Light green     #7FC97F     (127/255,201/255,127/255)
1  C (final state)   No rec         Brown           #994C00     (153/255,76/255,0)
2  O (final state)   Reconnected    Grey            #666666     (102/255,102/255,102/255)
3  C (final state)   Reconnected    Pink            #F0027F     (240/255,2/255,127/255)
4  Not connected to solar surface   Orange          #FF8000     (255/255,128/255,0)
"""


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


connection_colors=[(171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)]
connection_cmap=matplotlib.colors.ListedColormap(((171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)))


#Outdated, published previously
#Flip stark colors in 'Accent' colormap
Accent_colors=matplotlib.cm.get_cmap('Accent').colors
color_brown=(153/255,76/255,0) #Instead of Accent_colors[3]
color_orange=(255/255,128/255,0) #Instead of Accent_colors[3]
Accent_colors_new=[Accent_colors[0],color_brown,Accent_colors[7],Accent_colors[5],color_orange]
Accent_new=matplotlib.colors.ListedColormap((Accent_colors[0],color_brown,Accent_colors[7],Accent_colors[5],color_orange))


fig1=plt.figure("Connection legend",figsize=(20,15))


#Key
ax1_4=fig1.add_axes([0.13, 0.12, 0.15, 0.2])
color_plot1_4=plt.pcolormesh(np.array([0,1,2]),np.array([0,1,2]),np.array([[0.5,1.5,1.5],[2.5,3.5,3.5],[2.5,3.5,3.5]]),cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
ax1_4.set_xticks([])
ax1_4.set_yticks([])
plt.text(-2.05,0.6,"Connectivity",fontsize=42)
plt.text(-2.05,0.2,"retained",fontsize=42)
plt.text(-2.05,1.4,"Reconnected",fontsize=42)

plt.text(0.0,2.2,r"Open",fontsize=42)
#plt.text(1.2,2.2,r"C$\rightarrow$O",fontsize=20)
plt.text(1.1,2.2,r"Closed",fontsize=42)
#plt.text(3.2,2.2,r"C$\rightarrow$C",fontsize=20)
plt.text(-0.35,2.7,r"Field line end state",fontsize=42)


ax1_5=fig1.add_axes([0.45, 0.22, 0.075, 0.1])
color_plot1_5=plt.pcolormesh(np.array([0,1]),np.array([0,1]),np.array([[4.5,4.5],[4.5,4.5]]),cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
ax1_5.set_xticks([])
ax1_5.set_yticks([])
plt.text(-0.5,1.2,r"Unconnected",fontsize=42)

plt.savefig("Connectivity_Legend.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.2)



