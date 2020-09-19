"""
Split points along OCB into:
"regions" = high density points where surface drive will occur
"static"  = low density points which will no drive will occur
"""


filename="./PFLS/0044295/qslR1Fine.bin"
binary_outfile="./SurfaceLineData/SLine_StaticStart.bin"

num_per_region=200
region_phi=[[-14.0,-2.0]]		#Make sure regions do not overlap
region_theta=[[108.0,120.0]]

contour_bounds=[[-25.0,35.0],[90.0,120.0]]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from scipy.interpolate import interp1d


R,theta,phi,Q=parse_QSL_Rbinfile(filename)
Q_grid=np.log(abs(Q))
Q_sign=np.sign(Q)
Q_grid[np.isinf(Q_grid)]=0.0#np.nan		#Note: animate bug if nan in data

num_regions=min(len(region_phi),len(region_theta))

#Get contour where Q flips sign
fig2=plt.figure("Q map with all contours",figsize=(8,8))
ax2=fig2.gca()
plt.tight_layout()
color_plot2=plt.pcolormesh(phi,theta,Q_grid*Q_sign,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
contour1=plt.contour(phi,theta,Q_sign,[0.0],colors=["black"])


fig3=plt.figure("Raw contour",figsize=(8,8))
ax3=fig3.gca()
plt.tight_layout()
contour_found=False
for idx in range(len(contour1.allsegs[0])):
	if len(contour1.allsegs[0][idx])>10 and max(contour1.allsegs[0][idx][:,0])<contour_bounds[0][1]/180.0*np.pi and min(contour1.allsegs[0][idx][:,0])>contour_bounds[0][0]/180.0*np.pi and max(contour1.allsegs[0][idx][:,1])<contour_bounds[1][1]/180.0*np.pi and min(contour1.allsegs[0][idx][:,1])>contour_bounds[1][0]/180.0*np.pi:
		contour_found=True
		contour_raw=np.transpose(contour1.allsegs[0][idx])*180.0/np.pi
		plt.plot(contour_raw[0,:],contour_raw[1,:],color="black")
		plt.plot([contour_raw[0,0]],[contour_raw[1,0]],'o',color="orange")
		plt.plot([contour_raw[0,-1]],[contour_raw[1,-1]],'o',color="green")
if not contour_found:
	print("Warning! Contour not found")
	contour_raw=np.array([[],[]])

contour_raw=surfaceline_smooth_curve(contour_raw)
num_contour=len(contour_raw[0,:])

if num_regions>0:

	idx1=0
	region_idx0=surfaceline_inside_region(contour_raw[0,0],contour_raw[1,0],num_regions,region_phi,region_theta)
	region_idx1=region_idx0
	while idx1<num_contour-1 and region_idx0==region_idx1:
		idx1+=1
		region_idx1=surfaceline_inside_region(contour_raw[0,idx1],contour_raw[1,idx1],num_regions,region_phi,region_theta)

	#No regions in contour => only static closed contour
	if region_idx0==num_regions and idx1==num_contour-1:
		line_static=contour_raw
		surfaceline_print_cage(line_static,num_contour,'cage_static0')
		surfaceline_write_binary_cage(binary_outfile,[line_static],[])
		fig4=plt.figure("Final, single static contour",figsize=(8,8))
		plt.plot(line_static[0,:],line_static[1,:],color="orange")
		fig5=plt.figure("Final on Q map",figsize=(8,8))
		color_plot3=plt.pcolormesh(phi/np.pi*180.0,theta/np.pi*180.0,Q_grid*Q_sign,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
		plt.plot(line_static[0,:],line_static[1,:],color="orange")

	#Only one region in contour
	elif region_idx0!=num_regions and idx1==num_contour-1:
		x_param=np.linspace(1,num_per_region,num=num_contour)
		x_out=np.linspace(1,num_per_region,num=num_per_region)
		phi_func=interp1d(x_param,contour_raw[0,:],kind='linear')
		theta_func=interp1d(x_param,contour_raw[1,:],kind='linear')
		line_region=np.zeros((2,num_per_region))
		line_region[0,:]=phi_func(x_out)
		line_region[1,:]=theta_func(x_out)
		line_region=surfaceline_smooth_Npt(line_region,3)
		surfaceline_print_cage(line_region,num_per_region,'cage0')
		surfaceline_write_binary_cage(binary_outfile,[],[line_region])
		fig4=plt.figure("Final, single region",figsize=(8,8))
		plt.plot(line_region[0,:],line_region[1,:],color="black")
		fig5=plt.figure("Final on Q map",figsize=(8,8))
		color_plot3=plt.pcolormesh(phi/np.pi*180.0,theta/np.pi*180.0,Q_grid*Q_sign,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
		plt.plot(line_region[0,:],line_region[1,:],color="black")

	#One region + (other regions and/or static regions)
	else:
		line_region=[]
		line_static=[]
		contour_shifted=np.concatenate((contour_raw[:,idx1:-1],contour_raw[:,:idx1+1]),axis=1)
		region_idx0=surfaceline_inside_region(contour_shifted[0,0],contour_shifted[1,0],num_regions,region_phi,region_theta)
		idx0=0
		idx1=0
		while idx1<num_contour-1:
			idx1+=1
			region_idx1=surfaceline_inside_region(contour_shifted[0,idx1],contour_shifted[1,idx1],num_regions,region_phi,region_theta)
			if region_idx1!=region_idx0 or idx1==num_contour-1:
				#Static
				if region_idx0==num_regions:
					line_static.append(contour_shifted[:,idx0:idx1+1])
				#Region
				else:
					x_param=np.linspace(1,num_per_region,num=idx1+1-idx0)
					x_out=np.linspace(1,num_per_region,num=num_per_region)
					phi_func=interp1d(x_param,contour_shifted[0,idx0:idx1+1],kind='linear')
					theta_func=interp1d(x_param,contour_shifted[1,idx0:idx1+1],kind='linear')
					line_region_temp=np.zeros((2,num_per_region))
					line_region_temp[0,:]=phi_func(x_out)
					line_region_temp[1,:]=theta_func(x_out)
					line_region.append(surfaceline_smooth_Npt(line_region_temp,3))
				region_idx0=region_idx1
				idx0=idx1

		for idx in range(len(line_region)):
			surfaceline_print_cage(line_region[idx],num_per_region,'cage'+str(idx))
		for idx in range(len(line_static)):
			surfaceline_print_cage(line_static[idx],len(line_static[idx][0,:]),'cage_static'+str(idx))
		surfaceline_write_binary_cage(binary_outfile,line_static,line_region)
		fig4=plt.figure("Final, all regions (black)+statics (orange)",figsize=(8,8))
		for idx in range(len(line_region)):
			plt.plot(line_region[idx][0,:],line_region[idx][1,:],color="black")
		for idx in range(len(line_static)):
			plt.plot(line_static[idx][0,:],line_static[idx][1,:],color="orange")
		fig5=plt.figure("Final on Q map",figsize=(8,8))
		color_plot3=plt.pcolormesh(phi/np.pi*180.0,theta/np.pi*180.0,Q_grid*Q_sign,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
		for idx in range(len(line_region)):
			plt.plot(line_region[idx][0,:],line_region[idx][1,:],color="black")
		for idx in range(len(line_static)):
			plt.plot(line_static[idx][0,:],line_static[idx][1,:],color="orange")

#No regions defined => only static closed contour
else:
	line_static=contour_raw
	surfaceline_print_cage(line_static,num_contour,'cage_static0')
	surfaceline_write_binary_cage(binary_outfile,[line_static],[])
	fig4=plt.figure("Final, single static contour",figsize=(8,8))
	plt.plot(line_static[0,:],line_static[1,:],color="orange")
	fig5=plt.figure("Final on Q map",figsize=(8,8))
	color_plot3=plt.pcolormesh(phi/np.pi*180.0,theta/np.pi*180.0,Q_grid*Q_sign,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
	plt.plot(line_static[0,:],line_static[1,:],color="orange")

plt.show()




