import numpy as np
import matplotlib.pyplot as plt
import imageio as iio

def f(theta, title=' '):   #theta is the angles data file written with the write function from the Film class
	"""Function to plot angles. A black square corresponds to an angle of pi, a white square corresponds
	to an angle of 0. An angle of theta is the same as the angle -theta[2pi]."""
	mask=(theta>np.pi)  #We look for angles greater than pi
	theta[mask]-=np.pi #In this case, we compute the opposite angle [2pi]
	theta*=180/np.pi  #We convert the angle in degrees
	#We plot the data using a gray colorbar :
	fig=plt.figure()
	im=plt.imshow(theta,cmap='gray_r',vmin=0,vmax=180)
	c=fig.colorbar(im)
	c.set_label(r'$\theta$ [°]')
	plt.title(title)
	plt.show()

def gif(files) :   #files is the list containing all the data files that we want to put in the GIF
	"""Function to create a GIF with many film plots"""
	i=0
	for file in files :
		mask=(file>np.pi)
		file[mask]-=np.pi
		file*=180/np.pi
		plt.imshow(file,cmap='gray_r',vmin=0,vmax=180)
		plt.savefig(f'i'+str(i)+'.png')
		i+=1
	frames=np.stack([iio.imread(f'i{i}.png') for i in range(len(files))],axis=0)
	iio.mimwrite('aimantation_0.gif',frames)
