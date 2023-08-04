# These codes are developed by Panchajanya Dey, Indian Institute of Science, Bangalore
# as a part of the SRFP project prominence oscillation in solar corona funded by Indian Academy of Science (2023)
# this folder contains the codes are for downloading the necessary fits files, modifying them by multi-gaussian image processing
# and creating an animation that helps to identify the location of the prominence oscillation
# the next part is to determine the pixel coordinates and creating a 4-pixel wide slit and get the corresponding XT map
# which is also covered in this project
# finally handpicking the prominence oscillation from the XT map and fitting it to a damped harmonic oscillation is done
# the final result is the period, damping time and amplitude determination for a LAL prominence oscillation
# I thank Prof. Vaibhav Pant (ARIES) and his group member Upasna Baweja (ARIES) for helping me throughout the project
# and clearing my doubts
# for comments, suggestions and clarity regarding the code contact at the following email address
# panchajanyad@iisc.ac.in

from genericpath import isfile
import numpy as np
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps
from scipy import ndimage
from astropy.io import fits
import time
from scipy.interpolate import CubicSpline
import os

start = time.time()

path = r"E:\SRFP\Python and Jupyter\FITS Files"
target = r"E:\SRFP\Python and Jupyter\Data Cube"
fits_files = sorted(os.listdir(path))
#print(fits_files)

# creating a gaussian kernel
def gaussian_kernel(sigma):
    kernel_size = int(2 * np.ceil(2 * sigma) + 1)
    offset = int(kernel_size // 2)
    kernel = [[(1 / (2 * np.pi * sigma ** 2)) * np.exp(-((x - offset) ** 2 + (y - offset) ** 2) / (2 * sigma ** 2)) for x in range(kernel_size)] for y in range(kernel_size)]
    #normalise the kernel
    sum_kernel = sum(sum(kernel,[]))
    kernel = [[element / sum_kernel for element in row] for row in kernel]
    return kernel

def mgn(data,sigma,k):
    blur_data = ndimage.convolve(data, gaussian_kernel(sigma))
    sharp_data = data-blur_data
    squared_sharp = np.square(sharp_data)
    data_standard_deviation = np.sqrt(ndimage.convolve(squared_sharp,gaussian_kernel(sigma)))
    mask = data_standard_deviation == data_standard_deviation.min()
    data_standard_deviation[mask]=10**-6
    mgn_data_0 = sharp_data/data_standard_deviation
    mgn_data = np.arctan(k*mgn_data_0)
    return mgn_data

def global_gamma_trans(data,gamma):
    a0 = data.min()
    a1 = data.max()
    data_1 = ((data-a0)/(a1-a0))**(1/gamma)
    return data_1

def fits_to_cube(file_name):
    #print(file_name[:-5])
    with fits.open(path+'\\'+file_name) as f:
        primary_data = f[0].data
    data = primary_data.copy()
    
    xmin = 2000
    xmax = 3700
    ymin = 400
    ymax = 2000
    
    data_cropped = data[ymin:ymax,xmin:xmax]
    A = 0.7*global_gamma_trans(data_cropped,3.2)+(1-0.7)/5*(0.91*mgn(data_cropped,1.25,0.7)+0.976*mgn(data_cropped,2.5,0.7)+0.994*mgn(data_cropped,5,0.7)+0.998*mgn(data_cropped,10,0.7)+0.999*mgn(data_cropped,20,0.7))#+0.999*mgn(data_cropped,40,0.7)
    with open(target+'\\'+file_name[:-5]+'.txt','w') as file:
        for arrays in A:
            for element in arrays:
                file.write(str(element))
                file.write('\t')
            file.write('\n')

n = 0
'''for file_name in fits_files:
    if os.path.isfile(target+'\\'+file_name[:-5]+'.txt'):
        print(str(file_name)+' : '+ str(n)+' : '+str(time.time()-start))
        n += 1
print(str(n) + ' files already exists\nSkipping Existing Files')'''

for file_name in fits_files[n:]:
    if os.path.isfile(target+'\\'+file_name[:-5]+'.txt'):
        print(str(file_name)+' : '+ str(n)+' : '+str(time.time()-start))
        n += 1
    else:
        fits_to_cube(file_name)
        print(str(file_name)+' : '+ str(n)+' : '+str(time.time()-start))
        n += 1