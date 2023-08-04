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

# do not use this code as it is not complete!


import numpy as np
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps
from scipy import ndimage
import os
import time as tm
from scipy.interpolate import CubicSpline

path = r"E:\SRFP\Python and Jupyter\Data Cube"
data_files = sorted(os.listdir(path))

time = []
for f in data_files:    
    time.append((int(f[10])-6)*86400+int(f[12:14])*3600+int(f[14:16])*60+int(f[16:18]))
    
#############################################################################################
slit_name = 'new_slit_7'
slit_x = [2178,2226,2284,2350,2396,2464,2521,2566,2648]
slit_y = [1077,1090,1134,1169,1196,1223,1262,1300,1325]
x_data = np.array(slit_x)
y_data = np.array(slit_y)

# Create the cubic spline interpolation function
cubic_interpolation = CubicSpline(x_data, y_data)

cubic_interpolation = CubicSpline(x_data, y_data)
x_values = np.arange(x_data.min(),x_data.max(),1) #,dtype='int')
y_values = cubic_interpolation(x_values)
# y_values = np.array(y_values,dtype='int')
# here pixel values can be fraction as well
local_slope = []
for i in range(len(x_values)-1):
    local_slope.append((y_values[i+1]-y_values[i])/(x_values[i+1]-x_values[i]))
local_slope.append(local_slope[-1])
#local_normal = [-1/m for m in local_slope]
def x_values_add(n): # finds the pixel which is at n length far from the slit and is perpendicular to the slit
    x_arr = []
    for xx,mm in zip(x_values,local_slope):
        x_arr.append((xx-n*mm/(mm**2+1)**0.5))
    return x_arr
def y_values_add(n):
    y_arr = []
    for yy,mm in zip(y_values,local_slope):
        y_arr.append((yy+n/(mm**2+1)**0.5))
    return y_arr

############################################################################

# to extract 4 pixel width slit data
def slit_data(files):
    with open(path+'\\'+files) as f:
        #print(str(f.name[40:])+' : '+str(tm.time()-start))
        A = f.readlines()
    B0 = []
    for line in A:
        B0.append([])
        numbers = line.split()
        for number in numbers:
            B0[-1].append(float(number))
    # function that retrieves intensity value of a fractional pixel
    def pixel_value(x,y):
        xf = int(np.floor(x))
        yf = int(np.floor(y))
        xc = int(np.ceil(x))
        yc = int(np.ceil(y))
        if (xc==xf) and (yc==yf):
            return B0[yc-ymin][xc-xmin]
        elif (xc==xf):
            return B0[yc-ymin][xc-xmin]*(y-yf)/(yc-yf)+B0[yf-ymin][xc-xmin]*(yc-y)/(yc-yf)
        elif (yc==yf):
            return B0[yc-ymin][xc-xmin]*(x-xf)/(xc-xf)+B0[yc-ymin][xf-xmin]*(xc-x)/(xc-xf)
        else:
            return B0[yc-ymin][xc-xmin]*(y-yf)*(x-xf)/((yc-yf)*(xc-xf))+B0[yf-ymin][xc-xmin]*(yc-y)*(x-xf)/((yc-yf)*(xc-xf))+B0[yc-ymin][xf-xmin]*(y-yf)*(xc-x)/((yc-yf)*(xc-xf))+B0[yf-ymin][xf-xmin]*(yc-y)*(xc-x)/((yc-yf)*(xc-xf))

    B = []
    for i in range(len(x_values)):
        #print(x,y)
        ii = pixel_value(x_values[i],y_values[i])+pixel_value(x_values_add(1)[i],y_values_add(1)[i])+pixel_value(x_values_add(-1)[i],y_values_add(-1)[i])+pixel_value(x_values_add(2)[i],y_values_add(2)[i])
        B.append(ii/4)
    return B

################################################################################

Intensity = []
start = tm.time()
for files in data_files:
    print(files+' : '+str(tm.time()-start))
    Intensity.append(slit_data(files))
    
plt.imshow(np.transpose(Intensity),cmap = 'sdoaia171')
plt.savefig(slit_name+'xtmap.png')