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

# to improve speed and to avoid memory overflow in cache each video contains only 100 image frames

import os
import time
start = time.time()

# Iterate over the files in the folder
figure_files = []
for filename in os.listdir(r"E:\SRFP\Python and Jupyter\Images"):
    if os.path.isfile(os.path.join(r"E:\SRFP\Python and Jupyter\Images", filename)):
        figure_files.append(filename)

figure_filenames = sorted(figure_files)
#for i in figure_filenames:
#    print(i)


import matplotlib.pyplot as plt
import matplotlib.animation as animation

import os
import shutil
import subprocess

def clear_imagemagick_cache():
    # Get the cache directory for ImageMagick
    cache_dir = subprocess.check_output(['magick', 'convert', '-list', 'cache']).decode().strip()
    
    # Check if the cache directory exists
    if os.path.exists(cache_dir):
        # List all files in the cache directory and delete them
        for root, _, files in os.walk(cache_dir):
            for file in files:
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")
    
        print("ImageMagick cache cleared successfully.")
    else:
        print("ImageMagick cache directory not found.")


# Create a figure and axes for the animation

'''xmin = 2000
xmax = 3700
ymin = 400
ymax = 2000
dpi = 250
figx = (xmax-xmin)/dpi
figy = (ymax-ymin)/dpi
'''
dpi = 409.575
fig, ax = plt.subplots(figsize=(2785/dpi, 2621/dpi))

for i in range(35):
    if not os.path.isfile(r"E:\SRFP\Python and Jupyter\Movies\\"+"oscillation_movie"+str(i)+".mp4"):
        # Function to update the plot with each frame of the animation
        def update(frame):
            ax.clear()
            # Load the figure for the current frame
            img = plt.imread(r"E:\SRFP\Python and Jupyter\Images\\"+figure_filenames[frame+100*i])
            #image_shape = image.shape[:2]
            ax.imshow(img,)
            
            print(str(frame+100*i)+' : time : '+str(time.time()-start)+' : '+ str(figure_filenames[frame+100*i]))
            #ax.set_title(f'Frame {frame}')
            #print(img.shape[:2])
            ax.set_axis_off()
            
            return ax


        # Create the animation using the update function and figure filenames
        ani = animation.FuncAnimation(fig, update, frames=100, interval=100)

        # # Save the animation to a file (e.g., GIF)
        ani.save(r"E:\SRFP\Python and Jupyter\Movies\\"+"oscillation_movie"+str(i)+".mp4", writer='imagemagick',dpi = dpi)
        clear_imagemagick_cache()
