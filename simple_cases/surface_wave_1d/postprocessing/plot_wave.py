### PLOT WAVE 1D CASE ###

# import necessary modules
import numpy as np               
import matplotlib.pyplot as plt
import os 

# import cv2 to create a video out of png output images
import cv2
# check if OpenCV is installed
# print("OpenCV version:", cv2.__version__)

# To execute coode asynchronously
import asyncio

base_dir = "/Users/polselli/Code/funwave/FUNWAVE-TVD/simple_cases/surface_wave_1d/"
output_files_dir = rf"{base_dir}output_files_reg/" 
output_dir = rf"{base_dir}postprocessing/output/py/"

m = 1024
dx = 1.0
SLP = 0.05
Xslp = 800.0

# Bathymetry computation
x = np.asarray([float(xa)*dx for xa in range(m)])
dep = np.zeros((m,m))+10.0

for i in range(len(x)):
     if i < Xslp:
         dep[i,0]=dep[i,0]
     else:
         dep[i,0]=10-(x[i]-Xslp)*SLP

DEP = dep*-1

# Locate Wavemaker and sponge
wd = 10
x_wm = [250-wd,250+wd,250+wd,250-wd,250-wd]
x_sponge = [0,180,180,0,0]
yy = [-10,-10,10,10,-10]


# Figure size options
wid = 8     #width
length = 4  #length


# plot figure
fig = plt.figure(figsize=(wid,length),dpi=600)

# Specify the specific string to match
prefix = "eta_"

output_files_x_imgs = [of for of in sorted(os.listdir(output_files_dir), key=str.casefold) if of.startswith(prefix)]

imgs_paths = []

for index, of in enumerate(output_files_x_imgs): 
    eta = np.loadtxt(os.path.join(output_files_dir, of))

    ax = fig.add_subplot(1,1,1)

    fig.subplots_adjust(hspace=1,wspace=.5)
    plt.plot(x,DEP[:,0],'k',x,eta[0,:],'b',linewidth=2)
    plt.plot(x_wm,yy,'r')
    plt.plot(x_sponge,yy,'k')

    ##plt.hold(True)
    plt.grid()

    plt.text(x_wm[1],0.6,'Wavemaker',color='red', fontsize=12,style='oblique')
    plt.text(x_sponge[0]+20,0.6,'Sponge Layer',color='black',
            fontsize=12,style='oblique')

    ax.axis([0, 1024, -1, 1])
    plt.minorticks_on()
    plt.xlabel('X (m)')
    plt.ylabel(r'$\eta$'+' (m)')

    imgs_paths.append(output_dir + 'eta_1d_wave_' + of.lstrip(prefix) + '.png')
    # save figure
    fig.savefig(imgs_paths[index], dpi=fig.dpi) # save figure

#----------------------- VIDEO -------------------
# -- Build video out of the output pictures


# Define the codec and create a VideoWriter object
output_video = 'output_video.mp4'

async def save_video(par_imgs_paths, par_output_dir, par_output_video):
    # Define video parameters
    height, width = None, None  # Initialize dimensions
    out = None  # Initialize VideoWriter object

    # Check if the file exists
    file_path = par_output_dir + par_output_video
    if os.path.exists(file_path):
        # Delete the file
        os.remove(file_path)
        print(f"File {file_path} has been deleted.")
    else:
        print(f"File {file_path} does not exist.")

    try:
        for img_path in par_imgs_paths:
            frame = cv2.imread(img_path)
            if frame is None:
                print(f"Skipping image file: {img_path}")
                continue
            
            # Initialize VideoWriter lazily to get the frame size
            if height is None or width is None:
                height, width, _ = frame.shape
                size = (width, height)
                fourcc = cv2.VideoWriter_fourcc(*'MJPG')
                out = cv2.VideoWriter(par_output_dir + par_output_video, fourcc, 1.0, size)

            # Perform operations on frame
            frame = cv2.cvtColor(frame, cv2.COLOR_RGB2BGR)

            # Write frame to video
            out.write(frame)

        if out is not None:
            # Release the VideoWriter object
            out.release()
            print(f'Video saved as {par_output_video}')
        else:
            print("No frames to write, video not saved.")

    except Exception as e:
        print(f"Error occurred: {e}")

# Run the save_video coroutine asynchronously
asyncio.run(save_video(imgs_paths, output_dir, output_video))
