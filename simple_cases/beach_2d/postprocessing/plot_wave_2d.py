### PLOT WAVE for 2D CASE ###

# import necessary modules
import numpy as np               
import matplotlib.pyplot as plt
import os
# write your OWN PC folder path for fdir
# Remember that we use for Mac & Linux machines '/', while on windows '\'
fdir = r'C:\Users\User\Documents\USACE_WORK\Funwave_Seminar\results\beach_2d\work\output'

# upload eta file
eta=np.loadtxt(os.path.join(fdir,'eta_00001'))

# define plot location
n,m = np.shape(eta)   
dx = 2.0
dy = 2.0

x = np.asarray([float(xa)*dx for xa in range(m)])
y = np.asarray([float(ya)*dy for ya in range(n)])

# define sponge and wavemaker location
x_sponge = [0,100,100,0,0]
y_sponge = [0,0,y[len(y)-1],y[len(y)-1],0]

x_wavemaker = [150, 160, 160, 150, 150]
y_wavemaker = [0, 0, y[len(y)-1],y[len(y)-1],0]


nfile = [10, 50]    # range of eta files you want to plot
min = ['20','100']  # time  you want to plot

# figure size option 
wid=10    # width
length=5 # length


# plot figure
fig = plt.figure(figsize=(wid,length),dpi=200)

for num in range(len(nfile)):
    fnum= '%.5d' % nfile[num]
    eta = np.loadtxt(os.path.join(fdir,'eta_'+fnum))
    mask = np.loadtxt(os.path.join(fdir,'mask_'+fnum))

    eta_masked = np.ma.masked_where(mask==0,eta) # do nt plot where mask = 0

    ax = fig.add_subplot(1,len(nfile),num+1)
    fig.subplots_adjust(hspace=1,wspace=.25)
    plt.pcolor(x, y, eta_masked,cmap='coolwarm')

    title = 'Time = '+min[num]+ ' sec'
    plt.title(title)
    plt.hold(True)

    # plot sponge and wavemaker
    plt.plot(x_sponge,y_sponge,'g--',linewidth=3)
    plt.text(50,500,'Sponge',color='g',rotation=90)
    
    plt.plot(x_wavemaker,y_wavemaker,'k-',linewidth=3)
    plt.text(180,700,'Wavemaker',color='k',rotation=90)
    
    if num == 0:
        plt.ylabel('Y (m)')
        plt.xlabel('X (m)')
    else:
        plt.xlabel('X (m)')
        cbar=plt.colorbar()
        cbar.set_label(r'$\eta$'+' (m)', rotation=90)
        

# save figure  
fig.savefig('eta_2d_wave.png', dpi=fig.dpi)
