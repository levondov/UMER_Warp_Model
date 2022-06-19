import glob
import scipy.io as sio
from PIL import Image
import numpy as np
import os

files = sorted(glob.glob('*.tif'))

photos = []
zpos = []

# grab photos 1 at a time and their locations
for photo in files:
	photos.append(np.array(Image.open(photo)))
	zpos.append(int(photo[-11:-4]))

# grab calculated envelopes from equation
env_dt = np.genfromtxt('env.txt');
pos_dt = np.genfromtxt('pos.txt');

sio.savemat('data2.mat',{'photos': photos, 'zpos': zpos, 'envs': env_dt, 'pos': pos_dt})

os.system('rm *.tif')
