from scripts import testqueue as tqueue

import numpy as np
from time import clock


dx = 0.01
minX, maxX = 0.0, 10.0
minY, maxY = 0.0, 10.0

nx = int((maxX - minX)/dx)+1
ny = int((maxY - minY)/dx)+1

xcoords = np.linspace(minX, maxX, nx)
ycoords = np.linspace(minY, maxY, ny)

X, Y = np.meshgrid(xcoords, ycoords)

coords = np.vstack([X.ravel(), Y.ravel()])


height  = 1.0 - (np.sin(0.5*X*np.pi)**2 * np.sin(0.5*Y*np.pi)**2)
#height  += np.random.random(height.shape) * 0.001
height  += 0.25 * (X + 0.5*Y)

dem = height.reshape((ny,nx))
# dem = np.load('dem2fill.npy')
# dem = np.random.rand(10,10)*10.
# dem[:,:] = 5.
# dem[0,:] = 10.
# dem[-1,:] = 10.
# dem[:,0] = 10.
# dem[:,-1] = 10.
# dem[4:6,4:6] = 3

np.save('egg',dem)
eps = 1.e-6

t0 = clock()
filldem,labels = tqueue.pitfillingzhou(eps,dem)
fillpts = tqueue.pitfillingpts(labels.max(),4)
print('Filling zhou 2016 (%0.02f seconds)'% (clock() - t0))
np.save('fillPts',fillpts)
np.save('fillZhou',filldem)
np.save('labelZhou',labels)
#
# t0 = clock()
# filldem = tqueue.pitfillingwei(eps,dem)
# print('Filling wei 2018 (%0.02f seconds)'% (clock() - t0))
# print (filldem-dem).max()
# np.save('fillWei',filldem)
# t0 = clock()
# if eps>0.:
#     filldem,label = tqueue.pitfillingbarneseps(eps,dem)
# else:
#     filldem,label = tqueue.pitfillingbarnes(dem)
# print('Filling barnes 2014 (%0.02f seconds)'% (clock() - t0))
# #
# np.save('labelBarnes',label)
# print (filldem-dem).max()
# np.save('fillBarnes',filldem)

# np.save('label',labels)
