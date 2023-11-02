#!/usr/bin/env python

# code adapted from http://matplotlib.org/examples/animation/dynamic_image.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

def read_csv(filename):
    data = np.loadtxt(filename, delimiter=',')
    return data

fig = plt.figure()

index = 0
prefix = "wave"
filename = prefix+"_{0:03d}.csv".format(index)
data = read_csv(filename)
im = plt.imshow(data, cmap=plt.get_cmap('viridis'), animated=True)

def updatefig(*args):

    global index
    index = index+1

    filename = prefix+"_{0:03d}.csv".format(index)

    if not os.path.isfile(filename):
        index = 0
        filename = prefix+"_{0:03d}.csv".format(index)

    data = read_csv(filename)

    im.set_array(data)
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)

plt.show()
