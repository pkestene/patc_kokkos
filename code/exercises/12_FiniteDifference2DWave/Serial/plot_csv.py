#!/usr/bin/env python

# code adapted from http://matplotlib.org/examples/animation/dynamic_image.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_csv(filename):
    data = np.loadtxt(filename, delimiter=',')
    return data

fig, ax = plt.subplots()

index = 0
filename = "Serial_{0:03d}.csv".format(index)
data = read_csv(filename)
nx,ny = data.shape
line, = ax.plot(data[nx/2,:])

def animate(i):

    filename = "Serial_{0:03d}.csv".format(i)
    data = read_csv(filename)

    line.set_ydata(data[nx/2,:])
    return line,

def init():
    
    filename = "Serial_{0:03d}.csv".format(0)
    data = read_csv(filename)

    line.set_ydata(data[nx/2,:])
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(1, 99), init_func=init, interval=25, blit=True)

plt.show()
