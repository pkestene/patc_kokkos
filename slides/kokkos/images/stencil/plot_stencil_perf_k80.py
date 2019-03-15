import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,92,128,160,192,224,256,320,384,448,512,640,768])

# test_stencil_3d_flat
v0=np.array([67.3823,97.3598,108.883,114.836,125.416,143.3,173.187,181.124,181.019,179.88,179.355,179.011,178.958,178.639,178.281])

# test_stencil_3d_flat_vector without views
v1=np.array([33.53,72.5903,111.455,196.83,263.321,278.502,270.189,288.481,288.865,280.08,283.08,287.397,281.791,236.851,103.041])

# test_stencil_3d_flat_vector with    views
v2=np.array([36.4741,79.1084,111.399,202.249,342.694,347.427,332.803,346.803,346.371,364.279,361.609,363.749,365.302,362.064,202.509])

# test_stencil_3d_range
v3=np.array([85.8318,129.07,138.823,140.06,152.699,153.767,154.093,153.628,148.906,145.162,145.466,144.93,145.624,145.529,145.13])

# test_stencil_3d_range_vector
v4=np.array([26.8551,61.1739,92.5664,92.0654,124.113,141.287,136.345,138.627,144.92,138.407,137.408,136.513,135.789,132.401,110.693])

# test_stencil_3d_range_vector2
v5=np.array([35.3318,38.3573,65.0245,81.8839,128.156,75.7353,90.7364,105.525,120.237,99.8524,119.666,104.453,119.664,119.286,93.7895])

plt.plot(size,v0, label='# test_stencil_3d_flat')
plt.plot(size,v1, label='# test_stencil_3d_flat_vector without views')
plt.plot(size,v2, label='# test_stencil_3d_flat_vector with    views')
plt.plot(size,v3, label='# test_stencil_3d_range')
plt.plot(size,v4, label='# test_stencil_3d_range_vector')
plt.plot(size,v5, label='# test_stencil_3d_range_vector2')
plt.grid(True)
plt.title('3d Heat kernel performance')
plt.xlabel('N - linear size')
plt.ylabel(r'Bandwidth (GBytes/s)')
plt.legend()
plt.show()
