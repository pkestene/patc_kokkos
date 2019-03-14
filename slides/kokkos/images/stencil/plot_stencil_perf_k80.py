import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,92,128,160,192,224,256])

# test_stencil_3d_flat
v0=np.array([66.4126,108.116,127.477,135.095,140.738,141.193,141.341,141.743,157.292])

# test_stencil_3d_flat_vector without views
v1=np.array([29.6409,64.3669,103.893,181.612,261.148,255.709,255.378,280.793,265.961])

# test_stencil_3d_flat_vector with    views
v2=np.array([37.2533,84.1228,130.206,223.113,284.583,281.542,284.448,288.226,303.383])

# test_stencil_3d_range
v3=np.array([75.1473,124.079,134.099,137.047,150.487,151.771,152.334,158.308,163.668])

# test_stencil_3d_range_vector
v4=np.array([37.1099,87.3658,125.713,132.93,177.567,203.226,198.414,203.97,213.25])

# test_stencil_3d_range_vector2
v5=np.array([33.4572,44.0755,75.2012,92.249,146.407,85.488,102.245,118.502,134.375])

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
