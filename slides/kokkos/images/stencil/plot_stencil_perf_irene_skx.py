import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,92,128,160,192,224,256])

# test_stencil_3d_flat
v0=np.array([40.6474,55.5205,53.1085,56.1479,58.2191,58.4607,58.7187,57.4828,58.1501])

# test_stencil_3d_flat_vector without views
v1=np.array([194.174,313.152,371.928,381.576,398.154,391.202,384.016,318.343,265.672])

# test_stencil_3d_flat_vector with    views
v2=np.array([170.467,439.805,676.623,969.515,1163.12,1170.32,673.194,344.669,267.759])

# test_stencil_3d_range
v3=np.array([138.086,196.601,213.433,221.592,224.149,219.406,217.446,207.082,204.643])

# test_stencil_3d_range_vector
v4=np.array([179.146,395.297,558.926,1034.28,479.188,941.605,670.391,341.184,267.173])

# test_stencil_3d_range_vector2
v5=np.array([135.117,284.356,314.146,570.801,396.165,375.991,246.402,185.752,167.224])

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
