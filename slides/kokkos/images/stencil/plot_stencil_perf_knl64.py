import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,96,128,160,192,224,256,320,384,448,512,640])

# test_stencil_3d_flat
v0=np.array([13.7917,15.1331,15.5294,15.4053,15.4355,15.1985,15.2929,15.1704,15.2632,15.0697,15.1805,15.1022,15.1871,15.039])

# test_stencil_3d_flat_vector without views
v1=np.array([60.6125,111.605,143.075,167.752,173.991,171.932,170.261,170.549,168.96,177.532,175.474,175.671,177.891,178.773])

# test_stencil_3d_flat_vector with    views
v2=np.array([66.5161,214.16,351.369,621.421,779.558,686.778,578.014,629.346,567.316,629.876,620.867,638.653,606.513,634.713])

# test_stencil_3d_range
v3=np.array([40.1905,65.414,75.5865,81.3421,82.7137,82.076,82.1537,82.6232,80.3771,83.9939,83.0511,84.5146,82.5698,83.0552])

# test_stencil_3d_range_vector
v4=np.array([52.5611,123.417,238.28,529.885,799.486,412.364,452.157,575.584,575.375,539.866,622.613,571.542,609.59,638.393])

# test_stencil_3d_range_vector2
v5=np.array([52.0756,107.893,114.983,127.023,151.959,148.539,170.467,170.199,187.586,210.947,233.839,251.122,252.643,276.983])

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
