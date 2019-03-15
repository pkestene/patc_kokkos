import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,96,128,160,192,224,256,320,384,448,512,640])

# test_stencil_3d_flat
v0=np.array([17.0526,20.1192,20.8597,21.2713,21.5433,21.3717,21.4785,21.4326,21.5032,21.1176,21.4883,21.024,21.5393,20.9602])

# test_stencil_3d_flat_vector without views
v1=np.array([51.5658,105.832,143.236,175.288,188.928,187.224,190.218,185.766,177.068,194.432,182.32,193.811,185.365,191.734])

# test_stencil_3d_flat_vector with    views
v2=np.array([51.7569,165.561,313.014,617.077,831.131,743.104,623.098,696.911,624.087,659.44,643.795,661.474,636.55,655.334])

# test_stencil_3d_range
v3=np.array([33.9946,49.7425,70.354,77.1893,80.472,82.5614,83.2173,83.9458,82.0809,84.119,84.7935,85.4897,85.6791,86.8686])

# test_stencil_3d_range_vector
v4=np.array([30.9994,79.4466,141.164,313.251,525.169,451.435,534.583,647.969,612.814,423.384,494.973,589.604,638.997,552.043])

# test_stencil_3d_range_vector2
v5=np.array([38.7989,74.4964,84.8573,111.685,130.603,144.88,154.006,157.9,156.519,189.014,204.696,230.797,240.847,261.886])

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
