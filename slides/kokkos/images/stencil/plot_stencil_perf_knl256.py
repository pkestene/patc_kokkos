import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,96,128,160,192,224,256,320,384,448,512,640])

# test_stencil_3d_flat
v0=np.array([24.3625,32.806,35.2833,37.2498,37.7838,37.527,38.0637,37.9272,37.6389,38.0641,38.1907,38.2284,38.054,37.9101])

# test_stencil_3d_flat_vector without views
v1=np.array([55.4364,128.154,180.266,250.799,270.551,276.052,277.775,283.8,260.708,296.565,291.619,294.439,283.076,299.309])

# test_stencil_3d_flat_vector with    views
v2=np.array([50.3872,171.829,341.496,752.661,923.217,564.808,583.614,708.707,581.968,733.306,682.989,723.346,736.688,668.489])

# test_stencil_3d_range
v3=np.array([41.4695,64.1461,102.089,120.247,124.087,126.04,127.335,127.679,116.075,130.185,129.534,131.404,127.682,131.367])

# test_stencil_3d_range_vector
v4=np.array([27.5688,69.7463,123.154,268.025,440.796,482.432,595.822,744.896,644.248,748.336,750.909,768.409,733.876,628.894])

# test_stencil_3d_range_vector2
v5=np.array([36.3956,71.0347,84.4094,104.569,120.576,137.086,148.578,155.703,152.617,173.711,187.026,226.408,247.446,286.201])

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
