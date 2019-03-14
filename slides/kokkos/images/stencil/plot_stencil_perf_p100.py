import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,92,128,160,192,224,256])

# test_stencil_3d_flat
v0=np.array([116.199,389.134,541.676,595.897,621.277,627.336,629.644,630.534,634.056])

# test_stencil_3d_flat_vector without views
v1=np.array([79.4376,204.12,383.083,414.372,706.766,928.356,999.039,1075.1,975.581])

# test_stencil_3d_flat_vector with    views
v2=np.array([84.3015,215.201,394.736,415.256,704.972,966.197,1086.54,1122.35,901.766])

# test_stencil_3d_range
v3=np.array([153.624,416.464,539.891,576.55,634.021,644.714,649.318,651.65,652.956])

# test_stencil_3d_range_vector
v4=np.array([57.0523,139.46,259.677,279.792,485.862,664.186,859.097,1061.64,825.661])

# test_stencil_3d_range_vector2
v5=np.array([78.6841,53.1807,116.682,64.2638,67.939,68.1239,68.6165,68.0034,68.3372])

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
