import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([32,48,64,92,128,160,192,224,256,320,384,448,512,640,768])

# test_stencil_3d_flat
v0=np.array([56.3896,189.76,435.238,595.145,620.955,627.115,629.582,630.391,632.74,636.846,636.648,636.366,694.32,699.607,706.804])

# test_stencil_3d_flat_vector without views
v1=np.array([54.3867,180.913,203.892,401.841,708.21,944.497,1055.79,1121.36,1011.14,898.973,855.719,829.163,811.193,747.457,736.711])

# test_stencil_3d_flat_vector with    views
v2=np.array([54.4862,182.827,269.114,401.841,742.118,1017.67,1139.15,1168.95,948.047,964.942,1029.57,1057.37,1105.86,1004.2,1017.68])

# test_stencil_3d_range
v3=np.array([99.9939,342.602,563.508,634.679,702.87,715.909,721.909,725.045,726.569,727.969,728.464,728.692,728.759,728.656,728.683])

# test_stencil_3d_range_vector
v4=np.array([59.7411,148.366,143.092,270.472,502.529,696.09,900.563,1081.41,832.389,1056.22,991.174,962.012,1063.21,947.978,992.752])

# test_stencil_3d_range_vector2
v5=np.array([48.6967,35.4723,62.2566,57.7341,68.822,68.957,68.912,68.436,68.3321,67.9594,67.7381,67.3909,67.6223,67.3647,66.957])

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
