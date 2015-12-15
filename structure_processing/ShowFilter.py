import matplotlib.pyplot as plt
import scipy.ndimage.filters as f
import numpy as np

x = np.random.rand(100,100)
for i in range(10):
    x[(i*10):((i+1)*10)] = i
mask = np.ones(x.shape)
mask[40:60,40:60] = 0
x[mask==0] = 0

fig = plt.figure()
ax = fig.add_subplot(111)

blurred = f.gaussian_filter(x, 3)

blurred_mask = f.gaussian_filter(mask, 3)

crap = blurred/blurred_mask
crap[mask==0] = 0

cax = ax.matshow(mask)
fig.colorbar(cax)
plt.show()