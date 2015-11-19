import numpy as np
import matplotlib.pyplot as plt


samples = 400**3
temp_vals = np.random.rand(samples)
q75, q25 = np.percentile(temp_vals, [75 ,25])
iqr = q75 - q25
bins = (np.max(temp_vals) - np.min(temp_vals))/(2*iqr*samples**(-1/3.0))
hist, bins = np.histogram(temp_vals, bins)
hist = np.asarray(hist, dtype=np.float_)
hist /= float(np.sum(hist))
centers = (bins[:-1] + bins[1:])/2
handles = []
handle, = plt.plot(centers, hist)
handles.append(handle)
legends = []
legends.append("test")
plt.legend(handles, legends)
plt.show()