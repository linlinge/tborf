import numpy as np
from sklearn.neighbors.kde import KernelDensity

### create data ###
sample_count = 1000
n = 6
data = np.random.randn(sample_count, n)
data_norm = np.sqrt(np.sum(data*data, axis=1))
data = data/data_norm[:, None]   # Normalized data to be on unit sphere

print(data.shape)
print(data)

## estimate pdf using KDE with gaussian kernel
kde = KernelDensity(kernel='gaussian', bandwidth=0.8).fit(data)

log_p = kde.score_samples(data)  # returns log(p) of data sample
p = np.exp(log_p)                # estimate p of data sample
entropy = -np.sum(p*log_p)       # evaluate entropy
print(entropy)