import numpy as np
from scipy.stats import maxwell
import matplotlib.pyplot as plt

mean, var, skew, kurt = maxwell.stats(moments='mvsk')
fig, ax = plt.subplots(1, 1)


x = np.linspace(maxwell.ppf(0.01),
                maxwell.ppf(0.99), 100)
ax.plot(x, maxwell.pdf(x),
       'r-', lw=5, alpha=0.6, label='maxwell pdf')
rv = maxwell()
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
vals = maxwell.ppf([0.001, 0.5, 0.999])
np.allclose([0.001, 0.5, 0.999], maxwell.cdf(vals))
r = maxwell.rvs(size=1000)

ax.hist(r, density=True, bins='auto', histtype='stepfilled', alpha=0.2)
ax.set_xlim([x[0], x[-1]])
ax.legend(loc='best', frameon=False)
plt.show()

print(x)