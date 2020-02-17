import numpy as np

import matplotlib.pyplot as plt

    
data = np.loadtxt('aeroload.dat',dtype='float', delimiter=',')

plt.plot(data)
plt.title('Aero Load')
plt.show()
