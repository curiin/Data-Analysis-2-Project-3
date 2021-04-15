import matplotlib.pyplot as plt
import numpy as np

f = np.loadtxt('dec_lengths.txt')


plt.hist(f,bins=50)
plt.show()
