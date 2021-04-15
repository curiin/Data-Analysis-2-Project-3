import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt('dec_lengths.txt')


plt.hist(d,bins=50)
plt.show()
