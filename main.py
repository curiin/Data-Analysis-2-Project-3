import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt('dec_lengths.txt')


plt.hist(d,bins=50)
plt.show()


# from the website: mean lifetime of pion 2.6033*10^-8s
# from the script: average decay length of pion 4.188km
# velocity of beam: 4.188/(2.6*10^-8s) = 1.6*10^8 km/s
# from the website: mean lifetime of kaon 1.2380*10^-8s
# expected average decay length of kaon 1.922km
average_decay_length_keon = 1.922*10**3 #m

print(np.random.exponential(average_decay_length_keon))

