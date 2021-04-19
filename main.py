import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from scipy.constants import c
import scipy.optimize

#data = np.loadtxt('dec_lengths.txt')

#plt.hist(data,bins=50)
#plt.show()


pion_m_positive = 2.48807 * 10**(-34)   # kg, 139.57039*1e6 eV/c^2      from pdg
pion_m_neutral = 2.40618 * 10**(-34)    # kg, 134.9768*1e6   #eV/c^2    from pdg
kaon_m = 8.80059 * 10**(-34)            # kg, 493.677*1e6 eV/c^2        from pdg

pion_lifetime_positive = 2.6033*10**(-8) #seconds   from pdg
kaon_lifetime = 1.2380*10**(-8)          #seconds    from pdg

average_decay_length_pion_positive = 4.188*1e3 #meters  from script

# pion_v = pion_l / pion_t
# p = pion_m * pion_v = pion_m * pion_l / pion_t
# kaon_v = p / kaon_m = pion_m / kaon_m * pion_l / pion_t
# kaon_l = kaon_v * kaon_t = pion_m / kaon_m * kaon_t / pion_t * pion_l

p = pion_m_positive * average_decay_length_pion_positive / pion_lifetime_positive
average_decay_length_kaon = pion_m_positive / kaon_m * kaon_lifetime / pion_lifetime_positive * average_decay_length_pion_positive

pion_E_positive = sqrt((pion_m_positive*c**2)**2 + (p*c)**2)  # E = sqrt( (mc**2)**2 + (p*c)**2)
pion_E_neutral = sqrt((pion_m_neutral*c**2)**2 + (p*c)**2)

pion_fv_positive = np.array([pion_E_positive, 0, 0, p])  # four-vector parallel to z-axis
pion_fv_neutral = np.array([pion_E_neutral, 0, 0, p])


def random_isotropic_rotation(pion_fv):  # ~Marsaglia method
    fv = 1 * pion_fv
    while True:
        a = 2 * np.random.rand() - 1
        b = 2 * np.random.rand() - 1
        if a * a + b * b < 1:
            break
    fv[1] = fv[3] * 2 * a * sqrt(1 - a * a - b * b)
    fv[2] = fv[3] * 2 * b * sqrt(1 - a * a - b * b)
    fv[3] = fv[3] * (1 - 2 * (a * a + b * b))
    return fv


data_len = 2000
pions_positive = np.zeros((data_len, 4))
pions_neutral = np.zeros((data_len, 4))
velocity_pion_positive = np.zeros((data_len, 3))
velocity_pion_neutral = np.zeros((data_len, 3))
decay_position = np.zeros(data_len)

for i in range(data_len):
    v = (random_isotropic_rotation(pion_fv_positive))
    w = (random_isotropic_rotation(pion_fv_neutral))
    pions_positive[i] = v
    pions_neutral[i] = w
    velocity_pion_positive[i] = [x / pion_m_positive for x in v][1:]
    velocity_pion_neutral[i] = [x / pion_m_neutral for x in w][1:]
    s = np.random.exponential(average_decay_length_kaon)
    decay_position[i] = s


def number_of_detections(detector_position):
    results_positive = []
    distance = np.zeros(data_len)
    for k in range(data_len):
        distance[k] = detector_position-decay_position[k]
        if velocity_pion_positive[k][2] <= 0 and detector_position-decay_position[k]>0:
            results_positive.append(0)
        elif velocity_pion_positive[k][2] >= 0 and detector_position-decay_position[k]<0:
            results_positive.append(0)
        else:
            t = abs(distance[k]) / velocity_pion_positive[k][2]
            dx = t * velocity_pion_positive[k][0]
            dy = t * velocity_pion_positive[k][1]
            if dx**2+dy**2 > 4:
                results_positive.append(0)
            else:
                results_positive.append(1)
    count_positive = results_positive.count(1) # this is how many times the detector detects the positive pion

    results_neutral = []
    for k in range(data_len):
        if velocity_pion_neutral[k][2] <= 0:
            results_neutral.append(0)
        else:
            t = abs(distance[k])/velocity_pion_neutral[k][2]
            dx = t*velocity_pion_neutral[k][0]
            dy = t * velocity_pion_neutral[k][1]
            if dx**2+dy**2>4:
                results_neutral.append(0)
            else:
                results_neutral.append(1)
    count_neutral = results_neutral.count(1) # this is how many times the detector detects the neutral pion

    if count_positive != 0 and count_neutral != 0:
        indices_positive = [i for i, x in enumerate(results_positive) if x == 1]
        indices_neutral = [i for i, x in enumerate(results_neutral) if x == 1]
        success = [x for x in indices_positive if x in indices_neutral]
    else:
        success = []
    fails = data_len - len(success)
    return fails    # to maximise success, we minimise fails


optimal_distance = scipy.optimize.minimize(number_of_detections,x0=600)
print(optimal_distance["fun"])
print(optimal_distance["x"])


N = 100
distances = np.linspace(1, 10, N)
results = 1*distances
for i, d in enumerate(distances):
    results[i] = number_of_detections(d)
plt.plot(distances, results, marker = ".", color = "orange")
plt.show()
