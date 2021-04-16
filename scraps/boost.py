import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from scipy.constants import c


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


def boost(fv):
    beta = kaon_v/c
    gamma = 1 / sqrt(1 - beta * beta)
    LT = np.array([[gamma, 0, 0, beta * gamma * c],
                   [0, 1, 0, 0],
                   [0, 0, 1, 0],
                   [beta * gamma / c, 0, 0, gamma]])

    return LT @ fv


# assumed input: pion four-vector:

kaon_v = 0.4 * c # arbitrary
pion_v = 0.7 * c # arbitrary
pion_m = 139.57039 *1e6 # eV/c**2 from pdg

pion_p = pion_m * pion_v
pion_E = sqrt(pion_m * c * c * pion_m * c * c + pion_p * c * pion_p * c)
pion_fv = np.array([pion_E, 0, 0, pion_p])

data_len = 2000
pions = np.zeros((data_len, 4))

fig, ax = plt.subplots(1, 2)

for i in range(data_len):
    v = boost(random_isotropic_rotation(pion_fv))
    pions[i] = v

    ax[0].plot((0, v[3]), (0, v[2]), lw=0.2, marker=".", ms = 3, color="orange", alpha=0.5)
    ax[1].plot((0, v[3]), (0, v[1]), lw=0.2, marker=".", ms = 3, color="green", alpha=0.5)

ax[0].grid()
ax[1].grid()
ax[0].set_xlabel("p_z [eV]")
ax[0].set_ylabel("p_y [eV]")
ax[0].set_title("side_view")
ax[1].set_xlabel("p_z [eV]")
ax[1].set_ylabel("p_x [eV]")
ax[1].set_title("top_view")
ax[0].set_aspect("equal")
ax[1].set_aspect("equal")
fig.tight_layout()
plt.show()

# rotate pions after boosting instead of rotating kaons prior to boosting
