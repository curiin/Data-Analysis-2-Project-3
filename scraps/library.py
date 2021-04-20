import numpy as np
from scipy.constants import c
from math import sqrt


m_k = 4.93677 * 1e8 # eV/C**2
m_pp = 1.3957039 * 1e8
m_np = 1.349768 * 1e8

E_k = c*c * m_k # rest energies of the three particlesp
E_pp = c*c * m_pp
E_np = c*c * m_np

pion_p = sqrt((E_k**4 + E_pp**4 + E_np**4 - 2*E_pp*E_pp*E_np*E_np - 2*E_k*E_k*(E_pp*E_pp + E_np*E_np))
               /(4*(E_k*E_k + E_pp*E_pp + E_np*E_np)))/c

'''
momentum calculated based on the four-vector conservation,
pion+ momentum is same in value but opposite in direction to pion0 momentum
'''