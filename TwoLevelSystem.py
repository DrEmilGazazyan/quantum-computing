import numpy as np
from scipy.integrate import solve_ivp

# Определение функции
def system_of_equations(t, r, W, G, gam, delt):
    r = np.asarray(r, dtype=np.complex128)
    drdt = np.zeros_like(r, dtype=np.complex128)

    drdt[0] = 1j*W/2*(r[1]-r[2])+G*r[3]
    drdt[3] = -1j*W/2*(r[1]-r[2])-G*r[3]
    drdt[1] = -(gam+1j*delt)*r[1]-1j*W/2*(r[3]-r[0])
    drdt[2] = np.conj(r[1])
    return drdt

W = 1
G = 1
gam = G/2
delt = 0
r0 = [1+0j, 0+0j, 0+0j, 0+0j] # Начальные значения переменных
t_span = [0, 10] # Интервал интегрирования

# Решение системы уравнений
sol = solve_ivp(system_of_equations, t_span, r0, 'RK45', args=(W, G, gam, delt))