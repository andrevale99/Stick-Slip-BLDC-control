import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def drill_string_open_loop(t, y, params):
    """
    Define o sistema de equações diferenciais de malha aberta para a coluna de perfuração.
    y = [varphi_r, omega_r, varphi_b, omega_b]
    onde:
    - varphi_r: Ângulo no topo
    - omega_r: Velocidade angular no topo
    - varphi_b: Ângulo na broca (bottom)
    - omega_b: Velocidade angular na broca
    """
    Ir = params['Ir']
    Ib = params['Ib']
    C = params['C']
    K = params['K']
    cr = params['cr']
    cb = params['cb']
    
    varphi_r, omega_r, varphi_b, omega_b = y
    
    # Entrada de torque no topo (Tm) - Exemplo: degrau de torque em t = 1s
    Tm = params['Tm_func'](t)
    
    # Torques resistivos (atrito viscoso/estrutural no topo e na broca)
    Tr = cr * omega_r
    Tb = cb * omega_b
    
    # Equações de movimento (Dinâmica Torsional)
    # Ir * omega_dot_r + C*(omega_r - omega_b) + K*(varphi_r - varphi_b) = Tm - Tr
    d_varphi_r = omega_r
    d_omega_r = (Tm - Tr - C * (omega_r - omega_b) - K * (varphi_r - varphi_b)) / Ir
    
    # Ib * omega_dot_b - C*(omega_r - omega_b) - K*(varphi_r - varphi_b) = -Tb
    d_varphi_b = omega_b
    d_omega_b = (C * (omega_r - omega_b) + K * (varphi_r - varphi_b) - Tb) / Ib
    
    return [d_varphi_r, d_omega_r, d_varphi_b, d_omega_b]


# Parâmetros do sistema baseados na literatura de modelagem torsional[cite: 2]
params = {
    'Ir': 120.0,       # Momento de inércia no topo (kg.m^2)
    'Ib': 80.0,        # Momento de inércia na broca/BHA (kg.m^2)
    'C': 15.0,         # Amortecimento torsional da coluna (N.m.s/rad)
    'K': 800.0,        # Rigidez torsional da coluna (N.m/rad)
    'cr': 2.0,         # Coeficiente de atrito no topo
    'cb': 5.0,         # Coeficiente de atrito na broca
    'Tm_func': lambda t: 2000.0 if t >= 1.0 else 0.0  # Torque de entrada (Tm)
}
