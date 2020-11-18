import sympy as sym
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np


#Define symbols to included
M, m, ell, g, x1, x2, x3, x4, F = sym.symbols('M, m, ell, g, x1, x2, x3, x4, F')

#Define phi
phi = 4*m*ell*x4**2*sym.sin(x3) + 4*F - 3*m*g*sym.sin(x3)*sym.cos(x3)
phi /= 4*(M+m) - 3*m*sym.cos(x3)**2

#Derivatives of psi
dphi_F = phi.diff(F)
dphi_x3 = phi.diff(x3)
dphi_x4 = phi.diff(x4)

#Derivatives of phi taken equilibrium point
dphi_F_eq = dphi_F.subs([(F, 0), (x3, 0), (x4, 0)])
dphi_x3_eq = dphi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
dphi_x4_eq = dphi_x4.subs([(F, 0), (x3, 0), (x4, 0)])

#Pretty print phi
#sym.pprint(dphi_F_eq)
#sym.pprint(dphi_x3_eq)
#sym.pprint(dphi_x4_eq)

#Define psi
psi = -3*(m*ell*x4**2*sym.sin(x3)*sym.cos(x3) + F*sym.cos(x3) - (M + m)*g*sym.sin(x3)) / (4*(M + m) - 3*m*sym.cos(x3)**2) * ell

#Derivatives of psi
dpsi_F = psi.diff(F)
dpsi_x3 = psi.diff(x3)
dpsi_x4 = psi.diff(x4)

#Derivative of psi taken at equilibrium point
dpsi_F_eq = dpsi_F.subs([(F, 0), (x3, 0), (x4, 0)])
dpsi_x3_eq = dpsi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
dpsi_x4_eq = dpsi_x4.subs([(F, 0), (x3, 0), (x4, 0)])

#Pretty Print psi
#sym.pprint(dpsi_F_eq)
#sym.pprint(dpsi_x3_eq)
#sym.pprint(dpsi_x4_eq)

#Define the positive real constants
a = dphi_F_eq
b = -dphi_x3_eq
c = -dpsi_F_eq
d = dpsi_x3_eq

#Define parameter values
M_value = 0.3
m_value = 0.1
ell_value = 0.35
g_value = 9.81

#Define function to substitute values

def substitute_values (z):

    return float(z.subs([(M, M_value), (m, m_value), (ell, ell_value), (g, g_value)]))

a_value = substitute_values(a)
b_value = substitute_values(b)
c_value = substitute_values(c)
d_value = substitute_values(d)

G_theta_tf = ctrl.TransferFunction([-c_value], [1, 0, -d_value])

#Controller Design

#Define function that takes in Kp, Ki, Kd and returns a CTRL.TF for pid controller

def pid(Kp, Ki, Kd):
    G_pid = Kp
    # (1s + 0) / 1 ==> from s/1 = s
    G_pid += Kd * ctrl.TransferFunction([1, 0], [1]) #Kd.s/1
    G_pid += Ki * ctrl.TransferFunction([1], [1, 0]) #ki.1/s
    return G_pid

my_pid = -pid(Kp =100, Ki= 3 , Kd= 10)


G_closed_loop = ctrl.feedback(G_theta_tf, my_pid)

t_final = 1
num_points = 500
t_imp, theta_imp = ctrl.impulse_response(G_closed_loop, T=np.linspace(0, t_final, num_points))


theta_imp_deg = np.rad2deg(theta_imp) #convert to degrees
plt.plot(t_imp, theta_imp_deg)
plt.xlabel('Time (s)')
plt.ylabel('Angle (deg)')
plt.grid()
plt.show()