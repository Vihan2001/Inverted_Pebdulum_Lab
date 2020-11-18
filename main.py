import sympy as sym

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
sym.pprint(dphi_F_eq)
sym.pprint(dphi_x3_eq)
sym.pprint(dphi_x4_eq)

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
sym.pprint(dpsi_F_eq)
sym.pprint(dpsi_x3_eq)
sym.pprint(dpsi_x4_eq)

