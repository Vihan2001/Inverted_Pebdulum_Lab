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

    return z.subs([(M, M_value), (m, m_value), (ell, ell_value), (g, g_value)])

a_value = substitute_values(a)
b_value = substitute_values(b)
c_value = substitute_values(c)
d_value = substitute_values(d)

#Define Gtheta
s, t =sym.symbols('s, t')
a, b, c, d = sym.symbols('a, b, c, d', real = True, positive = True) #define values as real and positive
Gtheta_s = -c / (s**2 -d)
#sym.pprint(Gtheta_s) #print transfer function

# Impulse Response, F(t) is Dirac
F_s = 1
Theta_s = Gtheta_s * F_s
theta_t = sym.inverse_laplace_transform(Theta_s, s, t)
#sym.pprint(theta_t) #Print Impulse Response
#print(sym.latex(theta_t))

#Step Response, F(t) is Heaviside
F_s = 1 / s
Theta_s = Gtheta_s * F_s
theta_t = sym.inverse_laplace_transform(Theta_s, s, t)
#sym.pprint(theta_t) #Print Step Response
#print(sym.latex(theta_t))

#Frequency Response, F(t) = sin(wt)
w = sym.symbols('w', real = True)
F_s = w / (s**2 + w**2)
Theta_s = Gtheta_s * F_s
theta_t = sym.inverse_laplace_transform(Theta_s, s, t)
#sym.pprint(theta_t.simplify()) #Print Step Response and simplify result to remove imaginary component
#print(sym.latex(theta_t.simplify()))


#Define Gx
s, t =sym.symbols('s, t')
a, b, c, d = sym.symbols('a, b, c, d', real = True, positive = True) #define values as real and positive
Gx_s = (a)/(s**2) + (b*c) / (s**4 -s**2 * d)
#sym.pprint(Gx_s) #print transfer function


# Impulse Response, F(t) is Dirac
F_s = 1
X_s = Gx_s * F_s
x_t = sym.inverse_laplace_transform(X_s, s, t)
#sym.pprint(x_t.simplify()) #Print Impulse Response and simplify
#print(sym.latex(x_t.simplify()))

#Step Response, F(t) is Heaviside
F_s = 1 / s
X_s = Gx_s * F_s
x_t = sym.inverse_laplace_transform(X_s, s, t)
#sym.pprint(x_t.simplify()) #Print Step Response and simplify
#print(sym.latex(x_t.simplify()))


#Frequency Response, F(t) = sin(wt)
w = sym.symbols('w', real = True)
F_s = w / (s**2 + w**2)
X_s = Gx_s * F_s
x_t = sym.inverse_laplace_transform(X_s, s, t)
#sym.pprint(x_t.simplify()) #Print Step Response and simplify
print(sym.latex(x_t.simplify()))
