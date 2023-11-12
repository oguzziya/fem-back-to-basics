import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines

plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = "serif"

def gauss_quadrature(n, f):
  x = []
  w = []

  if n == 2:
    x = [-1, 1]
    w = [1, 1]

  elif n == 3:
    x = [-1., 0., 1.]
    w = [1./3., 4./3., 1./3.]
  
  val = 0

  for x_i, w_i in zip(x, w):
    val += f(x_i) * w_i 

  return val

def lagrange_polynomial(p, i, x):
  y = 1.

  ksi = lambda j : -1 + j*(2./p)

  for k in range(p+1):
    if k is not i-1:
      y *= (x - ksi(k)) / (ksi(i-1) - ksi(k))

  return y

def legendre_base(n, x):
  if n == 0:
    return 1
  elif n == 1:
    return x
  else:
    return ((2*n - 1)*x*legendre_base(n-1, x) - (n-1)*legendre_base(n-2, x))/n

def legendre_polynomial(i, x):
    if i == 1:
      return 0.5*(1 - x)
    elif i == 2:
      return 0.5*(1 + x)
    else:
      return (1./np.sqrt(4*i - 6)) * (legendre_base(i-1, x) - legendre_base(i-3, x))

class Lagrange:
  def __init__(self, p):
    self.name = "Lagrange"
    self.p = p
    self.shape_functions = [(lambda x, i=i: lagrange_polynomial(p, i+1, x)) for i in range(p+1)]

class Legendre:
  def __init__(self, p):
    self.name = "Legendre"
    self.p = p
    self.shape_functions = [(lambda x, i=i: legendre_polynomial(i+1, x)) for i in range(p+1)]

def plot_shape_functions(basis):
  ksi = np.linspace(-1, 1, 100)
  
  fig, ax = plt.subplots()
  ax.set_xlabel(r'$\xi$')
  ax.set_title(r"Order {} {} Basis Functions".format(basis.p, basis.name))
  
  ax.spines['left'].set_position(('data', 0))
  ax.spines['bottom'].set_position(('data', 0))
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)

  ax.grid()

  for i in range(basis.p+1):
    shape_function = basis.shape_functions[i]
    ax.plot(ksi, shape_function(ksi), label=r'$N_{}$'.format(i+1), color="black")

  labelLines(ax.get_lines(), align=False)

  ax.xaxis.set_label_coords(1, 0.3)

p = 5

plot_shape_functions(Legendre(p))
plot_shape_functions(Lagrange(p))

plt.show()