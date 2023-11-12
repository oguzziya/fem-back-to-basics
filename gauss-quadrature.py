def gauss_quadrature(n, f):
  x = []
  w = []

  if n == 2:
    x = [-1, 1]
    w = [1, 1]

  elif n == 3:
    x = [-1, 0, 1]
    w = [1./3., 4./3., 1./3.]
  
  val = 0

  for x_i, w_i in x, w:
    val += f(x_i) * w_i 

  return val