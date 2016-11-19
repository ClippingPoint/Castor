# 
# Transfer function described in balance.tex
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
import math

"""
Monkey patch dot notation
http://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
"""
class dotdict(dict):
  """dot.notation access to dictionary attributes"""
  __getattr__ = dict.get
  __setattr__ = dict.__setitem__
  __delattr__ = dict.__delitem__
"""
Simplify naming of math functions
"""
cos = math.cos
sin = math.sin
""" 
Transfer function coeffecitions described in balance.tex
"""
motor_params = {
  'Jm': 102.3,
  'kb': 0.0847,
  'kt': 0.084,
  'L': 2.18,
  'N': 45,
  'R': 4.887
}

bike_params = {
  'Cf': 86.53,
  'Cw': 1.08,
  'd': 4.5,
  'hf': 0.4110,
  'hr': 0.4096,
  'Ixx2': 0.0431,
  'Iyyf': 1.1467,
  'Iyzf': -0.012,
  'Izzf': 0.0530,
  'Iyyr': 4.7235,
  'Iyzr': 0.2204,
  'Izzr': 11.2341,
  'lf': -0.1073,
  'lr': 0.4,
  'mf': 13.9242,
  'mr': 105.3572,
  'rf': 0.205,
  'rr': 0.205,
  'lam': 27 #(degree) should we convert it to radius:w
}

"""
dotify dictionary
"""
motor_params = dotdict(motor_params)
bike_params = dotdict(bike_params)
"""
Term determination
http://stackoverflow.com/questions/20969773/exponentials-in-python-x-y-vs-math-powx-y
Using ** is faster without invoking function
"""
def pre_processor(motor_dict, bike_dict):
    m, b = motor_dict, bike_dict
    FyyPrime = b.mf * b.hf ** 2 +  b.Iyyr + b.Iyyf * cos(b.lam)**2 
    + b.Iyzf*sin(2*b.lam) + b.Izzf*sin(2*b.lam) + b.Izzf*sin(b.lam)**2
    Mxx = b.Iyyr + b.mr * b.he ** 2 + FyyPrime + b.mf * b.hf ** 2
    ## mt gt h
    # parameters not listed
    mt = b.mf + b.mr
    # Central point height
    ht =  (b.mr * b.hr + b.mf * b.hf)/mt
    # Gravitational constant
    g = 9.80035 
    Kxx = mb * ht * g
    # Mxphi
    # Flambdaf' + Tyz Cf/Cw
    Flambdaf = -m.Iyzf * cos(b.lam) - m.Izzf * sin(b.lam) - b.mf * b.hf * b.d

#    Mxp = 
    return Mxx
rtn = pre_processor(motor_params, bike_params)
print rtn
