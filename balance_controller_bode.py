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
    pt1 = b.mr * b.hr ** 2 + b.mf * b.hf ** 2 +  b.Iyyr + b.Iyyf * cos(b.lam)**2
    pt2 = b.Iyzf*sin(2*b.lam) + b.Izzf*sin(2*b.lam) + b.Izzf*sin(b.lam)**2
    return pt1 + pt2

rtn = pre_processor(motor_params, bike_params)
print rtn
