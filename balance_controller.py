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
# Maybe Ixx1 is front wheel moment of inertia?
# Measured about wheel axle
  'Ixx1': 0.0431,
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
Constitutional relationships
"""
def term_determination(motor_dict, bike_dict):
  m, b = motor_dict, bike_dict
  FyyPrime = b.Iyyf * cos(b.lam)**2 + b.Iyzf * sin(2 * b.lam) + b.Izzf * sin(b.lam) ** 2
  FyzPrime = b.Iyyf * sin(b.lam) * cos(b.lam) - b.Iyzf * cos(2 * b.lam) - b.Izzf * sin(b.lam) * cos(b.lam);
  FzzPrime = b.Iyyf * sin(b.lam)**2 - b.Iyzf * sin(2 * b.lam) + b.Izzf * cos(b.lam) ** 2
  mt = b.mr + b.mf
  ht = (b.mr * b.hr + b.mf * b.hf)/mt
  lt = (b.mr * b.lr + b.mf * b.lf)/mt
  Tyy = b.Iyyr + b.mr * b.hr**2 + FyyPrime + b.mf * b.hf ** 2
  Tzz = b.Izzr + b.mr * b.lr**2 + FzzPrime + b.mf * (b.Cw + b.lf) ** 2
  """
  check below reference equation A 11 for Tyz
  http://ruina.tam.cornell.edu/research/topics/bicycle_mechanics/*FinalBicyclePaperv45wAppendix.pdf
  """
  Tyz = b.Iyzr + b.mr * b.hr*b.lr + FyzPrime + b.mf * b.hf * (b.Cw + b.lf)

  FlamyP = -b.Iyzf * cos(b.lam) - b.Izzf * sin(b.lam) - b.mf * b.hf * b.d
  FlamzPP = -b.Iyzf * sin(b.lam) + b.Izzf * cos(b.lam) - b.mf * (b.Cw + b.lf) * b.d
  FlamlamP = b.Izzf + b.mf * b.d ** 2
  
  v = b.mf * b.d + mt * lt/b.Cw * b.Cf
  return dotdict({
      'FyyPrime': FyyPrime,
      'FyzPrime': FyzPrime,
      'FzzPrime': FzzPrime,
      'mt': mt,
      'ht': ht,
      'lt': lt,
      'Tyy': Tyy,
      'Tzz': Tzz,
      'Tyz': Tyz,
      'FlamyP': FlamyP,
      'FlamzPP': FlamzPP,
      'FlamlamP': FlamlamP,
      'v': v
  });

terms = term_determination(motor_params, bike_params)

"""
Term determination
http://stackoverflow.com/questions/20969773/exponentials-in-python-x-y-vs-math-powx-y
Using ** is faster without invoking function
"""
def pre_processor(terms, motor_dict, bike_dict, V):
    t, m, b = terms, motor_dict, bike_dict
    g = 9.80035 
# Roll Equation
    Mxx = t.Tyy
    Cxx = 0
    Kxx = -t.mt * g * t.ht
    Mxphi = t.FlamyP + t.Tyz * b.Cf/b.Cw
    Cxphi = V * (t.Tyz * cos(b.lam)/b.Cw - b.Cf/b.Cw * (b.Ixx1/b.rf + b.Ixx2/b.rr) - b.Ixx1/b.rf * cos(b.lam) - t.mt * t.ht * b.Cf / b.Cw)
    Kxphi = V ** 2 * (-cos(b.lam)/b.Cw*(b.Ixx1/b.rf + b.Ixx2/b.rr) - t.mt * t.ht * cos(b.lam)/b.Cw + g * t.v)
# Sterring equation
    Mphiphi = t.FlamlamP + 2 * t.FlamzPP * b.Cf/b.Cw + t.Tzz * (b.Cf/b.Cw) ** 2
    Cphiphi = V * (t.FlamzPP * cos(b.lam)/b.Cw + b.mf * b.d * b.Cf/b.Cw + t.Tzz * b.Cf/b.Cw ** 2 * cos(b.lam) + t.mt * t.lt * (b.Cf/b.Cw) ** 2)
    Kphiphi = V**2 * (b.Ixx1/b.rf * cos(b.lam)/b.Cw * sin(b.lam) + b.mf * b.d * cos(b.lam)/ b.Cw + t.mt * t.lt * b.Cf/b.Cw ** 2 * cos(b.lam)) - g*t.v*sin(b.lam)

    Mphix = t.FlamyP + t.Tyz * b.Cf/b.Cw
    Cphix = V * (b.Ixx1/b.rf*cos(b.lam) + (b.Ixx1/b.rf + b.Ixx2/b.rr)*b.Cf/b.Cw)
    Kphix = g * t.v

#   steering torque about the steering axis, applied to the front assembly by the sterring motor (control input)
#    Mphi
#   An external roll distrubance moment
#   

    return dotdict({
        'Mxx': Mxx,
        'Kxx': Kxx,
        'Mxphi': Mxphi,
        'Cxx': Cxx,
        'Cxphi': Cxphi,
        'Kxphi': Kxphi,

        'Mphiphi': Mphiphi,
        'Cphiphi': Cphiphi,
        'Kphiphi': Kphiphi,
        'Mphix': Mphix,
        'Cphix': Cphix,
        'Kphix': Kphix
      });
# Forward Velocity
V = 1
rtn = pre_processor(terms, motor_params, bike_params, V)
print rtn


