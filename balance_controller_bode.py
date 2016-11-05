# 
# Transfer function described in balance.tex
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

""" 
Transfer function coeffecitions described in balance.tex
"""
motor_params = {
  Jm: 102.3,
  kb: 0.0847,
  kt: 0.084,
  L: 2.18,
  N: 45,
  R: 4.887
}

bike_params = {
  Cf: 86.53,
  Cw: 1.08,
  d: 4.5,
  hf: 0.4110,
  hr: 0.4096,
  Ixx2: 0.0431,
  Iyyf: 1.1467,
  Iyzf: -0.012,
  Izzf: 0.0530,
  Iyyr: 4.7235,
  Iyzr: 0.2204,
  Izzr: 11.2341,
  lf: -0.1073,
  lr: 0.4,
  mf: 13.9242,
  mr: 105.3572,
  rf: 0.205,
  rr: 0.205,
  lam: 27 #(degree) should we convert it to radius:w
}

"""
Term determination
"""
def pre_processor(motor_dict, bike_dict):
    m, b = motor_dict, bike_dict
