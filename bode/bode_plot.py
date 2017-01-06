# Require matplotlib
# OSX installation guide
# http://matplotlib.org/users/installing.html#build-osx
# blog.analgomachine/2015/02/07/bode-plots-using-python/
# http://nbviewer.jupyter.org/url/aspicc.fs.cvut.cz/PVVR/CourseWare/IPython_Notebook_7.ipynb
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

# Coefficients in numerator of transfer function
num = [1]
# Coefficients in denominator of transfer function
# High order to low order, eg 1*s^2 + 0.1*s + 1
den = [1, 0.1, 1]

# Scan over zeta, a parameter for a second-order system
zetaRange = np.arange(0.1,1.1,0.1)

f1 = plt.figure()
for i in range(0,9):
    den = [1, 2*zetaRange[i], 1]
    # print den
    s1 = signal.lti(num, den)
    print s1
    # Specify our own frequency range: np.arange(0.1, 5, 0.01)
    w, mag, phase = signal.bode(s1, np.arange(0.1, 5, 0.01).tolist())
    plt.semilogx (w, mag, color="blue", linewidth="1")
plt.xlabel ("Frequency")
plt.ylabel ("Magnitude")
plt.savefig ("mag.png", dpi=300, format="png")

plt.figure()

for i in range(0,9):
    den = [1, zetaRange[i], 1]
    s1 = signal.lti(num, den)
    w, mag, phase = signal.bode(s1, np.arange(0.1, 10, 0.02).tolist())
    plt.semilogx (w, phase, color="red", linewidth="1.1")
plt.xlabel ("Frequency")
plt.ylabel ("Amplitude")
plt.savefig ("phase.png", dpi=300, format="png")
