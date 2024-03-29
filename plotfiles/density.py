import matplotlib.pyplot as plt
from numpy import loadtxt,reshape
import sys

re,im = loadtxt(sys.argv[1],unpack=True)
z = re + 1j * im
z = reshape(z,[513,513])

plt.imshow(abs(z)**2,cmap='CMRmap')
plt.show()
