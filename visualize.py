import matplotlib.pyplot as plt
from math import pi

num_frames = 10

def setBeautifulGrid ():
    plt.grid(which='major',
            color = 'k',
            linewidth = 1)
    plt.minorticks_on()
    plt.grid(which='minor',
            color = 'k',
            linestyle = ':')


T = list (map (float, input().split()))
X = list (map (float, input().split()))

mod = len(T) // num_frames

for i, t in enumerate (T):
    u = list (map (float, input().split()))
    if i % mod != 0: continue
    plt.figure ()
    setBeautifulGrid ()
    plt.plot (X, u)
    plt.show()