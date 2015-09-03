import numpy as np
import matplotlib.pyplot as plt


def plot_graph(values, marker, label):
   xx = [x[0] for x in values]
   val = [x[1] for x in values]
   plt.plot(xx, val, marker, label=label)

f = np.loadtxt("sol/error-35.txt")
plot_graph(f, "m-", "35")
f = np.loadtxt("sol/error-131.txt")
plot_graph(f, "r-", "131")
f = np.loadtxt("sol/error-373.txt")
plot_graph(f, "g-", "373")
f = np.loadtxt("sol/error-485.txt")
plot_graph(f, "b-", "485")


plt.legend(loc=3)
plt.grid()
plt.show()
