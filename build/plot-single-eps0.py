import numpy as np
import matplotlib.pyplot as plt


def plot_graph(values, position, marker, label):
   n = [x[0] for x in values]
   h = [1./x[0] for x in values]
   global minx, maxx
   minx = min(h)
   maxx = max(h)
   val = [x[position] for x in values]
   plt.plot(h, val, marker, label=label)

double = np.loadtxt("f-f-f-eps0.txt")
plot_graph(double, 1, "gx-", "$\epsilon = 0$, f-f-f - l2")
plot_graph(double, 2, "gx--", "$\epsilon = 0$, f-f-f - h1 semi")

double = np.loadtxt("f-f-d-eps0.txt")
plot_graph(double, 1, "r-", "$\epsilon = 0$, f-f-d - l2")
plot_graph(double, 2, "r--", "$\epsilon = 0$, f-f-d - h1 semi")

double = np.loadtxt("f-d-f-eps0.txt")
plot_graph(double, 1, "bo-", "$\epsilon = 0$, f-d-f - l2")
plot_graph(double, 2, "bo--", "$\epsilon = 0$, f-d-f - h1 semi")

double = np.loadtxt("f-d-d-eps0.txt")
plot_graph(double, 1, "m-", "$\epsilon = 0$, f-d-d - l2")
plot_graph(double, 2, "m--", "$\epsilon = 0$, f-d-d - h1 semi")

plt.legend(loc=3)
plt.xlabel("discretization step $h$")
plt.grid()
plt.xlim(maxx, minx)
plt.yscale('log')
plt.xscale('log')
plt.show()
