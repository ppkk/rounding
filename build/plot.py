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

double = np.loadtxt("double_1e6.txt")
plot_graph(double, 1, "g-", "$\epsilon = 10^6$, double - l2")
plot_graph(double, 2, "g--", "$\epsilon = 10^6$, double - h1 semi")

double = np.loadtxt("double_1e-6.txt")
plot_graph(double, 1, "m-", "$\epsilon = 10^{-6}$, double - l2")
plot_graph(double, 2, "m--", "$\epsilon = 10^{-6}$, double - h1 semi")

double = np.loadtxt("double_0.txt")
plot_graph(double, 1, "c-", "$\epsilon = 0$, double - l2")
plot_graph(double, 2, "c--", "$\epsilon = 0$, double - h1 semi")

   

double = np.loadtxt("double_1.txt")
plot_graph(double, 1, "r-", "$\epsilon = 1$, double - l2")
plot_graph(double, 2, "r--", "$\epsilon = 1$, double - h1 semi")
fl = np.loadtxt("float_1.txt")
plot_graph(fl, 1, "b-", "$\epsilon = 1$, single - l2")
plot_graph(fl, 2, "b--", "$\epsilon = 1$, single - h1 semi")

plt.legend(loc=3)
plt.xlabel("discretization step $h$")
plt.grid()
plt.xlim(maxx, minx)
plt.yscale('log')
plt.xscale('log')
plt.show()
