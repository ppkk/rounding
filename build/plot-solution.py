import numpy as np
import matplotlib.pyplot as plt


def plot_graph(values, marker, label):
   xx = [x[0] for x in values]
   val = [x[1] for x in values]
   plt.plot(xx, val, marker, label=label)

f = np.loadtxt("sol/solution-131.txt")
plot_graph(f, "g-", "131")
#f = np.loadtxt("sol/solution-3960.txt")
#plot_graph(f, "ro-", "3960")
#f = np.loadtxt("sol/solution-41993.txt")
#plot_graph(f, "b-", "41993")


plt.legend(loc=3)
plt.grid()
plt.show()
