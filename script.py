import numpy as np
import matplotlib.pyplot as plt

data = []
with open("data.txt", "r") as file:
    for line in file:
        data.append(float(line.strip()))
plt.xlabel("loops")
plt.ylabel("total energy")
plt.yscale('log')
plt.title("dt=0.0014")
plt.plot(data[2001:10000])
plt.grid()
plt.savefig("plot.png")
