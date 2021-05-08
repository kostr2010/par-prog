from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas
import sys

points = pandas.read_csv("res.csv")
fig1 = plt.figure(figsize=(16, 16))
ax1 = fig1.add_subplot(111, projection='3d')

x = points['x'].values
y = points['t'].values
z = points['u'].values

plt.xlabel('x')
plt.ylabel('t')

ax1.scatter(x, y, z, c='b', s=1)
plt.savefig("res.png", dpi=500)
