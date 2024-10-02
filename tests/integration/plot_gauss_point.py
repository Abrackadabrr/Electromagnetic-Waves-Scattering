import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import legendre as leg

nx = 20
ny = 20
x = leg.legroots([0] * nx + [1])
y = leg.legroots([0] * ny + [1])
x_grid, y_grid = np.meshgrid(x, y)
fig, ax = plt.subplots()
ax.scatter(x_grid, y_grid, color="red")
ax.grid()
ax.plot([-1] * 2, [-1, 1], color="blue")
ax.plot([1] * 2, [-1, 1], color="blue")
ax.plot([-1, 1], [-1] * 2, color="blue")
ax.plot([-1, 1], [1] * 2, color="blue")
ax.minorticks_on()
ax.set_aspect("equal")
plt.show()
