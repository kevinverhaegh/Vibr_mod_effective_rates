import matplotlib.pyplot as plt
import numpy as np

# create some sample data
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.tan(x)
c = x  # the values to map to line colors

# create a figure and axis object
fig, ax = plt.subplots()

# plot each line with a different color and the specified values
line1 = ax.plot(x, y1, c=c, cmap='viridis', label='sin(x)')
line2 = ax.plot(x, y2, c=c, cmap='viridis', label='cos(x)')
line3 = ax.plot(x, y3, c=c, cmap='viridis', label='tan(x)')

# add a colorbar to the plot
cbar = fig.colorbar(line1, ax=ax)
cbar.ax.set_ylabel('x')

# add a legend to the plot
ax.legend()

# show the plot
plt.show()