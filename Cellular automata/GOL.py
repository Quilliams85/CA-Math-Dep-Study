import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

N=60

diehard = [[0, 0, 0, 0, 0, 0, 1, 0],
           [1, 1, 0, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 0, 1, 1, 1]]

class Automaton():
    def __init__(self):
        self.size = N
        self.grid = np.random.randint(0, 2, size=(N, N))
        stamp = np.ndarray(diehard)
        self.grid[20:20+stamp.size, 20:20+stamp.size]

    def updateCell(self, X):
        nbrs_count = sum(np.roll(np.roll(X, i, 0), j, 1)
                     for i in (-1, 0, 1) for j in (-1, 0, 1)
                     if (i != 0 or j != 0))
        return (nbrs_count == 3) | (X & (nbrs_count == 2))

    def animate(self, i):
        self.grid = self.updateCell(self.grid)
        im.set_array(self.grid)  # Update the image data
        return [im]


# Initial setup
conway = Automaton()
array = conway.grid

# Plot setup
fig, ax = plt.subplots()
im = ax.imshow(array, cmap='binary')  # Display array using binary color map

# Create animation
ani = FuncAnimation(fig, conway.animate, frames=200, interval=20, blit=True)

plt.show()


