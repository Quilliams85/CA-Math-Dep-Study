import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy

N=60


class Automaton():
    def __init__(self):
        self.size = N
        self.grid = np.random.randint(0, 2, size=(N, N))

    def calculateKernel(self):
        K = np.asarray([[1,1,1], [1,0,1], [1,1,1]])
        return K

    def growth(self, U):
        return 0 + (U==3) - ((U<2)|(U>3))

    def updateCell(self, X):
        kernel = self.calculateKernel()
        U = scipy.signal.convolve2d(X, kernel, mode='same', boundary='wrap')
        X = np.clip(X + self.growth(U), 0, 1)
        return X

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
ani = FuncAnimation(fig, conway.animate, frames=200, interval=100, blit=True)

plt.show()