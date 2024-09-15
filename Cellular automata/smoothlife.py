import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy
from scipy.spatial import distance


N=400


size = 64;  mid = size // 2;  scale = 1;  cx, cy = 20, 20; R = 10
C = np.asarray([[0,0,0,0,0,0,0.1,0.14,0.1,0,0,0.03,0.03,0,0,0.3,0,0,0,0], [0,0,0,0,0,0.08,0.24,0.3,0.3,0.18,0.14,0.15,0.16,0.15,0.09,0.2,0,0,0,0], [0,0,0,0,0,0.15,0.34,0.44,0.46,0.38,0.18,0.14,0.11,0.13,0.19,0.18,0.45,0,0,0], [0,0,0,0,0.06,0.13,0.39,0.5,0.5,0.37,0.06,0,0,0,0.02,0.16,0.68,0,0,0], [0,0,0,0.11,0.17,0.17,0.33,0.4,0.38,0.28,0.14,0,0,0,0,0,0.18,0.42,0,0], [0,0,0.09,0.18,0.13,0.06,0.08,0.26,0.32,0.32,0.27,0,0,0,0,0,0,0.82,0,0], [0.27,0,0.16,0.12,0,0,0,0.25,0.38,0.44,0.45,0.34,0,0,0,0,0,0.22,0.17,0], [0,0.07,0.2,0.02,0,0,0,0.31,0.48,0.57,0.6,0.57,0,0,0,0,0,0,0.49,0], [0,0.59,0.19,0,0,0,0,0.2,0.57,0.69,0.76,0.76,0.49,0,0,0,0,0,0.36,0], [0,0.58,0.19,0,0,0,0,0,0.67,0.83,0.9,0.92,0.87,0.12,0,0,0,0,0.22,0.07], [0,0,0.46,0,0,0,0,0,0.7,0.93,1,1,1,0.61,0,0,0,0,0.18,0.11], [0,0,0.82,0,0,0,0,0,0.47,1,1,0.98,1,0.96,0.27,0,0,0,0.19,0.1], [0,0,0.46,0,0,0,0,0,0.25,1,1,0.84,0.92,0.97,0.54,0.14,0.04,0.1,0.21,0.05], [0,0,0,0.4,0,0,0,0,0.09,0.8,1,0.82,0.8,0.85,0.63,0.31,0.18,0.19,0.2,0.01], [0,0,0,0.36,0.1,0,0,0,0.05,0.54,0.86,0.79,0.74,0.72,0.6,0.39,0.28,0.24,0.13,0], [0,0,0,0.01,0.3,0.07,0,0,0.08,0.36,0.64,0.7,0.64,0.6,0.51,0.39,0.29,0.19,0.04,0], [0,0,0,0,0.1,0.24,0.14,0.1,0.15,0.29,0.45,0.53,0.52,0.46,0.4,0.31,0.21,0.08,0,0], [0,0,0,0,0,0.08,0.21,0.21,0.22,0.29,0.36,0.39,0.37,0.33,0.26,0.18,0.09,0,0,0], [0,0,0,0,0,0,0.03,0.13,0.19,0.22,0.24,0.24,0.23,0.18,0.13,0.05,0,0,0,0], [0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.09,0.07,0.05,0.01,0,0,0,0,0]])
C = scipy.ndimage.zoom(C, scale, order=0);  R *= scale


class Automaton():
    def __init__(self):
        self.size = N
        self.grid = np.zeros((self.size, self.size))
        self.last_grid = np.zeros((self.size, self.size))
        self.dt = 0.1
        self.bell = lambda x, m, s: np.exp(-((x-m)/s)**2 / 2)
        self.grid = self.generate_voronoi_noise(10)
        #self.grid = np.random.rand(self.size, self.size)
        self.kernel = self.calculateKernel(N, 10)

    def generate_voronoi_noise(self, num_points):
        size = self.size
        num_points = 100
        x, y = np.indices((size, size))
        grid = np.stack((x, y), axis=-1)
        points = np.random.rand(num_points, 2) * size  # Random points scaled to grid size
        distances = np.array([distance.cdist(grid.reshape(-1, 2), [p]) for p in points])
        min_distances = np.min(distances, axis=0).reshape(size, size)
        normalized_distances = (min_distances - min_distances.min()) / (min_distances.max() - min_distances.min())

        return normalized_distances


    def calculateKernel(self, size, radius, mu=0.5, sigma=0.15):
        y, x = np.ogrid[-size//2:size//2, -size//2:size//2]
        distance = np.sqrt(x**2 + y**2) / radius
        bell = np.exp(-((distance - mu) ** 2) / (2 * sigma ** 2))
        #kernel = ((distance < 1) * bell)
        kernel = (np.floor(distance)+1) * (distance<1) * (distance>0.5)
        kernel /= np.sum(kernel)
        fig, ax = plt.subplots()
        kn = ax.imshow(kernel, cmap='plasma')

        self.fK = np.fft.fft2(np.fft.fftshift(kernel / np.sum(kernel)))


        return kernel

    def growth(self, U):
        current = np.sum(U)
        total = U.size
        #return 2* np.exp(-1*(U**2)) - 1
        return (0 + ((U>=0.12)&(U<=0.15)) - ((U<0.12)|(U>0.15)))
        #return self.bell(U, 0.5, 0.15)*2-1
        #- 2*(current/total

    def updateCell(self, X):
        #U = scipy.signal.convolve2d(X, self.kernel, mode='same', boundary='wrap')
        U = np.real(np.fft.ifft2(self.fK * np.fft.fft2(X)))
        X = np.clip(X + self.dt * self.growth(U), 0, 1)
        return X

    def animate(self, i):
        self.last_grid = self.grid
        self.grid = self.updateCell(self.grid)
        im.set_array(self.grid)  # Update the image data
        return [im]


# Initial setup
conway = Automaton()
array = conway.grid

# Plot setup
fig, ax = plt.subplots()
im = ax.imshow(array, cmap='plasma')  # Display array using binary color map

# Create animation
ani = FuncAnimation(fig, conway.animate, frames=200, interval=1, blit=True)

plt.show()