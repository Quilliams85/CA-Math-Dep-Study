import numpy as np
import matplotlib.pyplot as plt
 

N = 50
max_iter = 100
iter = 16

def IX(x, y):
    return (x,y)
   # return x + y * N

def diffuse(b, x, x0, diff, dt):
    a = dt * diff * (N - 2) * (N - 2)
    lin_solve(b, x, x0, a, 1 + 6 * a)

def project(velocX, velocY, p, div):
    # Calculate the divergence of the velocity field and set initial pressure to zero
    for j in range(1, N - 1):
        for i in range(1, N - 1):
            div[IX(i, j)] = (
                -0.5 * (
                    velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] +
                    velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]
                )
            ) / N
            p[IX(i, j)] = 0

    # Apply boundary conditions to div and p
    set_bnd(0, div)
    set_bnd(0, p)
    
    # Solve the linear system to compute the pressure correction
    lin_solve(0, p, div, 1, 6)

    # Subtract the gradient of the pressure field from the velocity field to correct it
    for j in range(1, N - 1):
        for i in range(1, N - 1):
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N

    # Apply boundary conditions to the corrected velocity fields
    set_bnd(1, velocX)
    set_bnd(2, velocY)



def advect(b, d, d0, velocX, velocY, dt):
    dtx = dt * (N - 2)
    dty = dt * (N - 2)
    Nfloat = N - 2

    # Loop over the grid, skipping the boundaries
    for j in range(1, N - 1):
        jfloat = j
        for i in range(1, N - 1):
            ifloat = i

            # Compute the positions from which to interpolate
            tmp1 = dtx * velocX[IX(i, j)]
            tmp2 = dty * velocY[IX(i, j)]
            x = ifloat - tmp1
            y = jfloat - tmp2

            # Clamp x and y within the bounds of the grid
            x = max(0.5, min(Nfloat + 0.5, x))
            y = max(0.5, min(Nfloat + 0.5, y))

            # Calculate integer indices and interpolation weights
            i0 = np.floor(x)
            i1 = i0 + 1
            j0 = np.floor(y)
            j1 = j0 + 1

            s1 = x - i0
            s0 = 1.0 - s1
            t1 = y - j0
            t0 = 1.0 - t1

            # Convert indices to integers for array access
            i0 = int(i0)
            i1 = int(i1)
            j0 = int(j0)
            j1 = int(j1)

            # Perform bilinear interpolation to update the field value
            d[IX(i, j)] = (
                s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)])
            )

    # Apply boundary conditions to the updated field
    set_bnd(b, d)


def set_bnd(b, x):
    # Apply boundary conditions along the top and bottom edges
    for i in range(1, N - 1):
        x[IX(i, 0)] = -x[IX(i, 1)] if b == 2 else x[IX(i, 1)]  # Top edge
        x[IX(i, N - 1)] = -x[IX(i, N - 2)] if b == 2 else x[IX(i, N - 2)]  # Bottom edge
    
    # Apply boundary conditions along the left and right edges
    for j in range(1, N - 1):
        x[IX(0, j)] = -x[IX(1, j)] if b == 1 else x[IX(1, j)]  # Left edge
        x[IX(N - 1, j)] = -x[IX(N - 2, j)] if b == 1 else x[IX(N - 2, j)]  # Right edge

    # Set corner values as the average of adjacent edges
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)])  # Top-left corner
    x[IX(0, N - 1)] = 0.5 * (x[IX(1, N - 1)] + x[IX(0, N - 2)])  # Bottom-left corner
    x[IX(N - 1, 0)] = 0.5 * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)])  # Top-right corner
    x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)])  # Bottom-right corner


def lin_solve(b, x, x0, a, c):
    cRecip = 1.0 / c  # Calculate the reciprocal of c
    for t in range(iter):  # Repeat for the number of iterations
        for j in range(1, N - 1):  # Iterate over rows, skipping boundaries
            for i in range(1, N - 1):  # Iterate over columns, skipping boundaries
                # Update x based on neighboring values and x0
                x[IX(i, j)] = (
                    x0[IX(i, j)] +
                    a * (
                        x[IX(i + 1, j)] +
                        x[IX(i - 1, j)] +
                        x[IX(i, j + 1)] +
                        x[IX(i, j - 1)]
                    )
                ) * cRecip
        set_bnd(b, x)  # Apply boundary conditions




class fluid():
    def __init__(self):
        self.size = N
        self.diffusion = 0
        self.viscosity = 0.00001
        self.dt = 0.1
        self.step = 0

        self.Vx = np.zeros((N, N))
        self.Vy = np.zeros((N, N))

        self.Vy0 = np.zeros((N, N))
        self.Vx0 = np.zeros((N, N))

        self.s = np.zeros((N, N))
        self.density = np.zeros((N, N))

        plt.style.use("dark_background")
        plt.figure(figsize=(5, 5), dpi=160)


    def simStep(self):
        self.step += 1
        print(self.step)
        diffuse(1, self.Vx0, self.Vx, self.viscosity, self.dt)
        diffuse(2, self.Vy0, self.Vy, self.viscosity, self.dt)

        project(self.Vx0, self.Vy0, self.Vx, self.Vy)

        advect(1, self.Vx, self.Vx0, self.Vx0, self.Vy0, self.dt)
        advect(2, self.Vy, self.Vy0, self.Vx0, self.Vy0, self.dt)

        project(self.Vx, self.Vy, self.Vx0, self.Vy0)
        diffuse(0, self.s, self.density, self.diffusion, self.dt)
        advect(0, self.density, self.s, self.Vx, self.Vy, self.dt)

    def addDensity(self, x, y, amount):
        index = IX(x,y)
        self.density[index] =+ amount

    def addVelocity(self, x, y, direction, amount):
        self.Vy[(x,y)] = amount * np.sin(direction)
        self.Vx[(x,y)] = amount * np.cos(direction)

    def render(self):
        plt.contourf(self.density, levels = 100)
        plt.draw()
        plt.pause(0.0001)
        plt.clf()


newFluid = fluid()

newFluid.addDensity(range(10,20),range(10,20), 20)


for i in range(max_iter):
    newFluid.addVelocity(10, 10, np.pi/4, 1)
    newFluid.simStep()
    newFluid.render()
plt.show()
print(newFluid.Vx)