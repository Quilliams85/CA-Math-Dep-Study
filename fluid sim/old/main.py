import os
import numpy as np

class Fluid:
    def __init__(self):
        self.addGhostCells()
        self.setGhostCells()
        self.setGhostGradients()



    
    def addGhostCells(rho, vx, vy, P ):
        """
        Add ghost cells to the top and bottom
        rho      is matrix of cell densities
        vx       is matrix of cell x-velocity
        vy       is matrix of cell y-velocity
        P        is matrix of cell pressures
        """
        rho = np.hstack((rho[:,0:1], rho, rho[:,-1:]))
        vx  = np.hstack(( vx[:,0:1],  vx,  vx[:,-1:]))
        vy  = np.hstack(( vy[:,0:1],  vy,  vy[:,-1:]))
        P   = np.hstack((  P[:,0:1],   P,   P[:,-1:]))
        
        return rho, vx, vy, P
    
    
    def setGhostCells( rho, vx, vy, P ):
        """
        Set ghost cells at the top and bottom
        rho         is matrix of cell densities
        vx       is matrix of cell x-velocity
        vy       is matrix of cell y-velocity
        P        is matrix of cell pressures
        """
        
        # copy bottom row to ghosts, negate vy
        rho[:,0]  = rho[:,1]
        vx[:,0]   =  vx[:,1]
        vy[:,0]   = -vy[:,1]
        P[:,0]    =   P[:,1]
        
        # copy top row to ghosts, negate vy
        rho[:,-1] = rho[:,-2]
        vx[:,-1]  =  vx[:,-2]
        vy[:,-1]  = -vy[:,-2]
        P[:,-1]   =   P[:,-2]
        
        return rho, vx, vy, P
    

    def setGhostGradients( f_dx, f_dy ):
        """
        Set ghost cell y-gradients at the top and bottom to be reflections
        f_dx     is a matrix of derivative of f in the x-direction
        f_dy     is a matrix of derivative of f in the y-direction
        """
        # copy top and bottom row (y) gradients to ghosts, negate
        f_dy[:,0]  = -f_dy[:,1]  
        f_dy[:,-1] = -f_dy[:,-2] 
        
        return f_dx, f_dy
    

    def addSourceTerm( Mass, Momx, Momy, Energy, g, dt ):
        """
        Add gravitational source term to conservative variables
        Mass     is matrix of mass in cells
        Momx     is matrix of x-momentum in cells
        Momy     is matrix of y-momentum in cells
        Energy   is matrix of energy in cells
        g        is strength of gravity
        Y        is matrix of y positions of cells
        dt       is timestep to progress solution
        """

        Energy += dt * Momy * g
        Momy += dt * Mass * g

        return Mass, Momx, Momy, Energy  

def main():
    pass

main()