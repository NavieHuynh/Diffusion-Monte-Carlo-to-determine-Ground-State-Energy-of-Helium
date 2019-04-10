import numpy as np
def calc_magnitude(x1,x2,x3):
            return np.sqrt((x1**2  + x2**2 + x3**2))

def metropolis_walk(x):
    rmax =.5
    return x + rmax*(np.random.rand()- .5)       

class ElectronSystem:
    def __init__(self,position = [1,1,1,-1,-1,-1]):
        self.position = position
        self.calc_positionMagnitude()
        self.walk(metropolis_walk)

    def calc_positionMagnitude(self):
        self.r1 = calc_magnitude(self.position[0],self.position[1],self.position[2])
        self.r2 = calc_magnitude(self.position[3],self.position[4],self.position[5])
        self.r12 = calc_magnitude((self.position[3]-self.position[0]), (self.position[4]-self.position[1]) , (self.position[5]-self.position[2]))      
      
    def calc_local_energy(self):
        return ( -4 + (1/self.r12 ))

    def update_position(self):
        self.position = self.nextposition
        self.calc_positionMagnitude()

    def walk(self, walk_func = metropolis_walk):
        self.nextposition = list(map(walk_func , self.position))
        self.r1new = calc_magnitude(self.nextposition[0],self.nextposition[1],self.nextposition[2])
        self.r2new = calc_magnitude(self.nextposition[3],self.nextposition[4],self.nextposition[5])

    def calc_dmc_cutoff_value(self,reference_energy, tau):
        self.update_position()
        life_or_death_number = np.exp(tau*(reference_energy - self.calc_local_energy()) )
        return np.floor( life_or_death_number + np.random.rand())