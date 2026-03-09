import numpy as np


class StefanBoltzmannMasslessQuarks:
    def __init__(
        self, 
        number_of_colors: int = 3, 
        number_of_flavors: int = 3
    ):
        self.number_of_colors = number_of_colors
        self.number_of_flavors = number_of_flavors
        
    def entropy_density(
        self, 
        mu: float | np.ndarray, 
        T: float | np.ndarray
    ) -> float | np.ndarray:
        Nc = self.number_of_colors
        Nf = self.number_of_flavors
        pi2 = np.pi**2
        
        s = (Nc*Nf/(45))*( 7*pi2*T**3 + 15*T*mu**2 )
        
        return s
    
    def entropy_density_over_temp3(
        self, 
        mu: float | np.ndarray, 
        T: float | np.ndarray
    ) -> float | np.ndarray:
        return self.entropy_density(mu, T)/T**3
