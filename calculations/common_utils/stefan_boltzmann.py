import numpy as np


class StefanBoltzmannMasslessFlavorDegenerateQuarks:
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
        
        s = ( Nc*Nf/(45) )*( 7*pi2*T**3 + 15*T*mu**2 )
        
        return s
    
    def entropy_density_over_temp3(
        self, 
        mu: float | np.ndarray, 
        T: float | np.ndarray
    ) -> float | np.ndarray:
        return self.entropy_density(mu, T)/T**3


class StefanBoltzmannMasslessQuarks:
    def __init__(
        self, 
        number_of_colors: int = 3, 
        number_of_flavors: int = 3
    ):
        self.number_of_colors = number_of_colors
        self.number_of_flavors = number_of_flavors

    def verify_chemical_potential_input(
        self, 
        mu: list[float] | list[np.ndarray]
    ) -> None:
        if (len(mu)!=self.number_of_flavors):
            raise ValueError(f"Expected {self.number_of_flavors} chemical potentials, got {len(mu)}.")
        
        if all(isinstance(m, np.ndarray) for m in mu):
            lengths = []
            for m in mu:
                if (isinstance(m, np.ndarray)):
                    lengths.append(len(m))
            if (min(lengths)!=max(lengths)):
                raise ValueError("All np.ndarray in mu must have the same length.")
        elif any(isinstance(m, np.ndarray) for m in mu):
            if not all(isinstance(m, np.ndarray) for m in mu):
                raise ValueError("mu must be all scalars or all np.ndarray.")
        else:
            if not all(isinstance(m, float) for m in mu):
                raise ValueError("mu must be all scalars or all np.ndarray.")
  
    def entropy_density(
        self, 
        mu: list[float] | list[np.ndarray], 
        T: float | np.ndarray
    ) -> float | np.ndarray:
        Nc = self.number_of_colors
        Nf = self.number_of_flavors
        pi2 = np.pi**2
        
        self.verify_chemical_potential_input(mu)
            
        sum_over_mu2 = 0
        for i in range(0, len(mu)):
            sum_over_mu2 += mu[i]**2

        s = ( Nc/(45) )*( Nf*7*pi2*T**3 + 15*T*sum_over_mu2 )
        
        return s
    
    def entropy_density_over_temp3(
        self, 
        mu: list[float] | list[np.ndarray], 
        T: float | np.ndarray
    ) -> float | np.ndarray:
        return self.entropy_density(mu, T)/T**3
