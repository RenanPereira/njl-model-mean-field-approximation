import math


def energy(mass: float, momentum: float) -> float:
    energy = math.sqrt( mass**2 + momentum**2 )

    return energy


def fermi_distribution(energy: float, temperature: float) -> float:
    """
    Fermi-Dirac distribution: n = 1 / (exp(E/T) + 1)
    Safe for large E/T to avoid overflow.
    """
    if temperature <= 0.0:
        raise ValueError("Temperature must be positive and non zero!")

    x = energy/temperature

    if x > 700:
        return 0.0
    elif x < -700:
        return 1.0

    n_fermi = 1.0/( math.exp(x) + 1.0 )
    
    return n_fermi
