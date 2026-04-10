import math
import numpy as np
from scipy.integrate import quad

from common_utils.quark_relaxation_times_data import QuarkRelaxationTimesData
from common_utils.general_physics import energy, fermi_distribution
from common_utils.physical_constants import hbarc_gevfm
from common_utils.io_utils import save_columns_to_file


def simplified_shear_viscosity_integrand(
    momentum: float, 
    mass: float, 
    eta: float, 
    chemical_potential: float, 
    temperature: float
) -> float:
    """
    Shear viscosity integrand for a single particle species:

        integrand = ( p^6/E^2 ) * n_fermi * ( 1 - n_fermi )

    Arguments:
        momentum: particle momentum
        mass: particle mass
        eta: +1 for particle, -1 for antiparticle
        chemical_potential: mu
        temperature: T

    Returns:
        float: value of integrand
    """
    E = energy(mass, momentum)
    n_fermi = fermi_distribution(E - eta*chemical_potential, temperature)
    shear_integrand = ( (momentum**6)/(E**2) )*n_fermi*( 1.0 - n_fermi )

    return shear_integrand


def simplified_shear_viscosity_integral(
    mass: float, 
    eta: float, 
    chemical_potential: float, 
    temperature: float,
    epsabs: float,
    epsrel: float,
) -> float:
    arguments = (mass, eta, chemical_potential, temperature)
    
    result, error = quad(
        simplified_shear_viscosity_integrand,
        0.0,
        math.inf,
        args=arguments,
        epsabs=epsabs,
        epsrel=epsrel,
    )

    return result


def calculate_shear_viscosity(
    quark_rel_times_data: QuarkRelaxationTimesData, 
    epsabs: float,
    epsrel: float,
    number_of_colors: float = 3.0
) -> np.ndarray:
    shear_viscosity = [] 
    
    for index in range(quark_rel_times_data.size()):

        # get temperature for specific data index
        temperature = quark_rel_times_data.get_temperature()[index]
        if temperature <= 0.0:
            raise ValueError("Temperature must be positive and non zero!")

        # get quark quantities from data
        eff_mass_up_quark = quark_rel_times_data.get_up_quark_effective_mass()[index]
        eff_mass_down_quark = quark_rel_times_data.get_down_quark_effective_mass()[index]
        eff_mass_strange_quark = quark_rel_times_data.get_strange_quark_effective_mass()[index]

        eff_chem_pot_up_quark = quark_rel_times_data.get_up_quark_effective_chemical_potential()[index]
        eff_chem_pot_down_quark = quark_rel_times_data.get_down_quark_effective_chemical_potential()[index]
        eff_chem_pot_strange_quark = quark_rel_times_data.get_strange_quark_effective_chemical_potential()[index]

        # relaxation time is stored in the file in fermi units: convert to GeV^-1
        tau_up_quark = quark_rel_times_data.get_rel_time_up_quark()[index]/hbarc_gevfm
        tau_down_quark = quark_rel_times_data.get_rel_time_down_quark()[index]/hbarc_gevfm
        tau_strange_quark = quark_rel_times_data.get_rel_time_strange_quark()[index]/hbarc_gevfm

        tau_up_antiquark = quark_rel_times_data.get_rel_time_up_antiquark()[index]/hbarc_gevfm
        tau_down_antiquark = quark_rel_times_data.get_rel_time_down_antiquark()[index]/hbarc_gevfm
        tau_strange_antiquark = quark_rel_times_data.get_rel_time_strange_antiquark()[index]/hbarc_gevfm

        eta_integral_up_quark = simplified_shear_viscosity_integral(
            eff_mass_up_quark, 
            +1, 
            eff_chem_pot_up_quark, 
            temperature,
            epsabs,
            epsrel
        )

        eta_integral_down_quark = simplified_shear_viscosity_integral(
            eff_mass_down_quark, 
            +1, 
            eff_chem_pot_down_quark, 
            temperature,
            epsabs,
            epsrel
        )

        eta_integral_strange_quark = simplified_shear_viscosity_integral(
            eff_mass_strange_quark, 
            +1, 
            eff_chem_pot_strange_quark, 
            temperature,
            epsabs,
            epsrel
        )

        eta_integral_up_antiquark = simplified_shear_viscosity_integral(
            eff_mass_up_quark, 
            -1, 
            eff_chem_pot_up_quark, 
            temperature,
            epsabs,
            epsrel
        )

        eta_integral_down_antiquark = simplified_shear_viscosity_integral(
            eff_mass_down_quark, 
            -1, 
            eff_chem_pot_down_quark, 
            temperature,
            epsabs,
            epsrel
        )

        eta_integral_strange_antiquark = simplified_shear_viscosity_integral(
            eff_mass_strange_quark, 
            -1, 
            eff_chem_pot_strange_quark, 
            temperature,
            epsabs,
            epsrel
        )

        coefficient = ( (2.0*number_of_colors)/(15.0*temperature) )*( (4.0*math.pi)/( (2.0*math.pi)**3 ) )

        eta_up_quark = coefficient*tau_up_quark*eta_integral_up_quark
        eta_down_quark = coefficient*tau_down_quark*eta_integral_down_quark
        eta_strange_quark = coefficient*tau_strange_quark*eta_integral_strange_quark
        eta_up_antiquark = coefficient*tau_up_antiquark*eta_integral_up_antiquark
        eta_down_antiquark = coefficient*tau_down_antiquark*eta_integral_down_antiquark
        eta_strange_antiquark = coefficient*tau_strange_antiquark*eta_integral_strange_antiquark

        shear_viscosity.append(
            eta_up_quark + eta_down_quark + eta_strange_quark + eta_up_antiquark + eta_down_antiquark + eta_strange_antiquark
        )
    
    return np.array(shear_viscosity) 


class ShearViscosity:
    def __init__(
        self, 
        input_data_filepath: str, 
        output_data_filepath: str,
        integral_epsabs: float = 1.49e-8,
        integral_epsrel: float = 1.49e-8,
        number_of_colors: float = 3
    ):
        quark_rel_times_data = QuarkRelaxationTimesData(input_data_filepath)

        self.shear_viscosity = calculate_shear_viscosity(
            quark_rel_times_data, 
            integral_epsabs,
            integral_epsrel,
            number_of_colors
        )
        
        self.temperature = quark_rel_times_data.get_temperature()

        self._save_to_file(output_data_filepath)

    def _save_to_file(
        self, 
        output_data_filepath: str, 
    ) -> None:
        column_data = [
            self.temperature, 
            self.shear_viscosity
        ]
        column_labels = [
            "T[GeV]",
            "eta[GeV^3]",
        ]
        save_columns_to_file(
            output_data_filepath, 
            25, 
            15, 
            column_data, 
            column_labels
        )
