import math
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from common_utils.quark_relaxation_times_data import QuarkRelaxationTimesData
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData
from common_utils.general_physics import energy, fermi_distribution
from common_utils.physical_constants import hbarc_gevfm
from common_utils.io_utils import save_columns_to_file


def simplified_thermal_conductivity_integrand(
    momentum: float, 
    mass: float, 
    eta: float, 
    chemical_potential: float, 
    temperature: float,
    energy_exponent: float
) -> float:
    """
    Simplified thermal_conductivity integrand for a single particle species:

        integrand = ( p^4/E^p ) * n_fermi * ( 1 - n_fermi )
    
    with p an integer

    Arguments:
        momentum: particle momentum
        mass: particle mass
        eta: +1 for particle, -1 for antiparticle
        chemical_potential: mu
        temperature: T
        energy_exponent: exponent used in the denominator for the energy factor 

    Returns:
        float: value of integrand
    """
    E = energy(mass, momentum)
    n_fermi = fermi_distribution(E - eta*chemical_potential, temperature)
    
    integrand = ( (momentum**4)/(E**energy_exponent) )*n_fermi*( 1.0 - n_fermi )

    return integrand


def simplified_thermal_conductivity_integral(
    mass: float, 
    eta: float, 
    chemical_potential: float, 
    temperature: float,
    enthalpy_quark_dens_ratio: float,
    epsabs: float,
    epsrel: float,
) -> float:
    arguments = (mass, eta, chemical_potential, temperature, 2)
    result_2, error = quad(
        simplified_thermal_conductivity_integrand,
        0.0,
        math.inf,
        args=arguments,
        epsabs=epsabs,
        epsrel=epsrel,
    )
    
    arguments = (mass, eta, chemical_potential, temperature, 1)
    result_1, error = quad(
        simplified_thermal_conductivity_integrand,
        0.0,
        math.inf,
        args=arguments,
        epsabs=epsabs,
        epsrel=epsrel,
    )

    arguments = (mass, eta, chemical_potential, temperature, 0)
    result_0, error = quad(
        simplified_thermal_conductivity_integrand,
        0.0,
        math.inf,
        args=arguments,
        epsabs=epsabs,
        epsrel=epsrel,
    )
    
    x = eta*enthalpy_quark_dens_ratio
    
    result = result_0 - 2*x*result_1 + (x**2)*result_2

    return result


def calculate_enthalpy_quark_dens_ratio(
    thermodynamics_data: FixedChemPotTempData, 
) -> np.ndarray:
    pressure = thermodynamics_data.get_pressure()
    energy_dens = thermodynamics_data.get_energy_density()
    quark_density = thermodynamics_data.get_up_quark_density() + thermodynamics_data.get_down_quark_density() + thermodynamics_data.get_strange_quark_density()
    enthalpy = pressure + energy_dens
    enthalpy_quark_dens_ratio = enthalpy/quark_density
    
    return enthalpy_quark_dens_ratio


def interpolate_enthalpy_quark_dens_ratio(
    thermodynamics_data: FixedChemPotTempData, 
    thermo_interpolation_kind: str = 'linear'
) -> interp1d:
    
    temperature = []
    enthalpy_quark_dens_ratio = []
    for i in range(thermodynamics_data.size()):
        quark_density = thermodynamics_data.get_up_quark_density()[i] + thermodynamics_data.get_down_quark_density()[i] + thermodynamics_data.get_strange_quark_density()[i]
        
        if quark_density>0:
            pressure = thermodynamics_data.get_pressure()[i]
            energy_dens = thermodynamics_data.get_energy_density()[i]
            enthalpy = pressure + energy_dens
            
            temperature.append(thermodynamics_data.get_temperature()[i])
            enthalpy_quark_dens_ratio.append(enthalpy/quark_density)
    
    enthalpy_quark_dens_ratio_inter = interp1d(
        temperature, 
        enthalpy_quark_dens_ratio, 
        kind=thermo_interpolation_kind,
        bounds_error=True
    )
    
    return enthalpy_quark_dens_ratio_inter


def check_temperature_ranges(
    quark_rel_times_data: QuarkRelaxationTimesData, 
    thermodynamics_data: FixedChemPotTempData, 
) -> None:
    temp = thermodynamics_data.get_temperature()
    min_temp = np.min(temp)
    max_temp = np.max(temp)
    
    for temp in quark_rel_times_data.get_temperature():
        if temp<min_temp or temp>max_temp:
            raise ValueError("Temperature values in QuarkRelaxationTimesData outside the temperature values in thermodynamics data!")


def calculate_thermal_conductivity(
    quark_rel_times_data: QuarkRelaxationTimesData, 
    thermodynamics_data: FixedChemPotTempData, 
    epsabs: float,
    epsrel: float,
    number_of_colors: float = 3.0,
    thermo_interpolation_kind: str = 'linear'
) -> np.ndarray:    
    # first check if temperature of quark_rel_times_data are within the values of min and max values of thermodynamics_data
    check_temperature_ranges(quark_rel_times_data, thermodynamics_data)
    
    enthalpy_quark_dens_ratio_inter = interpolate_enthalpy_quark_dens_ratio(
        thermodynamics_data, thermo_interpolation_kind
    )
    
    thermal_conductivity = [] 
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
        
        ratio = enthalpy_quark_dens_ratio_inter(temperature)
        
        kappa_integral_up_quark = simplified_thermal_conductivity_integral(
            eff_mass_up_quark, 
            +1, 
            eff_chem_pot_up_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        kappa_integral_down_quark = simplified_thermal_conductivity_integral(
            eff_mass_down_quark, 
            +1, 
            eff_chem_pot_down_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        kappa_integral_strange_quark = simplified_thermal_conductivity_integral(
            eff_mass_strange_quark, 
            +1, 
            eff_chem_pot_strange_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )
        
        kappa_integral_up_antiquark = simplified_thermal_conductivity_integral(
            eff_mass_up_quark, 
            -1, 
            eff_chem_pot_up_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        kappa_integral_down_antiquark = simplified_thermal_conductivity_integral(
            eff_mass_down_quark, 
            -1, 
            eff_chem_pot_down_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        kappa_integral_strange_antiquark = simplified_thermal_conductivity_integral(
            eff_mass_strange_quark, 
            -1, 
            eff_chem_pot_strange_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        coefficient = ( (2.0*number_of_colors)/(3.0*temperature**2) )*( (4.0*math.pi)/( (2.0*math.pi)**3 ) )

        kappa_up_quark = coefficient*tau_up_quark*kappa_integral_up_quark
        kappa_down_quark = coefficient*tau_down_quark*kappa_integral_down_quark
        kappa_strange_quark = coefficient*tau_strange_quark*kappa_integral_strange_quark
        kappa_up_antiquark = coefficient*tau_up_antiquark*kappa_integral_up_antiquark
        kappa_down_antiquark = coefficient*tau_down_antiquark*kappa_integral_down_antiquark
        kappa_strange_antiquark = coefficient*tau_strange_antiquark*kappa_integral_strange_antiquark
        
        thermal_conductivity.append(
            kappa_up_quark + 
            kappa_down_quark + 
            kappa_strange_quark + 
            kappa_up_antiquark + 
            kappa_down_antiquark + 
            kappa_strange_antiquark
        )
    
    return np.array(thermal_conductivity) 


class ThermalConductivity:
    def __init__(
        self, 
        input_quark_rel_times_data_filepath: str, 
        input_thermodynamics_data_filepath: str,
        output_data_filepath: str,
        thermo_interpolation_kind: str = 'linear',
        integral_epsabs: float = 1.49e-8,
        integral_epsrel: float = 1.49e-8,
        number_of_colors: float = 3
    ):
        quark_rel_times_data = QuarkRelaxationTimesData(input_quark_rel_times_data_filepath)
        thermodynamics_data = FixedChemPotTempData(input_thermodynamics_data_filepath)
        
        self.thermal_conductivity = calculate_thermal_conductivity(
            quark_rel_times_data, 
            thermodynamics_data,
            integral_epsabs,
            integral_epsrel,
            number_of_colors,
            thermo_interpolation_kind
        )
        
        self.temperature = quark_rel_times_data.get_temperature()

        self._save_to_file(output_data_filepath)

    def _save_to_file(
        self, 
        output_data_filepath: str, 
    ) -> None:
        column_data = [
            self.temperature, 
            self.thermal_conductivity
        ]
        column_labels = [
            "T[GeV]",
            "kappa[GeV^2]"
        ]
        save_columns_to_file(
            output_data_filepath, 
            25, 
            15, 
            column_data, 
            column_labels
        )
