import math
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from common_utils.quark_relaxation_times_data import QuarkRelaxationTimesData
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData
from common_utils.thermal_conductivity import (
    simplified_thermal_conductivity_integrand, 
    check_temperature_ranges,
    calculate_enthalpy_quark_dens_ratio
)
from common_utils.physical_constants import (
    hbarc_gevfm, 
    up_quark_electric_charge, 
    down_quark_electric_charge, 
    strange_quark_electric_charge
)
from common_utils.io_utils import save_columns_to_file


def simplified_seebeck_sigmae_product_integral(
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
    
    result = result_1 - (eta*enthalpy_quark_dens_ratio)*result_2

    return result


def calculate_seebeck_sigmae_product(
    quark_rel_times_data: QuarkRelaxationTimesData, 
    thermodynamics_data: FixedChemPotTempData, 
    epsabs: float,
    epsrel: float,
    number_of_colors: float = 3.0,
    thermo_interpolation_kind: str = 'linear'
) -> np.ndarray:    
    # first check if temperature of quark_rel_times_data are within the values of min and max values of thermodynamics_data
    check_temperature_ranges(quark_rel_times_data, thermodynamics_data)
    
    # calculate and interpolate the enthalpy to total quark density ratio
    ratio = calculate_enthalpy_quark_dens_ratio(thermodynamics_data)
    enthalpy_quark_dens_ratio_inter = interp1d(
        thermodynamics_data.get_temperature(), 
        ratio, 
        kind=thermo_interpolation_kind,
        bounds_error=True
    )
    
    seebeck_sigmae_product = [] 
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
        
        seebeck_sigmae_integral_up_quark = simplified_seebeck_sigmae_product_integral(
            eff_mass_up_quark, 
            +1, 
            eff_chem_pot_up_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        seebeck_sigmae_integral_down_quark = simplified_seebeck_sigmae_product_integral(
            eff_mass_down_quark, 
            +1, 
            eff_chem_pot_down_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        seebeck_sigmae_integral_strange_quark = simplified_seebeck_sigmae_product_integral(
            eff_mass_strange_quark, 
            +1, 
            eff_chem_pot_strange_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )
        
        seebeck_sigmae_integral_up_antiquark = simplified_seebeck_sigmae_product_integral(
            eff_mass_up_quark, 
            -1, 
            eff_chem_pot_up_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        seebeck_sigmae_integral_down_antiquark = simplified_seebeck_sigmae_product_integral(
            eff_mass_down_quark, 
            -1, 
            eff_chem_pot_down_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        seebeck_sigmae_integral_strange_antiquark = simplified_seebeck_sigmae_product_integral(
            eff_mass_strange_quark, 
            -1, 
            eff_chem_pot_strange_quark, 
            temperature,
            ratio,
            epsabs,
            epsrel
        )

        coefficient = ( (2.0*number_of_colors)/(3.0*temperature**2) )*( (4.0*math.pi)/( (2.0*math.pi)**3 ) )

        seebeck_sigmae_up_quark = coefficient*(up_quark_electric_charge)*tau_up_quark*seebeck_sigmae_integral_up_quark
        seebeck_sigmae_down_quark = coefficient*(down_quark_electric_charge)*tau_down_quark*seebeck_sigmae_integral_down_quark
        seebeck_sigmae_strange_quark = coefficient*(strange_quark_electric_charge)*tau_strange_quark*seebeck_sigmae_integral_strange_quark
        seebeck_sigmae_up_antiquark = coefficient*(-up_quark_electric_charge)*tau_up_antiquark*seebeck_sigmae_integral_up_antiquark
        seebeck_sigmae_down_antiquark = coefficient*(-down_quark_electric_charge)*tau_down_antiquark*seebeck_sigmae_integral_down_antiquark
        seebeck_sigmae_strange_antiquark = coefficient*(-strange_quark_electric_charge)*tau_strange_antiquark*seebeck_sigmae_integral_strange_antiquark
        
        seebeck_sigmae_product.append(
            seebeck_sigmae_up_quark + 
            seebeck_sigmae_down_quark + 
            seebeck_sigmae_strange_quark + 
            seebeck_sigmae_up_antiquark + 
            seebeck_sigmae_down_antiquark + 
            seebeck_sigmae_strange_antiquark
        )
    
    return np.array(seebeck_sigmae_product) 


class SeebeckSigmaeProduct:
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
        
        self.seebeck_sigmae_product = calculate_seebeck_sigmae_product(
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
            self.seebeck_sigmae_product
        ]
        column_labels = [
            "T[GeV]",
            "seebeck_sigmae[GeV^2]"
        ]
        save_columns_to_file(
            output_data_filepath, 
            25, 
            15, 
            column_data, 
            column_labels
        )
