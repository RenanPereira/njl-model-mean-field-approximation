import numpy as np

from common_utils.integrated_cross_section_data import IntegratedCrossSectionData
from common_utils.physical_constants import hbarc_gevfm
from common_utils.io_utils import save_columns_to_file


class QuarkRelaxationTimes:
    def __init__(
        self, 
        path_input_data_folder: str, 
        parameter_set: str, 
        method: str, 
        physical_scenario: str,
        output_data_filepath: str
    ):
        print(f"Calculating quark relaxation times using integrated cross section files in found in path: {path_input_data_folder}")
        
        if (physical_scenario=="zero_chemical_potential_isospin_symmetric"):
            self._load_data_zero_mu_isospin_symmetric(path_input_data_folder, parameter_set, method)
        elif (physical_scenario=="finite_chemical_potential_isospin_symmetric"):
            self._load_data_finite_mu_isospin_symmetric(path_input_data_folder, parameter_set, method)
        else:
            raise ValueError(f"Unknown physical_scenario: {physical_scenario}")
        
        self._test_data_mismatch()
        self._calculate_quark_relaxation_times()
        self._save_to_file(output_data_filepath)

    def _load_data_zero_mu_isospin_symmetric(
        self, 
        path_data_folder: str, 
        parameter_set: str, 
        method: str
    ) -> None:
        # read data files for the specific scenario of zero chemical potential and isospin symmetry
        loader = IntegratedCrossSectionData.from_matching_files
        
        # 1) ud->ud
        self.data_udud = loader(path_data_folder, parameter_set, "UDUD", method)

        # 2) us->us
        self.data_usus = loader(path_data_folder, parameter_set, "USUS", method)

        # 3) ds->ds (equal to us->us by isospin symmetry)
        self.data_dsds = self.data_usus

        # 4) uu->uu
        self.data_uuuu = loader(path_data_folder, parameter_set, "UUUU", method)

        # 5) dd->dd    (equal to uu->uu by isospin symmetry)
        self.data_dddd = self.data_uuuu

        # 6) ss->ss
        self.data_ssss = loader(path_data_folder, parameter_set, "SSSS", method)

        # 7) udbar->udbar
        self.data_udbarudbar = loader(path_data_folder, parameter_set, "UDBarUDBar", method)

        # 8) dubar->dubar (equal to udbar-> udbar by isospin symmetry)
        self.data_dubardubar = self.data_udbarudbar

        # 9) usbar->usbar
        self.data_usbarusbar = loader(path_data_folder, parameter_set, "USBarUSBar", method)

        # 10) subar->subar (equal to usbar->usbar by charge symmetry at CP=0)
        self.data_subarsubar = self.data_usbarusbar

        # 11) dsbar->dsbar (equal to usbar->usbar by isospin symmetry)
        self.data_dsbardsbar = self.data_usbarusbar

        # 12) sdbar->sdbar (equal to dsbar->dsbar by charge symmetry at CP=0)
        self.data_sdbarsdbar = self.data_dsbardsbar

        # 13) uubar->uubar
        self.data_uubaruubar = loader(path_data_folder, parameter_set, "UUBarUUBar", method)

        # 14) uubar->ddbar
        self.data_uubarddbar = loader(path_data_folder, parameter_set, "UUBarDDBar", method)

        # 15) uubar-> ssbar
        self.data_uubarssbar = loader(path_data_folder, parameter_set, "UUBarSSBar", method)

        # 16) ddbar->uubar (equal to uubar->ddbar by isospin symmetry)
        self.data_ddbaruubar = self.data_uubarddbar

        # 17) ddbar->ddbar (equal to uubar->uubar by isospin symmetry)
        self.data_ddbarddbar = self.data_uubaruubar

        # 18) ddbar->ssbar (equal to uubar->ssbar by isospin symmetry)
        self.data_ddbarssbar = self.data_uubarssbar

        # 19) ssbar->uubar
        self.data_ssbaruubar = loader(path_data_folder, parameter_set, "SSBarUUBar", method)

        # 20) ssbar->ddbar (equal to ssbar->uubar by isospin symmetry)
        self.data_ssbarddbar = self.data_ssbaruubar

        # 21) ssbar->ssbar
        self.data_ssbarssbar = loader(path_data_folder, parameter_set, "SSBarSSBar", method)

        # 22) ubardbar->ubardbar (equal to ud->ud by charge symmetry at CP=0)
        self.data_ubardbarubardbar = self.data_udud

        # 23) ubarsbar-> ubarsbar (equal to us->us by charge symmetry at CP=0)
        self.data_ubarsbarubarsbar = self.data_usus

        # 24) dbarsbar->dbarsbar (equal to ds->ds by charge symmetry at CP=0)
        self.data_dbarsbardbarsbar = self.data_dsds

        # 25) ubarubar->ubarubar (equal to uu->uu by charge symmetry at CP=0)
        self.data_ubarubarubarubar = self.data_uuuu

        # 26) dbardbar->dbardbar (equal to dd->dd by charge symmetry at CP=0)
        self.data_dbardbardbardbar = self.data_dddd

        # 27) sbarsbar->sbarsbar (equal to ss->ss by charge symmetry at CP=0)
        self.data_sbarsbarsbarsbar = self.data_ssss

    def _load_data_finite_mu_isospin_symmetric(
        self, 
        path_data_folder: str, 
        parameter_set: str, 
        method: str
    ) -> None:
        # read data files for the specific scenario of finite chemical potential and isospin symmetry
        loader = IntegratedCrossSectionData.from_matching_files
        
        # 1) ud->ud
        self.data_udud = loader(path_data_folder, parameter_set, "UDUD", method)

        # 2) us->us
        self.data_usus = loader(path_data_folder, parameter_set, "USUS", method)

        # 3) ds->ds (equal to us->us by isospin symmetry)
        self.data_dsds = self.data_usus

        # 4) uu->uu
        self.data_uuuu = loader(path_data_folder, parameter_set, "UUUU", method)

        # 5) dd->dd    (equal to uu->uu by isospin symmetry)
        self.data_dddd = self.data_uuuu

        # 6) ss->ss
        self.data_ssss = loader(path_data_folder, parameter_set, "SSSS", method)

        # 7) udbar->udbar
        self.data_udbarudbar = loader(path_data_folder, parameter_set, "UDBarUDBar", method)

        # 8) dubar->dubar (equal to udbar-> udbar by isospin symmetry)
        self.data_dubardubar = self.data_udbarudbar

        # 9) usbar->usbar
        self.data_usbarusbar = loader(path_data_folder, parameter_set, "USBarUSBar", method)

        # 10) subar->subar
        self.data_subarsubar = loader(path_data_folder, parameter_set, "SUBarSUBar", method)

        # 11) dsbar->dsbar (equal to usbar->usbar by isospin symmetry)
        self.data_dsbardsbar = self.data_usbarusbar

        # 12) sdbar->sdbar (equal to subar->subar by isospin symmetry)
        self.data_sdbarsdbar = self.data_subarsubar

        # 13) uubar->uubar
        self.data_uubaruubar = loader(path_data_folder, parameter_set, "UUBarUUBar", method)

        # 14) uubar->ddbar
        self.data_uubarddbar = loader(path_data_folder, parameter_set, "UUBarDDBar", method)

        # 15) uubar-> ssbar
        self.data_uubarssbar = loader(path_data_folder, parameter_set, "UUBarSSBar", method)

        # 16) ddbar->uubar (equal to uubar->ddbar by isospin symmetry)
        self.data_ddbaruubar = self.data_uubarddbar

        # 17) ddbar->ddbar (equal to uubar->uubar by isospin symmetry)
        self.data_ddbarddbar = self.data_uubaruubar

        # 18) ddbar->ssbar (equal to uubar->ssbar by isospin symmetry)
        self.data_ddbarssbar = self.data_uubarssbar

        # 19) ssbar->uubar
        self.data_ssbaruubar = loader(path_data_folder, parameter_set, "SSBarUUBar", method)

        # 20) ssbar->ddbar (equal to ssbar->uubar by isospin symmetry)
        self.data_ssbarddbar = self.data_ssbaruubar

        # 21) ssbar->ssbar
        self.data_ssbarssbar = loader(path_data_folder, parameter_set, "SSBarSSBar", method)

        # 22) ubardbar->ubardbar
        self.data_ubardbarubardbar = loader(path_data_folder, parameter_set, "UBarDBarUBarDBar", method)

        # 23) ubarsbar-> ubarsbar
        self.data_ubarsbarubarsbar = loader(path_data_folder, parameter_set, "UBarSBarUBarSBar", method)

        # 24) dbarsbar->dbarsbar (equal to ubarsbar->ubarsbar by isospin symmetry)
        self.data_dbarsbardbarsbar = self.data_ubarsbarubarsbar

        # 25) ubarubar->ubarubar
        self.data_ubarubarubarubar = loader(path_data_folder, parameter_set, "UBarUBarUBarUBar", method)

        # 26) dbardbar->dbardbar (equal to ubarubar->ubarubar by isospin symmetry)
        self.data_dbardbardbardbar = self.data_ubarubarubarubar

        # 27) sbarsbar->sbarsbar
        self.data_sbarsbarsbarsbar = loader(path_data_folder, parameter_set, "SBarSBarSBarSBar", method)

    def _test_data_mismatch(
        self,
        equal_entries_rel_tolerance: float = 1e-8, 
        equal_entries_abs_tolerance: float = 1e-8
    ) -> None:
        # organize all data into list for easier testing
        data_list = [
            self.data_udud,
            self.data_usus,
            self.data_dsds,
            self.data_uuuu,
            self.data_dddd,
            self.data_ssss,
            self.data_udbarudbar,
            self.data_dubardubar,
            self.data_usbarusbar,
            self.data_subarsubar,
            self.data_dsbardsbar,
            self.data_sdbarsdbar,
            self.data_uubaruubar,
            self.data_uubarddbar,
            self.data_uubarssbar,
            self.data_ddbaruubar,
            self.data_ddbarddbar,
            self.data_ddbarssbar,
            self.data_ssbaruubar,
            self.data_ssbarddbar,
            self.data_ssbarssbar,
            self.data_ubardbarubardbar,
            self.data_ubarsbarubarsbar,
            self.data_dbarsbardbarsbar,
            self.data_ubarubarubarubar,
            self.data_dbardbardbardbar,
            self.data_sbarsbarsbarsbar,
        ]
        
        if len(data_list) < 2:
            return

        # check data if data lengths are all the same
        data_lengths = []
        for data in data_list:
            data_lengths.append(len(data.get_data()))

        if (min(data_lengths) != max(data_lengths)):
            raise ValueError("Data files with different sizes!")

        checks = [
            ("temperature", lambda data: data.get_temperature()),
            ("up quark number", lambda data: data.get_up_quark_number()),
            ("down quark number", lambda data: data.get_down_quark_number()),
            ("strange quark number", lambda data: data.get_strange_quark_number()),
            ("up antiquark number", lambda data: data.get_up_antiquark_number()),
            ("down antiquark number", lambda data: data.get_down_antiquark_number()),
            ("strange antiquark number", lambda data: data.get_strange_antiquark_number()),
            ("up quark effective mass", lambda data: data.get_up_quark_effective_mass()),
            ("down quark effective mass", lambda data: data.get_down_quark_effective_mass()),
            ("strange quark effective mass", lambda data: data.get_strange_quark_effective_mass()),
            ("up quark effective chemical potential", lambda data: data.get_up_quark_effective_chemical_potential()),
            ("down quark effective chemical potential", lambda data: data.get_down_quark_effective_chemical_potential()),
            ("strange quark effective chemical potential", lambda data: data.get_strange_quark_effective_chemical_potential()),
        ]
        
        for name, getter in checks:
            for i in range(len(data_list) - 1):
                if not np.allclose(
                    getter(data_list[i]),
                    getter(data_list[i + 1]),
                    rtol=equal_entries_rel_tolerance,
                    atol=equal_entries_abs_tolerance
                ):
                    raise ValueError(
                        f"Problem found! {name.capitalize()} is not "
                        f"the same in all data! Mismatch found between data indexes {i} and {i+1}"
                    )

    def _calculate_quark_relaxation_times(self) -> None:

        # Extract quark numbers from single data source
        nu = self.data_uuuu.get_up_quark_number()
        nd = self.data_uuuu.get_down_quark_number()
        ns = self.data_uuuu.get_strange_quark_number()

        nubar = self.data_uuuu.get_up_antiquark_number()
        ndbar = self.data_uuuu.get_down_antiquark_number()
        nsbar = self.data_uuuu.get_strange_antiquark_number()

        # Calculate the up quark relaxation time

        sigma_uuuu_nu = self.data_uuuu.get_integrated_cross_section()*nu
        sigma_udud_nd = self.data_udud.get_integrated_cross_section()*nd
        sigma_usus_ns = self.data_usus.get_integrated_cross_section()*ns

        sigma_uubaruubar_nubar = self.data_uubaruubar.get_integrated_cross_section()*nubar
        sigma_uubarddbar_nubar = self.data_uubarddbar.get_integrated_cross_section()*nubar
        sigma_uubarssbar_nubar = self.data_uubarssbar.get_integrated_cross_section()*nubar
        sigma_udbarudbar_ndbar = self.data_udbarudbar.get_integrated_cross_section()*ndbar
        sigma_usbarusbar_nsbar = self.data_usbarusbar.get_integrated_cross_section()*nsbar

        # relaxation time is in fermis
        self.tau_up_quark = hbarc_gevfm*1.0/(
            sigma_uuuu_nu +
            sigma_udud_nd +
            sigma_usus_ns +
            sigma_uubaruubar_nubar +
            sigma_uubarddbar_nubar +
            sigma_uubarssbar_nubar +
            sigma_udbarudbar_ndbar +
            sigma_usbarusbar_nsbar
        )

        # Calculate the down quark relaxation time

        sigma_dudu_nu = self.data_udud.get_integrated_cross_section()*nu
        sigma_dddd_nd = self.data_dddd.get_integrated_cross_section()*nd
        sigma_dsds_ns = self.data_dsds.get_integrated_cross_section()*ns

        sigma_dubardubar_nubar = self.data_dubardubar.get_integrated_cross_section()*nubar
        sigma_ddbaruubar_ndbar = self.data_ddbaruubar.get_integrated_cross_section()*ndbar
        sigma_ddbarddbar_ndbar = self.data_ddbarddbar.get_integrated_cross_section()*ndbar
        sigma_ddbarssbar_ndbar = self.data_ddbarssbar.get_integrated_cross_section()*ndbar
        sigma_dsbardsbar_nsbar = self.data_dsbardsbar.get_integrated_cross_section()*nsbar

        # relaxation time is in fermis
        self.tau_down_quark = hbarc_gevfm*1.0/(
            sigma_dudu_nu  +
            sigma_dddd_nd +
            sigma_dsds_ns +
            sigma_dubardubar_nubar +
            sigma_ddbaruubar_ndbar +
            sigma_ddbarddbar_ndbar +
            sigma_ddbarssbar_ndbar +
            sigma_dsbardsbar_nsbar
        )

        # Calculate the strange quark relaxation time

        sigma_susu_nu = self.data_usus.get_integrated_cross_section()*nu
        sigma_sdsd_nd = self.data_dsds.get_integrated_cross_section()*nd
        sigma_ssss_ns = self.data_ssss.get_integrated_cross_section()*ns

        sigma_subarsubar_nubar = self.data_subarsubar.get_integrated_cross_section()*nubar
        sigma_sdbarsdbar_ndbar = self.data_sdbarsdbar.get_integrated_cross_section()*ndbar
        sigma_ssbaruubar_nsbar = self.data_ssbaruubar.get_integrated_cross_section()*nsbar
        sigma_ssbarddbar_nsbar = self.data_ssbarddbar.get_integrated_cross_section()*nsbar
        sigma_ssbarssbar_nsbar = self.data_ssbarssbar.get_integrated_cross_section()*nsbar

        # relaxation time is in fermis
        self.tau_strange_quark = hbarc_gevfm*1.0/(
            sigma_susu_nu +
            sigma_sdsd_nd +
            sigma_ssss_ns +
            sigma_subarsubar_nubar +
            sigma_sdbarsdbar_ndbar +
            sigma_ssbaruubar_nsbar +
            sigma_ssbarddbar_nsbar +
            sigma_ssbarssbar_nsbar
        )

        # Calculate the up antiquark relaxation time

        sigma_ubaruubaru_nu = self.data_uubaruubar.get_integrated_cross_section()*nu
        sigma_ubarudbard_nu = self.data_uubarddbar.get_integrated_cross_section()*nu
        sigma_ubarusbars_nu = self.data_uubarssbar.get_integrated_cross_section()*nu
        sigma_ubardubard_nd = self.data_dubardubar.get_integrated_cross_section()*nd
        sigma_ubarsubars_ns = self.data_subarsubar.get_integrated_cross_section()*ns

        sigma_ubarubarubarubar_nubar = self.data_ubarubarubarubar.get_integrated_cross_section()*nubar
        sigma_ubardbarubardbar_ndbar = self.data_ubardbarubardbar.get_integrated_cross_section()*ndbar
        sigma_ubarsbarubarsbar_nsbar = self.data_ubarsbarubarsbar.get_integrated_cross_section()*nsbar

        # relaxation time is in fermis
        self.tau_up_antiquark = hbarc_gevfm*1.0/(
            sigma_ubaruubaru_nu +
            sigma_ubarudbard_nu +
            sigma_ubarusbars_nu +
            sigma_ubardubard_nd +
            sigma_ubarsubars_ns +
            sigma_ubarubarubarubar_nubar +
            sigma_ubardbarubardbar_ndbar +
            sigma_ubarsbarubarsbar_nsbar
        )

        # Calculate the down quark relaxation time

        sigma_dbarudbaru_nu = self.data_udbarudbar.get_integrated_cross_section()*nu
        sigma_dbardubaru_nd = self.data_ddbaruubar.get_integrated_cross_section()*nd
        sigma_dbarddbard_nd = self.data_ddbarddbar.get_integrated_cross_section()*nd
        sigma_dbardsbars_nd = self.data_ddbarssbar.get_integrated_cross_section()*nd
        sigma_dbarsdbars_ns = self.data_sdbarsdbar.get_integrated_cross_section()*ns

        sigma_dbarubardbarubar_nubar = self.data_ubardbarubardbar.get_integrated_cross_section()*nubar
        sigma_dbardbardbardbar_ndbar = self.data_dbardbardbardbar.get_integrated_cross_section()*ndbar
        sigma_dbarsbardbarsbar_nsbar = self.data_dbarsbardbarsbar.get_integrated_cross_section()*nsbar

        # relaxation time is in fermis
        self.tau_down_antiquark = hbarc_gevfm*1.0/(
            sigma_dbarudbaru_nu +
            sigma_dbardubaru_nd +
            sigma_dbarddbard_nd +
            sigma_dbardsbars_nd +
            sigma_dbarsdbars_ns + 
            sigma_dbarubardbarubar_nubar + 
            sigma_dbardbardbardbar_ndbar +
            sigma_dbarsbardbarsbar_nsbar
        )

        # Calculate the strange antiquark relaxation time

        sigma_sbarusbaru_nu = self.data_usbarusbar.get_integrated_cross_section()*nu
        sigma_sbardsbard_nd = self.data_dsbardsbar.get_integrated_cross_section()*nd
        sigma_sbarsubaru_ns = self.data_ssbaruubar.get_integrated_cross_section()*ns
        sigma_sbarsdbard_ns = self.data_ssbarddbar.get_integrated_cross_section()*ns
        sigma_sbarssbars_ns = self.data_ssbarssbar.get_integrated_cross_section()*ns

        sigma_sbarubarsbarubar_nubar = self.data_ubarsbarubarsbar.get_integrated_cross_section()*nubar
        sigma_sbardbarsbardbar_ndbar = self.data_dbarsbardbarsbar.get_integrated_cross_section()*ndbar
        sigma_sbarsbarsbarsbar_nsbar = self.data_sbarsbarsbarsbar.get_integrated_cross_section()*nsbar

        # relaxation time is in fermis
        self.tau_strange_antiquark = hbarc_gevfm*1.0/(
            sigma_sbarusbaru_nu +
            sigma_sbardsbard_nd +
            sigma_sbarsubaru_ns +
            sigma_sbarsdbard_ns +
            sigma_sbarssbars_ns +
            sigma_sbarubarsbarubar_nubar +
            sigma_sbardbarsbardbar_ndbar +
            sigma_sbarsbarsbarsbar_nsbar
        )

    def _save_to_file(
        self, 
        output_data_filepath: str
    ) -> None:
        temperature = self.data_uuuu.get_temperature()
        up_quark_effective_mass = self.data_uuuu.get_up_quark_effective_mass()
        down_quark_effective_mass = self.data_uuuu.get_down_quark_effective_mass()
        strange_quark_effective_mass = self.data_uuuu.get_strange_quark_effective_mass()
        up_quark_effective_chemical_potential = self.data_uuuu.get_up_quark_effective_chemical_potential()
        down_quark_effective_chemical_potential = self.data_uuuu.get_down_quark_effective_chemical_potential()
        strange_quark_effective_chemical_potential = self.data_uuuu.get_strange_quark_effective_chemical_potential()

        column_data = [
            temperature, 
            self.tau_up_quark, 
            self.tau_down_quark,
            self.tau_strange_quark,
            self.tau_up_antiquark,
            self.tau_down_antiquark,
            self.tau_strange_antiquark,
            up_quark_effective_mass,
            down_quark_effective_mass,
            strange_quark_effective_mass,
            up_quark_effective_chemical_potential,
            down_quark_effective_chemical_potential,
            strange_quark_effective_chemical_potential,
        ]
        
        column_labels = [
            "T[GeV]",
            "tauU[fm]",
            "tauD[fm]",
            "tauS[fm]",
            "tauUBar[fm]",
            "tauDBar[fm]",
            "tauSBar[fm]",
            "effMU[GeV]",
            "effMD[GeV]",
            "effMS[GeV]",
            "effCPU[GeV]",
            "effCPD[GeV]",
            "effCPS[GeV]",
        ]

        save_columns_to_file(
            output_data_filepath, 
            25, 
            15, 
            column_data, 
            column_labels
        )
