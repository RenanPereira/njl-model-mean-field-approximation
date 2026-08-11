import os
from common_utils.transport_coefficients.seebeck_sigmae_product import SeebeckSigmaeProduct


quark_rel_times_filepath_prefix = "su3_3d_cutoff_quark_relaxation_times/data/RelaxationTimes"
thermodynamics_filepath_prefix = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp"
output_filepath_prefix = "su3_3d_cutoff_transport_coefficients/data/SeebeckSigmaeProduct"

configs = [
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setA_COMPLETE_COV_CP0p318436.dat",
        "thermodynamics_filepath": f"{thermodynamics_filepath_prefix}_setA_TMin0p0_TMax0p5_CPU0p318436.dat",
        "output_filepath": f"{output_filepath_prefix}_setA_COMPLETE_COV_CP0p318436.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setB_COMPLETE_COV_CP0p231030.dat",
        "thermodynamics_filepath": f"{thermodynamics_filepath_prefix}_setB_TMin0p0_TMax0p5_CPU0p23103.dat",
        "output_filepath": f"{output_filepath_prefix}_setB_COMPLETE_COV_CP0p231030.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setC_COMPLETE_COV_CP0p164012.dat",
        "thermodynamics_filepath": f"{thermodynamics_filepath_prefix}_setC_TMin0p0_TMax0p5_CPU0p164012.dat",
        "output_filepath": f"{output_filepath_prefix}_setC_COMPLETE_COV_CP0p164012.dat"
    },
]

for config in configs:
    quark_rel_times_filepath = config["quark_rel_times_filepath"]
    thermodynamics_filepath = config["thermodynamics_filepath"]
    output_filepath = config["output_filepath"]
    
    if not os.path.exists(quark_rel_times_filepath):
        print(f"Skipping missing file: {quark_rel_times_filepath}")
        continue
    if not os.path.exists(thermodynamics_filepath):
        print(f"Skipping missing file: {thermodynamics_filepath}")
        continue

    print(f"Calculating the SeebeckSigmaeProduct using the files:")
    print(quark_rel_times_filepath)
    print(thermodynamics_filepath)
    SeebeckSigmaeProduct(quark_rel_times_filepath, thermodynamics_filepath, output_filepath)
    print()
