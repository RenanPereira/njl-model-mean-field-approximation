import os
from common_utils.transport_coefficients.shear_viscosity import ShearViscosity


quark_rel_times_filepath_prefix = "su3_3d_cutoff_quark_relaxation_times/data/RelaxationTimes"
output_filepath_prefix = "su3_3d_cutoff_transport_coefficients/data/ShearViscosity"

configs = [
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setA_COMPLETE_COV_CP0p0.dat",
        "output_filepath": f"{output_filepath_prefix}_setA_COMPLETE_COV_CP0p0.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setA_KLEVANSKY_CP0p0.dat",
        "output_filepath": f"{output_filepath_prefix}_setA_KLEVANSKY_CP0p0.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setA_ZHUANG_CP0p0.dat",
        "output_filepath": f"{output_filepath_prefix}_setA_ZHUANG_CP0p0.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setA_COMPLETE_COV_CP0p318436.dat",
        "output_filepath": f"{output_filepath_prefix}_setA_COMPLETE_COV_CP0p318436.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setB_COMPLETE_COV_CP0p0.dat",
        "output_filepath": f"{output_filepath_prefix}_setB_COMPLETE_COV_CP0p0.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setB_COMPLETE_COV_CP0p231030.dat",
        "output_filepath": f"{output_filepath_prefix}_setB_COMPLETE_COV_CP0p231030.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setC_COMPLETE_COV_CP0p0.dat",
        "output_filepath": f"{output_filepath_prefix}_setC_COMPLETE_COV_CP0p0.dat"
    },
    {
        "quark_rel_times_filepath": f"{quark_rel_times_filepath_prefix}_setC_COMPLETE_COV_CP0p164012.dat",
        "output_filepath": f"{output_filepath_prefix}_setC_COMPLETE_COV_CP0p164012.dat"
    },
]

for config in configs:
    quark_rel_times_filepath = config["quark_rel_times_filepath"]
    output_filepath = config["output_filepath"]
    if not os.path.exists(quark_rel_times_filepath):
        print(f"Skipping missing file: {quark_rel_times_filepath}")
        continue
            
    print(f"Calculating the Shear Viscosity using the file:")
    print(quark_rel_times_filepath)
    ShearViscosity(quark_rel_times_filepath, output_filepath)
    print()
