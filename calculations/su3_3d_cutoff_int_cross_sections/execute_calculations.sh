#!/bin/bash

# Parameter passed to the script
calculation="$1"

echo "Script that calculates the Integrated Cross Sections for the NJL model for different parameter sets \
as functions of temperature and chemical potential. Different methods are used to evaluate the integration of the cross sections."
echo ""

# Go to project root
cd ../.. || exit

# Clean previous build and re-build
make clean
make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_int_cross_sections/data/$calculation

cd calculations/su3_3d_cutoff_int_cross_sections/data/$calculation


if [ "$calculation" = "setA_CP0p000000" ]; then
    filename="integratedCrossSections_setA_TMin0p120000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setA_TMin0p205000_TMax0p220000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setA_CP0p000000_KLEVANSKY" ]; then
    filename="integratedCrossSections_setA_TMin0p120000_TMax0p300000.ini"
    sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = KLEVANSKY/g' $filename
    ./nambuJonaLasinioModel.out use-config-file $filename
    sed -i 's/approximationMethod = KLEVANSKY/approximationMethod = COMPLETE_COV/g' $filename

    filename="integratedCrossSections_setA_TMin0p205000_TMax0p220000.ini"
    sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = KLEVANSKY/g' $filename
    ./nambuJonaLasinioModel.out use-config-file $filename
    sed -i 's/approximationMethod = KLEVANSKY/approximationMethod = COMPLETE_COV/g' $filename

elif [ "$calculation" = "setA_CP0p000000_ZHUANG" ]; then
    filename="integratedCrossSections_setA_TMin0p120000_TMax0p300000.ini"
    sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = ZHUANG/g' $filename
    ./nambuJonaLasinioModel.out use-config-file $filename
    sed -i 's/approximationMethod = ZHUANG/approximationMethod = COMPLETE_COV/g' $filename

    filename="integratedCrossSections_setA_TMin0p205000_TMax0p220000.ini"
    sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = ZHUANG/g' $filename
    ./nambuJonaLasinioModel.out use-config-file $filename
    sed -i 's/approximationMethod = ZHUANG/approximationMethod = COMPLETE_COV/g' $filename

elif [ "$calculation" = "setA_CP0p050000" ]; then
    filename="integratedCrossSections_setA_TMin0p040000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setA_CP0p100000" ]; then
    filename="integratedCrossSections_setA_TMin0p040000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setA_CP0p200000" ]; then
    filename="integratedCrossSections_setA_TMin0p040000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setA_CP0p300000" ]; then
    filename="integratedCrossSections_setA_TMin0p040000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setA_CP0p318436" ]; then
    filename="integratedCrossSections_setA_TMin0p040000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setA_TMin0p062100_TMax0p072100.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setA_TMin0p067090_TMax0p068110.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setB_CP0p000000" ]; then
    filename="integratedCrossSections_setB_TMin0p120000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setB_TMin0p153000_TMax0p178000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setB_CP0p231030" ]; then
    filename="integratedCrossSections_setB_TMin0p070000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setB_TMin0p095100_TMax0p105100.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setC_CP0p000000" ]; then
    filename="integratedCrossSections_setC_TMin0p120000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setC_TMin0p135000_TMax0p160000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

elif [ "$calculation" = "setC_CP0p164012" ]; then
    filename="integratedCrossSections_setC_TMin0p084000_TMax0p300000.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

    filename="integratedCrossSections_setC_TMin0p108100_TMax0p118100.ini"
    ./nambuJonaLasinioModel.out use-config-file $filename

else
    echo "Error: unknown calculation '$calculation'"
fi

rm nambuJonaLasinioModel.out
