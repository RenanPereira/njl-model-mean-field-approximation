#ifndef UNITARYGROUP3DIMENSIONS_H
#define UNITARYGROUP3DIMENSIONS_H

#include <vector> 
#include <iostream>
#include "gsl_wrapper/ComplexSquareMatrixGSL.h"

using namespace std;

ComplexSquareMatrixGSL unitaryGroup3DGenerator0();

ComplexSquareMatrixGSL unitaryGroup3DGenerator1();

ComplexSquareMatrixGSL unitaryGroup3DGenerator2();

ComplexSquareMatrixGSL unitaryGroup3DGenerator3();

ComplexSquareMatrixGSL unitaryGroup3DGenerator4();

ComplexSquareMatrixGSL unitaryGroup3DGenerator5();

ComplexSquareMatrixGSL unitaryGroup3DGenerator6();

ComplexSquareMatrixGSL unitaryGroup3DGenerator7();

ComplexSquareMatrixGSL unitaryGroup3DGenerator8();

ComplexSquareMatrixGSL unitaryGroup3DGenerator(int );

gsl_complex unitaryGroup3DCalculateAntisymmetricStructureConstant(int , int , int );

gsl_complex unitaryGroup3DCalculateSymmetricStructureConstant(int , int , int );


#endif