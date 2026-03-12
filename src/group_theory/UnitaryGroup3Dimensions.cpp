#include "group_theory/UnitaryGroup3Dimensions.h"
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>


//Gell-Mann Matrices, Identity and Group Structure constants


//lambda0 = Sqrt[2/3] IdentityMatrix[3];
ComplexSquareMatrixGSL unitaryGroup3DGenerator0()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(sqrt(2.0/3.0), 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0,           0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0,           0.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0,           0.0));
    lambda.setValue(1, 1, gsl_complex_rect(sqrt(2.0/3.0), 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0,           0.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0,           0.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0,           0.0));
    lambda.setValue(2, 2, gsl_complex_rect(sqrt(2.0/3.0), 0.0));

    return lambda;
}


//lambda1 = {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator1()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(1.0, 0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(1, 0, gsl_complex_rect(1.0, 0.0));
    lambda.setValue(1, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda2 = {{0, -I, 0}, {I, 0, 0}, {0, 0, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator2()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(0.0,  0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0, -1.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0,  0.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0, 1.0));
    lambda.setValue(1, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda3 = {{1, 0, 0}, {0, -1, 0}, {0, 0, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator3()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(1.0, 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(1, 0, gsl_complex_rect( 0.0, 0.0));
    lambda.setValue(1, 1, gsl_complex_rect(-1.0, 0.0));
    lambda.setValue(1, 2, gsl_complex_rect( 0.0, 0.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda4 = {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator4()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 2, gsl_complex_rect(1.0, 0.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(2, 0, gsl_complex_rect(1.0, 0.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda5 = {{0, 0, -I}, {0, 0, 0}, {I, 0, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator5()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(0.0,  0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0,  0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0, -1.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0, 1.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda6 = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator6()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(1.0, 0.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 1, gsl_complex_rect(1.0, 0.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda7 = {{0, 0, 0}, {0, 0, -I}, {0, I, 0}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator7()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0, 0.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0,  0.0));
    lambda.setValue(1, 1, gsl_complex_rect(0.0,  0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0, -1.0));

    lambda.setValue(2, 0, gsl_complex_rect(0.0, 0.0));
    lambda.setValue(2, 1, gsl_complex_rect(0.0, 1.0));
    lambda.setValue(2, 2, gsl_complex_rect(0.0, 0.0));

    return lambda;
}


//lambda8 = Sqrt[1/3] {{1, 0, 0}, {0, 1, 0}, {0, 0, -2}};
ComplexSquareMatrixGSL unitaryGroup3DGenerator8()
{
    ComplexSquareMatrixGSL lambda(3);

    lambda.setValue(0, 0, gsl_complex_rect(1.0/sqrt(3.0), 0.0));
    lambda.setValue(0, 1, gsl_complex_rect(0.0,           0.0));
    lambda.setValue(0, 2, gsl_complex_rect(0.0,           0.0));

    lambda.setValue(1, 0, gsl_complex_rect(0.0,           0.0));
    lambda.setValue(1, 1, gsl_complex_rect(1.0/sqrt(3.0), 0.0));
    lambda.setValue(1, 2, gsl_complex_rect(0.0,           0.0));

    lambda.setValue(2, 0, gsl_complex_rect( 0.0,           0.0));
    lambda.setValue(2, 1, gsl_complex_rect( 0.0,           0.0));
    lambda.setValue(2, 2, gsl_complex_rect(-2.0/sqrt(3.0), 0.0));

    return lambda;
}


ComplexSquareMatrixGSL unitaryGroup3DGenerator(int n)
{
    if      ( n==0 ){ return unitaryGroup3DGenerator0(); }
    else if ( n==1 ){ return unitaryGroup3DGenerator1(); }
    else if ( n==2 ){ return unitaryGroup3DGenerator2(); }
    else if ( n==3 ){ return unitaryGroup3DGenerator3(); }
    else if ( n==4 ){ return unitaryGroup3DGenerator4(); }
    else if ( n==5 ){ return unitaryGroup3DGenerator5(); }
    else if ( n==6 ){ return unitaryGroup3DGenerator6(); }
    else if ( n==7 ){ return unitaryGroup3DGenerator7(); }
    else if ( n==8 ){ return unitaryGroup3DGenerator8(); }
    else
    {
        std::cout << "The requested index for the unitaryGroup3DGenerator does not exist! It must be in the set: [0,8]! Aborting!\n";
        abort();
        return ComplexSquareMatrixGSL(0);
    }
}


gsl_complex unitaryGroup3DCalculateAntisymmetricStructureConstant(int a, int b, int c)
{
    ComplexSquareMatrixGSL commutatorLambdaALambdaB = commutator(unitaryGroup3DGenerator(a), unitaryGroup3DGenerator(b));

    ComplexSquareMatrixGSL commutatorLambdaALambdaBTimesLambdaC = multiply(commutatorLambdaALambdaB, unitaryGroup3DGenerator(c));

    gsl_complex trace = commutatorLambdaALambdaBTimesLambdaC.trace();
    gsl_complex fabc = gsl_complex_mul(gsl_complex_rect(0.0, -1.0/4.0), trace);

    return fabc;
}


gsl_complex unitaryGroup3DCalculateSymmetricStructureConstant(int a, int b, int c)
{
    ComplexSquareMatrixGSL anticommutatorLambdaALambdaB = anticommutator(unitaryGroup3DGenerator(a), unitaryGroup3DGenerator(b));

    ComplexSquareMatrixGSL anticommutatorLambdaALambdaBTimesLambdaC = multiply(anticommutatorLambdaALambdaB, unitaryGroup3DGenerator(c));

    gsl_complex trace = anticommutatorLambdaALambdaBTimesLambdaC.trace();
    gsl_complex dabc = gsl_complex_mul(gsl_complex_rect(1.0/4.0, 0.0), trace);

    return dabc;
}


