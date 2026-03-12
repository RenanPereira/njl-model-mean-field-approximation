#ifndef COMPLEXSQUAREMATRIXGSL_H
#define COMPLEXSQUAREMATRIXGSL_H

#include <iostream>
#include <gsl/gsl_matrix.h>


class ComplexSquareMatrixGSL
{
private:
    int dimension = 0;
    gsl_matrix_complex* matrixGSLPtr;

public:
    ComplexSquareMatrixGSL(int dimensionAux)//constructor
    {
        dimension = dimensionAux;
        matrixGSLPtr = gsl_matrix_complex_alloc(dimension, dimension);
    };

    ComplexSquareMatrixGSL(gsl_matrix_complex* matrixGSLPtrAux)//constructor
    {
        dimension = matrixGSLPtrAux->size1;//read gsl_matrix dimension from struct
        matrixGSLPtr = gsl_matrix_complex_alloc(dimension, dimension);
        gsl_matrix_complex_memcpy(matrixGSLPtr, matrixGSLPtrAux);//copy matrix
    };

    ~ComplexSquareMatrixGSL()//destructor
    {
        gsl_matrix_complex_free(matrixGSLPtr);
    };

    ComplexSquareMatrixGSL(const ComplexSquareMatrixGSL &complexSquareMatrix)//copy constructor
    {
        dimension = complexSquareMatrix.dimension;
        matrixGSLPtr = gsl_matrix_complex_alloc(dimension, dimension);
        gsl_matrix_complex_memcpy(matrixGSLPtr, complexSquareMatrix.matrixGSLPtr);//copy matrix
    };

    ComplexSquareMatrixGSL& operator=(const ComplexSquareMatrixGSL& complexSquareMatrix)//copy assignment
    {
        if (this == &complexSquareMatrix){ return *this; }
        
        //free old matrix
        gsl_matrix_complex_free(matrixGSLPtr);

        dimension = complexSquareMatrix.dimension;
        matrixGSLPtr = gsl_matrix_complex_alloc(dimension, dimension);
        gsl_matrix_complex_memcpy(matrixGSLPtr, complexSquareMatrix.matrixGSLPtr);//copy matrix

        return *this;
    }

    int getDimension(){ return dimension; };
    gsl_matrix_complex* getMatrixGSLPtr(){ return matrixGSLPtr; };

    void setValue(int , int , gsl_complex );
    gsl_complex getValue(int , int );

    gsl_complex determinant();
    gsl_complex trace();

    ComplexSquareMatrixGSL inverse();
};

ComplexSquareMatrixGSL add(ComplexSquareMatrixGSL , ComplexSquareMatrixGSL );

ComplexSquareMatrixGSL subtract(ComplexSquareMatrixGSL , ComplexSquareMatrixGSL );

ComplexSquareMatrixGSL multiply(ComplexSquareMatrixGSL , ComplexSquareMatrixGSL );

ComplexSquareMatrixGSL kroneckerDelta(int );

ComplexSquareMatrixGSL commutator(ComplexSquareMatrixGSL , ComplexSquareMatrixGSL );

ComplexSquareMatrixGSL anticommutator(ComplexSquareMatrixGSL , ComplexSquareMatrixGSL );

#endif