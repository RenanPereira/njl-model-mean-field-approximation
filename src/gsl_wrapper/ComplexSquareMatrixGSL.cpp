#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include "gsl_wrapper/ComplexSquareMatrixGSL.h"


void ComplexSquareMatrixGSL::setValue(int a, int b, gsl_complex value)
{
    if ( (a>=0 && a<dimension) || (b>=0 && b<dimension) )
    {
        gsl_matrix_complex_set(matrixGSLPtr, a, b, value);//linha, coluna
    }
    else
    {
        cout << "Index set using setValue in ComplexSquareMatrixGSL out of range! Aborting!\n";
        abort();
    }
}


gsl_complex ComplexSquareMatrixGSL::getValue(int a, int b)
{	
	gsl_complex matrixElement;

	if ( (a>=0 && a<dimension) || (b>=0 && b<dimension) )
    {
        matrixElement = gsl_matrix_complex_get(matrixGSLPtr, a, b);
    }
    else
    {
        cout << "Index get using getValue in ComplexSquareMatrixGSL out of range! Aborting!\n";
        abort();
    }

	return matrixElement;
}


gsl_complex ComplexSquareMatrixGSL::determinant()
{   
    gsl_matrix_complex* matrixGSLPtrAux = gsl_matrix_complex_alloc(dimension, dimension);
    gsl_matrix_complex_memcpy(matrixGSLPtrAux, matrixGSLPtr);//copy gslMatrix

    int s;
    gsl_permutation* p = gsl_permutation_alloc(dimension);
    gsl_linalg_complex_LU_decomp(matrixGSLPtrAux, p, &s);   

    //calculate determinant using the LU decomposition
    gsl_complex determinant = gsl_linalg_complex_LU_det(matrixGSLPtrAux, s);

    //free memory
    gsl_matrix_complex_free(matrixGSLPtrAux);
    gsl_permutation_free(p);

    return determinant;
}


gsl_complex ComplexSquareMatrixGSL::trace()
{   
    gsl_complex trace = gsl_complex_rect(0, 0);

    for (int i = 0; i < dimension; ++i){ trace = gsl_complex_add(trace, getValue(i, i)); }

    return trace;
}


ComplexSquareMatrixGSL ComplexSquareMatrixGSL::inverse()
{
	gsl_matrix_complex* matrixGSLPtrAux = gsl_matrix_complex_alloc(dimension, dimension);
    gsl_matrix_complex_memcpy(matrixGSLPtrAux, matrixGSLPtr);//copy gslMatrix

    int s;
    gsl_permutation* p = gsl_permutation_alloc(dimension);
    gsl_matrix_complex* inverseMatrixAux = gsl_matrix_complex_alloc(dimension, dimension);

    //calculate the inverse matrix using the LU decomposition
    gsl_linalg_complex_LU_decomp(matrixGSLPtrAux, p, &s);   
    gsl_linalg_complex_LU_invert(matrixGSLPtrAux, p, inverseMatrixAux);

    ComplexSquareMatrixGSL inverseMatrix(inverseMatrixAux);

    //free memory
    gsl_matrix_complex_free(matrixGSLPtrAux);
    gsl_permutation_free(p);
    gsl_matrix_complex_free(inverseMatrixAux);

    return inverseMatrix;
}


ComplexSquareMatrixGSL add(ComplexSquareMatrixGSL matrixA, ComplexSquareMatrixGSL matrixB)
{   
    //check if dimensions of matrices A and B are the same
    int dimensionA = matrixA.getDimension();
    int dimensionB = matrixB.getDimension();
    if ( dimensionA!=dimensionB )
    {
        cout << "Trying to add ComplexSquareMatrixGSL of different dimensions! Aborting!\n"; 
        abort();
    }

    //create gsl matrix to store result of subtraction
    gsl_matrix_complex* matrixGSLPtrAux = gsl_matrix_complex_alloc(dimensionA, dimensionA);
    gsl_matrix_complex_memcpy(matrixGSLPtrAux, matrixA.getMatrixGSLPtr());//copy matrixA to matrix  

    //This function adds the elements of matrix b to the elements of matrix a. The result a(i, j) ← a(i, j) + b(i, j) is
    //stored in a and b remains unchanged. The two matrices must have the same dimensions.
    gsl_matrix_complex_add(matrixGSLPtrAux, matrixB.getMatrixGSLPtr());

    ComplexSquareMatrixGSL matrixAPlusMatrixB(matrixGSLPtrAux);

    //free memory
    gsl_matrix_complex_free(matrixGSLPtrAux);

    return matrixAPlusMatrixB;
}


ComplexSquareMatrixGSL subtract(ComplexSquareMatrixGSL matrixA, ComplexSquareMatrixGSL matrixB)
{ 	
	//check if dimensions of matrices A and B are the same
	int dimensionA = matrixA.getDimension();
	int dimensionB = matrixB.getDimension();
	if ( dimensionA!=dimensionB )
	{
		cout << "Trying to subtract ComplexSquareMatrixGSL of different dimensions! Aborting!\n"; 
		abort();
	}

	//create gsl matrix to store result of subtraction
	gsl_matrix_complex* matrixGSLPtrAux = gsl_matrix_complex_alloc(dimensionA, dimensionA);
    gsl_matrix_complex_memcpy(matrixGSLPtrAux, matrixA.getMatrixGSLPtr());//copy matrixA to matrix  

    //The function below subtracts the elements of matrix b from the elements of matrix a.
    //The result a(i, j) <- a(i, j) − b(i, j) is stored in a and b remains unchanged
    gsl_matrix_complex_sub(matrixGSLPtrAux, matrixB.getMatrixGSLPtr());

	ComplexSquareMatrixGSL matrixAMinusMatrixB(matrixGSLPtrAux);

	//free memory
    gsl_matrix_complex_free(matrixGSLPtrAux);

	return matrixAMinusMatrixB;
}


ComplexSquareMatrixGSL multiply(ComplexSquareMatrixGSL matrixA, ComplexSquareMatrixGSL matrixB)
{   
    //check if dimensions of matrices A and B are the same
    int dimensionA = matrixA.getDimension();
    int dimensionB = matrixB.getDimension();
    if ( dimensionA!=dimensionB )
    {
        cout << "Trying to multiply ComplexSquareMatrixGSL of different dimensions! Aborting!\n"; 
        abort();
    }

    ComplexSquareMatrixGSL matrixATimesMatrixB(dimensionA);
    for (int i = 0; i < dimensionA; ++i)
    {
        for (int j = 0; j < dimensionA; ++j)
        {
            gsl_complex Cij = gsl_complex_rect(0, 0);
            for (int k = 0; k < dimensionA; ++k)
            {   
                gsl_complex Aik = matrixA.getValue(i, k);
                gsl_complex Bkj = matrixB.getValue(k, j);
                Cij = gsl_complex_add(Cij, gsl_complex_mul(Aik, Bkj));
            }
            matrixATimesMatrixB.setValue(i , j, Cij);
        }
    }

    return matrixATimesMatrixB;
}


ComplexSquareMatrixGSL kroneckerDelta(int dimension)
{
    ComplexSquareMatrixGSL delta(dimension);

    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            if ( i==j ){ delta.setValue(i, j, gsl_complex_rect(1.0, 0.0)); }
            else       { delta.setValue(i, j, gsl_complex_rect(0.0, 0.0)); }
        }
    }

    return delta;
}


ComplexSquareMatrixGSL commutator(ComplexSquareMatrixGSL matrixA, ComplexSquareMatrixGSL matrixB)
{
    ComplexSquareMatrixGSL matrixATimesMatrixB = multiply(matrixA, matrixB);
    ComplexSquareMatrixGSL matrixBTimesMatrixA = multiply(matrixB, matrixA);

    ComplexSquareMatrixGSL matrixCommutatorAB = subtract(matrixATimesMatrixB, matrixBTimesMatrixA);

    return matrixCommutatorAB;
}


ComplexSquareMatrixGSL anticommutator(ComplexSquareMatrixGSL matrixA, ComplexSquareMatrixGSL matrixB)
{
    ComplexSquareMatrixGSL matrixATimesMatrixB = multiply(matrixA, matrixB);
    ComplexSquareMatrixGSL matrixBTimesMatrixA = multiply(matrixB, matrixA);

    ComplexSquareMatrixGSL matrixCommutatorAB = add(matrixATimesMatrixB, matrixBTimesMatrixA);

    return matrixCommutatorAB;
}


