//
// Curtis Glavin
// idealist@cs.bu.edu
//
// matrix.cpp
//
// The iterative inversion routine
// is based on code from the following
// website:
// http://rustam.uwp.edu/499/programs/matrices/invert.c
// whose author is not specified.
//
// This code may be used, modified, distributed, etc.
// as desired, so long as the entirety of any project
// that uses any or all of this code is available as
// open source and as freeware.
//
// No guarantees whatsoever (regarding the
// appropriateness of this code for a specific purpose,
// etc.) are made.
//
// In the short term, I would like to implement
// reduced row echelon form (rref()) and determinant
// methods.
//
// Please email suggestions, bugs, etc. to
// idealist@cs.bu.edu
//

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

#include "matrix.h"

#define SCALE 1.4

unsigned matrixClass::inverse_precision = DEFAULT_INVERSE_PRECISION;

//
// Constructor
// Dynamically allocates an array of
// m*n doubles and initializes each
// element to 0.
//
matrixClass::matrixClass(unsigned m, unsigned n)
{
        assert(m>0 && n>0);
        matrix = new double[m*n];
        assert(matrix != 0);
        rows = m;
        columns = n;
        for(unsigned counter=(m*n); counter-->0;)
                matrix[counter] = 0.0;
}

//
// Copy Constructor
// Dynamically allocates a new matrix
// with the same dimensions and data
// as the specified initializer matrix.
//
matrixClass::matrixClass(const matrixClass& initializer)
{
        unsigned initrows=0,initcolumns=0;
        initrows = initializer.rows;
        initcolumns = initializer.columns;
        assert(initrows > 0);
        assert(initcolumns > 0);
        assert(initializer.matrix != 0);
        matrix = new double[initrows*initcolumns];
        assert(matrix != 0);
        rows = initrows;
        columns = initcolumns;
        for(unsigned counter=(rows*columns); counter-->0;)
                matrix[counter] = initializer.matrix[counter];
}

//
// Destructor
//
matrixClass::~matrixClass()
{
        delete [] matrix;
}

//
// Transpose function; will not
// alter the original matrix but
// will return a matrix which
// is the transpose of the original.
//
const matrixClass matrixClass::transpose(void) const
{
        matrixClass transposed(columns, rows);
        for(unsigned i=0; i<rows; i++)
                for(unsigned j=0; j<columns; j++)
                        transposed.matrix[j*rows+i] = matrix[i*columns + j];
        return transposed;
}

//
// Overloaded function call operator
// which allows the caller to place
// or see a value at the specified
// row,column in the matrix.  The
// overloaded function operator method
// is itself overloaded to handle const
// and non-const objects.
//
const double& matrixClass::operator()(unsigned m, unsigned n) const
{
        assert(m>=0 && m<rows);
        assert(n>=0 && n<columns);
        return matrix[m*columns + n];
}
double& matrixClass::operator()(unsigned m, unsigned n)
{
        assert(m>=0 && m<rows);
        assert(n>=0 && n<columns);
        return matrix[m*columns + n];
}

//
// Overloaded assignment operator
// which allows the caller to assign
// one matrix to another.
// The return is const to prevent
// assignments such as:
// (matrix1 = matrix2) = matrix3;
// while allowing
// matrix1 = matrix2 = matrix3;
//
const matrixClass& matrixClass::operator=(const matrixClass& matrix1)
{
        //
        // This test is to avoid
        // self-assignment.
        //
        if(&matrix1 != this)
        {
                if(matrix1.rows != rows || matrix1.columns != columns)
                {
                        delete [] matrix;
                        rows = matrix1.rows;
                        columns = matrix1.columns;
                        matrix = new double[rows*columns];
                        assert(matrix != 0);
                }

                for(unsigned counter=(rows*columns); counter-->0;)
                        matrix[counter] = matrix1.matrix[counter];
        }
        return *this;
}

//
// Overloaded equivalence operators.
// Matrices are considered equivalent
// if and only if they are of the same
// dimension and all elements at
// corresponding locations in the 2
// matrices are equivalent.  True is
// returned if the matrices are
// equivalent, false is returned otherwise.
//
bool matrixClass::operator==(const matrixClass& matrix1) const
{
        if(matrix1.rows != rows || matrix1.columns != columns)
                return false;
        for(unsigned counter=(rows*columns); counter-->0;)
                if(matrix1.matrix[counter] != matrix[counter])
                        return false;
        return true;
}

//
// Overloaded multiplication operator
// for multiplying 2 matrices.  The
// method will terminate the program
// if the matrices' dimensions prevent
// multiplication.
//
const matrixClass matrixClass::operator*(const matrixClass& matrix1) const
{
        assert(columns == matrix1.rows);
        matrixClass product(rows, matrix1.columns);
        //
        // Each loop fills 1 row of
        // the product matrix.
        //
        for(unsigned i=0; i<rows; i++)
                //
                // Each loop fills 1 element
                // of the row.
                //
                for(unsigned j=0; j<matrix1.columns; j++)
                        //
                        // Each loop calculates one
                        // component of the sum.
                        for(unsigned k=0; k<columns; k++)
                                product(i,j) += ((*this)(i,k)) * matrix1(k,j);
        return product;
}

//
// Overloaded multiplication
// operator for scalar multiplication
// of the form:
//   matrix * scalar
//
const matrixClass matrixClass::operator*(const double& scalar) const
{
        matrixClass m(rows, columns);
        for(unsigned i=(rows*columns); i-->0;)
                 m.matrix[i] = matrix[i] * scalar;
        return m;
}

//
// Overloaded division - Added by Chris Jacobsen
// operator for scalar multiplication
// of the form:
//   matrix / scalar
//
const matrixClass matrixClass::operator/(const double& scalar) const
{
        matrixClass m(rows, columns);
        for(unsigned i=(rows*columns); i-->0;)
                 m.matrix[i] = matrix[i] / scalar;
        return m;
}

//
// Overloaded *= operator for
// matrix multiplication of the form:
// matrix *= matrix;
// The method will terminate the program
// if the matrices' dimensions prevent
// multiplication.
//
const matrixClass& matrixClass::operator*=(const matrixClass& matrix1)
{
        assert(columns == matrix1.rows);
        matrixClass product(rows, matrix1.columns);
        //
        // Each loop fills 1 row of
        // the product matrix.
        //
        for(unsigned i=0; i<rows; i++)
                //
                // Each loop fills 1 element
                // of the row.
                //
                for(unsigned j=0; j<matrix1.columns; j++)
                        //
                        // Each loop calculates one
                        // component of the sum.
                        for(unsigned k=0; k<columns; k++)
                                product(i,j) += ((*this)(i,k)) * matrix1(k,j);
        *this = product;
        return *this;
}

//
// Overloaded *= operator for
// matrix multiplication of the form:
// matrix *= scalar;
//
const matrixClass& matrixClass::operator*=(const double& scalar)
{
        for(unsigned i=(rows*columns); i-->0;)
                matrix[i] *= scalar;
        return *this;
}

//
// Overloaded addition operator
// for the addition of matrices
// of identical dimensions.
// The method will terminate the
// program if the matrices' dimensions
// prevent addition.
//
const matrixClass matrixClass::operator+(const matrixClass& matrix1) const
{
        assert(rows==matrix1.rows && columns==matrix1.columns);
        matrixClass sum(rows,columns);
        for(unsigned i=(rows*columns); i-->0;)
                sum.matrix[i]=matrix[i]+matrix1.matrix[i];
        return sum;
}

//
// Overloaded += operator for
// matrix addition of the form:
// matrix += matrix;
// The method will terminate the program
// if the matrices' dimensions prevent
// addition.
//
const matrixClass& matrixClass::operator+=(const matrixClass& matrix1)
{
        assert(rows==matrix1.rows && columns==matrix1.columns);
        for(unsigned i=(rows*columns); i-->0;)
                matrix[i]+=matrix1.matrix[i];
        return *this;
}

//
// Overloaded subtraction operator
// for the subtraction of matrices
// of identical dimensions.
// The method will terminate the
// program if the matrices' dimensions
// prevent subtraction.
//
const matrixClass matrixClass::operator-(const matrixClass& matrix1) const
{
        assert(rows==matrix1.rows && columns==matrix1.columns);
        matrixClass difference(rows,columns);
        for(unsigned i=(rows*columns); i-->0;)
                difference.matrix[i]=matrix[i]-matrix1.matrix[i];
        return difference;
}

//
// Overloaded -= operator for
// matrix subtraction of the form:
// matrix -= matrix;
// The method will terminate the program
// if the matrices' dimensions prevent
// subtraction.
//
const matrixClass& matrixClass::operator-=(const matrixClass& matrix1)
{
        assert(rows==matrix1.rows && columns==matrix1.columns);
        for(unsigned i=(rows*columns); i-->0;)
                matrix[i]-=matrix1.matrix[i];
        return *this;
}

//
// The iterative inversion routine
// is based on code from the following
// website:
// http://rustam.uwp.edu/499/programs/matrices/invert.c
// whose author is not specified.
//
// Overloaded ^ operator for raising
// a square matrix to an integer power.
// Examples for n*n matrix m:
//
// m^-1 yields inverse of m if it exists
// m^0 yields n*n identity matrix
// m^1 yields m
// m^-i is the product of i inverses of m
// m^3 is the product (m*m*m)
// etc.
//
// The method will terminate the program
// if the matrix is not square (i.e., the
// matrix must have as many rows as it has
// columns).  This is an iterative inversion
// algorithm and its precision is based on
// the value of the class' (static) variable
// inverse_precision.  If the matrix is
// square but singular, the inversion method
// often generates pseudo-meaningless results
// depending on numerical precision. Inverses
// should be checked to ensure that
// matrix * (matrix^-1) = Identity
//
const matrixClass matrixClass::operator^(const int& exponent) const
{
        assert(rows == columns);
        assert(rows > 0);
        matrixClass result(rows,columns);
        unsigned i;
        //
        // This will generate an identity
        // matrix of the same dimension as
        // the original/this matrix.
        //
        if(exponent == 0)
        {
                for(i=0; i<rows; i++)
                        for(unsigned j=0; j<columns; j++)
                        {
                                //
                                // If currently on matrix diagonal.
                                //
                                if(i==j)
                                        result.matrix[i*columns + j] = 1;
                                else
                                        result.matrix[i*columns + j] = 0;
                        }
        }
        else
        {
                if(exponent > 0)
                {
                        result = (*this);
                        for(i=1; i<static_cast<unsigned>(exponent); i++)
                                result *= (*this);
                }
                //
                // Exponent must be less than 0,
                // attempt to find inverse.
                //
                else
                {
                        matrixClass inverse(rows,columns), transposed;
                        double normal=0;

                        for(i=(rows*columns); i-->0;)
                                normal += matrix[i] * matrix[i];

                        normal = SCALE * sqrt(normal);

                        for(i=(rows*columns); i-->0;)
                                inverse.matrix[i]= matrix[i]/normal;

                        transposed = inverse.transpose();

                        for(i=0; i< inverse_precision; i++)
                                transposed = (transposed*2)-((transposed * inverse) * transposed);

                        inverse = (1/normal) * transposed;
                        result = inverse;

                        //
                        // Now find the product of n inverses,
                        // where n is the absolute value of the
                        // exponent.
                        //
                        for(i=(exponent*-1); i>1; i--)
                                result *= inverse;
                }
        }
        return result;
}

//
// Overloaded ^= operator for raising
// a square matrix to an integer power.
// See the description of the overloaded
// ^ operator above for details.
//
const matrixClass& matrixClass::operator^=(const int& exponent)
{
        (*this) = (*this)^exponent;
        return (*this);
}

//
// Methods to view and adjust the
// precision of the inversion method;
// the larger the value passed to this
// method is, the more attempts the
// inversion method will make to
// converge on the inverse and the more
// precise its results will be.  There
// is a performance penatly for
// this additional accuracy.  The
// precision variable is static and
// changing it affects all matrixClass
// objects.
// 
unsigned matrixClass::currentInversePrecision() const
{
        return inverse_precision;
}
//
// Function returns true if successful,
// false otherwise.  Calling with no
// parameters resets precision to its
// default value.  As explained above,
// after each change, the precision
// value is fixed until the next call to
// this method or until the end of
// the program, affecting all matrixClass
// objects.
//
bool matrixClass::setInversePrecision(unsigned precision)
{
        assert(precision > 0);
        inverse_precision = precision;
        return true;
}

//
// This friend function overloads
// the << operator so it can be used
// to display the matrix on the
// screen using standard ostreams.
//
ostream& operator<<(ostream& output, const matrixClass& matrix1)
{
        unsigned lastElement = matrix1.rows*matrix1.columns;
        output << endl;
        for(unsigned i=0; i < lastElement; i++)
        {
                output << setw(18) << matrix1.matrix[i];
                //
                // Begins the next row on a new line.
                //
                if((i+1) % matrix1.columns == 0)
                        output << endl;
        }
        output << endl;
        return output;
}

//
// This friend function overloads
// the >> operator so it can be used
// to input matrix data using standard
// istreams.
//
istream& operator>>(istream& input, matrixClass& matrix1)
{
        unsigned lastElement = matrix1.rows*matrix1.columns;
        for(unsigned i=0; i < lastElement; i++)
                input >> matrix1.matrix[i];
        return input;
}

//
// This friend function overloads
// the multiplication operator
// for scalar multiplication of the
// form:
//   scalar * matrix
//
const matrixClass operator*(const double& scalar, const matrixClass& matrix1)
{
        matrixClass product(matrix1.rows, matrix1.columns);
        for(unsigned i=(product.rows*product.columns); i-->0;)
                product.matrix[i] = matrix1.matrix[i] * scalar;
        return product;
}

//
// This friend function overloads
// the negation operator
// for negation of the form:
// matrix = -matrix1;
//
const matrixClass operator-(const matrixClass& matrix1)
{
        matrixClass negative(matrix1.rows, matrix1.columns);
        for(unsigned i=(negative.rows*negative.columns); i-->0;)
                negative.matrix[i] = matrix1.matrix[i] * -1;
        return negative;
}