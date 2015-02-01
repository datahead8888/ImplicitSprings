//
// Curtis Glavin
// idealist@cs.bu.edu
//
// matrix.h
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

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

using namespace std;

//
// The following can be used to
// change the default dimensions
// of matrices.  These values
// are only used by the
// constructor when dimensions
// are not specified; they do
// not affect the rest of the
// class in any way.
//
#define DEFAULT_ROWS 3
#define DEFAULT_COLUMNS 3

//
// The following can be used
// to change the default matrix
// inverse method's precision.
// The method produces more
// accurate results but is slower
// for larger values.
//
#define DEFAULT_INVERSE_PRECISION 20

class matrixClass {

        //
        // This friend function overloads
        // the << operator so it can be used
        // to display the matrix on the
        // screen using standard ostreams.
        //
        friend ostream& operator<<(ostream&, const matrixClass&);

        //
        // This friend function overloads
        // the >> operator so it can be used
        // to input matrix data using standard
        // istreams.
        //
        friend istream& operator>>(istream&, matrixClass&);

        //
        // This friend function overloads
        // the multiplication operator
        // for scalar multiplication of the
        // form:
        //   scalar * matrix
        //
        friend const matrixClass operator*(const double&, const matrixClass&);


        //
        // This friend function overloads
        // the negation operator
        // for negation of the form:
        // matrix = -matrix1;
        //
        friend const matrixClass operator-(const matrixClass&);

public:

        //
        // Constructor
        // Dynamically allocates a matrix of
        // m*n doubles and initializes each
        // element to 0.
        //
        matrixClass(unsigned m = DEFAULT_ROWS, unsigned n = DEFAULT_COLUMNS);

        //
        // Copy Constructor
        // Dynamically allocates a new matrix
        // with the same dimensions and data
        // as the specified initializer matrix.
        //
        matrixClass(const matrixClass&);

        //
        // Destructor
        //
        ~matrixClass();

        //
        // Transpose function; will not
        // alter the original matrix but
        // will return a matrix which
        // is the transpose of the original.
        //
        const matrixClass transpose(void) const;

        //
        // Overloaded function call operator
        // which allows the caller to place
        // or see a value at the specified
        // row,column in the matrix.  The
        // overloaded function operator method
        // is itself overloaded to handle const
        // and non-const cases.
        //
        const double& operator()(unsigned m, unsigned n) const;
        double& operator()(unsigned m, unsigned n);

        //
        // Overloaded assignment operator
        // which allows the caller to assign
        // one matrix to another.
        //
        const matrixClass& operator=(const matrixClass&);

        //
        // Overloaded equivalence operators.
        // Matrices are considered equivalent
        // if and only if they are of the same
        // dimension and all elements at
        // corresponding locations in the 2
        // matrices are equivalent.
        //
        bool operator==(const matrixClass&) const;
        bool operator!=(const matrixClass& matrix1) const
                {return !(*this == matrix1);}

        //
        // Overloaded multiplication operator
        // for multiplying 2 matrices.  The
        // method will terminate the program
        // if the matrices' dimensions prevent
        // multiplication.
        //
        const matrixClass operator*(const matrixClass&) const;

        //
        // Overloaded multiplication
        // operator for scalar multiplication
        // of the form:
        //   matrix * scalar
        //
        const matrixClass operator*(const double&) const;

		//
        // Overloaded multiplication
        // operator for scalar division
        // of the form:
        //   matrix * scalar
        //
        const matrixClass operator/(const double&) const;

        //
        // Overloaded *= operator for
        // matrix multiplication of the form:
        // matrix *= matrix;
        // The method will terminate the program
        // if the matrices' dimensions prevent
        // multiplication.
        //
        const matrixClass& operator*=(const matrixClass&);

        //
        // Overloaded *= operator for
        // matrix multiplication of the form:
        // matrix *= scalar;
        //
        const matrixClass& operator*=(const double&);

        //
        // Overloaded addition operator
        // for the addition of matrices
        // of identical dimensions.
        // The method will terminate the
        // program if the matrices' dimensions
        // prevent addition.
        //
        const matrixClass operator+(const matrixClass&) const;

        //
        // Overloaded += operator for
        // matrix addition of the form:
        // matrix += matrix;
        // The method will terminate the program
        // if the matrices' dimensions prevent
        // addition.
        //
        const matrixClass& operator+=(const matrixClass&);

        //
        // Overloaded subtraction operator
        // for the subtraction of matrices
        // of identical dimensions.
        // The method will terminate the
        // program if the matrices' dimensions
        // prevent subtraction.
        //
        const matrixClass operator-(const matrixClass&) const;

        //
        // Overloaded -= operator for
        // matrix subtraction of the form:
        // matrix -= matrix;
        // The method will terminate the program
        // if the matrices' dimensions prevent
        // subtraction.
        //
        const matrixClass& operator-=(const matrixClass&);

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
        const matrixClass operator^(const int&) const;

        //
        // Overloaded ^= operator for raising
        // a square matrix to an integer power.
        // See the description of the overloaded
        // ^ operator above for details.
        //
        const matrixClass& operator^=(const int&);

        //
        // Methods (defined inline) to
        // return the number of rows
        // or columns to the user.
        //
        unsigned numberOfRows() const {return rows;}
        unsigned numberOfColumns() const {return columns;}

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
        unsigned currentInversePrecision() const;
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
        bool setInversePrecision(unsigned precision = DEFAULT_INVERSE_PRECISION);

//private: //good times
public:

        unsigned rows;
        unsigned columns;
        double* matrix;
        static unsigned inverse_precision;
};

#endif