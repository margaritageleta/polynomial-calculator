//===================================================================//
//===================================================================//
/* The Polynomial class stores a polynomial and performs many
 mathematical operations with it.
 A basic polynomial is of the following form: p(x) = a0 + a1x + a2x^2
 + ... + adx^d
where /a/'s are the coefficients of the polynomial and /d/ is the degree
 of the polynomial.
Note that the /a/'s subindex /i/ corresponds to the exponent of the /x/
 that is multiplying that coefficient.
The Polynomial class is aimed to do the following functions which will
 be all implemented in \methods/:
 - Create a polynomial.
 - Consult a polynomial.
 - Modify any term (monomial).
 - Evaluate it on a certain x.
 - Basic arithmetic operations: addition, subtraction, multiplication
    and division.
 - Get the maximum common divisor of two polynomials.
 - Represent a polynomial as (1) a string, (2) a vector of coefficients.
 */
//===================================================================//
//===================================================================//
#ifndef Polynomial_hh
#define Polynomial_hh
//--------------------//
#include <complex>    //
#include <vector>     //
#include <string>     //
//--------------------//
using namespace std;
typedef vector<double> Coefficients; // Coefficients of a Polynomial.
typedef vector<complex<double>> iCoefficients; // Complex coefficients.

class Polycomplex; // Complex polynomial...
//===================================================================//
//===================================================================//
class Polynomial { // Real polynomial...
    
public:
//===================================================================//
    // INICIALIZATION FUNCTIONS ://
    /*---------------------------------------------------------------*/
    // Trivial empty constructor. To be able to clear the map.
    Polynomial () {}
    // Constructor from vector.
    Polynomial (const vector<double>& coeff);
//===================================================================//
    // VISUALIZATION FUNCTIONS ://
    /*---------------------------------------------------------------*/
    // Write the elements of the vector of coefficients.
    void write () const;
    /*---------------------------------------------------------------*/
    // Get the degree of the polynomial.
    int get_degree () const;
    /*---------------------------------------------------------------*/
    // Get the size of the vector of coefficients.
    int get_size () const;
    /*---------------------------------------------------------------*/
    // Get the coefficient of a given degree.
    double get_coeff (int deg) const;
    /*---------------------------------------------------------------*/
    // Transforms a Polynomial to a poly-string.
    void toString () const; // with cout's.
    /*---------------------------------------------------------------*/
    // Transforms a Polynomial to a poly-string.
    string toString2 () const; // returning a string.
//===================================================================//
    // MODIFICATION FUNCTIONS ://
    /*---------------------------------------------------------------*/
    // Modify the coefficient of a given degree.
    Polynomial& mod_coeff (int deg, double value);
//===================================================================//
    // ARITHMETIC FUNCTIONS ://
    /*---------------------------------------------------------------*/
    // Perform the addition of two polynomials.
    Polynomial operator+ (const Polynomial& poly) const;
    /*---------------------------------------------------------------*/
    // Perform the substraction of two polynomials.
    Polynomial operator- (const Polynomial& poly) const;
    /*---------------------------------------------------------------*/
    // Perform the multiplication of two polynomials.
    Polynomial operator* (const Polynomial& poly) const;
    /*---------------------------------------------------------------*/
    // Perform the division of two polynomials.
    void divide (const Polynomial& dividend, const Polynomial& divisor,
                 Polynomial& quotient, Polynomial& remainder);
    /*---------------------------------------------------------------*/
    // Return the quotient of a division.
    Polynomial operator/ (const Polynomial& divisor);
    /*---------------------------------------------------------------*/
    // Perform the multiplication of two polynomials. FOURIER.
    Polynomial mul_fft (const Polynomial& poly);
//===================================================================//
    // OTHER UTILITIES ://
    /*---------------------------------------------------------------*/
    // Evaluate the polynomial on a certain x.
    double operator() (double x) const;
    /*---------------------------------------------------------------*/
    // Get the maximum common divisor of two polynomials.
    Polynomial gcd (const Polynomial& poly) const;
    /*---------------------------------------------------------------*/
    // Convert to complex Polynomial.
    Polycomplex toComplex ();
    /*---------------------------------------------------------------*/
    // Updates the polynomial's degree so that the leading
    // coefficient is not zero.
    Polynomial& normalize ();

private:
//===================================================================//
    // INICIALIZATION OF PRIVATE FIELDS ://
    vector<double> coeff; // this vector stores each coeff. per degree.
    // string spoly; // this is a poly-string.
//===================================================================//
    // PRIVATE METHODS ://
    /*---------------------------------------------------------------*/
    // Get the maximum common divisor of two polynomials.
    Polynomial operator% (const Polynomial& divisor);
    /*---------------------------------------------------------------*/
    // Convert to monic Polynomial.
    Polynomial mono ();
    /*---------------------------------------------------------------*/
    // Calculate an appropriate power for /omega/ to execute FFT.
    int power_2 (int n);
    /*---------------------------------------------------------------*/
    // Clean the vector of coefficients of floating-point errors.
    Polynomial& clean();
};
//===================================================================//
//===================================================================//
class Polycomplex { // The class of the Complex Polynomial...
public:
    
    // INICIALIZATION FUNCTIONS ://
    /*---------------------------------------------------------------*/
    // Trivial empty constructor. To be able to clear the map.
    Polycomplex () {}
    // Constructor from vector.
    Polycomplex (const iCoefficients& icoeff);
//===================================================================//
    // Get the size of the vector of coefficients.
    int get_size () const;
    /*---------------------------------------------------------------*/
    // Get the coefficient of a given degree.
    complex<double> get_coeff (int deg) const;
    /*---------------------------------------------------------------*/
    // Perform the FFT, obtaining the point-value representation. FOURIER.
    Polycomplex FFT ();
    /*---------------------------------------------------------------*/
    // Perform the iFFT. INVERSE FOURIER.
    Polycomplex iFFT ();
    /*---------------------------------------------------------------*/
    // Convert a Polycomplex into Polynomial.
    Polynomial toReal ();
    /*---------------------------------------------------------------*/
    // Perform the multiplication of two polycomplexes.
    Polycomplex complex_mul (const Polycomplex& pcomplex) const;
    
private:
//===================================================================//
    // INICIALIZATION OF PRIVATE FIELDS ://
    // Coefficients with Real and Complex (Imaginary) parts.
    iCoefficients icoeff;
};
//===================================================================//
#endif


















