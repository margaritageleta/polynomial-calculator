#include "Polynomial.hh"
//--------------------//
#include <stdlib.h>   //
#include <iostream>   //
#include <sstream>    //
#include <complex>    //
#include <cmath>      //
//--------------------//
using namespace std;
typedef vector<double> Coefficients;
typedef vector<complex<double>> iCoefficients;
//===================================================================//
//===================================================================//
// Constructor from vector.
Polynomial::Polynomial (const vector<double>& coeff)
: coeff(coeff)
{
    /* \DEBUGGING TOOL/
    int siz = coeff.size();
    if (siz != 0){
        for (int i = 0; i < siz; ++i) cout << coeff[i] << " ";
        cout << endl;
    }
    else cout << "ERROR!!!" << endl;
    */
}
//===================================================================//
//===================================================================//
// Write the elements of the vector of coefficients.
void Polynomial::write () const {
    int siz = int(coeff.size()); // Size of the vector of coefficients.
    if (siz != 0){
        for (int i = 0; i < siz; ++i) cout << coeff[i] << " ";
    }
    else cout << "(the polynomial is empty)";
}
//===================================================================//
//===================================================================//
/* OBS: if the vector of coeffs is empty, the degree will be -1, showing
 that the Polynomial is empty. */
// Get the degree of the polynomial.
int Polynomial::get_degree () const {
    return int(coeff.size()) - 1;
}
//===================================================================//
//===================================================================//
// Get the size of the vector of coefficients.
int Polynomial::get_size () const {
    return int(coeff.size());
}
//===================================================================//
//===================================================================//
// Get the coefficient of a given degree.
double Polynomial::get_coeff (int deg) const {
    return coeff[deg];
}
//===================================================================//
//===================================================================//
// Modify the coefficient of a given degree.
Polynomial& Polynomial::mod_coeff (int deg, double value) {
    coeff[deg] = value;
    return *this;
}
//===================================================================//
//===================================================================//
// Evaluate the polynomial on a certain x.
double Polynomial::operator() (double x) const {
    double eval = 0;
    // IN USE: Horner's rule.
    for (int i = get_degree(); i >= 0; --i) eval = eval*x + get_coeff(i);
    return eval;
}
//===================================================================//
//===================================================================//
// Update the polynomial's degree so that the leading
// coefficient is not zero.
Polynomial& Polynomial::normalize (){
    //===============================================================//
    // INVARIANT: while there are still zeroes ....000, do work.
    //              if not, return.
    //===============================================================//
    for (int i = get_degree(); i >= 0; --i) {
        if (get_coeff(i) == 0) coeff.pop_back();
        else return *this;
    }
    return *this;
}
//===================================================================//
//===================================================================//
// Perform the addition of two polynomials.
Polynomial Polynomial::operator+ (const Polynomial& poly) const {
    // /max_degree/ is the maximum possible degree of the result.
    int max_degree = max(get_degree(), poly.get_degree());
    Polynomial result;
    for (int i = 0; i <= max_degree; ++i) {
        // The element of the result-polynomial of degree /i/:
        double sum_element;
        //===========================================================//
        if (i > get_degree()) sum_element = poly.get_coeff(i);
        else if (i > poly.get_degree()) sum_element = get_coeff(i);
        else sum_element = get_coeff(i) + poly.get_coeff(i);
        //===========================================================//
        result.coeff.push_back(sum_element);
    }
    //===============================================================//
    // Normalize the result in case there are unnecessary zeroes.
    result.normalize();
    //===============================================================//
    return result;
}
//===================================================================//
//===================================================================//
// Perform the substraction of two polynomials.
Polynomial Polynomial::operator- (const Polynomial& poly) const {
    Coefficients coeffs; // Vector of coefficients.
    for (int i = 0; i <= poly.get_degree(); ++i) {
        // Change the sign of each element of the second Polynomial.
        coeffs.push_back(poly.get_coeff(i) * -1);
    }
    Polynomial p (coeffs); // Polynomial with changed elements.
    //===============================================================//
    // Now sum:
    Polynomial result = *this + p;
    return result;
}
//===================================================================//
//===================================================================//
// Perform the multiplication of two polynomials. NOT FOURIER.
Polynomial Polynomial::operator* (const Polynomial& poly) const {
    // The degree of the resulting Polynomial.
    int degree = get_degree() + poly.get_degree();
    Coefficients coeffs (degree + 1); // Vector of coefficients.
    //===============================================================//
    for (int i = 0; i <= get_degree(); ++i) {
        for (int j = 0; j <= poly.get_degree(); ++j) {
            coeffs[i + j] += get_coeff(i)*poly.get_coeff(j);
        }
    }
    //===============================================================//
    Polynomial result (coeffs);
    result.normalize();
    return result;
}
//===================================================================//
//===================================================================//
// Perform the division of two polynomials.
void Polynomial::divide (const Polynomial& dividend,
                         const Polynomial& divisor,
                         Polynomial& quotient,
                         Polynomial& remainder) {
    remainder = dividend;
    int size_dividend = dividend.get_size();
    int size_divisor = divisor.get_size();
    //===============================================================//
    // Case when degree of dividend (*this) < degree of divisor (poly).
    if (size_dividend < size_divisor) {
        quotient = Polynomial ();
        return;
    }
    //===============================================================//
    // Other cases:
    Coefficients q (size_dividend - size_divisor + 1, 0);
    quotient = Polynomial (q);
    remainder = dividend;
    // Index of the leading coefficient of R.
    int lead_coeff_R = remainder.get_size() - 1;
    //===============================================================//
    // Every iteration performs a step of the division.
    // INVARIANT: dividend = divisor * quotient + remainder.
    for (int lead_coeff_Q = quotient.get_size() - 1; lead_coeff_Q >= 0;
         --lead_coeff_Q) {
        // New coefficient of the quotient vector:
        quotient.coeff[lead_coeff_Q] = remainder.get_coeff(lead_coeff_R)
        / divisor.get_coeff(size_divisor - 1);
        // Guarantee an /exact/ zero at the leading coefficient:
        remainder.coeff[lead_coeff_R] = 0;
        for (int k = size_divisor - 2; k >= 0; --k) remainder.coeff[k + lead_coeff_Q]
            -= quotient.get_coeff(lead_coeff_Q) * divisor.get_coeff(k);
        // Index to the next leading coefficient of the quotient vector:
        --lead_coeff_R;
    }
    //===============================================================//
    // Normalize both polynomials in case there are unnecessary zeroes.
    quotient.normalize();
    remainder.normalize();
}
//===================================================================//
//===================================================================//
// Return the quotient of a division.
Polynomial Polynomial::operator/ (const Polynomial& divisor) {
    Polynomial quotient, remainder;
     // If divisor is empty.
    if (divisor.get_degree() < 0) return Polynomial ();
    divide(*this, divisor, quotient, remainder); // Divides.
    return quotient;
}
//===================================================================//
//===================================================================//
// Return the remainder of a division.
Polynomial Polynomial::operator% (const Polynomial& divisor) {
    Polynomial quotient, remainder;
    divide(*this, divisor, quotient, remainder); // Divides.
    return remainder;
}
//===================================================================//
//===================================================================//
// Convert to monic Polynomial.
Polynomial Polynomial::mono () {
    Polynomial p = *this;
    double c = p.coeff[p.get_degree()];
    for (int j = 0; j < p.get_size(); ++j) p.coeff[j] /= c;
    return p;
}
//===================================================================//
//===================================================================//
// Get the maximum common divisor of two polynomials.
Polynomial Polynomial::gcd (const Polynomial& poly) const {
    Polynomial dividend = *this;
    Polynomial divisor = poly;
    while (divisor.get_size() > 0) {
        Polynomial remainder = dividend % divisor;
        dividend = divisor;
        divisor = remainder;
    }
    dividend.mono(); // Convert to monic Polynomial.
    return dividend;
}
//===================================================================//
//===================================================================//
// It is the function that given a vector of coefficients prints
// the polynomial like: -2 3 0 -4 5 --> -2x^0 + 3x^1 - 4x^3 + 5x^4.
void Polynomial::toString () const {
    bool first = true; // If the coefficient is the first one.
    Coefficients coefficients = coeff; // Vector of coefficients.
    for (int j = 0; j < int(coeff.size()); ++j) {
        //===========================================================//
        if (get_coeff(j) == 0) continue; // Do nothing.
        //===========================================================//
        else if (get_coeff(j) < 0) { // Case for a negative coefficient.
            if (first) { // If it si first, print sign.
                cout << "-";
                first = false;
            }
            else cout << " - ";
            cout << get_coeff(j)*-1 << "x^" << j;
        }
        //===========================================================//
        else { // Case for a postive coefficient.
            if (first) first = false; // If it is first.
            else cout << " + ";
            cout << get_coeff(j) << "x^" << j;
        }
        //===========================================================//
    }
}
//===================================================================//
//===================================================================//
// The same function as toString() but it actually returns the string
// of the polynomial.
string Polynomial::toString2 () const {
    bool first = true; // If the coefficient is the first one.
    Coefficients coefficients = coeff; // Vector of coefficients.
    string poly_string = "";
    for (int j = 0; j < int(coeff.size()); ++j) {
        //===========================================================//
        if (get_coeff(j) == 0) continue; // Do nothing.
        //===========================================================//
        else if (get_coeff(j) < 0) { // Case for a negative coefficient.
            if (first) { // If it is first, print sign.
                poly_string.append("-");
                first = false;
            }
            else poly_string.append(" - ");
            poly_string.append(to_string(get_coeff(j)*-1));
            poly_string.append("x^");
            poly_string.append(to_string(j));
        }
        //===========================================================//
        else { // Case for a positive coefficient.
            if (first) first = false; // If it is first.
            else poly_string.append(" + ");
            poly_string.append(to_string(get_coeff(j)));
            poly_string.append("x^");
            poly_string.append(to_string(j));
        }
        //===========================================================//
    }
    return poly_string; // The string.
}
//===================================================================//
//===================================================================//
//             FUNCTIONS RELATED TO POLYCOMPLEX and FFT:
//-------------------------------------------------------------------//
// Constructor from ivector (vector<complex<double>>).
Polycomplex::Polycomplex (const iCoefficients& icoeff)
: icoeff(icoeff) {}
//===================================================================//
//===================================================================//
// Converts a Polynomial into Polycomplex;
Polycomplex Polynomial::toComplex () {
    // Real coefficients:
    Coefficients real_vector = coeff;
    // Complementary imaginary coefficients (zero at the moment):
    Coefficients im_vector (real_vector.size(), 0);
    /* -----------------STRUCTURE OF A POLYCOMPLEX-------------------//
     v: {|(R,I)|(R,I)|(R,I)|(R,I)|(R,I)|...|(R,I)|(R,I)|(R,I)|}
        where R is an element of /real_vector/,
        and I is an element of /im_vector/.
    //-------------------------------------------------------------- */
    iCoefficients complex_vector;
    for (int i = 0; i < real_vector.size(); ++i) {
        // /i_num/ is an element of the complex_vector of coefficients.
        // i_num = (R,I)
        complex<double> i_num(real_vector[i], im_vector[i]);
        complex_vector.push_back(i_num);
    }
    return Polycomplex (complex_vector);
}
//===================================================================//
//===================================================================//
// Get the size of the vector of coefficients.
int Polycomplex::get_size () const {
    return int(icoeff.size());
}
//===================================================================//
//===================================================================//
// Get the coefficient of a given degree.
complex<double> Polycomplex::get_coeff (int deg) const {
    return icoeff[deg];
}
//===================================================================//
//===================================================================//
// Calculate an appropriate power for /omega/ to execute FFT.
int Polynomial::power_2 (int n) {
    // Such that n is a power of 2.
    return pow(2, ceil(log(n)/log(2)));
}
//===================================================================//
//===================================================================//
// Perform the FFT, obtaining the point-value representation. FOURIER.
Polycomplex Polycomplex::FFT () {
    Polycomplex complex_vector = *this;
    int n = int (complex_vector.icoeff.size());
    //===============================================================//
    // If input (vector) contains just one element.
    if (n == 1) return Polycomplex (iCoefficients (1, complex_vector.get_coeff(0)));
    //===============================================================//
    // Store /n/ complex n-th roots of unity.
    iCoefficients omega (n);
    Polycomplex w(omega); // w stands for /omega/.
    for (int i = 0; i < n; i++) {
        double a = 2 * M_PI * i/n; // a stands for /alpha/
        // FORMULA: e^(ai) = cos(a) + i*sin(a)
        // where a = /alpha/ = (2k_pi_)/n
        // for k = 0, 1, ..., n - 1.
        w.icoeff[i] = complex<double>(cos(a), sin(a));
    }
    //===============================================================//
    Polycomplex Ae(iCoefficients (n/2)), Ao(iCoefficients (n/2));
    for (int i = 0; i < n/2; i++) {
        // Even-indexed coefficients.
        Ae.icoeff[i] = complex_vector.get_coeff(i * 2);
        // Odd-indexed coefficients.
        Ao.icoeff[i] = complex_vector.get_coeff(i * 2 + 1);
    }
    //===============================================================//
    // Recursive call for even-indexed coefficients.
    Polycomplex Re = Ae.FFT();
    // Recursive call for odd-indexed coefficients.
    Polycomplex Ro = Ao.FFT();
    // For storing values of R0, R1, R2, ..., Rn-1.
    iCoefficients icoeffs (n);
    Polycomplex R(icoeffs);
    //===============================================================//
    for (int k = 0; k < n/2; k++) { // Computing...
        R.icoeff[k] = Re.get_coeff(k) + w.get_coeff(k) * Ro.get_coeff(k);
        R.icoeff[k + n/2] = Re.get_coeff(k) - w.get_coeff(k) * Ro.get_coeff(k);
    }
    return R;
}
//===================================================================//
//===================================================================//
// Perform the iFFT. INVERSE FOURIER.
Polycomplex Polycomplex::iFFT () {
    Polycomplex complex_vector = *this;
    Polycomplex conjugate;
    //===============================================================//
    // /conj/ stands for conjugate.
    for (const auto &elem : complex_vector.icoeff) conjugate.icoeff.push_back(conj(elem));
    //===============================================================//
    Polycomplex vecfft = conjugate.FFT(); // FFT of the conjugate values.
    Polycomplex result;
    //===============================================================//
    /* OBS: as complex is composed of doubles, we use double(complex_vector.
     get_size()) in order to be able to multiply it with conj(elem). */
    for (const auto &elem : vecfft.icoeff) {
        result.icoeff.push_back(conj(elem)/double(complex_vector.get_size()));
    }
    //===============================================================//
    return result;
}
//===================================================================//
//===================================================================//
// Convert a Polycomplex into Polynomial.
Polynomial Polycomplex::toReal () {
    Polycomplex ivector = *this;
    // Real coefficients:
    Coefficients real_vector;
    //===============================================================//
    for (int i = 0; i < ivector.get_size(); ++i) {
        // /i_num/ is an element of the complex_vector of coefficients.
        // i_num = (R,I)
        double num = ivector.icoeff[i].real(); // We take /R/.
        real_vector.push_back(num); // And put it into the /real_vector/.
    }
    //===============================================================//
    return Polynomial (real_vector);
}
//===================================================================//
//===================================================================//
// Perform the multiplication of two polycomplexes.
Polycomplex Polycomplex::complex_mul (const Polycomplex& pcomplex) const {
    Polycomplex p = *this;
    iCoefficients complex_vector;
    //===============================================================//
    // Say R is the result vector, and A and B the operands, respectively,
    // Multiply the complex numbers A[j] * B[j] to obtain R[j].
    for (int j = 0; j < p.get_size(); ++j) {
        complex<double> i_num(p.icoeff[j]*pcomplex.icoeff[j]);
        complex_vector.push_back(i_num);
    }
    //===============================================================//
    Polycomplex result (complex_vector);
    return result;
}
//===================================================================//
//===================================================================//
// Clean the vector of coefficients of floating-point errors.
Polynomial& Polynomial::clean() {
    // IN USE: /epsilon/ = ±0.000001.
    for (int i = 0; i != get_size(); ++i) if (coeff[i] <= 0.000001 and
                                              coeff[i] >= -0.000001) coeff[i] = 0;
    //===============================================================//
    Polynomial result = *this;
    return result.normalize();
}
//===================================================================//
//===================================================================//
// Perform the multiplication of two polynomials. FOURIER.
Polynomial Polynomial::mul_fft (const Polynomial& poly) {
    Polynomial A = *this;
    Polynomial B = poly;
    //===============================================================//
    // Case when one of the Polynomials is empty:
    if (A.get_degree() < 0 or B.get_degree() < 0) {
        return Polynomial ();
    }
    //===============================================================//
    // The degree of the resulting Polynomial.
    int degree = power_2(A.get_degree() + B.get_degree());
    //===============================================================//
    // Pad zeroes:
    while (A.get_degree() != degree - 1) A.coeff.push_back(0);
    while (B.get_degree() != degree - 1) B.coeff.push_back(0);
    //===============================================================//
    // Perform FFT on complex A and B:
    Polycomplex iA = (A.toComplex()).FFT();
    Polycomplex iB = (B.toComplex()).FFT();
    //===============================================================//
    // Multiply complex A and B:
    Polycomplex iResult = iA.complex_mul(iB);
    //===============================================================//
    // Perform iFFT and convert Polycomplex to real Polynomial:
    Polynomial result = (iResult.iFFT()).toReal();
    //===============================================================//
    // Delete very small coefficients:
    result.clean();
    return result; // :)
}
//===================================================================//
//===================================================================//
