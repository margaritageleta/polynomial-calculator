#include "Polynomial.hh"
//--------------------//
#include <stdlib.h>   //
#include <iostream>   //
#include <sstream>    //
#include <complex>    //
#include <cassert>    //
#include <string>     //
#include <cmath>      //
#include <stack>      //
#include <map>        //
//--------------------//
using namespace std;
typedef vector<double> Coefficients;

// CONVERSION FUNCTIONS ://
//*******************************************************************//
//===================================================================//
// Extract a polynomial term from a string.
Coefficients extract_term(string poly_term, vector<double>& coefficients, bool& neg_coeff) {
    // Get the coefficient.
    string coeff_str = "";
    int i;
    for (i = 0; poly_term[i] != 'x'; i++)
        coeff_str.push_back(poly_term[i]);
    double coeff = atof(coeff_str.c_str()); // atof: str to double. c_str: pointer to string.
    if (neg_coeff) coeff *= -1;
    //===============================================================//
    // Get the power (skipping 2 characters for x and ^).
    string pow_str = "";
    for (i = i + 2; i != int(poly_term.size()); i++)
        pow_str.push_back(poly_term[i]);
    int power = atoi(pow_str.c_str()); // atoi: str to int.
    coefficients[power] = coeff;
    return coefficients;
}
//===================================================================//
//===================================================================//
// Convert an input string to a vector of coeffs.
vector<double> toVector(string& poly_string, int degree) {
    istringstream iss(poly_string);
    Coefficients coefficients(degree, 0);
    bool neg_coeff = false;
    //===============================================================//
    string poly_term;
    while (iss >> poly_term) {
        if (poly_term == "+") { // If positive.
            neg_coeff = false;
            continue;
        }
        else if (poly_term == "-") { // If negative.
            neg_coeff = true;
            continue;
        }
        // Extract term and convert to double.
        else extract_term(poly_term, coefficients, neg_coeff);
    }
    return coefficients;
}
//===================================================================//
//===================================================================//
// Given a poly-string, find the maximum degree.
int find_degree_string(const string& str) {
    stack<char> temp;
    for (int j = int(str.size()) - 1; str[j] != '^'; --j) {
        temp.push(str[j]);
    }
    //===============================================================//
    string deg = "";
    while (not temp.empty()) {
        deg.push_back(temp.top());
        temp.pop();
    }
    int degree = atoi(deg.c_str());
    return degree;
}
//===================================================================//
//******************************************************************************************************//
//===================================================================//
// Help function, especially for the user, with love.
void help () {
    cout << "POLYNOMIAL CALCULATOR -> LIST OF FUNCTIONS:" << endl;
    cout << "==========================================" << endl;
    cout << "INICIALIZATION FUNCTIONS: new, vector, delete." << endl;
    cout << "VISUALIZATION FUNCTIONS: see, print, coeff, degree." << endl;
    cout << "MODIFICATION FUNCTIONS: mod." << endl;
    cout << "ARITHMETIC FUNCTIONS: +, -, *, **, /." << endl;
    cout << "OTHER UTILITIES: help, eval, gcd." << endl;
}
//===================================================================//
//===================================================================//
void help_more () {
    string f; // Name of the requested function.
    cin >> f;
    cout << "POLYNOMIAL CALCULATOR -> LIST OF FUNCTIONS -> "  << f << ":" << endl;
    cout << "==========================================" << endl;
    if (f == "help") cout << "Print \"help\" in order to see the functions of the calculator." << endl;
    else if (f == "new") {
        cout << "Print \"new (name)\" to create a polynomial. " << endl;
        cout << endl;
        cout << "You must introduce the polynomial coefficients with their respective degree with the following form:" << endl;
        cout << endl;
        cout << "new P(x) (coefficient)(x)(^)(degree)()(operand if continued). Repeat." << endl;
        cout << endl;
        cout << "EXAMPLE:" << endl;
        cout << "new P(x) -4x^0 - 5x^2 + 98x^8 - 0.1413x^10" << endl;
        cout << endl;
        cout << "OBS: Note that the degrees are introduced in ascendent order," << endl;
        cout << "otherwise, it will not work correctly." << endl;
    }
    else if (f == "vector") {
        cout << "Print \"vector (name) (degree)\" to create a polynomial. " << endl;
        cout << endl;
        cout << "You must introduce the polynomial coefficients with the following form:" << endl;
        cout << endl;
        cout << "vector P(x) (degree) (coefficient 0) ... (coefficient degree)" << endl;
        cout << endl;
        cout << "EXAMPLE:" << endl;
        cout << "vector P(x) 10 -4 0 -5 0 0 0 0 0 98 0 -0.1413" << endl;
        cout << endl;
        cout << "OBS: Note that the coefficients are introduced in ascendent order of degrees," << endl;
        cout << "otherwise, it will not work correctly." << endl;
    }
    else if (f == "see") cout << "Print \"see (name)\" in order to see the coefficients of the polynomial." << endl;
    else if (f == "coeff") cout << "Print \"coeff (name) (degree)\" in order to see a coefficient of a given degree." << endl;
    else if (f == "mod") cout << "Print \"mod (name) (degree) (value)\" in order to modify the coefficient of a given degree." << endl;
    else if (f == "print") cout << "Print \"print (name)\" in order to print the polynomial of a given name." << endl;
    else if (f == "degree") {
        cout << "Print \"degree (name)\" in order to see the degree of the polynomial." << endl;
        cout << endl;
        cout << "OBS: Note that if degree is -1 it means that the Polynomial is empty." << endl;
    }
    else if (f == "eval") cout << "Print \"eval (name) (value)\" in order to evaluate the polynomial on a given value." << endl;
    else if (f == "+") {
        cout << "Print \"+ (name_1) (name_2) (name_3)\" in order to get the sum of name_1 polynomial" << endl;
        cout << "and name_2 polynomial and save the result in name_3 polynomial." << endl;
    }
    else if (f == "-") {
        cout << "Print \"- (name_1) (name_2) (name_3)\" in order to subtract name_2 polynomial" << endl;
        cout << "from name_1 polynomial and save the result in name_3 polynomial." << endl;
    }
    else if (f == "*") {
        cout << "Print \"* (name_1) (name_2) (name_3)\" in order to multiply name_1 polynomial" << endl;
        cout << "and name_2 polynomial and save the result in name_3 polynomial." << endl;
    }
    else if (f == "**") {
        cout << "Print \"** (name_1) (name_2) (name_3)\" in order to multiply name_1 polynomial" << endl;
        cout << "and name_2 polynomial and save the result in name_3 polynomial." << endl;
        cout << endl;
        cout << "OBS: Note that ** uses the Fast Fourier Transform to perform the opertion." << endl;
    }
    else if (f == "/") {
        cout << "Print \"/ (name_1) (name_2) (name_3)\" in order to divide name_1 polynomial" << endl;
        cout << "and name_2 polynomial and save the result in name_3 polynomial." << endl;
    }
    else if (f == "gcd") {
        cout << "Print \"gcd (name_1) (name_2) (name_3)\" in order to get the Greatest Common" << endl;
        cout << "Divisor of name_1 polynomialand name_2 polynomial and save the ith in name_3 polynomial." << endl;
    }
    else if (f == "delete") cout << "Print \"delete (name)\" in order to delete the polynomial." << endl;
    else cout << "There is no such a function." << endl;
}
//===================================================================//
//===================================================================//
// Define a polynomial. Example:
/*  Input: P(x) 1.89x^0 + 4x^2 - 0.1212x^8
    Output: 1.89 0 4 0 0 0 0 0 0.1212 (Not seen by user) */
void def_poly (map<string, Polynomial>& dictionary) {
    string def; // Its name.
    cin >> def;
    //===============================================================//
    cin.ignore(); // Clean the cin entry in order to execute getline.
    string str;
    getline(cin, str);
    int deg = find_degree_string(str); // Its degree.
    Coefficients coefficients = toVector(str, deg + 1);
    //===============================================================//
    Polynomial p (coefficients);
    p.normalize();
    dictionary.insert({def, p});
}
//===================================================================//
//===================================================================//
// Define a polynomial using its vector of coefficients. Example:
/*  Input: P(x) 8 1.89 0 4 0 0 0 0 0 0.1212
 Output: 1.89x^0 + 4x^2 - 0.1212x^8 (Not seen by user) */
void def_vec (map<string, Polynomial>& dictionary){
    string def;                                    // Name.
    int degree;                                    // Its degree.
    cin >> def >> degree;
    Coefficients coefficients;
    //===============================================================//
    for (int i = 0; i <= degree; ++i) {
        double coeff;
        cin >> coeff;
        coefficients.push_back(coeff);
    }
    //===============================================================//
    Polynomial p (coefficients);
    p.normalize();
    dictionary.insert({def, p});
}
//===================================================================//
//===================================================================//
// Delete a polynomial from the dictionary.
void del (map<string, Polynomial>& dictionary){
    string def;
    cin >> def;
    if (dictionary.count(def)) {
        cout << "Are you sure? (y/n)" << endl;
        char c;
        cin >> c;
        if (c == 'y') {
            cout << def << " was deleted." << endl;
            dictionary.erase(def);
        }
        else return;
    }
    else cout << "Polynomial " << def << " is not in the dictionary." << endl;
}
//===================================================================//
//===================================================================//
// Print the polynomial. Example:
/*  Input: P(x)
 Output: 1.89x^0 + 4x^2 - 0.1212x^8 */
void print_poly (map<string, Polynomial>& dictionary){
    string def;                                     // Name.
    cin >> def;
    //===============================================================//
    /* We can also store the result in a string like this:
    string s = dictionary[def].toString2();
    cout << s << endl;
    That way allows us to use the resulting string in other
    functions or to store it in a dictionary of strings. If needed.
     
    The composition of toString2() and toVector() is the identity.*/
    //===============================================================//
    // Or do it with the cout's: :
    dictionary[def].toString();
    cout << endl;
}
//===================================================================//
//===================================================================//
// See the polynomial's coefficients. Example:
/*  Input: P(x)
 Output: 1.89 0 4 0 0 0 0 0 0.1212 */
void see (map<string, Polynomial>& dictionary) {
    string def;                                     // Name.
    cin >> def;
    cout << "The coefficients of the polynomial " << def << " are ";
    dictionary[def].write();
    cout << "." << endl;
}
//===================================================================//
//===================================================================//
// See a particular coefficient. Example:
/*  Input: P(x) 8
 Output: 0.1212 */
void coeff (map<string, Polynomial>& dictionary) {
    string def;                                     // Name.
    int deg;                                        // Its degree.
    cin >> def >> deg;
    if (deg <= dictionary[def].get_degree()) {
        cout << "The coefficient of degree " << deg << " in polynomial " << def << " is ";
        cout << dictionary[def].get_coeff(deg) << endl;
    }
    else cout << "WARNING: the polynomial's degree is " << dictionary[def].get_degree() << "!" << endl;
}
//===================================================================//
//===================================================================//
// Get the polynomial's degree. Example:
/*  Input: P(x)
 Output: 8 */
void degree (map<string, Polynomial>& dictionary){
    string def;                                     // Name.
    cin >> def;
    cout << "The degree of the polynomial " << def << " is ";
    cout << dictionary[def].get_degree() << endl;
}
//===================================================================//
//===================================================================//
// Modify the coefficient of a given degree. Example:
/*  Input: P(x) 2 -3
 Output: 1.89x^0 - 3x^2 - 0.1212x^8 (Not seen by user) */
void mod (map<string, Polynomial>& dictionary){
    string def;                                     // Name.
    int deg;                                        // Degree.
    double value;                                   // New value.
    cin >> def >> deg >> value;
    double old = dictionary[def].get_coeff(deg); // old value.
    dictionary[def].mod_coeff(deg, value);
    cout << old << " has been changed to " << dictionary[def].get_coeff(deg) << endl;
}
//===================================================================//
//===================================================================//
// Evaluate the polynomial on a certain x. Example:
/*  Input: P(x) 3
 Output: 1.89*3^0 + 4*3^2 - 0.1212*3^8 = 833.083 */
void eval (map<string, Polynomial>& dictionary){
    string def;                                     // Name.
    double x;                                       // Value of /x/.
    cin >> def >> x;
    cout << "The polynomial " << def << " evaluated on " << x << " gives the result ";
    cout << dictionary[def](x) << endl;
    // Could implement a solver by steps, if had time;
}
//===================================================================//
//===================================================================//
// Perform the addition of two polynomials.
void sum (map<string, Polynomial>& dictionary){
    string def1, def2, def3;
    cin >> def1 >> def2 >> def3;
    cout << def1 << " + " << def2 << " = " << def3 << endl;
    cout << "(";
    dictionary[def1].toString();
    cout << ") + (";
    dictionary[def2].toString();
    cout << ") = (";
    dictionary.insert({def3, dictionary[def1] + dictionary[def2]});
    dictionary[def3].toString();
    cout << ")" << endl;
}
//===================================================================//
//===================================================================//
// Perform the substraction of two polynomials.
void subtract (map<string, Polynomial>& dictionary){
    string def1, def2, def3;
    cin >> def1 >> def2 >> def3;
    cout << def1 << " - " << def2 << " = " << def3 << endl;
    cout << "(";
    dictionary[def1].toString();
    cout << ") - (";
    dictionary[def2].toString();
    cout << ") = (";
    dictionary.insert({def3, dictionary[def1] - dictionary[def2]});
    dictionary[def3].toString();
    cout << ")" << endl;
}
//===================================================================//
//===================================================================//
// Perform the multiplication of two polynomials.
void multiply (map<string, Polynomial>& dictionary) {
    string def1, def2, def3;
    cin >> def1 >> def2 >> def3;
    cout << def1 << " * " << def2 << " = " << def3 << endl;
    cout << "(";
    dictionary[def1].toString();
    cout << ") * (";
    dictionary[def2].toString();
    cout << ") = (";
    dictionary.insert({def3, dictionary[def1] * dictionary[def2]});
    dictionary[def3].toString();
    cout << ")" << endl;
}
//===================================================================//
//===================================================================//
// Perform the multiplication of two polynomials using FOURIER!
void fast_multiply (map<string, Polynomial>& dictionary) {
    string def1, def2, def3;
    cin >> def1 >> def2 >> def3;
    cout << def1 << " * " << def2 << " = " << def3 << endl;
    cout << "(";
    dictionary[def1].toString();
    cout << ") * (";
    dictionary[def2].toString();
    cout << ") = (";
    dictionary.insert({def3, dictionary[def1].mul_fft(dictionary[def2])});
    dictionary[def3].toString();
    cout << ")" << endl;
}
//===================================================================//
//===================================================================//
// Perform the division of two polynomials.
void divide (map<string, Polynomial>& dictionary) {
    string def1, def2, def3;
    cin >> def1 >> def2 >> def3;
    cout << def1 << " / " << def2 << " = " << def3 << endl;
    cout << "(";
    dictionary[def1].toString();
    cout << ") / (";
    dictionary[def2].toString();
    cout << ") = (";
    dictionary.insert({def3, dictionary[def1] / dictionary[def2]});
    dictionary[def3].toString();
    cout << ")" << endl;
}
//===================================================================//
//===================================================================//
// Compute the gcd of two polynomials.
void gcd (map<string, Polynomial>& dictionary){
    string def1, def2, def3;
    cin >> def1 >> def2 >> def3;
    cout << "The gcd of ";
    dictionary[def1].toString();
    cout << " and ";
    dictionary[def2].toString();
    cout << " is:" << endl;
    dictionary.insert({def3, dictionary[def1].gcd(dictionary[def2])});
    dictionary[def3].toString();
    cout << endl;
}
//===================================================================//
//*******************************************************************//
int main () {
    // \Inicialization of the polynomial calculator/
    cout << "*******************************************************************************" << endl;
    cout << "POLYNOMIAL CALCULATOR" << endl;
    cout << "Print \"help\" to see the list of the functions." << endl;
    cout << "Print \"help: (function name)\" to see more details about a requested function." << endl;
    cout << "N.B. Use a distinct name for every Polynomial!" << endl;
    map<string, Polynomial> dictionary;
    cout << "*******************************************************************************" << endl;
    //===============================================================//
    string action;
    while (cin >> action) {                                          //
        if (action == "help")           help();                      //
        else if (action == "help:")     help_more();                 //
        /*-----------------------------------------------------------*/
        else if (action == "new")       def_poly (dictionary);       //
        else if (action == "vector")    def_vec (dictionary);        //
        /*-----------------------------------------------------------*/
        else if (action == "see")       see (dictionary);            //
        else if (action == "print")     print_poly (dictionary);     //
        else if (action == "coeff")     coeff (dictionary);          //
        else if (action == "degree")    degree (dictionary);         //
        /*-----------------------------------------------------------*/
        else if (action == "mod")       mod (dictionary);            //
        /*-----------------------------------------------------------*/
        else if (action == "+")         sum (dictionary);            //
        else if (action == "-")         subtract (dictionary);      //
        else if (action == "*")         multiply (dictionary);       //
        else if (action == "**")        fast_multiply (dictionary);  //
        else if (action == "/")         divide (dictionary);         //
        /*-----------------------------------------------------------*/
        else if (action == "eval")      eval (dictionary);           //
        else if (action == "gcd")       gcd (dictionary);            //
        /*-----------------------------------------------------------*/
        else if (action == "delete")    del (dictionary);            //
        /*-----------------------------------------------------------*/
        else assert(false);
        cout << "*******************************************************************************" << endl;
    }
}
