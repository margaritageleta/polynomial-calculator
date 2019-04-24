# POLYNOMIAL CALCULATOR  

## INCLUDED
————————
*	All source files (Polynomial.hh, Polynomial.cc).
*	This README.
*	Makefile.
* 	Input example (with some extreme cases).

## INTRODUCTION
————————————
The Polynomial Calculator program is aimed to perform many mathematical operations with polynomials. It is a useful tool to compute tediously long computations on polynomials and it does implement many algorithms like Long Polynomial Division, Fast Fourier Transform for polynomial product, Horner’s method, among others, all of which simplify the dealing with polynomials. 

## FEATURES
————————
The Polynomial Calculator is based on two main classes: the Polynomial class and the Polycomplex class. The Polycomplex class is an exclusive and unique class created especially for this Polynomial Calculator.

## REQUIREMENTS
————————————
Source files included the following libraries:
*	stdlib.h		*	string		*	cmath		*	map		*	cassert
*	iostream		*	vector		*	stack		* 	complex		* 	sstream

## IMPLEMENTED CLASSES
———————————————————
The Polynomial class implements the following functions:
 - Create a polynomial.
 - Consult a polynomial.
 - Modify any term (monomial).
 - Evaluate it on a certain x.
 - Basic arithmetic operations: addition, subtraction, multiplication and division.
 - Get the maximum common divisor of two polynomials.
 - Represent a polynomial as (1) a string, (2) a vector of coefficients.

As said above, the Polycomplex class is an exclusive class, which deals with complex polynomials, i.e. polynomials with complex coefficients, which have a real part and an imaginary part. Basically it does implement all the functions that are related to FFT, iFFT (inverse Fast Fourier Transform) and point-value polynomial multiplication.

There are also methods which convert a Polynomial into a Polycomplex and viceversa.

## COMPILE & RUN
—————————————
Notation: ‘#’ is the prompt of the shell.
*	Compile with Makefile: 			# make
*	Run program: 				# ./main.exe
*	Run program with input example:		# ./main.exe < jocsproves.txt
* 	Exit with ctrl + C

## HOW TO USE?
———————————
Excited to get started? First of all, at any moment you can use the help() function in order to find out how all the functions work! Just like that:
#### help
You will instantly get a list with all the utilities of the calculator. You will see that all of them are sorted in four main groups, which are, namely:
 - INICIALIZATION FUNCTIONS: new, vector, delete
 - VISUALIZATION FUNCTIONS: see, print, coeff, degree
 - MODIFICATION FUNCTIONS: mod
 - ARITHMETIC FUNCTIONS: +, -, *, **, /
 - OTHER UTILITIES: help, eval, gcd
To inspect a particular function use ‘help:’ and type the function name after the function call. Example:
#### help: vector
Inicialization functions are the ones which can create and destroy objects. Which objects? The polynomials of course! ‘new’ is used to create a polynomial using the habitual mathematical notation in ascending order of degrees:
#### new P(x) 1x^0 - 2x^2 + 0.25x^3 + 45x^9
‘vector’ is used to create a polynomial using its coefficients in ascending order of degrees, but here we must type the polynomial’s degree before: in the example we have used 9 as the maximum degree, followed by the coefficients.
#### vector P(x) 9 1 0 -2 0.25 0 0 0 0 0 45
There I should warn you that every polynomial should have a unique name, in order to not to get confused. If you do not use distinct names, it can be a mess and the results can be wrong. Please pay attention when naming a new polynomial, check that the name is unique. What if you wanted to use a name of a polynomial which you do not use anymore? You can also destroy polynomials with ‘delete’. When you use delete upon a polynomial, it will immediately deleted from the data base (in our case, our dictionary):
#### delete P(x)
‘delete’ is an interactive function, and before deleting straight away the polynomial it will ask you whether you REALLY mean it, the following text will be displayed:
#### Are you sure? (y/n)
So you just need to print ‘y’ in affirmative case, and ’n’ otherwise.

Visualization functions are the ones that are related to get some information about a particular polynomial and consulting. ‘see’ prints all the coefficients of a polynomial.
#### see Q(x)
‘print’ outputs that polynomial in the beautiful mathematical notation.
#### print Q(x)
‘coeff’ shows the coefficient of a given degree.
#### coeff P(x) 3
‘degree’ returns the degree of the polynomial.
#### degree P(x)

There is only one modification function whose function is modify. Just as simple as that:
#### mod Q(x) 3 2
So what it is doing in this example is change the value of the coefficient of degree 3 to value = 2. So it useful when you have typed a long polynomial and made a mistake on one of coefficients. You do not have to the delete all the polynomial, just use mod!

Arithmetic functions always play with two polynomials and save the result in another one. So when using +, -, *, ** or /, you will have to type three names, two of existing polynomials and one of the resulting polynomial. If you use a name of a non-existing polynomial, you will be performing operations with an empty polynomial. To be honest, empty polynomials are not really useful… Well, let’s see a couple of examples:
#### + P(x) Q(x) R(x)
#### / P(x) Q(x) R(x)
#### ** P(x) Q(x) R(x)
Incidentally, have you noticed **? What is **? * is the brute force multiplication of the nasty O(n^2) running time. ** uses the point-value multiplication of polynomials and the Fast Fourier Transform (Polycomplex and iFFT). It is more efficient than * as evaluation and interpolation take O(blown) time and multiplication just O(n). For really huge polynomials, use **.

Other utilities of the calculator? You can compute the greatest common divisor of two polynomials with ‘gcd’:
#### gcd P(x) Q(x) R(x)
Another greatest feature: you can evaluate any value on any polynomial with ‘eval’:
#### eval P(x) -0.0025
And of course, use help().

Enjoy!

UPC C++ <3
