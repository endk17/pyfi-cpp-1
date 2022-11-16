#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <Python.h>

// when C compilers don't define PI 
#ifndef PI 
#define PI 3.141592653589793238462643
#endif
 
#ifndef PI 
const double PI=3.141592653589793238462643;
#endif

// normal distribution function
double n(double z) {
    return (1.0/sqrt(2.0*PI))*exp(-0.5*z*z);
}

// cumulative normal
double N(double z) {
    if (z &gt; 6.0) { return 1.0; }; // this guards against overflow
    if (z &lt; -6.0) { return 0.0; };
    
    double b1 =  0.31938153;
    double b2 = -0.356563782;
    double b3 =  1.781477937;
    double b4 = -1.821255978;
    double b5 =  1.330274429;
    double p  =  0.2316419;
    double c2 =  0.3989423;
    
    double a = fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z &lt; 0.0 ) n = 1.0 - n;
    return n;
}

// black scholes call
double _bs_call(double S, double K, double r, double t, double sigma) {
    double time_sqrt = sqrt(t);
    double d1 = (log(S/K)+r*t)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*t)*N(d2);
}

// Module objects for python
static PyObject *
bs_call(PyObject *self, PyObject *args)
{
    double S, K, r, t, sigma;
    if (!PyArg_ParseTuple(args, "ddddd", &amp;S, &amp;K, &amp;r, &amp;t, &amp;sigma))
        return NULL;
    return Py_BuildValue("d", _bs_call(S, K, r, t, sigma));
}