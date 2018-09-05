#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "em_i2em.h"
#include <vector>


std::vector<double> I2EM(double fr, double sig, double l, double theta_d, double el, double ei, double sp){
    /**
        Emissivity from I2EM model
        
        Ulaby, F.T. and Long D.G.(2014), Microwave Radar and Radiometric Remote Sensing, The University of Michigan Press
        
        fr - frequency [ghz], double

        sig - standard deviation of the surface height variation (rms) [m], double

        l - correlation length [m], double

        theta_d - incidence angle [deg], double

        el - real relative permittivity (complex dielectric constant) of the surface, double
        
        ei - imaginary relative permittivity (complex dielectric constant) of the surface, double
    
        sp - type of surface correlation function (1) exponential (2) Gaussian, integer
     */

    cdouble er;
    er = cdouble(el, ei);

    std::vector<double> e(2);

    Calc_emissivity_I2EMmodel(fr, sig, l, theta_d, er, sp, &e[0], &e[1]);

    return e;
};

std::vector<double> test_I2EM(){
    // set frequency in ghz
    double fr = 3.0;

    double el = 11.3;
    double ei = 1.5;

    // set correlation length [m]
    double l = 0.10; // 10 cm

    // set standard deviation of the surface height variation (rms) [m]
    double sig = 0.0025; // .25 cm

    // incidence angle [deg]
    double theta_d = 30.0;

    // type of surface correlation function (1) exponential (2) Gaussian
    int sp = 2;

    vector<double> e;

    e = I2EM(fr, sig, l, theta_d, el, ei, sp);
    
    return e;
}


namespace py = pybind11;

PYBIND11_MODULE(pyi2em, m) {
    m.doc() = "I2EM library";
    m.def("I2EM", &I2EM, R"pbdoc(
    
        Emissivity from I2EM model

        Ulaby, F.T. and Long D.G.(2014), Microwave Radar and Radiometric Remote Sensing, The University of Michigan Press
        
        usage: I2EM(fr, sig, l, theta_d, el, ei, sp)

        fr - frequency [ghz], double

        sig - standard deviation of the surface height variation (rms) [m], double

        l - correlation length [m], double

        theta_d - incidence angle [deg], double

        el - real relative permittivity (complex dielectric constant) of the surface, double
        
        ei - imaginary relative permittivity (complex dielectric constant) of the surface, double
    
        sp - type of surface correlation function (1) exponential (2) Gaussian, integer

    )pbdoc");

    m.def("test_I2EM", &test_I2EM, "test I2EM - pybind11");
}