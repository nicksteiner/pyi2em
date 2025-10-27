#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "em_i2em.h"
#include <vector>
#include <limits>

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

std::vector<double> I2EM_Bistatic(double freq_ghz, double rmsheight, double correl_length, 
                                double thi, double ths, double phs, double el, double ei, 
                                int correl_func, double xcoeff){
    /**
        Bistatic scattering from I2EM model
        
        Ulaby, F.T. and Long D.G.(2014), Microwave Radar and Radiometric Remote Sensing, The University of Michigan Press
        
        freq_ghz - frequency [GHz], double
        
        rmsheight - standard deviation of the surface height variation (rms) [m], double
        
        correl_length - correlation length [m], double
        
        thi - incident angle [degrees], double
        
        ths - scattered angle [degrees], double
        
        phs - scattered azimuth angle [degrees], double
        
        el - real relative permittivity (complex dielectric constant) of the surface, double
        
        ei - imaginary relative permittivity (complex dielectric constant) of the surface, double
        
        correl_func - type of surface correlation function:
                      (1) CORREL_EXPONENTIAL
                      (2) CORREL_GAUSSIAN
                      (3) CORREL_X_POWER
                      (4) CORREL_X_EXPONENTIAL, integer
        
        xcoeff - correlation coefficient for exponential power spectrum (used with CORREL_X_EXPONENTIAL), double
        
        Returns: [sigma0_hh, sigma0_vv] - backscatter coefficients in dB
     */

    cdouble er(el, ei);
    
    std::vector<double> sigma0(2);
    
    // Convert angles from degrees to radians
    double thi_rad = thi * DEG2RAD;
    double ths_rad = ths * DEG2RAD;
    double phs_rad = phs * DEG2RAD;
    
    I2EM_Bistat_model((float)freq_ghz, (float)rmsheight, (float)correl_length, 
                      thi_rad, ths_rad, phs_rad, er, 
                      correl_func, (float)xcoeff,
                      &sigma0[0], &sigma0[1]);
    
    return sigma0;
}

// Code 10.1: I2EM Backscattering from Single-Scale Random Surface
// Returns [sigma0_hh_dB, sigma0_vv_dB, sigma0_hv_dB or NaN]
std::vector<double> I2EM_Backscatter_Single(double freq_ghz, double rmsheight, double correl_length,
                                           double theta_deg, double el, double ei,
                                           int correl_func, double xcoeff, bool include_hv)
{
    cdouble er(el, ei);

    // Backscatter geometry: ths = thi, phs = 180 deg
    double thi_rad = theta_deg * DEG2RAD;
    double ths_rad = thi_rad;
    double phs_rad = 180.0 * DEG2RAD;

    std::vector<double> out(3);

    // HH and VV via core I2EM bistatic (in backscatter configuration)
    I2EM_Bistat_model((float)freq_ghz, (float)rmsheight, (float)correl_length,
                      thi_rad, ths_rad, phs_rad, er,
                      correl_func, (float)xcoeff,
                      &out[0] /*sigma0_hh_dB*/, &out[1] /*sigma0_vv_dB*/);

    // Optionally compute cross-pol HV using IEMX model
    if (include_hv) {
        double sigma0_hv_dB = std::numeric_limits<double>::quiet_NaN();
        IEMX_model((float)freq_ghz, (float)rmsheight, (float)correl_length,
                   thi_rad, er, correl_func, (float)xcoeff, 0, &sigma0_hv_dB);
        out[2] = sigma0_hv_dB;
    } else {
        out[2] = std::numeric_limits<double>::quiet_NaN();
    }

    return out;
}

// Code 10.2: I2EM Backscattering from Multi-Scale Random Surface (placeholder)
// Not implemented in core C++ sources; raise a runtime error with guidance.
std::vector<double> I2EM_Backscatter_Multi(double /*freq_ghz*/, double /*rms_small*/, double /*L_small*/, double /*rms_large*/, double /*L_large*/,
                                          double /*theta_deg*/, double /*el*/, double /*ei*/, int /*correl_func_small*/, int /*correl_func_large*/)
{
    throw std::runtime_error("I2EM multi-scale random surface backscatter is not implemented in the core library. "
                             "Consider using the periodic-surface mode (Code 10.4) as an approximation for deterministic slopes, "
                             "or open an issue to add a true two-scale model.");
}

// Code 10.4: I2EM Backscattering from Periodic Sinusoidal Surface
// Arguments: theta_0_deg, phi_0_deg, el, ei, freq_ghz, rmsheight_m, correl_length_m, correl_func, Gmm_m, A_m
// Returns: [sigma0_vv_dB, sigma0_hh_dB, sigma0_vh_dB]
std::vector<double> I2EM_Backscatter_Periodic(double theta_0_deg, double phi_0_deg,
                                              double el, double ei,
                                              double freq_ghz, double rmsheight, double correl_length,
                                              int correl_func,
                                              double Gmm, double A)
{
    cdouble eps(el, ei);
    double s_vv, s_hh, s_vh;
    I2EM_Periodic(theta_0_deg, phi_0_deg, eps,
                  freq_ghz, rmsheight, correl_length, correl_func,
                  Gmm, A,
                  &s_vv, &s_hh, &s_vh);
    return {s_vv, s_hh, s_vh};
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

    
    m.def("I2EM_Bistat", &I2EM_Bistat, R"pbdoc(
    
        Bistatic scattering from I2EM model

        Ulaby, F.T. and Long D.G.(2014), Microwave Radar and Radiometric Remote Sensing, The University of Michigan Press
        
        usage: I2EM_Bistat(freq_ghz, rmsheight, correl_length, thi, ths, phs, el, ei, correl_func, xcoeff)

        freq_ghz - frequency [GHz], double

        rmsheight - standard deviation of the surface height variation (rms) [m], double

        correl_length - correlation length [m], double

        thi - incident angle [degrees], double

        ths - scattered angle [degrees], double

        phs - scattered azimuth angle [degrees], double

        el - real relative permittivity (complex dielectric constant) of the surface, double
        
        ei - imaginary relative permittivity (complex dielectric constant) of the surface, double
    
        correl_func - type of surface correlation function:
                      (1) CORREL_EXPONENTIAL
                      (2) CORREL_GAUSSIAN
                      (3) CORREL_X_POWER
                      (4) CORREL_X_EXPONENTIAL, integer

        xcoeff - correlation coefficient for exponential power spectrum, double

        Returns: [sigma0_hh, sigma0_vv] - backscatter coefficients in dB

    )pbdoc",
    py::arg("freq_ghz"),
    py::arg("rmsheight"),
    py::arg("correl_length"),
    py::arg("thi"),
    py::arg("ths"),
    py::arg("phs"),
    py::arg("el"),
    py::arg("ei"),
    py::arg("correl_func") = CORREL_GAUSSIAN,
    py::arg("xcoeff") = 1.0);
    
    // Export correlation function constants
    m.attr("CORREL_EXPONENTIAL") = py::int_(CORREL_EXPONENTIAL);
    m.attr("CORREL_GAUSSIAN") = py::int_(CORREL_GAUSSIAN);
    m.attr("CORREL_X_POWER") = py::int_(CORREL_X_POWER);
    m.attr("CORREL_X_EXPONENTIAL") = py::int_(CORREL_X_EXPONENTIAL);

    // Single-scale backscatter
    m.def("I2EM_Backscatter_Single", &I2EM_Backscatter_Single, R"pbdoc(

        I2EM Backscattering from Single-Scale Random Surface (Code 10.1)

        usage: I2EM_Backscatter_Single(freq_ghz, rmsheight, correl_length,
                                       theta_deg, el, ei,
                                       correl_func=CORREL_GAUSSIAN, xcoeff=1.0, include_hv=False)

        Returns [sigma0_hh_dB, sigma0_vv_dB, sigma0_hv_dB or NaN]

        Notes:
        - Backscatter geometry is enforced internally (ths=thi, phs=180°)
        - When include_hv=True, the HV channel is computed using IEMX_model
    )pbdoc",
        py::arg("freq_ghz"),
        py::arg("rmsheight"),
        py::arg("correl_length"),
        py::arg("theta_deg"),
        py::arg("el"),
        py::arg("ei"),
        py::arg("correl_func") = CORREL_GAUSSIAN,
        py::arg("xcoeff") = 1.0,
        py::arg("include_hv") = false);

    // Multi-scale backscatter 
    m.def("I2EM_Backscatter_Multi", &I2EM_Backscatter_Multi, R"pbdoc(

        I2EM Backscattering from Multi-Scale Random Surface (Code 10.2)

        Not implemented in the current C++ core. This function will raise a runtime error.
        Please open an issue if you need this capability.
    )pbdoc",
        py::arg("freq_ghz"),
        py::arg("rms_small"),
        py::arg("L_small"),
        py::arg("rms_large"),
        py::arg("L_large"),
        py::arg("theta_deg"),
        py::arg("el"),
        py::arg("ei"),
        py::arg("correl_func_small"),
        py::arg("correl_func_large"));

    // Backscatter from periodic sinusoidal surface
    m.def("I2EM_Backscatter_Periodic", &I2EM_Backscatter_Periodic, R"pbdoc(

        I2EM Backscattering from Periodic Sinusoidal Surface (Code 10.4)

        usage: I2EM_Backscatter_Periodic(theta_0_deg, phi_0_deg, el, ei,
                                         freq_ghz, rmsheight, correl_length,
                                         correl_func, Gmm, A)

        Returns [sigma0_vv_dB, sigma0_hh_dB, sigma0_vh_dB]

        Units:
        - Angles in degrees
        - Lengths in meters (rmsheight, correl_length, Gmm period, A amplitude)
        - Frequency in GHz
    )pbdoc",
        py::arg("theta_0_deg"),
        py::arg("phi_0_deg"),
        py::arg("el"),
        py::arg("ei"),
        py::arg("freq_ghz"),
        py::arg("rmsheight"),
        py::arg("correl_length"),
        py::arg("correl_func") = CORREL_GAUSSIAN,
        py::arg("Gmm"),
        py::arg("A"));
}