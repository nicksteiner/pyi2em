#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <complex>
#include <stdexcept>
#include <cmath>
#include "em_i2em.h"

namespace py = pybind11;
using namespace pybind11::literals;

// Map human-readable correl func to legacy ints
static int correl_to_int(const std::string& s){
  if(s=="exponential") return CORREL_EXPONENTIAL;
  if(s=="gaussian") return CORREL_GAUSSIAN;
  if(s=="x_power") return CORREL_X_POWER;
  if(s=="x_exponential") return CORREL_X_EXPONENTIAL;
  throw std::invalid_argument("Unknown correlation function: "+s);
}

static void validate_common(double freq_ghz, double rms_m, double L_m, double theta_deg){
  if(freq_ghz <= 0) throw std::invalid_argument("freq_ghz must be > 0");
  if(rms_m < 0) throw std::invalid_argument("rms_height_m must be >= 0");
  if(L_m <= 0) throw std::invalid_argument("corr_length_m must be > 0");
  if(theta_deg < 0 || theta_deg > 89.999) throw std::invalid_argument("theta_deg must be in [0, 90)");
}

PYBIND11_MODULE(pyi2em, m) {
  m.doc() = "Python bindings for the I2EM surface scattering/emissivity model";

  m.def("emissivity",
    [](double freq_ghz, double rms_height_m, double corr_length_m, double theta_deg,
       std::complex<double> er_complex, const std::string& correl, bool return_linear){
        validate_common(freq_ghz, rms_height_m, corr_length_m, theta_deg);
        int sp = correl_to_int(correl);

        double eh=0.0, ev=0.0;
        // Core expects meters for s and L and does internal cm conversion
        cdouble er(er_complex.real(), er_complex.imag());
        Calc_emissivity_I2EMmodel(freq_ghz, rms_height_m, corr_length_m, theta_deg, er, sp, &eh, &ev);

        // emissivity is already linear 0..1
        (void)return_linear; // kept for API symmetry
        return py::make_tuple(eh, ev);
    },
    py::arg("freq_ghz"), py::arg("rms_height_m"), py::arg("corr_length_m"),
    py::arg("theta_deg"), py::arg("er_complex"), py::arg("correl")="exponential",
    py::arg("return_linear")=true,
    "Rough-surface emissivity (eh, ev). Units: GHz, m, m, deg, complex(er).");

  m.def("sigma0_backscatter",
    [](double freq_ghz, double rms_height_m, double corr_length_m, py::object theta_deg,
       std::complex<double> er_complex, const std::string& correl, bool include_hv, bool return_db){
        // Accept scalar or array for theta
        auto to_array = [](py::object obj)->py::array_t<double>{
          if (py::isinstance<py::float_>(obj) || py::isinstance<py::int_>(obj)) {
            auto arr = py::array_t<double>({1});
            *arr.mutable_data() = obj.cast<double>();
            return arr;
          }
          return py::cast<py::array_t<double>>(obj);
        };
        py::array_t<double> thetas = to_array(theta_deg);
        auto r = thetas.request();
        auto out_hh = py::array_t<double>(r.size);
        auto out_vv = py::array_t<double>(r.size);
        auto out_hv = py::array_t<double>(r.size);

        int sp = correl_to_int(correl);
        cdouble er(er_complex.real(), er_complex.imag());

        const double deg2rad = 3.14159265358979323846/180.0;
        for(ssize_t i=0;i<r.size;i++){
          double th = thetas.data()[i];
          validate_common(freq_ghz, rms_height_m, corr_length_m, th);

          double hh=0.0, vv=0.0, hv=0.0;
          // Use bistatic with backscatter geometry for HH/VV
          I2EM_Bistat_model((float)freq_ghz, (float)rms_height_m, (float)corr_length_m,
                            th*deg2rad, th*deg2rad, 180.0*deg2rad, er, sp, 1.0f, &hh, &vv);
          if(include_hv){
            IEMX_model((float)freq_ghz, (float)rms_height_m, (float)corr_length_m,
                       th*deg2rad, er, sp, 1.0f, /*auto_select*/1, &hv);
          }
          if(!return_db){
            // convert dB to linear if requested
            hh = std::pow(10.0, hh/10.0);
            vv = std::pow(10.0, vv/10.0);
            hv = std::pow(10.0, hv/10.0);
          }
          out_hh.mutable_data()[i]=hh;
          out_vv.mutable_data()[i]=vv;
          out_hv.mutable_data()[i]=hv;
        }
        return py::dict("hh"_a=out_hh, "vv"_a=out_vv, "hv"_a=out_hv);
    },
    py::arg("freq_ghz"), py::arg("rms_height_m"), py::arg("corr_length_m"),
    py::arg("theta_deg"), py::arg("er_complex"), py::arg("correl")="exponential",
    py::arg("include_hv")=true, py::arg("return_db")=true,
    "Monostatic sigma0 at backscatter. Accepts scalar or array for theta_deg.");

  m.def("sigma0_bistatic",
    [](double freq_ghz, double rms_height_m, double corr_length_m,
       double thi_deg, double ths_deg, double phs_deg,
       std::complex<double> er_complex, const std::string& correl,
       double xcoeff, bool return_db){
        validate_common(freq_ghz, rms_height_m, corr_length_m, thi_deg);
        if(ths_deg<0 || ths_deg>89.999) throw std::invalid_argument("ths_deg must be in [0,90)");
        int sp = correl_to_int(correl);

        double hh=0.0, vv=0.0;
        cdouble er(er_complex.real(), er_complex.imag());
        const double deg2rad = 3.14159265358979323846/180.0;
        I2EM_Bistat_model((float)freq_ghz, (float)rms_height_m, (float)corr_length_m,
                          thi_deg*deg2rad, ths_deg*deg2rad, phs_deg*deg2rad,
                          er, sp, (float)xcoeff, &hh, &vv);
        if(!return_db){
          hh = std::pow(10.0, hh/10.0);
          vv = std::pow(10.0, vv/10.0);
        }
        return py::dict("hh"_a=hh, "vv"_a=vv);
    },
    py::arg("freq_ghz"), py::arg("rms_height_m"), py::arg("corr_length_m"),
    py::arg("thi_deg"), py::arg("ths_deg"), py::arg("phs_deg"),
    py::arg("er_complex"), py::arg("correl")="exponential",
    py::arg("xcoeff")=1.0, py::arg("return_db")=true,
    "Bistatic sigma0 (HH/VV). Angles in degrees.");
}
