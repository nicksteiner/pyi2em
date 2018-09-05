#ifndef em_i2em_h
#define em_i2em_h

const double PI=3.1415927;
const double DEG2RAD=0.0174533;

#define CORREL_EXPONENTIAL 1
#define CORREL_GAUSSIAN    2
#define CORREL_X_EXPONENTIAL 4
#define CORREL_X_POWER 3

#include <complex>
#include <array>

using namespace std;

typedef complex<double> cdouble;
typedef std::array<double, 2> twoarr;

void parse_data(char *data, float *real_eps1, float *imag_eps1, 
                float *real_eps2, float *imag_eps2, 
                float *freq_ghz, float *rmsheight,  
                float *correl_length, float *albedo, float *extinction, 
                float *depth);

void calc_model(float real_eps1, float imag_eps1, 
                float real_eps2, float imag_eps2, 
                float freq_ghz, float rmsheight,  
                float correl_length, float albedo, 
                float extinction, float depth,
                float *x, float *y1, float *y2, float *y3);

void I2EM_Bistat_model(float freq_ghz, float rmsheight, float correl_length, 
                       double thi, double ths, double phs, cdouble er, 
                       int correl_func, float xcoeff,
                       double *sigma0_hh, double *sigma0_vv);

void IEMX_model(float freq_ghz, float rmsheight, float correl_length,  
                double thi, cdouble er, 
                int correl_func, float xcoeff, int auto_select, double *sigma0_hv);


void roughness_spectrum(int correl_func, float xx, double wvnb, double sig, 
                        double L, int Ts, double *wn, double *rss);

double x_exponential_spectrum(double z, void *params);

void do_rav_integration(double cs,double s,cdouble er,
                        double s2,double sigx, double sigy, double xxx,   
                        cdouble *Rav);

int Rav_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval);


void do_rah_integration(double cs,double s,cdouble er,
                        double s2,double sigx, double sigy, double xxx,   
                        cdouble *Rah);

int Rah_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval);

void Fppupdn_is_calculations(int ud, int is, cdouble Rvi, cdouble Rhi,
                             cdouble er,double k, double kz, double ksz,
                             double s, double cs, double ss, double css, double cf,
                             double cfs, double sfs, cdouble *vv, cdouble *hh);

void integrate_xpol(int sp, double xx, double ks2, double cs, double s, double kl2, 
               double L, cdouble er, double rss, cdouble rvh, 
               int n_spec, double *svh);

int xpol_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval);

void spectrm1(int sp, double xx, double kl2, double L, double rx, double ry, 
              double s, int np, double *wn);

void spectrm2(int sp, double xx, double kl2, double L, double rx, double ry, 
              double s, int np, double *wn);


void I2EM_Periodic(double theta_0, double phi_0, cdouble eps, 
                   double f, double s, double l, int sp, 
                   double Gmm, double A, 
                   double *sigma_0_vv, double *sigma_0_hh, double *sigma_0_vh);

void do_e_integration(double Zx,double theta_0, double phi_0,double f,
                        double s, double l, cdouble er,int sp,
                        double Gmm, double A,
                        double *ev, double *eh);

void Calc_scatt_coeff(double Zx,double Zy,double th,double ph,cdouble eps, 
                      double f, double s, double l, int sp,
                      double *sig_ss_vv, double *sig_ss_hh, double *sig_ss_vh);


double dot3( double *a, double *b);

double norm3( double *a);

void cross3( double *a, double *b, double *c);

double trapz(int n, double *x, double *y);

void  vectors2(double Zx,double Zy,double th, double ph,
               double *v, double *h,double *v_pr,double *h_pr,double *thpr);

int sig_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval);

void do_sig_integration(double Zx,double theta_0, double phi_0,double f,
                        double s, double l, cdouble eps,int sp,
                        double Gmm, double A,
                        double *sig_ss_vv, double *sig_ss_hh, double *sig_ss_vh);


void Emissivity_I2EM_Periodic(double theta_0, double phi_0, cdouble eps, 
                   double f, double s, double l, int sp, 
                   double Gmm, double A, 
                   double *ev, double *eh);


int e_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval);

void Calc_scatt_emiss_coeff(double Zx,double Zy,double th,double ph,cdouble eps, 
                            double f, double s, double l, int sp,
                            double *ev, double *eh);


void Calc_emissivity_I2EMmodel(double fr, double sig, double L, 
                               double theta_d, cdouble er, 
                               int sp, double *eh, double *ev);

double Emissivity_I2EMmodel_VPOL(double fr, double sig, double L, 
                             double theta_d, double el, double ei, int sp);

std::string Emissivity_I2EMmodel(double fr, double sig, double L, 
                             double theta_d, double el, double ei, int sp);

int emsv_integrand(unsigned ndim, const double *x, void *params, 
                   unsigned fdim, double *fval);

void roughness_spectrum_12_2(int correl_func, double kl2, double L, 
                        double wvnb, int np, double *wn);



void ZRTemission_veg(cdouble er, double s, double l,
                     double freq_ghz, double thi, double albedo, 
                     double extinction, double depth, 
                     double *ev, double *eh);

void ReflTransm_PlanarBoundary(cdouble eps1, cdouble eps2, 
                               double theta1d, 
                               cdouble *rhoh, cdouble *rhov,
                               double *gammah, double *gammav,
                               cdouble *tauh, cdouble *tauv,
                               double *Th, double *Tv);




void ZRTemission_DUB(cdouble eps2,cdouble eps3, 
                     double theta_i, double f, double s, 
                     double l, double a, double d, double kappa_e,
                     double *e_v, double *e_h);


int add(int i, int j);

#endif