#include <stdio.h>
#include <string.h>
#include <complex> //changed to C++;
#include <array>

#include "cubature.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include "em_i2em.h"

using namespace std;

//*********************************************************************************
void parse_data(char *data, float *real_eps1, float *imag_eps1,
                float *real_eps2, float *imag_eps2,
                float *freq_ghz,  float *rmsheight,  
                float *correl_length, 
                float *albedo, float *extinction, float *depth)
{
  /* expect: density=%20%2015.85&pressure=%20%201000&temp=%20%2027 */
  char *start;
  float r;
  int i;
  // remove the "%20"'s if there.
  i=0;
  while(i<strlen(data))
    {
      if( data[i]=='%' && data[i+1]=='2' && data[i+2]=='0') 
        {
          data[i]=' ';
          data[i+1]=' ';
          data[i+2]=' ';
        }//endif
      i++;
    }//endwhile
  
  //mod10_1.cgi?real_eps1=27&imag_eps1=0&freq=1&correl_func=1&rmsheight=0&correl_length=1&xcoeff=1

  start = strstr(data,"real_eps1=");
  sscanf(start+10,"%f",real_eps1);

  start = strstr(data,"imag_eps1=");
  sscanf(start+10,"%f",imag_eps1);


  start = strstr(data,"real_eps2=");
  sscanf(start+10,"%f",real_eps2);

  start = strstr(data,"imag_eps2=");
  sscanf(start+10,"%f",imag_eps2);


  start = strstr(data,"freq=");
  sscanf(start+5,"%f",freq_ghz);

  start = strstr(data,"rmsheight=");
  sscanf(start+10,"%f",rmsheight);

  start = strstr(data,"correl_length=");
  sscanf(start+14,"%f",correl_length);

  start = strstr(data,"albedo=");
  sscanf(start+7,"%f",albedo);

  start = strstr(data,"extinction=");
  sscanf(start+11,"%f",extinction);

  start = strstr(data,"depth=");
  sscanf(start+6,"%f",depth);


}

//*********************************************************************************
void calc_model(float real_eps1, float imag_eps1, 
                float real_eps2, float imag_eps2, 
                float freq, float rmsheight,  
                float correl_length, 
                float falbedo, float fextinction, float fdepth,
                float *x, float *y1, float *y2, float *y3)
{
  int i;

  double thi, ths, phs;

  double ev, eh;

  double freq_ghz, s,l;

  double albedo, extinction ,depth;

  // ask adib: -imag or +imag part?
  //er = real_eps1 + imag_eps1*I; change to C++
  //er2 = real_eps2 + imag_eps2*I;
  
  cdouble er(real_eps1, imag_eps1);
  cdouble er2(real_eps2, imag_eps2);

  freq_ghz = freq;
  s=rmsheight;
  l=correl_length;

  albedo = falbedo;
  extinction = fextinction;
  depth= fdepth;

  for(i=0;i<45;i++)
    {
      thi = 2.*i+0.001;

      ZRTemission_DUB(er,er2,thi,freq_ghz, s, l, albedo, depth,
                      extinction, &ev, &eh);

      x[i] = i*2;
      y1[i] = eh;
      y2[i] = ev;

    }
  
}


//*********************************************************************************
//*********************************************************************************
//*********************************************************************************
void I2EM_Bistat_model(float freq_ghz, float rmsheight, float correl_length, 
                       double thi, double ths, double phs, cdouble er, 
                       int correl_func, float xcoeff,
                       double *sigma0_hh, double *sigma0_vv)
{
  double error;

  double sig, L;
  double mu_r;

  double k, ks, kl, ks2;
  double kx, ky, kz;
  double ksx, ksy, ksz;
  double theta, phi, thetas, phis;

  double cs, s, s2;
  double sf, cf;
  double ss, css;
  double cfs, sfs;

  cdouble rt, Rvi, Rhi;
  cdouble Rvt, Rht;
  double wvnb;
  int Ts;
  double scale;

  double wn[1000];
  double rss;

  cdouble Rv0, Rh0;
  cdouble Ft;
  double a0, a1;
  cdouble b1;
  cdouble tmp1;
  cdouble St, St0, Tf;

  double sigx, sigy, xxx;

  cdouble Rav, Rah;
  cdouble fvv, fhh;
  cdouble Fvvupi, Fhhupi;
  cdouble Fvvups, Fhhups;
  cdouble Fvvdni, Fhhdni;
  cdouble Fvvdns, Fhhdns;


  double qi, qs;
  cdouble Ivv[1000], Ihh[1000];

  double ct, cts, rslp, ctorslp, ctsorslp, shadf, shadfs, ShdwS;

  double ssv, ssh;

  int i, n;

  double temp1;
  double sigmavv, sigmahh;

  sig = rmsheight * 100.; // change from m to cm scale
  L = correl_length * 100.;
  mu_r = .01; // relative permeability



  k = 2*PI *freq_ghz/30.; // wavenumber in free space. Speed of light is in cm/sec
  theta = thi; // radians
  phi = 0.;
  thetas = ths; // radians
  phis = phs;

  ks = k * sig; // roughness parameter
  kl = k * L;  

  ks2 = ks * ks; 

  cs = cos(theta+ 0.01);
  s =  sin(theta+ 0.01);

  sf = sin(phi);
  cf = cos(phi);

  ss =  sin(thetas);
  css = cos(thetas);

  cfs = cos(phis);
  sfs = sin(phis);

  s2 = s * s;

  kx = k * s *cf;
  ky = k * s *sf;
  kz = k * cs;

  ksx = k * ss *cfs;
  ksy = k * ss *sfs;
  ksz = k * css;

  //-- reflection coefficients
  rt = sqrt(er - s2); // changed from sqrt
  Rvi = (er * cs - rt) / (er * cs + rt);
  Rhi = (cs - rt) / (cs + rt);


  wvnb = k * sqrt((ss * cfs - s * cf) * (ss * cfs - s * cf) + (ss * sfs - s * sf) * (ss * sfs - s * sf));

  Ts = 1;

  error=1.0e8;
  scale=1.;
  while (error > 1.0e-8) {
    Ts++;
    scale /= Ts;
    error = pow(ks2*(cs + css)*(cs + css) ,Ts ) * scale;   // was: / factorial(Ts); 
  } // end while



  //---------------- calculating roughness spectrum -----------
  // returns: wn[] and rss
  roughness_spectrum(correl_func,xcoeff, wvnb, sig, L, Ts, wn, &rss);


  //----------- compute R- transition ------------

  Rv0 = (sqrt(er)-1.) /(sqrt(er)+1.); // changed from sqrt
  Rh0 = -Rv0;

  Ft = 8.* Rv0*Rv0 *ss*(cs + sqrt(er - s2)) / (cs * sqrt(er - s2));  // changed from sqrt
  a1 = 0.0;
  b1 = 0.0;
  scale =1.;
  for(i=0; i<Ts; i++) {
      n= i + 1;
      scale /= n;
      a0 = pow(ks * cs, 2. * n) *scale;  /// was:./factorial(n);
      a1 += a0 * wn[i];
      tmp1 = abs(Ft * 0.5 + pow(2. ,n+1.) * Rv0 / cs *exp(-(ks * ks * cs * cs)));
      b1 += a0 * tmp1 * tmp1 * wn[i];
  }// endfor
  
  St = 0.25 * (abs(Ft)) * (abs(Ft)) * a1 / b1;  // changed from abs
  tmp1 = abs(1. + 8. * Rv0 / (cs * Ft));  // changed from abs
  St0 = 1./ (tmp1 * tmp1);
  Tf = 1. - St / St0;


  //----------- compute average reflection coefficients ------------
  //-- these coefficients account for slope effects, especially near the
  //brewster angle. They are not important if the slope is small.

  sigx = 1.1*sig/L;
  sigy = sigx;
  xxx = 3.*sigx;

  // integrations:
  do_rav_integration(cs,s,er,s2,sigx, sigy, xxx, &Rav);
  //Rav = dblquad(@(Zx, Zy)Rav_integration(Zx, Zy, cs,s,er,s2,sigx, sigy),-xxx,xxx, -xxx, xxx );

  do_rah_integration(cs,s,er,s2,sigx, sigy, xxx, &Rah);
  //Rah = dblquad(@(Zx, Zy)Rah_integration(Zx, Zy, cs,s,er,s2,sigx, sigy),-xxx,xxx, -xxx, xxx );

  Rav = Rav /(2.*PI * sigx * sigy);
  Rah = Rah /(2.*PI * sigx * sigy);



  //-- select proper reflection coefficients

  if (fabs(thi - ths) < 1.e-3 && fabs(phs -180.*DEG2RAD) < 1.e-3) { // %i.e. operating in backscatter mode
    Rvt = Rvi + (Rv0 - Rvi)*Tf;
    Rht = Rhi + (Rh0 - Rhi)*Tf;
  } else  {
    // in this case, it is the bistatic configuration and average R is used
    //     Rvt = Rav + (Rv0 - Rav) .* Tf;
    //     Rht = Rah + (Rh0 - Rah) .* Tf;
    Rvt = Rav;
    Rht = Rah;
  }// endif
  

  fvv =  2.* Rvt *(s * ss - (1. + cs * css) * cfs)/(cs + css);
  fhh = -2.* Rht *(s * ss - (1. + cs * css) * cfs)/(cs + css);

  //------- Calculate the Fppup(dn) i(s) coefficients ----
  Fppupdn_is_calculations(+1, 1, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs, &Fvvupi, &Fhhupi);
  Fppupdn_is_calculations(+1, 2, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs, &Fvvups, &Fhhups);
  Fppupdn_is_calculations(-1, 1, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs, &Fvvdni, &Fhhdni);
  Fppupdn_is_calculations(-1, 2, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs, &Fvvdns, &Fhhdns);


  qi = k * cs;
  qs = k * css;

  //----- calculating  Ivv and Ihh ----

  for(i=0;i<Ts;i++){
    n = i+1;
    Ivv[i] = pow(kz + ksz,n) * fvv * exp(-sig*sig * kz * ksz) + 
      0.25*(Fvvupi *pow(ksz-qi,n-1.) *exp(-sig*sig *(qi*qi - qi*(ksz-kz)))+ 
            Fvvdni *pow(ksz+qi,n-1.) *exp(-sig*sig *(qi*qi + qi*(ksz-kz)))+ 
            Fvvups *pow(kz+qs ,n-1.) *exp(-sig*sig *(qs*qs - qs*(ksz-kz)))+ 
            Fvvdns *pow(kz-qs ,n-1.) *exp(-sig*sig *(qs*qs + qs*(ksz-kz))));

    
    Ihh[i] = pow(kz + ksz,n) * fhh * exp(-sig*sig * kz * ksz) +
      0.25*(Fhhupi *pow(ksz-qi,n-1.) *exp(-sig*sig *(qi*qi - qi*(ksz-kz)))+ 
            Fhhdni *pow(ksz+qi,n-1.) *exp(-sig*sig *(qi*qi + qi*(ksz-kz)))+ 
            Fhhups *pow(kz+qs ,n-1.) *exp(-sig*sig *(qs*qs - qs*(ksz-kz)))+ 
            Fhhdns *pow(kz-qs ,n-1.) *exp(-sig*sig *(qs*qs + qs*(ksz-kz))));
  }//endfor
  


  //-- Shadowing function calculations
  if (fabs(thi - ths) < 1.e-3 && fabs(phs -180.*DEG2RAD) < 1.e-3) { // %i.e. operating in backscatter mode
    ct = 1./tan(theta);
    cts = 1./tan(thetas);
    rslp = rss;
    ctorslp = ct / sqrt(2) /rslp;
    ctsorslp = cts / sqrt(2) /rslp;
    shadf = 0.5 *(exp(-ctorslp*ctorslp) / sqrt(PI)/ctorslp - gsl_sf_erfc(ctorslp));
    shadfs = 0.5 *(exp(-ctsorslp*ctsorslp) / sqrt(PI)/ctsorslp - gsl_sf_erfc(ctsorslp));
    ShdwS = 1./(1. + shadf + shadfs); 
  } else {
    ShdwS = 1.;
  }// endif
  

  //------- calculate the values of sigma_nought --------------

  sigmavv = 0.;
  sigmahh = 0.; 

  scale =1.;
  for(i=0;i<Ts;i++) {
      n=i+1;
      scale = scale/n;  // instead of 1/factorial....
      a0 = wn[i] *scale *pow(sig,2.*n);
    
      sigmavv += abs(Ivv[i])*abs(Ivv[i]) *a0; // changed from abs
      sigmahh += abs(Ihh[i])*abs(Ihh[i]) *a0;

  }// endfor
  
  temp1 = ShdwS * 0.5*k*k * exp(-sig*sig *(kz*kz +ksz*ksz));
  sigmavv *= temp1;
  sigmahh *= temp1;


  ssv = 10. * log10(fabs(sigmavv)); 
  ssh = 10. * log10(fabs(sigmahh)); 

  // return values in dB:
  *sigma0_vv = ssv;
  *sigma0_hh = ssh;

}


//*********************************************************************************************
void Fppupdn_is_calculations(int ud, int is, cdouble Rvi, cdouble Rhi,
                             cdouble er,double k, double kz, double ksz,
                             double s, double cs, double ss, double css, double cf,
                             double cfs, double sfs, cdouble *vv, cdouble *hh)
{
  double Gqi;
  cdouble Gqti;
  double qi;

  double Gqs;
  cdouble Gqts;
  double qs;

  cdouble c11, c21, c31, c41, c51;
  cdouble c12, c22, c32, c42, c52;

  cdouble temp1, temp2;
  cdouble q, invq;
  cdouble qt, invqt;

  switch(is) {
  case 1:
    Gqi = ud * kz;
    Gqti = ud *k * sqrt(er-s*s); // changed from sqrt
    qi = ud * kz;
    
    c11 = k * cfs *(ksz - qi);

    temp1 = k*k *s*cf*(ss *cfs - s * cf);
    temp2 = k*k *cf * s *ss *sfs*sfs;
    //c21 = cs *(cfs *(k*k *s*cf*(ss *cfs - s * cf) + Gqi*(k * css - qi)) 
    //           + k*k *cf * s *ss *sfs*sfs);
    c21 = cs *(cfs *(temp1 + Gqi*(k * css - qi)) + temp2);


    //c22 = cs *(cfs *(k*k *s*cf*(ss *cfs - s * cf) + Gqti*(k * css - qi))
    //    + k*k *cf * s *ss *sfs*sfs);
    c22 = cs *(cfs *(temp1 + Gqti*(k * css - qi)) + temp2);



    temp1 = s*cf*cfs*(k*css-qi);
    temp2 = (cfs*(ss*cfs -s*cf)+ ss *sfs*sfs);
    //c31 = k*s*(s*cf*cfs*(k*css-qi) - Gqi*(cfs*(ss*cfs -s*cf)+ ss *sfs*sfs));
    c31 = k*s*(temp1 - Gqi*temp2);

    //c32 = k*s*(s*cf*cfs*(k*css-qi) - Gqti*(cfs*(ss*cfs -s*cf)- ss *sfs*sfs));
    c32 = k*s*(temp1 - Gqti*temp2);


    c41 = k *cs*(cfs*css*(k*css - qi) + k *ss*(ss*cfs-s*cf));

    temp1 = (cfs *css*(qi-k*css) - k *ss*(ss*cfs-s*cf));
    //c51 = Gqi*(cfs *css*(qi-k*css) - k *ss*(ss*cfs-s*cf));
    c51 = Gqi*temp1;
  
    //c52 = Gqti*(cfs *css*(qi-k*css) - k *ss*(ss*cfs-s*cf));    
    c52 = Gqti*temp1;

  
    //c12 = k * cfs *(ksz - qi);
    c12 = c11;
    //c42 = k *cs*(cfs*css*(k*css - qi) + k *ss*(ss*cfs-s*cf));
    c42 = c41;
    break;

    case 2:
    Gqs = ud * ksz;
    Gqts = ud *k * sqrt(er-ss*ss); // changed from sqrt
    qs = ud * ksz;
    
    c11 = k * cfs *(kz + qs);
    temp1 = k*s*(ss *cfs-s*cf);
    temp2 = k*s*ss*sfs*sfs;
    //c21 = Gqs  *(cfs*(cs*(k*cs+qs)-k*s*(ss *cfs-s*cf))-k*s*ss*sfs*sfs);
    //c22 = Gqts *(cfs*(cs*(kz+qs)  -k*s*(ss *cfs-s*cf))-k*s*ss*sfs*sfs);
    c21 = Gqs  *(cfs*(cs*(k*cs+qs)-temp1)-temp2);
    c22 = Gqts *(cfs*(cs*(kz+qs)  -temp1)-temp2);

    c31 = k *ss*(k*cs*(ss*cfs - s*cf)+ s*(kz+qs));

    c41 = k*css*(cfs*(cs*(kz+qs)-k*s*(ss*cfs-s*cf))-k*s*ss*sfs*sfs);

    temp1 = k*k *ss *(ss*cfs -s*cf);
    temp2=cfs*(kz+qs);
    //c51 = -css *(k*k *ss *(ss*cfs -s*cf)+ Gqs *cfs*(kz+qs));
    //c52 = -css *(k*k *ss *(ss*cfs -s*cf)+ Gqts*cfs*(kz+qs));
    c51 = -css *(temp1+ Gqs *temp2);
    c52 = -css *(temp1+ Gqts*temp2);
         
    //c12 = k * cfs *(kz + qs);
    c12 = c11;
    //c32 = k *ss*(k*cs*(ss*cfs - s*cf)+ s*(kz+qs));
    c32 = c31;
    //c42 = k*css*(cfs*(cs*(kz+qs)-k*s*(ss*cfs-s*cf))-k*s*ss*sfs*sfs);
    c42 = c41;
    break;
    }// end-switch
  

  q=kz;
  qt = k * sqrt(er - s*s); // changed from sqrt

  invq = 1./q;
  invqt= 1./qt;

  *vv = (1. + Rvi) *(-(1. - Rvi) * c11 * invq +    (1.+ Rvi) * c12      * invqt) + 
        (1. - Rvi) * ((1. - Rvi) * c21 * invq -    (1.+ Rvi) * c22      * invqt) + 
        (1. + Rvi) * ((1. - Rvi) * c31 * invq -    (1.+ Rvi) * c32 / er * invqt) + 
        (1. - Rvi) * ((1. + Rvi) * c41 * invq - er*(1.- Rvi) * c42      * invqt) + 
        (1. + Rvi) * ((1. + Rvi) * c51 * invq -    (1.- Rvi) * c52      * invqt);
  
  *hh = (1. + Rhi) *( (1.-Rhi) *c11 *invq - er*(1.+Rhi) *c12 *invqt) - 
    (1. - Rhi) *( (1.-Rhi) *c21 *invq - (1.+Rhi)    *c22 *invqt) - 
    (1. + Rhi) *( (1.-Rhi) *c31 *invq - (1.+Rhi)    *c32 *invqt) - 
    (1. - Rhi) *( (1.+Rhi) *c41 *invq - (1. - Rhi)  *c42 *invqt) - 
    (1. + Rhi) *( (1.+Rhi) *c51 *invq - (1.-Rhi)    *c52 *invqt);
  
}


//*********************************************************************************************
void do_rav_integration(double cs,double s,cdouble er,
                        double s2,double sigx, double sigy, double xxx,   
                        cdouble *Rav)
{

  double params[7];
  double xmin[2], xmax[2];
  double val[2], err[2];

  params[0] = cs;
  params[1] = s;
  params[2] = real(er);
  params[3] = imag(er);
  params[4] = s2;
  params[5] = sigx;
  params[6] = sigy;

  // integration ranges:
  xmin[0]=-xxx;
  xmin[1]=-xxx;
  xmax[0]=xxx;
  xmax[1]=xxx;

  pcubature(2, Rav_integrand, params,
            2, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);

  //*Rav = val[0] + val[1]*I;
  cdouble T (val[0], val[1]);
  *Rav = T;

}


//*********************************************************************************************
int Rav_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval)
{
  double A;
  cdouble B;
  double CC;

  cdouble Rv;
  double pd;

  double cs,s, s2, sigx, sigy;
  //cdouble er;
  double Zx, Zy;

  double Zx_squared, Zy_squared;
  cdouble temp1;

  cs = ((double *)params)[0];
  s  = ((double *)params)[1];
  //er = ((double *)params)[2] + (((double *)params)[3])*I;
  cdouble er (((double *)params)[2], (((double *)params)[3]));
  s2 = ((double *)params)[4];
  sigx = ((double *)params)[5];
  sigy = ((double *)params)[6];

  Zx = x[0];
  Zy = x[1];

  Zx_squared = Zx*Zx;
  Zy_squared = Zy*Zy;

  A = cs + Zx * s;
  B = er * (1 + Zx_squared + Zy_squared);
  CC = s2 - 2.*Zx*s*cs + Zx_squared * cs*cs + Zy_squared;

  temp1 = sqrt(B-CC);
  Rv = (er*A - temp1)/(er*A + temp1); 

  pd = exp(-Zx_squared/(2.*sigx*sigx) -Zy_squared/(2*sigy*sigy));

  temp1 = Rv * pd;
  fval[0] = real(temp1);
  fval[1] = imag(temp1);

  return 0; //success
}


//*********************************************************************************************
void do_rah_integration(double cs,double s,cdouble er,
                        double s2,double sigx, double sigy, double xxx,   
                        cdouble *Rah)
{

  double params[7];
  double xmin[2], xmax[2];
  double val[2], err[2];

  params[0] = cs;
  params[1] = s;
  params[2] = real(er);
  params[3] = imag(er);
  params[4] = s2;
  params[5] = sigx;
  params[6] = sigy;

  // integration ranges:
  xmin[0]=-xxx;
  xmin[1]=-xxx;
  xmax[0]=xxx;
  xmax[1]=xxx;

  pcubature(2, Rah_integrand, params,
            2, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);

  cdouble T (val[0], val[1]);
  //*Rah = val[0] + val[1]*I;
  *Rah = T;

}


//*********************************************************************************************
int Rah_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval)
{
  double A;
  cdouble B;
  double CC;

  cdouble Rh;
  double pd;

  double cs,s, s2, sigx, sigy;
  //cdouble er;
  double Zx, Zy;

  double Zx_squared, Zy_squared;
  cdouble temp1;

  cs = ((double *)params)[0];
  s  = ((double *)params)[1];
  cdouble er (((double *)params)[2],  (((double *)params)[3]));
  //er = ((double *)params)[2] + (((double *)params)[3])*I;
  s2 = ((double *)params)[4];
  sigx = ((double *)params)[5];
  sigy = ((double *)params)[6];

  Zx = x[0];
  Zy = x[1];

  Zx_squared = Zx*Zx;
  Zy_squared = Zy*Zy;

  A = cs + Zx * s;
  B = er * (1. + Zx_squared + Zy_squared);
  CC = s2 - 2.*Zx*s*cs + Zx_squared * cs*cs + Zy_squared;

  temp1 = sqrt(B-CC);sqrt(B-CC);
  Rh = (A - temp1)/(A + temp1); 

  pd = exp(-Zx_squared/(2.*sigx*sigx) -Zy_squared/(2*sigy*sigy));

  temp1 = Rh * pd;
  fval[0] = real(temp1);
  fval[1] = imag(temp1);

  return 0; //success
}



//*********************************************************************************************
void roughness_spectrum(int correl_func, float xx, double wvnb, double sig, 
                        double L, int Ts, double *wn, double *rss)
{
  int i, n;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;

  double the_params[4];

  the_params[0]=wvnb;
  the_params[1]=L;
  the_params[3]=xx;


  switch(correl_func)
    {
    case CORREL_EXPONENTIAL:
      for (i=0;i<Ts;i++){
        n=i+1;
        wn[i] = L*L/ (n*n) * pow(1.+(wvnb*L/n)*(wvnb*L/n), -1.5);
      }// endfor
      *rss= sig/L;
      break;


    case CORREL_GAUSSIAN:
      for (i=0;i<Ts;i++){
        n=i+1;
        wn[i] =  L*L/(2*n) * exp( -(wvnb*L)*(wvnb*L)/(4*n) );
      }// endfor
      *rss = sqrt(2.) * sig/L;
      break;


    case CORREL_X_POWER:
      for (i=0;i<Ts;i++){
        n=i+1;
        if (fabs(wvnb) < 1.0e-50) {
          wn[i] = L*L /(3 * n - 2.);
        } else {
          wn[i] = L*L * pow(wvnb*L, -1 + xx*n) * gsl_sf_bessel_Knu(1.-xx*n, wvnb*L)/
            (pow(2.,xx*n -1) * gsl_sf_gamma(xx*n)); 
        }// endif
      }// endfor
      *rss  = sqrt(xx *2.)*sig/L;
      break;

    case CORREL_X_EXPONENTIAL:
      F.function = &x_exponential_spectrum;
      F.params = the_params;

      for (i=0;i<Ts;i++){
        n=i+1;
        the_params[2]= n;

        gsl_integration_qag (&F, 0.0, 9.0,   0, 1e-3, 1000, 1, 
                             w, &result, &error);

        wn[i] = L*L/pow(n,2./xx)*result;
      }
      *rss  =  sig/ L;
      break;
    }// end-switch
}


//*********************************************************************************************
double x_exponential_spectrum(double z, void *params)
{
  double wvnb, L, n, xx;

  wvnb = ((double *)params)[0];
  L    = ((double *)params)[1];
  n    = ((double *)params)[2];
  xx   = ((double *)params)[3];
  
  
  return exp(pow(-fabs(z),xx)) * gsl_sf_bessel_j0(z*wvnb*L/pow(n,1./xx) ) *z;
}


//*********************************************************************************
//*********************************************************************************
//*********************************************************************************
//*********************************************************************************
void IEMX_model(float freq_ghz, float rmsheight, float correl_length,  
                double thi, cdouble er, 
                int correl_func, float xcoeff, int auto_select, double *sigma0_vh)
{

  double fr, sig, L, theta;
  int sp;
  double xx;
  double error;
  double k, ks, kl, ks2, kl2;
  double cs, s, s2;

  cdouble rt, rv, rh, rvh;
  double sig_l;
  double rss;
  double scale;
  double ct, farg, gamma, Shdw; 
  double svh;
  int n_spec;
  
  fr = freq_ghz;
  sig = rmsheight*100.; //change to cm scale
  L   = correl_length*100.;//change to cm scale
  theta = thi; // in radians
  sp = correl_func;
  xx = xcoeff;

  //- fr: frequency in GHz
  //- sig: rms height of surface in cm
  //- L: correlation length of surface in cm
  //- theta_d: incidence angle in degrees
  //- er: relative permittivity
  //- sp: type of surface correlation function

  error = 1.0e8;

  k = 2*PI *fr/30.; // wavenumber in free space. Speed of light is in cm/sec

  ks = k * sig; // roughness parameter
  kl = k * L;  

  ks2 = ks * ks; 
  kl2 = kl *kl;

  cs = cos(theta);
  s = sin(theta+ 0.001);
  s2 = s*s;

  //-- calculation of reflection coefficients
  rt = sqrt(er - s2);

  rv = (er *cs - rt) /(er*cs +rt);
  rh = (cs - rt)/(cs + rt);

  rvh = 0.5*(rv - rh);

  //-- rms slope values
  sig_l = sig/L;
  switch(sp) {
  case CORREL_X_EXPONENTIAL:
    // not allowed....
    // silently use exponential
  case CORREL_EXPONENTIAL:
    rss = sig_l;
    break;
  case CORREL_GAUSSIAN:
    rss = sig_l * sqrt(2.);
    break;
  case CORREL_X_POWER:    //-- 1.5-power spectra correl func  (??)
    rss = sig_l * sqrt(2.*xx);
    break;
  }// end switch


  //--- Selecting number of spectral components of the surface roughness
  if (auto_select == 0 ){
    n_spec = 15; // number of terms to include in the surface roughness spectra
  } else {
   n_spec = 1; 
   scale =1.;
   while (error > 1.0e-8) {
    n_spec++;
    scale = scale/n_spec;
    error = pow(ks2 *4.*cs*cs , n_spec) * scale;  // was: ./ factorial(n_spec); 
   }// endwhile
  }// endif


  //-- calculating shadow consideration in single scat (Smith, 1967)

  ct = 1./tan(theta+ 0.001);
  farg = ct /sqrt(2.)/rss;
  gamma = 0.5 *(exp(-farg*farg) / 1.772 / farg - gsl_sf_erfc(farg));
  Shdw = 1. / (1. + gamma);

  //-- calculating multiple scattering contribution
  //------ a double integration function
  integrate_xpol(sp,xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec, &svh);
  //svh = dblquad(@(r,phi)xpol_integralfunc(r, phi, sp,xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec), 0.1, 1, 0, pi);

  //printf("svh=%f\n",svh);

  *sigma0_vh = 10.*log10(fabs(svh* Shdw));

}

//*********************************************************************************
void integrate_xpol(int sp, double xx, double ks2, double cs, double s, double kl2, 
               double L, cdouble er, double rss, cdouble rvh, 
               int n_spec, double *svh)
{
  double params[13];
  double xmin[2], xmax[2];
  double val[1], err[1];

  params[0] = sp;
  params[1] = xx;
  params[2] = ks2;
  params[3] = cs;
  params[4] = s;
  params[5] = kl2;
  params[6] = L;
  params[7] = real(er);
  params[8] = imag(er);
  params[9] = rss;
  params[10] = real(rvh);
  params[11] = imag(rvh);
  params[12] = n_spec+0.1;

  // integration ranges:
  xmin[0]=0.1;
  xmax[0]=1.0;

  xmin[1]=0.0;
  xmax[1]=PI;

  pcubature(1, xpol_integrand, params,
            2, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);

  *svh = val[0];
}

//*********************************************************************************************
int xpol_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval)
{

  double sp, xx, ks2, cs, s, kl2, L, rss;
  //cdouble er, rvh;
  int n_spec;

  double r, phi;


  double cs2, r2;
  double sf, csf, rx, ry;
  cdouble rp, rm;
  double q;
  cdouble qt;
  cdouble a,b,c,d;
  double B3;
  cdouble fvh1, fvh2, Fvh;
  double au, fsh, sha;
  double wn[1000], wm[1000];
  double acc;
  double vhmnsum;

  int i,j;
  int m, n;
  double scalen, scalem;
  double VH;





  sp = ((double *)params)[0];
  xx  = ((double *)params)[1];
  ks2  = ((double *)params)[2];
  cs  = ((double *)params)[3];
  s  = ((double *)params)[4];
  kl2  = ((double *)params)[5];
  L  = ((double *)params)[6];
  //er = ((double *)params)[7] + (((double *)params)[8])*I;
  cdouble er (((double *)params)[2],  (((double *)params)[3]));
  rss  = ((double *)params)[9];
  //rvh = ((double *)params)[10] + (((double *)params)[11])*I;
  cdouble rvh (((double *)params)[10],  (((double *)params)[11]));
  n_spec = (int) (((double *)params)[12]);



  r   = x[0];
  phi = x[1];



  cs2 = cs*cs;
  r2 = r*r;

  //nr = length(r);   nr=1 for my case....

  sf = sin(phi);
  csf = cos(phi);
  rx = r * csf;
  ry = r * sf;

  //printf("r=%lf, phi=%lf, ks2=%lf, cs=%lf, csf=%lf\n",r,phi,ks2,cs,csf);


  //-- calculation of the field coefficients
  rp = 1. + rvh;
  rm = 1. - rvh; 

  q = sqrt(1.0001 - r2);
  qt = sqrt(er - r2);

  a = rp /q;
  b = rm /q;
  c = rp /qt;
  d = rm /qt;

  //--calculate cross-pol coefficient
  B3 = rx * ry /cs;
  fvh1 = (b-c)*(1.- 3.*rvh) - (b - c/er) * rp; 
  fvh2 = (a-d)*(1.+ 3.*rvh) - (a - d*er) * rm;
  Fvh = abs( (fvh1 + fvh2) *B3)*abs( (fvh1 + fvh2) *B3) ;

  //printf("rx=%lf, ry=%lf, cs=%lf\n",rx,ry,cs);
  //printf("B3=%f, fvh1=%f+%fI, fvh2=%f+%fI, Fvh=%f\n",B3, real(fvh1),imag(fvh1),
  //       real(fvh2),imag(fvh2),Fvh);


  //-- calculate shadowing func for multiple scattering 
  au = q /r/1.414/rss;
  fsh = (0.2821/au) *exp(-au*au) -0.5*(1.- gsl_sf_erf(au));
  sha = 1./(1 + fsh); 

  //-- calculate expressions for the surface spectra
  spectrm1(sp, xx, kl2, L, rx, ry, s, n_spec, wn);
  spectrm2(sp, xx, kl2, L, rx, ry, s, n_spec, wm);


  //--compute VH scattering coefficient
  acc = exp(-2.* ks2 *cs2) /(16. * PI);
  vhmnsum = 0.0;
  scalen=1.;
  for(i=0;i<n_spec;i++) {
    n=i+1;
    scalen /= (1.*n);
    scalem=1.;
    for(j=0;j<n_spec;j++) {
      m=j+1;
      scalem /= (1.*m);
      vhmnsum += wn[i]*wm[j]*pow(ks2*cs2,1.*n+m)*scalen*scalem; // was /factorials()
      //printf("vhmnsum=%lf\n",vhmnsum);
    }// endfor: j
  }// endfor: i

  VH = real(4. * acc * Fvh * vhmnsum *r *sha);

  fval[0]= VH;

  //printf("VH=%lf, acc=%lf\n",VH, acc);
  //printf("Fvh=%lf, vhmnsum=%lf\n",Fvh, vhmnsum);
  //printf("r= %lf, sha=%lf, n_spec=%d\n",r, sha, n_spec);

  return 0;

}


//*********************************************************************************
void spectrm1(int sp, double xx, double kl2, double L, double rx, double ry, 
              double s, int np, double *wn)
{
  int i,n;
  switch(sp)  {
    case CORREL_X_EXPONENTIAL: // not allowed, silently use exp
    case CORREL_EXPONENTIAL:
      for(i=0;i<np;i++) {
          n=i+1;
          wn[i] = n* kl2 /pow(n*n + kl2 *((rx-s)*(rx-s)+ry*ry),1.5);
      }//endfor
      break;
    case CORREL_GAUSSIAN:
      for(i=0;i<np;i++) {
          n=i+1;
          wn[i] = 0.5 * kl2 /n * exp(-kl2*((rx-s)*(rx-s) + ry*ry)/(4.*n)) ;
      }//endfor
      break;
    case CORREL_X_POWER:    //-- 1.5-power spectra correl func  (??)
      for(i=0;i<np;i++) {
          n=i+1;
          wn[i] = kl2 /( pow(2.,xx*n-1.) *gsl_sf_gamma(xx*n))* 
                   pow( ( (rx-s)*(rx-s) + ry*ry)*L,xx*n-1.) * 
                   gsl_sf_bessel_Knu(-xx*n+1, L*((rx-s)*(rx-s) + ry*ry));
      }//endfor
      break;
  }// endswitch

}


//*********************************************************************************
void spectrm2(int sp, double xx, double kl2, double L, double rx, double ry, 
              double s, int np, double *wn)
{
  int i,n;
  switch(sp)  {
    case CORREL_X_EXPONENTIAL: // not allowed, silently use exp
    case CORREL_EXPONENTIAL:
      for(i=0;i<np;i++) {
          n=i+1;
          wn[i] = n* kl2 /pow(n*n + kl2 *((rx+s)*(rx+s)+ry*ry),1.5);
      }//endfor
      break;
    case CORREL_GAUSSIAN:
      for(i=0;i<np;i++) {
          n=i+1;
          wn[i] = 0.5 * kl2 /n * exp(-kl2*((rx+s)*(rx+s) + ry*ry)/(4.*n));
      }//endfor
      break;
    case CORREL_X_POWER:    //-- 1.5-power spectra correl func  (??)
      for(i=0;i<np;i++) {
          n=i+1;
          wn[i] = kl2 /( pow(2.,xx*n-1.) *gsl_sf_gamma(xx*n))* 
                   pow( ( (rx+s)*(rx+s) + ry*ry)*L,xx*n-1.) * 
                   gsl_sf_bessel_Knu(-xx*n+1, L*((rx+s)*(rx+s) + ry*ry));
      }//endfor
      break;
  }// endswitch

}


//**********************************************************************************
void I2EM_Periodic(double theta_0, double phi_0, cdouble eps, 
                   double f, double s, double l, int sp, 
                   double Gmm, double A, 
                   double *sigma_0_vv, double *sigma_0_hh, double *sigma_0_vh)
{
  double Zx;
  double sig_ss_vv, sig_ss_hh, sig_ss_vh;
  // both angles in degrees.

  Zx = 0.0; //partial derivative of periodic surface along x

  do_sig_integration(Zx,theta_0,phi_0,f,s,l,eps,sp,Gmm,A,
                     &sig_ss_vv,&sig_ss_hh,&sig_ss_vh);


  *sigma_0_vv = 10.*log10(fabs(sig_ss_vv / Gmm));
  *sigma_0_hh = 10.*log10(fabs(sig_ss_hh / Gmm));
  *sigma_0_vh = 10.*log10(fabs(sig_ss_vh / Gmm));

}


//*********************************************************************************************
void do_sig_integration(double Zx,double theta_0, double phi_0,double f,
                        double s, double l, cdouble er,int sp,
                        double Gmm, double A,
                        double *sig_ss_vv, double *sig_ss_hh, double *sig_ss_vh)
{

  double params[11];
  double xmin[1], xmax[1];
  double val[3], err[3];

  params[0] = Zx;
  params[1] = theta_0;
  params[2] = phi_0;
  params[3] = f;
  params[4] = s;
  params[5] = l;
  params[6] = real(er);
  params[7] = imag(er);
  params[8] = sp+0.1;
  params[9] = Gmm;
  params[10] = A;

  // integration ranges:
  xmin[0]= 0.0;
  xmax[0]= Gmm;

  pcubature(3, sig_integrand, params,
            1, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);

  *sig_ss_vv = val[0];
  *sig_ss_hh = val[1];
  *sig_ss_vh = val[2];
  
}

//*****************************************************************************************
int sig_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval)
{
  double Zx, Zy, theta_0, phi_0, f, s, l, epsr, epsi, dsp;
  //cdouble eps;
  double y;
  double alpha, sec_alpha;
  double vv, hh, vh;
  double Gmm, A;
  int sp;

  Zx = ((double *)params)[0];
  theta_0 = ((double *)params)[1];
  phi_0 = ((double *)params)[2];
  f = ((double *)params)[3];
  s = ((double *)params)[4];
  l = ((double *)params)[5];
  epsr = ((double *)params)[6];
  epsi = ((double *)params)[7];
  dsp =  ((double *)params)[8];
  Gmm = ((double *)params)[9];
  A = ((double *)params)[10];

  //eps = epsr + epsi*I;
  cdouble eps(epsr, epsi);
  sp = (int)dsp;
  
  y = x[0];

  //printf("y=%lf, Gmm=%lf, A=%lf\n", y, Gmm, A);

  Zy = 2.*PI/Gmm *A * sin(2.*PI/Gmm * y);
  alpha = atan(Zy);
  sec_alpha = 1./cos(alpha);
  //printf("in sig_integrand: Zy=%lf, alpha=%lf, sec_alpha=%lf\n",Zy, alpha, sec_alpha);
    
  Calc_scatt_coeff(Zx,Zy, theta_0, phi_0, eps, f, s, l, sp, &vv, &hh, &vh);
     
  fval[0] = vv *sec_alpha;
  fval[1] = hh *sec_alpha;
  fval[2] = vh *sec_alpha;


  //printf("in sig_integrand: vv=%lf\n",vv);

  return 0; // success
}


//*****************************************************************************************
void Calc_scatt_coeff(double Zx,double Zy,double th,double ph,cdouble eps, 
                      double f, double s, double l, int sp,
                      double *sig_ss_vv, double *sig_ss_hh, double *sig_ss_vh)
{
  
  //-Calculates the scattering coefficients of the differential area at a
  //-given point on the periodic surface due to rough surface scattering from
  //-the small scale roughness. 

  double th_rad, ph_rad;
  double v[3], h[3], v_pr[3], h_pr[3], thpr;
  double phs;
  double temp1;
  double sig_s_vv, sig_s_hh, sig_s_hv;

  //-calculate the polarization vectors in both local and global coordinates
  th_rad = th*DEG2RAD; //present angle in radian
  ph_rad = ph*DEG2RAD;

  vectors2(Zx,Zy,th_rad, ph_rad, v,h,v_pr,h_pr,&thpr);

  //-calculate local backscattering coefficients due to small scale roughness

  //-calculate co-pol using I2EM_Bistat_model. Model is run in backscater
  phs = 180.*DEG2RAD;
  I2EM_Bistat_model((float)f, (float)s, (float)l, thpr, thpr, phs, eps, sp, (float)1.0,
                    &sig_s_vv, &sig_s_hh);


  //-- Calculation of x-pol can be done as shown below. However, it will
  //slow down the entire code considerably. It have been commented out here.
  //However, the user may choose to use it by removing the comments.
  //--calculate x-pol using IEMX_model. 
  //IEMX_model((float)f, (float)s, (float)l, thpr, eps, sp,(float)1.0, 0, &sig_s_hv);
  sig_s_hv=0.0;
 
 
  sig_s_vv = pow(10.,0.1*sig_s_vv); // transform from dB power to linear power
  sig_s_hh = pow(10.,0.1*sig_s_hh);
  sig_s_hv = pow(10.,0.1*sig_s_hv);


  //-- The following expressions keep only the strong terms.
  *sig_ss_vv = pow(dot3(v, v_pr),4.) * sig_s_vv;
  
  *sig_ss_hh = pow(dot3(h, h_pr),4.) * sig_s_hh ;

  temp1 = dot3(v, v_pr) * dot3(h_pr, h) + dot3(v, h_pr) * dot3(v_pr, h);
  *sig_ss_vh = temp1*temp1 * sig_s_hv;

}

//************************************************************************************
double dot3( double *a, double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
//************************************************************************************
double norm3( double *a)
{
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
//************************************************************************************
void cross3( double *a, double *b, double *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}



//************************************************************************************
void  vectors2(double Zx,double Zy,double th, double ph,
               double *v, double *h,double *v_pr,double *h_pr,double *thpr)
{

  // Calculates the local coordinates and polarization vectors in global
  // coordinate system. 

  double snt, cst, snp, csp;
  double d0;
  double n[3], ni[3], x_pr[3], y_pr[3], z_pr[3], tmp[3];
  double norm_n, norm_ni, norm_tmp;
  register int i;
  double costhpr, sinthpr;

  snt = sin(th);
  cst = cos(th);
  snp = sin(ph);
  csp = cos(ph);

  d0 = 1./sqrt(1. + Zx*Zx + Zy*Zy);

  //-- calculate polarization vectors, unit normal vector, and incident wave
  // directions in global coordinates

  n[0] = -Zx*d0; // unit normal vector at point on periodic surface
  n[1] = -Zy*d0;
  n[2] = d0;

  norm_n = norm3(n);
  for(i=0;i<3;i++) n[i] /= norm_n;

  
  h[0] = -snp;
  h[1] = csp;
  h[2] = 0.;

  v[0] = -cst*csp;
  v[1] = -cst*snp;
  v[2] = -snt;

  ni[0] = snt * csp; // incident wave direction
  ni[1] = snt * snp;
  ni[2] = -cst;

  norm_ni = norm3(ni);
  for(i=0;i<3;i++) ni[i] /= norm_ni;


  costhpr = -dot3(ni, n);
  *thpr = acos(costhpr);
  sinthpr = sin(*thpr);

  //--calculate local coordinate axes
  for(i=0;i<3;i++) z_pr[i] = n[i];



  if(fabs(fabs(dot3(n,ni))-1.)<1.e-3) {
    for(i=0;i<3;i++) y_pr[i] = h[i];
    //printf("vectors2: y_pr: copied h\n");
  } else  {
    cross3(n, ni, tmp);
    norm_tmp = norm3(tmp);
    for(i=0;i<3;i++) y_pr[i] = tmp[i]/norm_tmp;
    //printf("vectors2: y_pr: calculated\n");
  }

  cross3(y_pr, z_pr, x_pr);

  for(i=0;i<3;i++) h_pr[i] = y_pr[i];

  cross3(h_pr, ni,tmp);
  norm_tmp = norm3(tmp);
  for(i=0;i<3;i++) v_pr[i] = tmp[i]/norm_tmp;

}










//**********************************************************************************
void Emissivity_I2EM_Periodic(double theta_0, double phi_0, cdouble eps, 
                   double f, double s, double l, int sp, 
                   double Gmm, double A, 
                   double *ev, double *eh)
{
  double Zx;
  double evv, ehh;
  // both angles in degrees.

  Zx = 0.0; //partial derivative of periodic surface along x

  do_e_integration(Zx,theta_0,phi_0,f,s,l,eps,sp,Gmm,A, &evv, &ehh);


  *ev =  evv/ Gmm;
  *eh = ehh / Gmm;

}



//*********************************************************************************************
void do_e_integration(double Zx,double theta_0, double phi_0,double f,
                        double s, double l, cdouble er,int sp,
                        double Gmm, double A,
                        double *ev, double *eh)
{

  double params[11];
  double xmin[1], xmax[1];
  double val[3], err[3];

  params[0] = Zx;
  params[1] = theta_0;
  params[2] = phi_0;
  params[3] = f;
  params[4] = s;
  params[5] = l;
  params[6] = real(er);
  params[7] = imag(er);
  params[8] = sp+0.1;
  params[9] = Gmm;
  params[10] = A;

  // integration ranges:
  xmin[0]= 0.0;
  xmax[0]= Gmm;

  pcubature(2, e_integrand, params,
            1, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);

  *ev = val[0];
  *eh = val[1];
  
}

//*****************************************************************************************
int e_integrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval)
{
  double Zx, Zy, theta_0, phi_0, f, s, l, epsr, epsi, dsp;
  //cdouble eps;
  double y;
  double alpha, sec_alpha;
  double vv, hh, vh;
  double Gmm, A;
  int sp;

  Zx = ((double *)params)[0];
  theta_0 = ((double *)params)[1];
  phi_0 = ((double *)params)[2];
  f = ((double *)params)[3];
  s = ((double *)params)[4];
  l = ((double *)params)[5];
  epsr = ((double *)params)[6];
  epsi = ((double *)params)[7];
  dsp =  ((double *)params)[8];
  Gmm = ((double *)params)[9];
  A = ((double *)params)[10];

  cdouble eps(epsr, epsi);
  sp = (int)dsp;
  
  y = x[0];

  //printf("y=%lf, Gmm=%lf, A=%lf\n", y, Gmm, A);

  Zy = 2.*PI/Gmm *A * sin(2.*PI/Gmm * y);
  alpha = atan(Zy);
  sec_alpha = 1./cos(alpha);
  //printf("in sig_integrand: Zy=%lf, alpha=%lf, sec_alpha=%lf\n",Zy, alpha, sec_alpha);
    
  Calc_scatt_emiss_coeff(Zx,Zy, theta_0, phi_0, eps, f, s, l, sp, &vv, &hh);
     
  fval[0] = vv *sec_alpha;
  fval[1] = hh *sec_alpha;


  //printf("in sig_integrand: vv=%lf\n",vv);

  return 0; // success
}


//*****************************************************************************************
void Calc_scatt_emiss_coeff(double Zx,double Zy,double th,double ph,cdouble eps, 
                            double f, double s, double l, int sp,
                            double *ev_ss, double *eh_ss)
{
  
  //-Calculates the scattering coefficients of the differential area at a
  //-given point on the periodic surface due to rough surface scattering from
  //-the small scale roughness. 

  double th_rad, ph_rad;
  double v[3], h[3], v_pr[3], h_pr[3], thpr;
  double phs;
  double temp1;
  double ev_s, eh_s;
  double thpr_deg;

  //-calculate the polarization vectors in both local and global coordinates
  th_rad = th*DEG2RAD; //present angle in radian
  ph_rad = ph*DEG2RAD;

  vectors2(Zx,Zy,th_rad, ph_rad, v,h,v_pr,h_pr,&thpr);

  //-calculate local backscattering coefficients due to small scale roughness

  //-calculate local emissivity due to small-scale roughness using I2EM
  thpr_deg = thpr*180./PI; //need this angle in degrees
  Calc_emissivity_I2EMmodel((float)f, (float)s, (float)l, thpr_deg, eps, sp, &eh_s, &ev_s);


  *ev_ss = dot3(v, h_pr)*dot3(v, h_pr)* eh_s + dot3(v, v_pr)*dot3(v, v_pr)* ev_s ;
  *eh_ss = dot3(h, h_pr)*dot3(h, h_pr)* eh_s + dot3(h, v_pr)*dot3(h, v_pr)* ev_s ;

}




//************************************************************************
void Calc_emissivity_I2EMmodel(double fr, double sig, double L, double theta_d, cdouble er, int sp, double *eh, double *ev)
{
  
  //-- This code calculates the emission from rough surfaces using the I2EM
  //model

  double error;
  double k;
  double theta;
  double ks, kl, ks2, kl2, cs, s, s2;
  cdouble sq, rv0, rh0, rv, rh;
  int pol;
  double refv, refh;

  double params[16];
  double xmin[2], xmax[2];
  double val[2], err[2];


  //printf("Calc_emissivity_I2EMmodel: fr=%lf, sig=%lf, L=%lf, theta_d=%lf, er=%lf+I%lf, sp=%d\n",
  //       fr, sig, L, theta_d, real(er), imag(er), sp);

  sig = sig * 100.; // transform to cm
  L = L* 100.; //transform to cm.


  error = 1.0e8;

  k = 2*PI *fr/30.; // wavenumber in free space. Speed of light is in cm/sec
  theta = theta_d *DEG2RAD;

  ks = k * sig; // roughness parameter
  kl = k * L;  

  ks2 = ks * ks; 
  kl2 = kl * kl;

  cs = cos(theta);
  s = sin(theta);

  s2 = s*s;

  //-- calculation of reflection coefficints
  sq = sqrt(er - s2);

  rv0 = (sqrt(er) - 1.) /(sqrt(er) + 1.);
  rh0 = -rv0;

  rv = (er *cs - sq)/(er*cs +sq);
  rh = (cs - sq)/(cs + sq);

  //printf("Calc_emissivity_I2EMmodel: sq=%lf+I%lf, rv0=%lf+I%lf, rv=%lf+I%lf\n",
  //       real(sq),imag(sq),real(rv0),imag(rv0),real(rv),imag(rv));



  pol = 1;
  params[0] = pol;
  params[1] = sp+0.1;
  params[2] = k;
  params[3] = L;
  params[4] = ks;
  params[5] = cs;
  params[6] = s;
  params[7] = kl;
  params[8] = real(er);
  params[9] = imag(er);
  params[10] = real(sq);
  params[11] = imag(sq);
  params[12] = real(rv);
  params[13] = imag(rv);
  params[14] = real(rh);
  params[15] = imag(rh);


  //refv = dblquad(@(ths,phs)emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs,s, kl, er, sq,rv,rh), 0.0, pi/2, 0, pi);

  // integration ranges:
  xmin[0]=0.;
  xmin[1]=0.;
  xmax[0]=0.5*PI;
  xmax[1]=PI;

  pcubature(2, emsv_integrand, params,
            2, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);
  refv = val[0];
  refh = val[1];

  //printf("Calc_emissivity_I2EMmodel: refv=%lf, refh=%lf\n",refv,refh);



  //pol = 2;
  //params[0] = pol;

  //refh = dblquad(@(ths,phs)emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs,s, kl, er, sq,rv,rh), 0.0, pi/2, 0, pi);

  //pcubature(2, emsv_integrand, params,
  //          2, xmin, xmax, 1000, 0., 1.e-3,
  //          ERROR_L2, val, err);
  //refh = val[0];

  *ev = 1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv);
  *eh = 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh);


  //printf("Calc_emissivity_I2EMmodel: *ev=%lf\n", *ev);

}


double Emissivity_I2EMmodel_VPOL(double fr, double sig, double L, double theta_d, double el, double ei, int sp)
{
  
  //-- This code calculates the emission from rough surfaces using the I2EM
  //model
  cdouble er(el, ei);
  double error;
  double k;
  double theta;
  double ks, kl, ks2, kl2, cs, s, s2;
  cdouble sq, rv0, rh0, rv, rh;
  int pol;
  double refv, refh;

  double params[16];
  double xmin[2], xmax[2];
  double val[2], err[2];


  //printf("Calc_emissivity_I2EMmodel: fr=%lf, sig=%lf, L=%lf, theta_d=%lf, er=%lf+I%lf, sp=%d\n",
  //       fr, sig, L, theta_d, real(er), imag(er), sp);

  sig = sig * 100.; // transform to cm
  L = L* 100.; //transform to cm.


  error = 1.0e8;

  k = 2*PI *fr/30.; // wavenumber in free space. Speed of light is in cm/sec
  theta = theta_d *DEG2RAD;

  ks = k * sig; // roughness parameter
  kl = k * L;  

  ks2 = ks * ks; 
  kl2 = kl * kl;

  cs = cos(theta);
  s = sin(theta);

  s2 = s*s;

  //-- calculation of reflection coefficints
  sq = sqrt(er - s2);

  rv0 = (sqrt(er) - 1.) /(sqrt(er) + 1.);
  rh0 = -rv0;

  rv = (er *cs - sq)/(er*cs +sq);
  rh = (cs - sq)/(cs + sq);

  //printf("Calc_emissivity_I2EMmodel: sq=%lf+I%lf, rv0=%lf+I%lf, rv=%lf+I%lf\n",
  //       real(sq),imag(sq),real(rv0),imag(rv0),real(rv),imag(rv));



  pol = 1;
  params[0] = pol;
  params[1] = sp+0.1;
  params[2] = k;
  params[3] = L;
  params[4] = ks;
  params[5] = cs;
  params[6] = s;
  params[7] = kl;
  params[8] = real(er);
  params[9] = imag(er);
  params[10] = real(sq);
  params[11] = imag(sq);
  params[12] = real(rv);
  params[13] = imag(rv);
  params[14] = real(rh);
  params[15] = imag(rh);


  //refv = dblquad(@(ths,phs)emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs,s, kl, er, sq,rv,rh), 0.0, pi/2, 0, pi);

  // integration ranges:
  xmin[0]=0.;
  xmin[1]=0.;
  xmax[0]=0.5*PI;
  xmax[1]=PI;

  pcubature(2, emsv_integrand, params,
            2, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);
  refv = val[0];
  refh = val[1];

  //printf("Calc_emissivity_I2EMmodel: refv=%lf, refh=%lf\n",refv,refh);



  //pol = 2;
  //params[0] = pol;

  //refh = dblquad(@(ths,phs)emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs,s, kl, er, sq,rv,rh), 0.0, pi/2, 0, pi);

  //pcubature(2, emsv_integrand, params,
  //          2, xmin, xmax, 1000, 0., 1.e-3,
  //          ERROR_L2, val, err);
  //refh = val[0];

  //array<double, 2> ER {1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv), 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh)};
  //return cdouble(1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv), 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh));
  //return ER;

  //*ev = 1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv);
  //*eh = 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh);


  //printf("Calc_emissivity_I2EMmodel: *ev=%lf\n", *ev);


  //double *eh, double *ev
  return 1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv);
}

std::string Emissivity_I2EMmodel(double fr, double sig, double L, double theta_d, double el, double ei, int sp)
{
  
  //-- This code calculates the emission from rough surfaces using the I2EM
  //model
  cdouble er(el, ei);
  double error;
  double k;
  double theta;
  double ks, kl, ks2, kl2, cs, s, s2;
  cdouble sq, rv0, rh0, rv, rh;
  int pol;
  double refv, refh;

  double params[16];
  double xmin[2], xmax[2];
  double val[2], err[2];


  //printf("Calc_emissivity_I2EMmodel: fr=%lf, sig=%lf, L=%lf, theta_d=%lf, er=%lf+I%lf, sp=%d\n",
  //       fr, sig, L, theta_d, real(er), imag(er), sp);

  sig = sig * 100.; // transform to cm
  L = L* 100.; //transform to cm.


  error = 1.0e8;

  k = 2*PI *fr/30.; // wavenumber in free space. Speed of light is in cm/sec
  theta = theta_d *DEG2RAD;

  ks = k * sig; // roughness parameter
  kl = k * L;  

  ks2 = ks * ks; 
  kl2 = kl * kl;

  cs = cos(theta);
  s = sin(theta);

  s2 = s*s;

  //-- calculation of reflection coefficints
  sq = sqrt(er - s2);

  rv0 = (sqrt(er) - 1.) /(sqrt(er) + 1.);
  rh0 = -rv0;

  rv = (er *cs - sq)/(er*cs +sq);
  rh = (cs - sq)/(cs + sq);

  //printf("Calc_emissivity_I2EMmodel: sq=%lf+I%lf, rv0=%lf+I%lf, rv=%lf+I%lf\n",
  //       real(sq),imag(sq),real(rv0),imag(rv0),real(rv),imag(rv));



  pol = 1;
  params[0] = pol;
  params[1] = sp+0.1;
  params[2] = k;
  params[3] = L;
  params[4] = ks;
  params[5] = cs;
  params[6] = s;
  params[7] = kl;
  params[8] = real(er);
  params[9] = imag(er);
  params[10] = real(sq);
  params[11] = imag(sq);
  params[12] = real(rv);
  params[13] = imag(rv);
  params[14] = real(rh);
  params[15] = imag(rh);


  //refv = dblquad(@(ths,phs)emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs,s, kl, er, sq,rv,rh), 0.0, pi/2, 0, pi);

  // integration ranges:
  xmin[0]=0.;
  xmin[1]=0.;
  xmax[0]=0.5*PI;
  xmax[1]=PI;

  pcubature(2, emsv_integrand, params,
            2, xmin, xmax, 1000, 0., 1.e-3,
            ERROR_L2, val, err);

  refv = val[0];
  refh = val[1];

  //printf("Calc_emissivity_I2EMmodel: refv=%lf, refh=%lf\n",refv,refh);



  //pol = 2;
  //params[0] = pol;

  //refh = dblquad(@(ths,phs)emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs,s, kl, er, sq,rv,rh), 0.0, pi/2, 0, pi);

  //pcubature(2, emsv_integrand, params,
  //          2, xmin, xmax, 1000, 0., 1.e-3,
  //          ERROR_L2, val, err);
  //refh = val[0];

  //array<double, 2> ER {1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv), 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh)};
  //return cdouble(1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv), 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh));
  //return ER;

  //*ev = 1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv);
  //*eh = 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh);


  //printf("Calc_emissivity_I2EMmodel: *ev=%lf\n", *ev);

  //double *eh, double *ev


  char buff[100];
  snprintf(buff, sizeof(buff), "%f,%f", 1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv), 1.- refh - exp(-ks2 * cs*cs)* abs(rh)*abs(rh));
  std::string buffAsStdStr = buff;
  return buffAsStdStr;
  
  //return 1.- refv - exp(-ks2 * cs*cs)* abs(rv)*abs(rv);
}

//*************************************************************************************
int emsv_integrand(unsigned ndim, const double *x, void *params, 
                   unsigned fdim, double *fval)
{
  double ths, phs;

  int pol, sp;
  double k, L, ks, cs, s, kl;
  //cdouble sq, rv, rh;
  //cdouble er;

  int n_spec;
  double error;
  double cs2, ks2, s2, kl2;
  double css, css2, ss, ss2, sf, sf2, csf;
  cdouble sqs, rc, tv, th, tv2, th2;
  double wvnb, thsmin, costhsmin;
  int nspec;
  double scale;

  double wn[1000];
  double ff;
  cdouble fvv, fhh, fvh, fhv;
  double cm1;
  cdouble T, cm2;
  double ex, de;
  cdouble Fvv, Fhv, Fvvs, Fhvs;
  double svv, shh;
  cdouble Ivv, Ihh, Ihv, Ivh;
  double ref;
  cdouble Fhh, Fvh, Fhhs, Fvhs;
  cdouble hh, vv, hv, vh;
  int i, n;
  double wnn;



  ths = x[0];
  phs = x[1];


  pol = (int) (((double *)params)[0]);
  sp = (int) (((double *)params)[1]);
  k = ((double *)params)[2];
  L = ((double *)params)[3];
  ks = ((double *)params)[4];
  cs = ((double *)params)[5];
  s = ((double *)params)[6];
  kl = ((double *)params)[7];
  //er = ((double *)params)[8] + ((double *)params)[9]*I;
  cdouble er(((double *)params)[8], ((double *)params)[9]);
  //sq = ((double *)params)[10] + ((double *)params)[11]*I;
  cdouble sq(((double *)params)[10], ((double *)params)[11]);
  //rv = ((double *)params)[12]+ ((double *)params)[13]*I;
  cdouble rv(((double *)params)[12], ((double *)params)[13]);
  //rh = ((double *)params)[14]+ ((double *)params)[15]*I;

  cdouble rh(((double *)params)[14], ((double *)params)[15]);

  //printf("emsv_integrand: pol=%d, sp=%d, k=%lf, L=%lf, ks=%lf, cs=%lf, s=%lf, kl=%lf\n",
  //       pol,sp,k,L,ks,cs,s,kl);
  //printf("emsv_integrand: er=%lf+I%lf, sq=%lf+I%lf, rv=%lf+I%lf, rh%lf+I%lf\n",
  //       real(er), imag(er), real(sq), imag(sq), real(rv), imag(rv), real(rh), imag(rh));


  error = 1.0e3;
  
  cs2 = cs*cs;
  ks2 = ks*ks;
  s2 = s*s;
  kl2 = kl* kl;

  css = cos(ths);
  css2 = css * css;

  ss = sin(ths);
  ss2 = ss * ss;

  sf = sin(phs);
  sf2 = sf *sf;
  csf= cos(phs);


  sqs = sqrt(er - ss2);

    // trace fvvs nan:
  //printf("emsv_integrand: sqs=%lf+I%lf, \n",real(sqs),imag(sqs));


  rc = (rv - rh) / 2.;
  tv = 1. + rv;
  th = 1. + rh; 
  tv2 = tv * tv;
  th2 = th * th;

  //-- calc coefficients for surface correlation spectra
  wvnb = k * sqrt(s2 - 2.*s*ss*csf + ss2); 

  //printf("emsv_integrand: wvnb=%lf\n",wvnb);

  //nr = length(ths);
  //thsmin = min(ths); % selected as the smallest angle in ths vector
  //costhsmin = cos(thsmin); 
  thsmin = ths;
  costhsmin = cos(thsmin);

  //-- calculate number of spectral components needed 
  n_spec = 1; 
  scale =1.;
  while (error > 1.0e-3) {
    n_spec++;
    scale /= n_spec;
    //%   error = (ks2 .*(cs + css).^2 ).^n_spec ./ factorial(n_spec); 
    //%---- in this case we will use the smallest ths to determine the number of
    //%spectral components to use. It might be more than needed for other angles
    //%but this is fine. This option is used to simplify calculations.
    error = pow(ks2 *(cs + costhsmin)*(cs + costhsmin),n_spec) *scale; //   was: / factorial(n_spec); 
  }// endwhile


  //%-- calculate expressions for the surface spectra
  roughness_spectrum_12_2(sp, kl2, L, wvnb, n_spec, wn);



  //%-- calculate fpq!

  ff = 2.*(s* ss - (1. + cs * css) * csf)/(cs + css);

  fvv =  rv * ff;
  fhh = -rh * ff;

  fvh = -2.* rc * sf;
  fhv =  2.* rc * sf;

  //printf("emsv_integrand: fvh=%lf+I%lf\n",real(fvh),imag(fvh));


  //-- calculate Fpq and Fpqs -----
  cm1 = s *(ss - s *csf) /(cs2 * css);
  T = (sq *(cs + sq) + cs*(er*cs + sq)) / (er *cs *(cs+sq) + sq *(er*cs +sq));
  cm2 = css * sq /cs /sqs - 1.;
    // trace fvvs nan:
  //printf("emsv_integrand: cm2=%lf+I%lf,  css=%lf, sq=%lf+I%lf, cs=%lf, sqs=\n",
  //       real(cm2),imag(cm2),css, real(sq),imag(sq), cs,  real(sqs),imag(sqs));

  ex = exp(-ks2 *cs *css);
  de = 0.5 *exp(-ks2 *(cs2 + css2)); 




  // pol ==1:
    Fvv = (er - 1.) *s2 *tv2 *cm1 / (er*er);
    Fhv = (T *s *s -1.+ cs/css + (er *T *cs *css *(er *T - s *s) - 
            sq* sq) / (T *er *sq * css)) *(1. - rc *rc) *sf;
 
    Fvvs = -cm2 *sq *tv2 *(csf - s *ss) /cs2 /er - cm2 * sqs *tv2 *csf /er
      -(css *sq /cs /sqs /er -1.) * ss * tv2 * (s - ss *csf) /cs;
    Fhvs = -(ss*ss /T -1. + css /cs + (cs *css *(1.- ss*ss*T) - 
            T *T*sqs*sqs) /(T *sqs*cs))*(1.-rc *rc)*sf;

    //printf("emsv_integrand: Fvv=%lf+I%lf, Fvvs=%lf+I%lf\n",real(Fvv),imag(Fvv),
    //       real(Fvvs),imag(Fvvs));

    // trace fvvs nan:
    //printf("emsv_integrand: cm2=%lf+I%lf,, sq=%lf+I%lf,, tv2=%lf+I%lf,, csf=%lf, s=%lf, ss=%lf, cs2=%lf, er=%lf+I%lf, sqs=%lf+I%lf,, css=%lf, cs=%lf\n",
    //       real(cm2),imag(cm2),real(sq),imag(sq),real(tv2),imag(tv2),
    //       csf,s,ss,cs2,real(er),imag(er),real(sqs),imag(sqs),css,cs);


    //-- calculate the bistatic field coefficients ---

    scale =1.;
    ref =0.;
    for(i=0;i<n_spec;i++){
      n=i+1;
      scale /= (1.*n);
      Ivv = fvv *ex *pow(ks *(cs + css),n) + 0.5*(Fvv *pow(ks *css,n) + Fvvs *pow(ks*cs,n) );
      Ihv = fhv *ex *pow(ks *(cs + css),n) + 0.5*(Fhv *pow(ks *css,n) + Fhvs *pow(ks*cs,n) );

      wnn = wn[i] * scale;    // was: /factorial(n);
      vv = wnn * abs(Ivv)*abs(Ivv);
      hv = wnn * abs(Ihv)*abs(Ihv);
      svv = real(de *(vv + hv) * ss /4./PI/cs);

      //printf("emsv_integrand: i=%d, svv=%lf\n",i,svv);

      ref += svv;
    }//endfor
    
    fval[0] = ref;
    //printf("emsv_integrand: fval[0]=%lf\n",fval[0]);



    // pol ==2: 
    Fhh = -(er -1.) * th2 * cm1;
    Fvh = (s *s /T -1.+ cs/css + (cs *css *(1.- s *s *T) - 
           T*T*sq*sq)/ (T*sq*css))*(1. - rc*rc)*sf;

    Fhhs = cm2 * sq *th2 *(csf - s *ss) /cs2 + cm2 * sqs * th2 * csf + 
           cm2 *ss *th2 *(s - ss*csf)/cs;
    Fvhs = -(T * ss * ss - 1. + css/cs + (er *T *cs *css *(er*T - ss *ss)
           - sqs * sqs) /(T * er *sqs *cs)) *(1.- rc *rc)*sf;

    scale=1.;
    ref=0.;
    for(i=0;i<n_spec;i++){
      n=i+1;
      scale /= (1.*n);
      Ihh = fhh *ex *pow(ks *(cs + css),n) + 0.5*(Fhh *pow(ks *css,n) + Fhhs *pow(ks*cs,n) );
      Ivh = fvh *ex *pow(ks *(cs + css),n) + 0.5*(Fvh *pow(ks *css,n) + Fvhs *pow(ks*cs,n) );
  
      wnn = wn[i]*scale;
      hh = wnn * abs(Ihh)*abs(Ihh);
      vh = wnn * abs(Ivh)*abs(Ivh);
      shh = real(de *(hh + vh) * ss /4./PI/cs);
      ref += shh;
    }// endfor
  
  
  fval[1]= ref;

  return 0;
}
  











//*********************************************************************************************
void roughness_spectrum_12_2(int correl_func, double kl2, double L, 
                        double wvnb, int np, double *wn)
{
  int i, n;

  switch(correl_func)
    {
    case CORREL_EXPONENTIAL:
      for (i=0;i<np;i++){
        n=i+1;
        wn[i] = n*kl2/pow(n*n + wvnb*L*wvnb*L,1.5);
      }// endfor
      break;


    case CORREL_GAUSSIAN:
      for (i=0;i<np;i++){
        n=i+1;
        wn[i] =  0.5*kl2/n* exp( -(wvnb*L)*(wvnb*L)/(4.*n) );
      }// endfor
      break;

    }
}






//************************************************************
void ReflTransm_PlanarBoundary(cdouble eps1, cdouble eps2, 
                               double theta1d, 
                               cdouble *rhoh, cdouble *rhov,
                               double *gammah, double *gammav,
                               cdouble *tauh, cdouble *tauv,
                               double *Th, double *Tv)
{
  double theta1;
  double theta, h, v;
  cdouble sin_theta2;
  cdouble cos_theta2;

  theta1 = theta1d*PI/180.;
    
  sin_theta2 = sqrt(eps1)/sqrt(eps2)*sin(theta1);
  cos_theta2 = sqrt(1. - sin_theta2*sin_theta2);
    
  *rhoh = (sqrt(eps1)*cos(theta1)-sqrt(eps2)*cos_theta2) / (sqrt(eps1)*cos(theta1) + sqrt(eps2)*cos_theta2);
  *rhov = (sqrt(eps1)*cos_theta2 -sqrt(eps2)*cos(theta1))/ (sqrt(eps1)*cos_theta2  + sqrt(eps2)*cos(theta1));

  *tauh = 1. + *rhoh;
  *tauv = (1. + *rhov)*(cos(theta1)/cos_theta2);
       

  *gammah = abs(*rhoh)*abs(*rhoh);
  *gammav = abs(*rhov)*abs(*rhov);
        
  *Th = 1.- *gammah;
  *Tv = 1.- *gammav;

}

//********************************************************************
void ZRTemission_veg(cdouble er, double s, double l,
                     double freq_ghz, double thi, double albedo, 
                     double extinction, double depth, 
                     double *e_v, double *e_h)
{
  double k, thetar;
  cdouble rhoh, rhov, tauh, tauv;
  double gammah, gammav, Th, Tv;
  //cdouble one;
  double gamma_coh_v, gamma_coh_h;
  double temp1;
  double e_h_soil, e_v_soil;
  double veg_att;
  double kappa_e, a;

  kappa_e  = extinction;
  a = albedo;
  k = 2.*PI*freq_ghz/0.3; // wave number
  thetar = thi * PI/180.;

  //--call function to calculate reflectivity of planar boundary!
  cdouble one(1.0, 0.0);
  ReflTransm_PlanarBoundary(one, er, thi, 
                            &rhoh, &rhov, &gammah, &gammav,&tauh, &tauv, &Th, &Tv);



  //--calculate coherent reflectivity
  temp1 = 2.*k*s* cos(thetar);
  temp1 *= temp1;
  gamma_coh_v = gammav * exp(-temp1);
  gamma_coh_h = gammah * exp(-temp1);


  // emissivity of soil. Use I2EM emissivity code for rough surface
  Calc_emissivity_I2EMmodel(freq_ghz, s, l, thi, er, 1, &e_h_soil, &e_v_soil);


  //attentuation through vegetation layer
  veg_att = exp(-kappa_e*depth /cos(thetar));

  //total emissivity of the vegetation layer over rough soil
  *e_v = (1.+ gamma_coh_v * veg_att) *(1.-a)*(1.- veg_att) + veg_att* e_v_soil;

  *e_h = (1.+ gamma_coh_h * veg_att) *(1.-a)*(1.- veg_att) + veg_att* e_h_soil;

}





//*************************************************************************
void ZRTemission_DUB(cdouble eps2,cdouble eps3, 
                     double theta_i, double f, double s, 
                     double l, double a, double d, double kappa_e,
                     double *e_v, double *e_h) 
{
  
  double theta1, theta2;
  cdouble n1, n2;
  cdouble rhoh, rhov;
  double gammah12, gammav12;
  cdouble tauh, tauv;
  double Th, Tv;
  double eh23, ev23;
  double gammav23, gammah23;
  double trans;

  //eps1 = 1.0 + 0.0*I;
  cdouble eps1(1.0, 0.0);
  theta1 = theta_i*PI/180.; //transform to radians
  n2 = sqrt(eps2); //calc index of refraction
  n1 = cdouble (1.0, 0.0);

  theta2 = asin(abs(n1/n2) * sin(theta1)); // incidence angle in medium 2 (refracted)

  //printf("DUB: theta1=%lf, theta2=%lf\n",theta1, theta2); 
  //-- calculate reflectivies at the two interfaces

  //- 12 interface
  ReflTransm_PlanarBoundary(eps1, eps2, theta_i, 
                            &rhoh, &rhov, &gammah12, &gammav12,
                            &tauh, &tauv, &Th, &Tv);

  //printf("DUB: gammah12=%lf\n",gammah12);

  //- 23 interface - use I2EM model 
  Calc_emissivity_I2EMmodel(f,s,l,theta2*180./PI, eps3/eps2, 1, &eh23, &ev23);
  gammav23 = 1.- ev23;
  gammah23 = 1.- eh23;

  //printf("DUB: eh23=%lf\n",eh23);


  trans = exp(-kappa_e * d /cos(theta2)); //extinction coefficient inside medium.

  //printf("DUB: trans=%lf\n",trans);

  *e_v = ( (1. - gammav12)/(1. - gammav12*gammav23*trans*trans))* 
         ( (1.+gammav23*trans)*(1.-a)*(1. - trans)+ (1. - gammav23)*trans);

  *e_h = ( (1. - gammah12)/(1. - gammah12*gammah23*trans*trans))*
         ( (1.+gammah23*trans)*(1.-a)*(1. - trans)+ (1. - gammah23)*trans);

  //printf("DUB: *e_h=%lf\n",*e_h);

}