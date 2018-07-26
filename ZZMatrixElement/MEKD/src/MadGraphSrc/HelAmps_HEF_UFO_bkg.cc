//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.3, 2012-11-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "HelAmps_HEF_UFO_bkg.h"

namespace MG5_HEF_UFO_bkg
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[im]; 
      fo[3] = ip * nsf * sqm[im]; 
      fo[4] = im * nsf * sqm[ - ip]; 
      fo[5] = ip * sqm[ - ip]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void FFV43_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP11; 
  complex<double> TMP10; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP11 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP10 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (OM3 * P3[0] * (+2. * (TMP11) - TMP10) + (F2[4] *
      F1[2] + F2[5] * F1[3] - 2. * (F2[2] * F1[4] + F2[3] * F1[5])));
  V3[3] = denom * - cI * (OM3 * P3[1] * (+2. * (TMP11) - TMP10) + (-2. * (F2[3]
      * F1[4] + F2[2] * F1[5]) - F2[5] * F1[2] - F2[4] * F1[3]));
  V3[4] = denom * - cI * (OM3 * P3[2] * (+2. * (TMP11) - TMP10) + (-cI * (F2[5]
      * F1[2]) + cI * (F2[4] * F1[3]) - 2. * cI * (F2[3] * F1[4]) + 2. * cI *
      (F2[2] * F1[5])));
  V3[5] = denom * - cI * (OM3 * P3[3] * (+2. * (TMP11) - TMP10) + (F2[5] *
      F1[3] + 2. * (F2[3] * F1[5]) - 2. * (F2[2] * F1[4]) - F2[4] * F1[2]));
}


void FFV43_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP13; 
  TMP13 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP12 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * (-cI * (TMP12) + 2. * cI * (TMP13)); 
}


void FFV44_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 2. * (V3[2] - V3[5]) + 2. * (F1[5]
      * (+cI * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * - 2. * (V3[3] + cI * (V3[4])) + 2. *
      (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * 2. * (V3[2] - V3[5]) + (P2[1] * - 2. *
      (V3[3] + cI * (V3[4])) + (P2[2] * 2. * (+cI * (V3[3]) - V3[4]) + 2. *
      (P2[3] * (V3[2] - V3[5]))))) + (F1[5] * (P2[0] * 2. * (+cI * (V3[4]) -
      V3[3]) + (P2[1] * 2. * (V3[2] + V3[5]) + (P2[2] * - 2. * (+cI * (V3[2] +
      V3[5])) + 2. * (P2[3] * (+cI * (V3[4]) - V3[3]))))) + M2 * (F1[2] *
      (V3[2] + V3[5]) + F1[3] * (V3[3] - cI * (V3[4])))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 2. * (V3[3] + cI * (V3[4])) + (P2[1]
      * 2. * (V3[2] - V3[5]) + (P2[2] * 2. * (-cI * (V3[5]) + cI * (V3[2])) +
      2. * (P2[3] * (V3[3] + cI * (V3[4])))))) + (F1[5] * (P2[0] * 2. * (V3[2]
      + V3[5]) + (P2[1] * 2. * (+cI * (V3[4]) - V3[3]) + (P2[2] * - 2. * (V3[4]
      + cI * (V3[3])) - 2. * (P2[3] * (V3[2] + V3[5]))))) + M2 * (F1[2] *
      (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] - V3[5]))));
}


void FFV42_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))));
  F2[4] = denom * cI * M2 * (F1[2] * (V3[2] + V3[5]) + F1[3] * (V3[3] - cI *
      (V3[4])));
  F2[5] = denom * cI * M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] -
      V3[5]));
}

void FFV42_43_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV42_2(F1, V3, COUP1, M2, W2, F2); 
  FFV43_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV42_44_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV42_2(F1, V3, COUP1, M2, W2, F2); 
  FFV44_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV42_45_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV42_2(F1, V3, COUP1, M2, W2, F2); 
  FFV45_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV44_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP13; 
  TMP13 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP12 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * (-cI * (TMP12) - 2. * cI * (TMP13)); 
}


void FFV42_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  TMP12 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - cI * TMP12; 
}

void FFV42_45_0(complex<double> F1[], complex<double> F2[], complex<double>
    V3[], complex<double> COUP1, complex<double> COUP2, complex<double> &
    vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV42_0(F1, F2, V3, COUP1, vertex); 
  FFV45_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV42_43_0(complex<double> F1[], complex<double> F2[], complex<double>
    V3[], complex<double> COUP1, complex<double> COUP2, complex<double> &
    vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV42_0(F1, F2, V3, COUP1, vertex); 
  FFV43_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV42_44_0(complex<double> F1[], complex<double> F2[], complex<double>
    V3[], complex<double> COUP1, complex<double> COUP2, complex<double> &
    vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV42_0(F1, F2, V3, COUP1, vertex); 
  FFV44_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV45_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP11; 
  complex<double> TMP10; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP11 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP10 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (OM3 * - P3[0] * (TMP10 + 4. * (TMP11)) + (F2[4] *
      F1[2] + F2[5] * F1[3] + 4. * (F2[2] * F1[4] + F2[3] * F1[5])));
  V3[3] = denom * - cI * (OM3 * - P3[1] * (TMP10 + 4. * (TMP11)) + (+4. *
      (F2[3] * F1[4] + F2[2] * F1[5]) - F2[5] * F1[2] - F2[4] * F1[3]));
  V3[4] = denom * - cI * (OM3 * - P3[2] * (TMP10 + 4. * (TMP11)) + (-4. * cI *
      (F2[2] * F1[5]) - cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) + 4. * cI *
      (F2[3] * F1[4])));
  V3[5] = denom * - cI * (OM3 * - P3[3] * (TMP10 + 4. * (TMP11)) + (F2[5] *
      F1[3] + 4. * (F2[2] * F1[4]) - 4. * (F2[3] * F1[5]) - F2[4] * F1[2]));
}


void FFV41_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP14; 
  TMP14 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * - cI * TMP14; 
}


void FFV41_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * - 1. * (V3[3] + cI * (V3[4])) +
      F1[5] * (V3[2] + V3[5]))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * (V3[2] - V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) + P2[3] *
      (V3[2] - V3[5])))) + (F1[5] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] *
      (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[2] * (V3[2] + V3[5]) + F1[3] * (V3[3] -
      cI * (V3[4])))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void FFV43_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 2. * (V3[5] - V3[2]) + 2. * (F1[5]
      * (V3[3] - cI * (V3[4]))))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * 2. * (V3[3] + cI * (V3[4])) - 2. *
      (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * 2. * (V3[5] - V3[2]) + (P2[1] * 2. *
      (V3[3] + cI * (V3[4])) + (P2[2] * 2. * (V3[4] - cI * (V3[3])) + 2. *
      (P2[3] * (V3[5] - V3[2]))))) + (F1[5] * (P2[0] * 2. * (V3[3] - cI *
      (V3[4])) + (P2[1] * - 2. * (V3[2] + V3[5]) + (P2[2] * 2. * (+cI * (V3[2]
      + V3[5])) + 2. * (P2[3] * (V3[3] - cI * (V3[4])))))) + M2 * (F1[2] *
      (V3[2] + V3[5]) + F1[3] * (V3[3] - cI * (V3[4])))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * 2. * (V3[3] + cI * (V3[4])) + (P2[1] *
      2. * (V3[5] - V3[2]) + (P2[2] * 2. * (-cI * (V3[2]) + cI * (V3[5])) - 2.
      * (P2[3] * (V3[3] + cI * (V3[4])))))) + (F1[5] * (P2[0] * - 2. * (V3[2] +
      V3[5]) + (P2[1] * 2. * (V3[3] - cI * (V3[4])) + (P2[2] * 2. * (V3[4] + cI
      * (V3[3])) + 2. * (P2[3] * (V3[2] + V3[5]))))) + M2 * (F1[2] * (V3[3] +
      cI * (V3[4])) + F1[3] * (V3[2] - V3[5]))));
}


void FFV44_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP11; 
  complex<double> TMP10; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP11 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP10 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (OM3 * - P3[0] * (TMP10 + 2. * (TMP11)) + (F2[4] *
      F1[2] + F2[5] * F1[3] + 2. * (F2[2] * F1[4] + F2[3] * F1[5])));
  V3[3] = denom * - cI * (OM3 * - P3[1] * (TMP10 + 2. * (TMP11)) + (+2. *
      (F2[3] * F1[4] + F2[2] * F1[5]) - F2[5] * F1[2] - F2[4] * F1[3]));
  V3[4] = denom * - cI * (OM3 * - P3[2] * (TMP10 + 2. * (TMP11)) + (-cI *
      (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - 2. * cI * (F2[2] * F1[5]) + 2. *
      cI * (F2[3] * F1[4])));
  V3[5] = denom * - cI * (OM3 * - P3[3] * (TMP10 + 2. * (TMP11)) + (F2[5] *
      F1[3] + 2. * (F2[2] * F1[4]) - 2. * (F2[3] * F1[5]) - F2[4] * F1[2]));
}


void FFV42_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP10; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP10 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP10); 
  V3[3] = denom * - cI * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 *
      TMP10);
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2]
      * OM3 * TMP10);
  V3[5] = denom * - cI * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP10); 
}

void FFV42_43_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P3[4];
//   double OM3;
  int i; 
  complex<double> Vtmp[6]; 
  FFV42_3(F1, F2, COUP1, M3, W3, V3); 
  FFV43_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV42_44_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P3[4]; 
//   double OM3;
  int i; 
  complex<double> Vtmp[6]; 
  FFV42_3(F1, F2, COUP1, M3, W3, V3); 
  FFV44_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV42_45_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P3[4]; 
//   double OM3; 
  int i; 
  complex<double> Vtmp[6]; 
  FFV42_3(F1, F2, COUP1, M3, W3, V3); 
  FFV45_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFV44_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * 2. * (V3[5] - V3[2]) + (P1[1] * 2. *
      (V3[3] - cI * (V3[4])) + (P1[2] * 2. * (V3[4] + cI * (V3[3])) + 2. *
      (P1[3] * (V3[5] - V3[2]))))) + (F2[3] * (P1[0] * 2. * (V3[3] + cI *
      (V3[4])) + (P1[1] * - 2. * (V3[2] + V3[5]) + (P1[2] * - 2. * (+cI *
      (V3[2] + V3[5])) + 2. * (P1[3] * (V3[3] + cI * (V3[4])))))) + M1 * (F2[4]
      * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4])))));
  F1[3] = denom * cI * (F2[2] * (P1[0] * 2. * (V3[3] - cI * (V3[4])) + (P1[1] *
      2. * (V3[5] - V3[2]) + (P1[2] * 2. * (-cI * (V3[5]) + cI * (V3[2])) + 2.
      * (P1[3] * (+cI * (V3[4]) - V3[3]))))) + (F2[3] * (P1[0] * - 2. * (V3[2]
      + V3[5]) + (P1[1] * 2. * (V3[3] + cI * (V3[4])) + (P1[2] * 2. * (V3[4] -
      cI * (V3[3])) + 2. * (P1[3] * (V3[2] + V3[5]))))) + M1 * (F2[4] * (V3[3]
      - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 2. * (V3[5] - V3[2]) + 2. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + M1 * (F2[2] * 2. * (V3[3] - cI * (V3[4])) - 2. *
      (F2[3] * (V3[2] + V3[5])))));
}


void FFV45_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP13; 
  TMP13 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP12 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * (-cI * (TMP12) - 4. * cI * (TMP13)); 
}


void FFV41_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (F2[2] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) + P1[3] *
      (V3[2] - V3[5])))) + (F2[3] * (P1[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[4] * - 1. * (V3[2] + V3[5]) - F2[5]
      * (V3[3] + cI * (V3[4])))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + M1 * (F2[2] * (V3[3] - cI * (V3[4])) - F2[3] *
      (V3[2] + V3[5]))));
}


void FFV42_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])));
  F1[3] = denom * cI * M1 * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] -
      V3[5]));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))));
}

void FFV42_44_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV42_1(F2, V3, COUP1, M1, W1, F1); 
  FFV44_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV45_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 4. * (V3[2] - V3[5]) + 4. * (F1[5]
      * (+cI * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * - 4. * (V3[3] + cI * (V3[4])) + 4. *
      (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * 4. * (V3[2] - V3[5]) + (P2[1] * - 4. *
      (V3[3] + cI * (V3[4])) + (P2[2] * 4. * (+cI * (V3[3]) - V3[4]) + 4. *
      (P2[3] * (V3[2] - V3[5]))))) + (F1[5] * (P2[0] * 4. * (+cI * (V3[4]) -
      V3[3]) + (P2[1] * 4. * (V3[2] + V3[5]) + (P2[2] * - 4. * (+cI * (V3[2] +
      V3[5])) + 4. * (P2[3] * (+cI * (V3[4]) - V3[3]))))) + M2 * (F1[2] *
      (V3[2] + V3[5]) + F1[3] * (V3[3] - cI * (V3[4])))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 4. * (V3[3] + cI * (V3[4])) + (P2[1]
      * 4. * (V3[2] - V3[5]) + (P2[2] * 4. * (-cI * (V3[5]) + cI * (V3[2])) +
      4. * (P2[3] * (V3[3] + cI * (V3[4])))))) + (F1[5] * (P2[0] * 4. * (V3[2]
      + V3[5]) + (P2[1] * 4. * (+cI * (V3[4]) - V3[3]) + (P2[2] * - 4. * (V3[4]
      + cI * (V3[3])) - 4. * (P2[3] * (V3[2] + V3[5]))))) + M2 * (F1[2] *
      (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] - V3[5]))));
}


void FFV41_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP15; 
  double P3[4]; 
  double OM3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP15 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3]
      * F1[5] - P3[0] * OM3 * TMP15);
  V3[3] = denom * - cI * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] * F1[2] - F2[4]
      * F1[3] - P3[1] * OM3 * TMP15);
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4] *
      F1[3] + F2[3] * F1[4]) - P3[2] * OM3 * TMP15);
  V3[5] = denom * - cI * (F2[5] * F1[3] + F2[2] * F1[4] - F2[4] * F1[2] - F2[3]
      * F1[5] - P3[3] * OM3 * TMP15);
}


}  // end namespace $(namespace)s_HEF_UF

