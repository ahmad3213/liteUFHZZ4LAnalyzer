//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "HelAmps_HEF_MEKD2_1.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
using namespace std; 

namespace MG5_HEF_MEKD2_1 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
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
      fo[2] = im * sqm[abs(ip)]; 
      fo[3] = ip * nsf * sqm[abs(ip)]; 
      fo[4] = im * nsf * sqm[abs(im)]; 
      fo[5] = ip * sqm[abs(im)]; 
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


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = pow(0.5, 0.5); 
  sqs = pow(0.5/3, 0.5); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], pow(pt2 + p[3] * p[3], 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(fabs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      std::cerr <<  "Invalid helicity in txxxxx.\n"; 
      std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
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

void FFV7_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
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
  F2[4] = denom * - 2. * cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + (+1./2. * (M2 * (F1[3] * (+cI * (V3[4]) - V3[3]) + 2. *
      (F1[2] * - 1./2. * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI *
      (V3[4])) + (P2[1] * - 1. * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] +
      V3[5])) + P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * - 2. * cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./2. * (M2 * (F1[3] * (V3[5] - V3[2]) + 2.
      * (F1[2] * - 1./2. * (V3[3] + cI * (V3[4]))))) + F1[5] * (P2[0] * - 1. *
      (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI
      * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void FFV2_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
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
  F2[3] = denom * - cI * (F1[2] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * - cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void VVT12_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP19; 
  complex<double> TMP31; 
  complex<double> TMP16; 
  complex<double> TMP15; 
  complex<double> TMP32; 
  complex<double> TMP18; 
  complex<double> TMP9; 
  complex<double> TMP13; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP12 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP19 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP8 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP15 = (P2[0] * - 1. * (T3[6] * V1[3] + T3[10] * V1[4] + T3[14] * V1[5] -
      T3[2] * V1[2]) + (P2[1] * (T3[7] * V1[3] + T3[11] * V1[4] + T3[15] *
      V1[5] - T3[3] * V1[2]) + (P2[2] * (T3[8] * V1[3] + T3[12] * V1[4] +
      T3[16] * V1[5] - T3[4] * V1[2]) + P2[3] * (T3[9] * V1[3] + T3[13] * V1[4]
      + T3[17] * V1[5] - T3[5] * V1[2]))));
  TMP32 = (V1[2] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (V1[3] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (V1[4] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + V1[5] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP31 = (V1[2] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (V1[3] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (V1[4] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + V1[5] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP16 = (P2[0] * - 1. * (T3[3] * V1[3] + T3[4] * V1[4] + T3[5] * V1[5] -
      T3[2] * V1[2]) + (P2[1] * (T3[7] * V1[3] + T3[8] * V1[4] + T3[9] * V1[5]
      - T3[6] * V1[2]) + (P2[2] * (T3[11] * V1[3] + T3[12] * V1[4] + T3[13] *
      V1[5] - T3[10] * V1[2]) + P2[3] * (T3[15] * V1[3] + T3[16] * V1[4] +
      T3[17] * V1[5] - T3[14] * V1[2]))));
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP13 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP18 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  vertex = COUP * (TMP12 * - 1. * (+cI * (TMP10 + TMP11)) + (TMP13 * - 1. *
      (+cI * (TMP31 + TMP32)) + (TMP8 * (+cI * (TMP15 + TMP16)) + TMP9 * (+cI *
      (TMP18 + TMP19)))));
}


void FFV5_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  complex<double> TMP20; 
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
  TMP20 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP20); 
  V3[3] = denom * - cI * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 *
      TMP20);
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2]
      * OM3 * TMP20);
  V3[5] = denom * - cI * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP20); 
}

void FFV5_7_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P3[4]; 
//   double OM3; 
  int i; 
  complex<double> Vtmp[6]; 
  FFV5_3(F1, F2, COUP1, M3, W3, V3); 
  FFV7_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void VVT11_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP31; 
  complex<double> TMP32; 
  TMP32 = (V1[2] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (V1[3] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (V1[4] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + V1[5] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP31 = (V1[2] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (V1[3] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (V1[4] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + V1[5] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  vertex = COUP * - 1. * (+cI * (TMP31 + TMP32)); 
}


void FFV5_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
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
  F1[3] = denom * - cI * M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5]
      - V3[2]));
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

void FFV5_7_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV5_1(F2, V3, COUP1, M1, W1, F1); 
  FFV7_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void VVT7_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP12 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  vertex = COUP * - TMP12 * (+cI * (TMP10 + TMP11)); 
}


void FFT2_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP5; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP6; 
  complex<double> TMP4; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP5 = (P2[0] * (F1[2] * (F2[4] * (T3[2] + T3[14]) + F2[5] * (T3[6] + cI *
      (T3[10]))) + (F1[3] * (F2[4] * (T3[6] - cI * (T3[10])) + F2[5] * (T3[2] -
      T3[14])) + (F1[4] * (F2[2] * (T3[2] - T3[14]) - F2[3] * (T3[6] + cI *
      (T3[10]))) + F1[5] * (F2[2] * (+cI * (T3[10]) - T3[6]) + F2[3] * (T3[2] +
      T3[14]))))) + (P2[1] * (F1[2] * (F2[4] * - 1. * (T3[3] + T3[15]) - F2[5]
      * (T3[7] + cI * (T3[11]))) + (F1[3] * (F2[4] * (+cI * (T3[11]) - T3[7]) +
      F2[5] * (T3[15] - T3[3])) + (F1[4] * (F2[2] * (T3[15] - T3[3]) + F2[3] *
      (T3[7] + cI * (T3[11]))) + F1[5] * (F2[2] * (T3[7] - cI * (T3[11])) -
      F2[3] * (T3[3] + T3[15]))))) + (P2[2] * (F1[2] * (F2[4] * - 1. * (T3[4] +
      T3[16]) - F2[5] * (T3[8] + cI * (T3[12]))) + (F1[3] * (F2[4] * (+cI *
      (T3[12]) - T3[8]) + F2[5] * (T3[16] - T3[4])) + (F1[4] * (F2[2] * (T3[16]
      - T3[4]) + F2[3] * (T3[8] + cI * (T3[12]))) + F1[5] * (F2[2] * (T3[8] -
      cI * (T3[12])) - F2[3] * (T3[4] + T3[16]))))) + P2[3] * (F1[2] * (F2[4] *
      - 1. * (T3[5] + T3[17]) - F2[5] * (T3[9] + cI * (T3[13]))) + (F1[3] *
      (F2[4] * (+cI * (T3[13]) - T3[9]) + F2[5] * (T3[17] - T3[5])) + (F1[4] *
      (F2[2] * (T3[17] - T3[5]) + F2[3] * (T3[9] + cI * (T3[13]))) + F1[5] *
      (F2[2] * (T3[9] - cI * (T3[13])) - F2[3] * (T3[5] + T3[17]))))))));
  TMP4 = (P1[0] * (F1[2] * (F2[4] * (T3[2] + T3[14]) + F2[5] * (T3[6] + cI *
      (T3[10]))) + (F1[3] * (F2[4] * (T3[6] - cI * (T3[10])) + F2[5] * (T3[2] -
      T3[14])) + (F1[4] * (F2[2] * (T3[2] - T3[14]) - F2[3] * (T3[6] + cI *
      (T3[10]))) + F1[5] * (F2[2] * (+cI * (T3[10]) - T3[6]) + F2[3] * (T3[2] +
      T3[14]))))) + (P1[1] * (F1[2] * (F2[4] * - 1. * (T3[3] + T3[15]) - F2[5]
      * (T3[7] + cI * (T3[11]))) + (F1[3] * (F2[4] * (+cI * (T3[11]) - T3[7]) +
      F2[5] * (T3[15] - T3[3])) + (F1[4] * (F2[2] * (T3[15] - T3[3]) + F2[3] *
      (T3[7] + cI * (T3[11]))) + F1[5] * (F2[2] * (T3[7] - cI * (T3[11])) -
      F2[3] * (T3[3] + T3[15]))))) + (P1[2] * (F1[2] * (F2[4] * - 1. * (T3[4] +
      T3[16]) - F2[5] * (T3[8] + cI * (T3[12]))) + (F1[3] * (F2[4] * (+cI *
      (T3[12]) - T3[8]) + F2[5] * (T3[16] - T3[4])) + (F1[4] * (F2[2] * (T3[16]
      - T3[4]) + F2[3] * (T3[8] + cI * (T3[12]))) + F1[5] * (F2[2] * (T3[8] -
      cI * (T3[12])) - F2[3] * (T3[4] + T3[16]))))) + P1[3] * (F1[2] * (F2[4] *
      - 1. * (T3[5] + T3[17]) - F2[5] * (T3[9] + cI * (T3[13]))) + (F1[3] *
      (F2[4] * (+cI * (T3[13]) - T3[9]) + F2[5] * (T3[17] - T3[5])) + (F1[4] *
      (F2[2] * (T3[17] - T3[5]) + F2[3] * (T3[9] + cI * (T3[13]))) + F1[5] *
      (F2[2] * (T3[9] - cI * (T3[13])) - F2[3] * (T3[5] + T3[17]))))))));
  TMP7 = (P2[0] * (F1[2] * (F2[4] * (T3[2] + T3[5]) + F2[5] * (T3[3] + cI *
      (T3[4]))) + (F1[3] * (F2[4] * (T3[3] - cI * (T3[4])) + F2[5] * (T3[2] -
      T3[5])) + (F1[4] * (F2[2] * (T3[2] - T3[5]) - F2[3] * (T3[3] + cI *
      (T3[4]))) + F1[5] * (F2[2] * (+cI * (T3[4]) - T3[3]) + F2[3] * (T3[2] +
      T3[5]))))) + (P2[1] * (F1[2] * (F2[4] * - 1. * (T3[6] + T3[9]) - F2[5] *
      (T3[7] + cI * (T3[8]))) + (F1[3] * (F2[4] * (+cI * (T3[8]) - T3[7]) +
      F2[5] * (T3[9] - T3[6])) + (F1[4] * (F2[2] * (T3[9] - T3[6]) + F2[3] *
      (T3[7] + cI * (T3[8]))) + F1[5] * (F2[2] * (T3[7] - cI * (T3[8])) - F2[3]
      * (T3[6] + T3[9]))))) + (P2[2] * (F1[2] * (F2[4] * - 1. * (T3[10] +
      T3[13]) - F2[5] * (T3[11] + cI * (T3[12]))) + (F1[3] * (F2[4] * (+cI *
      (T3[12]) - T3[11]) + F2[5] * (T3[13] - T3[10])) + (F1[4] * (F2[2] *
      (T3[13] - T3[10]) + F2[3] * (T3[11] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (T3[11] - cI * (T3[12])) - F2[3] * (T3[10] + T3[13]))))) + P2[3] * (F1[2]
      * (F2[4] * - 1. * (T3[14] + T3[17]) - F2[5] * (T3[15] + cI * (T3[16]))) +
      (F1[3] * (F2[4] * (+cI * (T3[16]) - T3[15]) + F2[5] * (T3[17] - T3[14]))
      + (F1[4] * (F2[2] * (T3[17] - T3[14]) + F2[3] * (T3[15] + cI * (T3[16])))
      + F1[5] * (F2[2] * (T3[15] - cI * (T3[16])) - F2[3] * (T3[14] +
      T3[17]))))))));
  TMP6 = (P1[0] * (F1[2] * (F2[4] * (T3[2] + T3[5]) + F2[5] * (T3[3] + cI *
      (T3[4]))) + (F1[3] * (F2[4] * (T3[3] - cI * (T3[4])) + F2[5] * (T3[2] -
      T3[5])) + (F1[4] * (F2[2] * (T3[2] - T3[5]) - F2[3] * (T3[3] + cI *
      (T3[4]))) + F1[5] * (F2[2] * (+cI * (T3[4]) - T3[3]) + F2[3] * (T3[2] +
      T3[5]))))) + (P1[1] * (F1[2] * (F2[4] * - 1. * (T3[6] + T3[9]) - F2[5] *
      (T3[7] + cI * (T3[8]))) + (F1[3] * (F2[4] * (+cI * (T3[8]) - T3[7]) +
      F2[5] * (T3[9] - T3[6])) + (F1[4] * (F2[2] * (T3[9] - T3[6]) + F2[3] *
      (T3[7] + cI * (T3[8]))) + F1[5] * (F2[2] * (T3[7] - cI * (T3[8])) - F2[3]
      * (T3[6] + T3[9]))))) + (P1[2] * (F1[2] * (F2[4] * - 1. * (T3[10] +
      T3[13]) - F2[5] * (T3[11] + cI * (T3[12]))) + (F1[3] * (F2[4] * (+cI *
      (T3[12]) - T3[11]) + F2[5] * (T3[13] - T3[10])) + (F1[4] * (F2[2] *
      (T3[13] - T3[10]) + F2[3] * (T3[11] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (T3[11] - cI * (T3[12])) - F2[3] * (T3[10] + T3[13]))))) + P1[3] * (F1[2]
      * (F2[4] * - 1. * (T3[14] + T3[17]) - F2[5] * (T3[15] + cI * (T3[16]))) +
      (F1[3] * (F2[4] * (+cI * (T3[16]) - T3[15]) + F2[5] * (T3[17] - T3[14]))
      + (F1[4] * (F2[2] * (T3[17] - T3[14]) + F2[3] * (T3[15] + cI * (T3[16])))
      + F1[5] * (F2[2] * (T3[15] - cI * (T3[16])) - F2[3] * (T3[14] +
      T3[17]))))))));
  vertex = COUP * (-cI * (TMP4 + TMP6) + cI * (TMP5 + TMP7)); 
}


void VVT13_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP15; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP31; 
  double P3[4]; 
  complex<double> TMP16; 
  complex<double> TMP33; 
  complex<double> TMP17; 
  complex<double> TMP14; 
  complex<double> TMP32; 
  complex<double> TMP34; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP33 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP19 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP18 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP14 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP15 = (P2[0] * - 1. * (T3[6] * V1[3] + T3[10] * V1[4] + T3[14] * V1[5] -
      T3[2] * V1[2]) + (P2[1] * (T3[7] * V1[3] + T3[11] * V1[4] + T3[15] *
      V1[5] - T3[3] * V1[2]) + (P2[2] * (T3[8] * V1[3] + T3[12] * V1[4] +
      T3[16] * V1[5] - T3[4] * V1[2]) + P2[3] * (T3[9] * V1[3] + T3[13] * V1[4]
      + T3[17] * V1[5] - T3[5] * V1[2]))));
  TMP32 = (V1[2] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (V1[3] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (V1[4] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + V1[5] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP31 = (V1[2] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (V1[3] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (V1[4] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + V1[5] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP16 = (P2[0] * - 1. * (T3[3] * V1[3] + T3[4] * V1[4] + T3[5] * V1[5] -
      T3[2] * V1[2]) + (P2[1] * (T3[7] * V1[3] + T3[8] * V1[4] + T3[9] * V1[5]
      - T3[6] * V1[2]) + (P2[2] * (T3[11] * V1[3] + T3[12] * V1[4] + T3[13] *
      V1[5] - T3[10] * V1[2]) + P2[3] * (T3[15] * V1[3] + T3[16] * V1[4] +
      T3[17] * V1[5] - T3[14] * V1[2]))));
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP17 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP34 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  vertex = COUP * (TMP14 * (TMP17 * - 1. * (+cI * (TMP10 + TMP11)) + TMP33 *
      (+cI * (TMP15 + TMP16))) + TMP34 * (TMP17 * (+cI * (TMP18 + TMP19)) -
      TMP33 * (+cI * (TMP31 + TMP32))));
}


void VVT8_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP9; 
  complex<double> TMP13; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP8 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP13 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP12 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  vertex = COUP * (TMP12 * TMP13 * (+cI * (TMP10 + TMP11)) - TMP8 * TMP9 * (+cI
      * (TMP10 + TMP11)));
}


void VVT3_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP39; 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP40; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP39 = -1. * (P1[0] * (P3[0] * (T3[14] * (V2[4] * V1[3] - V2[3] * V1[4]) +
      (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) +
      (T3[14] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[5] * V1[2] -
      V2[2] * V1[5]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) +
      (T3[14] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[6] * (V2[2] * V1[5] -
      V2[5] * V1[2]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      (T3[6] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[2] * V1[3] -
      V2[3] * V1[2])))))) + (P1[1] * (P3[0] * (T3[11] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (T3[15] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[7] * (V2[4] *
      V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[11] * (V2[2] * V1[5] - V2[5] *
      V1[2]) + (T3[15] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[3] * (V2[5] *
      V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[15] * (V2[2] * V1[3] - V2[3] *
      V1[2]) + (T3[3] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[7] * (V2[5] *
      V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[11] * (V2[3] * V1[2] - V2[2] *
      V1[3]) + (T3[3] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[7] * (V2[2] *
      V1[4] - V2[4] * V1[2])))))) + (P1[2] * (P3[0] * (T3[12] * (V2[5] * V1[3]
      - V2[3] * V1[5]) + (T3[16] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[8] *
      (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] -
      V2[5] * V1[2]) + (T3[16] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] *
      (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[16] * (V2[2] * V1[3] -
      V2[3] * V1[2]) + (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[8] *
      (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[12] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[4] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[8] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))) + P1[3] * (P3[0] * (T3[13] * (V2[5]
      * V1[3] - V2[3] * V1[5]) + (T3[17] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      T3[9] * (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[13] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + (T3[17] * (V2[4] * V1[2] - V2[2] * V1[4]) +
      T3[5] * (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]) + (T3[5] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[9]
      * (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[13] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[5] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[9] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))))));
  TMP38 = -1. * (P2[0] * (P3[0] * (T3[3] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[5] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) + (T3[4] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[5] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) + (T3[3] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[5] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] * (V2[2] * V1[3] - V2[3] *
      V1[2])))))) + (P2[1] * (P3[0] * (T3[7] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[8] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[9] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[8] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[9] * (V2[4] * V1[2] - V2[2]
      * V1[4]))) + (P3[2] * (T3[6] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[7] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[9] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P3[3] * (T3[6] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[7] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[8] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P2[2] * (P3[0] * (T3[11] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[12] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[13] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] - V2[5] * V1[2]) +
      (T3[13] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[2] * (T3[11] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[13] * (V2[2] * V1[3] - V2[3] * V1[2]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + P3[3] * (T3[11] * (V2[2] * V1[4] - V2[4] * V1[2]) +
      (T3[12] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[10] * (V2[4] * V1[3] -
      V2[3] * V1[4])))))) + P2[3] * (P3[0] * (T3[15] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (T3[16] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[17] * (V2[3] *
      V1[4] - V2[4] * V1[3]))) + (P3[1] * (T3[14] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (T3[16] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[17] * (V2[4] *
      V1[2] - V2[2] * V1[4]))) + (P3[2] * (T3[14] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (T3[15] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + P3[3] * (T3[14] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (T3[15] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[16] * (V2[3] *
      V1[2] - V2[2] * V1[3])))))))));
  TMP37 = -1. * (P1[0] * (P3[0] * (T3[3] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[5] * (V2[4] * V1[3] - V2[3]
      * V1[4]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) + (T3[4] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[5] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) + (T3[3] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + T3[5] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) + (T3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] * (V2[2] * V1[3] - V2[3] *
      V1[2])))))) + (P1[1] * (P3[0] * (T3[7] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[8] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[9] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) +
      (T3[8] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[9] * (V2[4] * V1[2] - V2[2]
      * V1[4]))) + (P3[2] * (T3[6] * (V2[3] * V1[5] - V2[5] * V1[3]) + (T3[7] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + T3[9] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P3[3] * (T3[6] * (V2[4] * V1[3] - V2[3] * V1[4]) + (T3[7] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + T3[8] * (V2[3] * V1[2] - V2[2] *
      V1[3])))))) + (P1[2] * (P3[0] * (T3[11] * (V2[4] * V1[5] - V2[5] * V1[4])
      + (T3[12] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[13] * (V2[3] * V1[4] -
      V2[4] * V1[3]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] - V2[5] * V1[2]) +
      (T3[13] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[5] * V1[4] -
      V2[4] * V1[5]))) + (P3[2] * (T3[11] * (V2[5] * V1[2] - V2[2] * V1[5]) +
      (T3[13] * (V2[2] * V1[3] - V2[3] * V1[2]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + P3[3] * (T3[11] * (V2[2] * V1[4] - V2[4] * V1[2]) +
      (T3[12] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[10] * (V2[4] * V1[3] -
      V2[3] * V1[4])))))) + P1[3] * (P3[0] * (T3[15] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (T3[16] * (V2[5] * V1[3] - V2[3] * V1[5]) + T3[17] * (V2[3] *
      V1[4] - V2[4] * V1[3]))) + (P3[1] * (T3[14] * (V2[5] * V1[4] - V2[4] *
      V1[5]) + (T3[16] * (V2[2] * V1[5] - V2[5] * V1[2]) + T3[17] * (V2[4] *
      V1[2] - V2[2] * V1[4]))) + (P3[2] * (T3[14] * (V2[3] * V1[5] - V2[5] *
      V1[3]) + (T3[15] * (V2[5] * V1[2] - V2[2] * V1[5]) + T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + P3[3] * (T3[14] * (V2[4] * V1[3] - V2[3] *
      V1[4]) + (T3[15] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[16] * (V2[3] *
      V1[2] - V2[2] * V1[3])))))))));
  TMP40 = -1. * (P2[0] * (P3[0] * (T3[14] * (V2[4] * V1[3] - V2[3] * V1[4]) +
      (T3[6] * (V2[5] * V1[4] - V2[4] * V1[5]) + T3[10] * (V2[3] * V1[5] -
      V2[5] * V1[3]))) + (P3[1] * (T3[2] * (V2[4] * V1[5] - V2[5] * V1[4]) +
      (T3[14] * (V2[2] * V1[4] - V2[4] * V1[2]) + T3[10] * (V2[5] * V1[2] -
      V2[2] * V1[5]))) + (P3[2] * (T3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) +
      (T3[14] * (V2[3] * V1[2] - V2[2] * V1[3]) + T3[6] * (V2[2] * V1[5] -
      V2[5] * V1[2]))) + P3[3] * (T3[2] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      (T3[6] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[10] * (V2[2] * V1[3] -
      V2[3] * V1[2])))))) + (P2[1] * (P3[0] * (T3[11] * (V2[5] * V1[3] - V2[3]
      * V1[5]) + (T3[15] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[7] * (V2[4] *
      V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[11] * (V2[2] * V1[5] - V2[5] *
      V1[2]) + (T3[15] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[3] * (V2[5] *
      V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[15] * (V2[2] * V1[3] - V2[3] *
      V1[2]) + (T3[3] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[7] * (V2[5] *
      V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[11] * (V2[3] * V1[2] - V2[2] *
      V1[3]) + (T3[3] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[7] * (V2[2] *
      V1[4] - V2[4] * V1[2])))))) + (P2[2] * (P3[0] * (T3[12] * (V2[5] * V1[3]
      - V2[3] * V1[5]) + (T3[16] * (V2[3] * V1[4] - V2[4] * V1[3]) + T3[8] *
      (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[12] * (V2[2] * V1[5] -
      V2[5] * V1[2]) + (T3[16] * (V2[4] * V1[2] - V2[2] * V1[4]) + T3[4] *
      (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[16] * (V2[2] * V1[3] -
      V2[3] * V1[2]) + (T3[4] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[8] *
      (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[12] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[4] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[8] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))) + P2[3] * (P3[0] * (T3[13] * (V2[5]
      * V1[3] - V2[3] * V1[5]) + (T3[17] * (V2[3] * V1[4] - V2[4] * V1[3]) +
      T3[9] * (V2[4] * V1[5] - V2[5] * V1[4]))) + (P3[1] * (T3[13] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + (T3[17] * (V2[4] * V1[2] - V2[2] * V1[4]) +
      T3[5] * (V2[5] * V1[4] - V2[4] * V1[5]))) + (P3[2] * (T3[17] * (V2[2] *
      V1[3] - V2[3] * V1[2]) + (T3[5] * (V2[3] * V1[5] - V2[5] * V1[3]) + T3[9]
      * (V2[5] * V1[2] - V2[2] * V1[5]))) + P3[3] * (T3[13] * (V2[3] * V1[2] -
      V2[2] * V1[3]) + (T3[5] * (V2[4] * V1[3] - V2[3] * V1[4]) + T3[9] *
      (V2[2] * V1[4] - V2[4] * V1[2])))))))));
  vertex = COUP * (-cI * (TMP37 + TMP39) + cI * (TMP38 + TMP40)); 
}


void FFV7_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
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
  F1[2] = denom * 2. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3]
      - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] -
      V3[2])))) + (+1./2. * (M1 * (F2[5] * (V3[3] + cI * (V3[4])) + 2. * (F2[4]
      * 1./2. * (V3[2] + V3[5])))) + F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) +
      (P1[1] * - 1. * (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * 2. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (+1./2. * (M1 * (F2[5] * (V3[2] - V3[5]) + 2. *
      (F2[4] * 1./2. * (V3[3] - cI * (V3[4]))))) + F2[3] * (P1[0] * - 1. *
      (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI
      * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 2. * (V3[5] - V3[2]) + 2. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * 2. * (+cI * (V3[4]) - V3[3]) + 2. * (F2[3] *
      (V3[2] + V3[5])))));
}


void FFT5_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP47; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP47 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  vertex = COUP * - TMP47 * (+cI * (TMP10 + TMP11)); 
}


void FFV5_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
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
  F2[4] = denom * - cI * M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]));
  F2[5] = denom * cI * M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] -
      V3[5]));
}

void FFV5_7_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV5_2(F1, V3, COUP1, M2, W2, F2); 
  FFV7_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV7_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  complex<double> TMP20; 
  complex<double> TMP21; 
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
  TMP20 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  TMP21 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 2. * cI * (OM3 * - 1./2. * P3[0] * (TMP20 + 2. * (TMP21)) +
      (+1./2. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * - 2. * cI * (OM3 * - 1./2. * P3[1] * (TMP20 + 2. * (TMP21)) +
      (-1./2. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * 2. * cI * (OM3 * 1./2. * P3[2] * (TMP20 + 2. * (TMP21)) +
      (+1./2. * cI * (F2[5] * F1[2]) - 1./2. * cI * (F2[4] * F1[3]) - cI *
      (F2[3] * F1[4]) + cI * (F2[2] * F1[5])));
  V3[5] = denom * 2. * cI * (OM3 * 1./2. * P3[3] * (TMP20 + 2. * (TMP21)) +
      (+1./2. * (F2[4] * F1[2]) - 1./2. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void FFV2_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
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
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * - 1. * (V3[2] +
      V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
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
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}


void FFT3_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP3; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP1 = (P2[0] * (F1[2] * (F2[4] * - 1. * (T3[2] + T3[14]) - F2[5] * (T3[6] +
      cI * (T3[10]))) + (F1[3] * (F2[4] * (+cI * (T3[10]) - T3[6]) + F2[5] *
      (T3[14] - T3[2])) + (F1[4] * (F2[2] * (T3[2] - T3[14]) - F2[3] * (T3[6] +
      cI * (T3[10]))) + F1[5] * (F2[2] * (+cI * (T3[10]) - T3[6]) + F2[3] *
      (T3[2] + T3[14]))))) + (P2[1] * (F1[2] * (F2[4] * (T3[3] + T3[15]) +
      F2[5] * (T3[7] + cI * (T3[11]))) + (F1[3] * (F2[4] * (T3[7] - cI *
      (T3[11])) + F2[5] * (T3[3] - T3[15])) + (F1[4] * (F2[2] * (T3[15] -
      T3[3]) + F2[3] * (T3[7] + cI * (T3[11]))) + F1[5] * (F2[2] * (T3[7] - cI
      * (T3[11])) - F2[3] * (T3[3] + T3[15]))))) + (P2[2] * (F1[2] * (F2[4] *
      (T3[4] + T3[16]) + F2[5] * (T3[8] + cI * (T3[12]))) + (F1[3] * (F2[4] *
      (T3[8] - cI * (T3[12])) + F2[5] * (T3[4] - T3[16])) + (F1[4] * (F2[2] *
      (T3[16] - T3[4]) + F2[3] * (T3[8] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (T3[8] - cI * (T3[12])) - F2[3] * (T3[4] + T3[16]))))) + P2[3] * (F1[2] *
      (F2[4] * (T3[5] + T3[17]) + F2[5] * (T3[9] + cI * (T3[13]))) + (F1[3] *
      (F2[4] * (T3[9] - cI * (T3[13])) + F2[5] * (T3[5] - T3[17])) + (F1[4] *
      (F2[2] * (T3[17] - T3[5]) + F2[3] * (T3[9] + cI * (T3[13]))) + F1[5] *
      (F2[2] * (T3[9] - cI * (T3[13])) - F2[3] * (T3[5] + T3[17]))))))));
  TMP0 = (P1[0] * (F1[2] * (F2[4] * - 1. * (T3[2] + T3[14]) - F2[5] * (T3[6] +
      cI * (T3[10]))) + (F1[3] * (F2[4] * (+cI * (T3[10]) - T3[6]) + F2[5] *
      (T3[14] - T3[2])) + (F1[4] * (F2[2] * (T3[2] - T3[14]) - F2[3] * (T3[6] +
      cI * (T3[10]))) + F1[5] * (F2[2] * (+cI * (T3[10]) - T3[6]) + F2[3] *
      (T3[2] + T3[14]))))) + (P1[1] * (F1[2] * (F2[4] * (T3[3] + T3[15]) +
      F2[5] * (T3[7] + cI * (T3[11]))) + (F1[3] * (F2[4] * (T3[7] - cI *
      (T3[11])) + F2[5] * (T3[3] - T3[15])) + (F1[4] * (F2[2] * (T3[15] -
      T3[3]) + F2[3] * (T3[7] + cI * (T3[11]))) + F1[5] * (F2[2] * (T3[7] - cI
      * (T3[11])) - F2[3] * (T3[3] + T3[15]))))) + (P1[2] * (F1[2] * (F2[4] *
      (T3[4] + T3[16]) + F2[5] * (T3[8] + cI * (T3[12]))) + (F1[3] * (F2[4] *
      (T3[8] - cI * (T3[12])) + F2[5] * (T3[4] - T3[16])) + (F1[4] * (F2[2] *
      (T3[16] - T3[4]) + F2[3] * (T3[8] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (T3[8] - cI * (T3[12])) - F2[3] * (T3[4] + T3[16]))))) + P1[3] * (F1[2] *
      (F2[4] * (T3[5] + T3[17]) + F2[5] * (T3[9] + cI * (T3[13]))) + (F1[3] *
      (F2[4] * (T3[9] - cI * (T3[13])) + F2[5] * (T3[5] - T3[17])) + (F1[4] *
      (F2[2] * (T3[17] - T3[5]) + F2[3] * (T3[9] + cI * (T3[13]))) + F1[5] *
      (F2[2] * (T3[9] - cI * (T3[13])) - F2[3] * (T3[5] + T3[17]))))))));
  TMP3 = (P2[0] * (F1[2] * (F2[4] * - 1. * (T3[2] + T3[5]) - F2[5] * (T3[3] +
      cI * (T3[4]))) + (F1[3] * (F2[4] * (+cI * (T3[4]) - T3[3]) + F2[5] *
      (T3[5] - T3[2])) + (F1[4] * (F2[2] * (T3[2] - T3[5]) - F2[3] * (T3[3] +
      cI * (T3[4]))) + F1[5] * (F2[2] * (+cI * (T3[4]) - T3[3]) + F2[3] *
      (T3[2] + T3[5]))))) + (P2[1] * (F1[2] * (F2[4] * (T3[6] + T3[9]) + F2[5]
      * (T3[7] + cI * (T3[8]))) + (F1[3] * (F2[4] * (T3[7] - cI * (T3[8])) +
      F2[5] * (T3[6] - T3[9])) + (F1[4] * (F2[2] * (T3[9] - T3[6]) + F2[3] *
      (T3[7] + cI * (T3[8]))) + F1[5] * (F2[2] * (T3[7] - cI * (T3[8])) - F2[3]
      * (T3[6] + T3[9]))))) + (P2[2] * (F1[2] * (F2[4] * (T3[10] + T3[13]) +
      F2[5] * (T3[11] + cI * (T3[12]))) + (F1[3] * (F2[4] * (T3[11] - cI *
      (T3[12])) + F2[5] * (T3[10] - T3[13])) + (F1[4] * (F2[2] * (T3[13] -
      T3[10]) + F2[3] * (T3[11] + cI * (T3[12]))) + F1[5] * (F2[2] * (T3[11] -
      cI * (T3[12])) - F2[3] * (T3[10] + T3[13]))))) + P2[3] * (F1[2] * (F2[4]
      * (T3[14] + T3[17]) + F2[5] * (T3[15] + cI * (T3[16]))) + (F1[3] * (F2[4]
      * (T3[15] - cI * (T3[16])) + F2[5] * (T3[14] - T3[17])) + (F1[4] * (F2[2]
      * (T3[17] - T3[14]) + F2[3] * (T3[15] + cI * (T3[16]))) + F1[5] * (F2[2]
      * (T3[15] - cI * (T3[16])) - F2[3] * (T3[14] + T3[17]))))))));
  TMP2 = (P1[0] * (F1[2] * (F2[4] * - 1. * (T3[2] + T3[5]) - F2[5] * (T3[3] +
      cI * (T3[4]))) + (F1[3] * (F2[4] * (+cI * (T3[4]) - T3[3]) + F2[5] *
      (T3[5] - T3[2])) + (F1[4] * (F2[2] * (T3[2] - T3[5]) - F2[3] * (T3[3] +
      cI * (T3[4]))) + F1[5] * (F2[2] * (+cI * (T3[4]) - T3[3]) + F2[3] *
      (T3[2] + T3[5]))))) + (P1[1] * (F1[2] * (F2[4] * (T3[6] + T3[9]) + F2[5]
      * (T3[7] + cI * (T3[8]))) + (F1[3] * (F2[4] * (T3[7] - cI * (T3[8])) +
      F2[5] * (T3[6] - T3[9])) + (F1[4] * (F2[2] * (T3[9] - T3[6]) + F2[3] *
      (T3[7] + cI * (T3[8]))) + F1[5] * (F2[2] * (T3[7] - cI * (T3[8])) - F2[3]
      * (T3[6] + T3[9]))))) + (P1[2] * (F1[2] * (F2[4] * (T3[10] + T3[13]) +
      F2[5] * (T3[11] + cI * (T3[12]))) + (F1[3] * (F2[4] * (T3[11] - cI *
      (T3[12])) + F2[5] * (T3[10] - T3[13])) + (F1[4] * (F2[2] * (T3[13] -
      T3[10]) + F2[3] * (T3[11] + cI * (T3[12]))) + F1[5] * (F2[2] * (T3[11] -
      cI * (T3[12])) - F2[3] * (T3[10] + T3[13]))))) + P1[3] * (F1[2] * (F2[4]
      * (T3[14] + T3[17]) + F2[5] * (T3[15] + cI * (T3[16]))) + (F1[3] * (F2[4]
      * (T3[15] - cI * (T3[16])) + F2[5] * (T3[14] - T3[17])) + (F1[4] * (F2[2]
      * (T3[17] - T3[14]) + F2[3] * (T3[15] + cI * (T3[16]))) + F1[5] * (F2[2]
      * (T3[15] - cI * (T3[16])) - F2[3] * (T3[14] + T3[17]))))))));
  vertex = COUP * (-cI * (TMP0 + TMP2) + cI * (TMP1 + TMP3)); 
}


void VVT10_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP10; 
  double P3[4]; 
  complex<double> TMP14; 
  complex<double> TMP43; 
  complex<double> TMP9; 
  complex<double> TMP11; 
  double P2[4]; 
  complex<double> TMP46; 
  complex<double> TMP33; 
  complex<double> TMP34; 
  complex<double> TMP42; 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP45; 
  complex<double> TMP41; 
  complex<double> TMP13; 
  complex<double> TMP44; 
  complex<double> TMP17; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP42 = (P2[0] * (P2[1] * - 1. * (T3[6] + T3[3]) + (P2[2] * - 1. * (T3[10] +
      T3[4]) + (P2[3] * - 1. * (T3[14] + T3[5]) + T3[2] * P2[0]))) + (P2[1] *
      (P2[2] * (T3[11] + T3[8]) + (P2[3] * (T3[15] + T3[9]) + T3[7] * P2[1])) +
      (P2[2] * (P2[3] * (T3[16] + T3[13]) + T3[12] * P2[2]) + T3[17] * P2[3] *
      P2[3])));
  TMP43 = (P1[0] * - 1. * (T3[6] * V1[3] + T3[10] * V1[4] + T3[14] * V1[5] -
      T3[2] * V1[2]) + (P1[1] * (T3[7] * V1[3] + T3[11] * V1[4] + T3[15] *
      V1[5] - T3[3] * V1[2]) + (P1[2] * (T3[8] * V1[3] + T3[12] * V1[4] +
      T3[16] * V1[5] - T3[4] * V1[2]) + P1[3] * (T3[9] * V1[3] + T3[13] * V1[4]
      + T3[17] * V1[5] - T3[5] * V1[2]))));
  TMP41 = (P1[0] * (P1[1] * - 1. * (T3[6] + T3[3]) + (P1[2] * - 1. * (T3[10] +
      T3[4]) + (P1[3] * - 1. * (T3[14] + T3[5]) + T3[2] * P1[0]))) + (P1[1] *
      (P1[2] * (T3[11] + T3[8]) + (P1[3] * (T3[15] + T3[9]) + T3[7] * P1[1])) +
      (P1[2] * (P1[3] * (T3[16] + T3[13]) + T3[12] * P1[2]) + T3[17] * P1[3] *
      P1[3])));
  TMP46 = (P2[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P2[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P2[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P2[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP44 = (P1[0] * - 1. * (T3[3] * V1[3] + T3[4] * V1[4] + T3[5] * V1[5] -
      T3[2] * V1[2]) + (P1[1] * (T3[7] * V1[3] + T3[8] * V1[4] + T3[9] * V1[5]
      - T3[6] * V1[2]) + (P1[2] * (T3[11] * V1[3] + T3[12] * V1[4] + T3[13] *
      V1[5] - T3[10] * V1[2]) + P1[3] * (T3[15] * V1[3] + T3[16] * V1[4] +
      T3[17] * V1[5] - T3[14] * V1[2]))));
  TMP45 = (P2[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P2[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P2[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P2[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP8 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP33 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP14 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP17 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP12 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP13 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP34 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  vertex = COUP * (TMP13 * (TMP12 * (-cI * (TMP10 + TMP11) + cI * (TMP41 +
      TMP42)) + (+1./2. * (TMP17 * (+cI * (TMP45 + TMP46))) + TMP14 * 1./2. *
      (+cI * (TMP43 + TMP44)))) + (TMP8 * (TMP9 * (-cI * (TMP41 + TMP42) + cI *
      (TMP10 + TMP11)) + (TMP34 * - 1./2. * (+cI * (TMP43 + TMP44)) - cI *
      (TMP17 * TMP42))) + (TMP33 * (TMP9 * - 1./2. * (+cI * (TMP45 + TMP46)) +
      cI * (TMP12 * TMP42)) + TMP41 * (-cI * (TMP9 * TMP14) + cI * (TMP12 *
      TMP34)))));
}

void VVT10_11_12_13_2_3_6_7_8_9_0(complex<double> V1[], complex<double> V2[],
    complex<double> T3[], complex<double> COUP1, complex<double> COUP2,
    complex<double> COUP3, complex<double> COUP4, complex<double> COUP5,
    complex<double> COUP6, complex<double> COUP7, complex<double> COUP8,
    complex<double> COUP9, complex<double> COUP10, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
  complex<double> tmp; 
//   double P2[4]; 
//   double P1[4]; 
  VVT10_0(V1, V2, T3, COUP1, vertex); 
  VVT11_0(V1, V2, T3, COUP2, tmp); 
  vertex = vertex + tmp; 
  VVT12_0(V1, V2, T3, COUP3, tmp); 
  vertex = vertex + tmp; 
  VVT13_0(V1, V2, T3, COUP4, tmp); 
  vertex = vertex + tmp; 
  VVT2_0(V1, V2, T3, COUP5, tmp); 
  vertex = vertex + tmp; 
  VVT3_0(V1, V2, T3, COUP6, tmp); 
  vertex = vertex + tmp; 
  VVT6_0(V1, V2, T3, COUP7, tmp); 
  vertex = vertex + tmp; 
  VVT7_0(V1, V2, T3, COUP8, tmp); 
  vertex = vertex + tmp; 
  VVT8_0(V1, V2, T3, COUP9, tmp); 
  vertex = vertex + tmp; 
  VVT9_0(V1, V2, T3, COUP10, tmp); 
  vertex = vertex + tmp; 
}

void VVT9_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP17; 
  double P3[4]; 
  complex<double> TMP16; 
  complex<double> TMP15; 
  complex<double> TMP14; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP19 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP18 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP15 = (P2[0] * - 1. * (T3[6] * V1[3] + T3[10] * V1[4] + T3[14] * V1[5] -
      T3[2] * V1[2]) + (P2[1] * (T3[7] * V1[3] + T3[11] * V1[4] + T3[15] *
      V1[5] - T3[3] * V1[2]) + (P2[2] * (T3[8] * V1[3] + T3[12] * V1[4] +
      T3[16] * V1[5] - T3[4] * V1[2]) + P2[3] * (T3[9] * V1[3] + T3[13] * V1[4]
      + T3[17] * V1[5] - T3[5] * V1[2]))));
  TMP14 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP17 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP16 = (P2[0] * - 1. * (T3[3] * V1[3] + T3[4] * V1[4] + T3[5] * V1[5] -
      T3[2] * V1[2]) + (P2[1] * (T3[7] * V1[3] + T3[8] * V1[4] + T3[9] * V1[5]
      - T3[6] * V1[2]) + (P2[2] * (T3[11] * V1[3] + T3[12] * V1[4] + T3[13] *
      V1[5] - T3[10] * V1[2]) + P2[3] * (T3[15] * V1[3] + T3[16] * V1[4] +
      T3[17] * V1[5] - T3[14] * V1[2]))));
  vertex = COUP * - 1. * (TMP14 * (+cI * (TMP15 + TMP16)) + TMP17 * (+cI *
      (TMP18 + TMP19)));
}


void VVT6_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP17; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP26; 
  complex<double> TMP14; 
  complex<double> TMP28; 
  complex<double> TMP27; 
  complex<double> TMP29; 
  complex<double> TMP24; 
  complex<double> TMP25; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  TMP24 = -1. * (P1[0] * (P2[0] * (P3[1] * (T3[5] * V1[4] - T3[4] * V1[5]) +
      (P3[2] * (T3[3] * V1[5] - T3[5] * V1[3]) + P3[3] * (T3[4] * V1[3] - T3[3]
      * V1[4]))) + (P2[1] * (P3[0] * (T3[4] * V1[5] - T3[5] * V1[4]) + (P3[2] *
      (T3[5] * V1[2] - T3[2] * V1[5]) + P3[3] * (T3[2] * V1[4] - T3[4] *
      V1[2]))) + (P2[2] * (P3[0] * (T3[5] * V1[3] - T3[3] * V1[5]) + (P3[1] *
      (T3[2] * V1[5] - T3[5] * V1[2]) + P3[3] * (T3[3] * V1[2] - T3[2] *
      V1[3]))) + P2[3] * (P3[0] * (T3[3] * V1[4] - T3[4] * V1[3]) + (P3[1] *
      (T3[4] * V1[2] - T3[2] * V1[4]) + P3[2] * (T3[2] * V1[3] - T3[3] *
      V1[2])))))) + (P1[1] * (P2[0] * (P3[1] * (T3[8] * V1[5] - T3[9] * V1[4])
      + (P3[2] * (T3[9] * V1[3] - T3[7] * V1[5]) + P3[3] * (T3[7] * V1[4] -
      T3[8] * V1[3]))) + (P2[1] * (P3[0] * (T3[9] * V1[4] - T3[8] * V1[5]) +
      (P3[2] * (T3[6] * V1[5] - T3[9] * V1[2]) + P3[3] * (T3[8] * V1[2] - T3[6]
      * V1[4]))) + (P2[2] * (P3[0] * (T3[7] * V1[5] - T3[9] * V1[3]) + (P3[1] *
      (T3[9] * V1[2] - T3[6] * V1[5]) + P3[3] * (T3[6] * V1[3] - T3[7] *
      V1[2]))) + P2[3] * (P3[0] * (T3[8] * V1[3] - T3[7] * V1[4]) + (P3[1] *
      (T3[6] * V1[4] - T3[8] * V1[2]) + P3[2] * (T3[7] * V1[2] - T3[6] *
      V1[3])))))) + (P1[2] * (P2[0] * (P3[1] * (T3[12] * V1[5] - T3[13] *
      V1[4]) + (P3[2] * (T3[13] * V1[3] - T3[11] * V1[5]) + P3[3] * (T3[11] *
      V1[4] - T3[12] * V1[3]))) + (P2[1] * (P3[0] * (T3[13] * V1[4] - T3[12] *
      V1[5]) + (P3[2] * (T3[10] * V1[5] - T3[13] * V1[2]) + P3[3] * (T3[12] *
      V1[2] - T3[10] * V1[4]))) + (P2[2] * (P3[0] * (T3[11] * V1[5] - T3[13] *
      V1[3]) + (P3[1] * (T3[13] * V1[2] - T3[10] * V1[5]) + P3[3] * (T3[10] *
      V1[3] - T3[11] * V1[2]))) + P2[3] * (P3[0] * (T3[12] * V1[3] - T3[11] *
      V1[4]) + (P3[1] * (T3[10] * V1[4] - T3[12] * V1[2]) + P3[2] * (T3[11] *
      V1[2] - T3[10] * V1[3])))))) + P1[3] * (P2[0] * (P3[1] * (T3[16] * V1[5]
      - T3[17] * V1[4]) + (P3[2] * (T3[17] * V1[3] - T3[15] * V1[5]) + P3[3] *
      (T3[15] * V1[4] - T3[16] * V1[3]))) + (P2[1] * (P3[0] * (T3[17] * V1[4] -
      T3[16] * V1[5]) + (P3[2] * (T3[14] * V1[5] - T3[17] * V1[2]) + P3[3] *
      (T3[16] * V1[2] - T3[14] * V1[4]))) + (P2[2] * (P3[0] * (T3[15] * V1[5] -
      T3[17] * V1[3]) + (P3[1] * (T3[17] * V1[2] - T3[14] * V1[5]) + P3[3] *
      (T3[14] * V1[3] - T3[15] * V1[2]))) + P2[3] * (P3[0] * (T3[16] * V1[3] -
      T3[15] * V1[4]) + (P3[1] * (T3[14] * V1[4] - T3[16] * V1[2]) + P3[2] *
      (T3[15] * V1[2] - T3[14] * V1[3])))))))));
  TMP25 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[5] * V2[4] - T3[4] * V2[5]) +
      (P3[2] * (T3[3] * V2[5] - T3[5] * V2[3]) + P3[3] * (T3[4] * V2[3] - T3[3]
      * V2[4]))) + (P1[1] * (P3[0] * (T3[4] * V2[5] - T3[5] * V2[4]) + (P3[2] *
      (T3[5] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] - T3[4] *
      V2[2]))) + (P1[2] * (P3[0] * (T3[5] * V2[3] - T3[3] * V2[5]) + (P3[1] *
      (T3[2] * V2[5] - T3[5] * V2[2]) + P3[3] * (T3[3] * V2[2] - T3[2] *
      V2[3]))) + P1[3] * (P3[0] * (T3[3] * V2[4] - T3[4] * V2[3]) + (P3[1] *
      (T3[4] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] - T3[3] *
      V2[2])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[8] * V2[5] - T3[9] * V2[4])
      + (P3[2] * (T3[9] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] * V2[4] -
      T3[8] * V2[3]))) + (P1[1] * (P3[0] * (T3[9] * V2[4] - T3[8] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[9] * V2[2]) + P3[3] * (T3[8] * V2[2] - T3[6]
      * V2[4]))) + (P1[2] * (P3[0] * (T3[7] * V2[5] - T3[9] * V2[3]) + (P3[1] *
      (T3[9] * V2[2] - T3[6] * V2[5]) + P3[3] * (T3[6] * V2[3] - T3[7] *
      V2[2]))) + P1[3] * (P3[0] * (T3[8] * V2[3] - T3[7] * V2[4]) + (P3[1] *
      (T3[6] * V2[4] - T3[8] * V2[2]) + P3[2] * (T3[7] * V2[2] - T3[6] *
      V2[3])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[12] * V2[5] - T3[13] *
      V2[4]) + (P3[2] * (T3[13] * V2[3] - T3[11] * V2[5]) + P3[3] * (T3[11] *
      V2[4] - T3[12] * V2[3]))) + (P1[1] * (P3[0] * (T3[13] * V2[4] - T3[12] *
      V2[5]) + (P3[2] * (T3[10] * V2[5] - T3[13] * V2[2]) + P3[3] * (T3[12] *
      V2[2] - T3[10] * V2[4]))) + (P1[2] * (P3[0] * (T3[11] * V2[5] - T3[13] *
      V2[3]) + (P3[1] * (T3[13] * V2[2] - T3[10] * V2[5]) + P3[3] * (T3[10] *
      V2[3] - T3[11] * V2[2]))) + P1[3] * (P3[0] * (T3[12] * V2[3] - T3[11] *
      V2[4]) + (P3[1] * (T3[10] * V2[4] - T3[12] * V2[2]) + P3[2] * (T3[11] *
      V2[2] - T3[10] * V2[3])))))) + P2[3] * (P1[0] * (P3[1] * (T3[16] * V2[5]
      - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[15] * V2[5]) + P3[3] *
      (T3[15] * V2[4] - T3[16] * V2[3]))) + (P1[1] * (P3[0] * (T3[17] * V2[4] -
      T3[16] * V2[5]) + (P3[2] * (T3[14] * V2[5] - T3[17] * V2[2]) + P3[3] *
      (T3[16] * V2[2] - T3[14] * V2[4]))) + (P1[2] * (P3[0] * (T3[15] * V2[5] -
      T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[14] * V2[5]) + P3[3] *
      (T3[14] * V2[3] - T3[15] * V2[2]))) + P1[3] * (P3[0] * (T3[16] * V2[3] -
      T3[15] * V2[4]) + (P3[1] * (T3[14] * V2[4] - T3[16] * V2[2]) + P3[2] *
      (T3[15] * V2[2] - T3[14] * V2[3])))))))));
  TMP26 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[5] * V1[4] - T3[4] * V1[5]) +
      (P3[2] * (T3[3] * V1[5] - T3[5] * V1[3]) + P3[3] * (T3[4] * V1[3] - T3[3]
      * V1[4]))) + (P1[1] * (P3[0] * (T3[4] * V1[5] - T3[5] * V1[4]) + (P3[2] *
      (T3[5] * V1[2] - T3[2] * V1[5]) + P3[3] * (T3[2] * V1[4] - T3[4] *
      V1[2]))) + (P1[2] * (P3[0] * (T3[5] * V1[3] - T3[3] * V1[5]) + (P3[1] *
      (T3[2] * V1[5] - T3[5] * V1[2]) + P3[3] * (T3[3] * V1[2] - T3[2] *
      V1[3]))) + P1[3] * (P3[0] * (T3[3] * V1[4] - T3[4] * V1[3]) + (P3[1] *
      (T3[4] * V1[2] - T3[2] * V1[4]) + P3[2] * (T3[2] * V1[3] - T3[3] *
      V1[2])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[8] * V1[5] - T3[9] * V1[4])
      + (P3[2] * (T3[9] * V1[3] - T3[7] * V1[5]) + P3[3] * (T3[7] * V1[4] -
      T3[8] * V1[3]))) + (P1[1] * (P3[0] * (T3[9] * V1[4] - T3[8] * V1[5]) +
      (P3[2] * (T3[6] * V1[5] - T3[9] * V1[2]) + P3[3] * (T3[8] * V1[2] - T3[6]
      * V1[4]))) + (P1[2] * (P3[0] * (T3[7] * V1[5] - T3[9] * V1[3]) + (P3[1] *
      (T3[9] * V1[2] - T3[6] * V1[5]) + P3[3] * (T3[6] * V1[3] - T3[7] *
      V1[2]))) + P1[3] * (P3[0] * (T3[8] * V1[3] - T3[7] * V1[4]) + (P3[1] *
      (T3[6] * V1[4] - T3[8] * V1[2]) + P3[2] * (T3[7] * V1[2] - T3[6] *
      V1[3])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[12] * V1[5] - T3[13] *
      V1[4]) + (P3[2] * (T3[13] * V1[3] - T3[11] * V1[5]) + P3[3] * (T3[11] *
      V1[4] - T3[12] * V1[3]))) + (P1[1] * (P3[0] * (T3[13] * V1[4] - T3[12] *
      V1[5]) + (P3[2] * (T3[10] * V1[5] - T3[13] * V1[2]) + P3[3] * (T3[12] *
      V1[2] - T3[10] * V1[4]))) + (P1[2] * (P3[0] * (T3[11] * V1[5] - T3[13] *
      V1[3]) + (P3[1] * (T3[13] * V1[2] - T3[10] * V1[5]) + P3[3] * (T3[10] *
      V1[3] - T3[11] * V1[2]))) + P1[3] * (P3[0] * (T3[12] * V1[3] - T3[11] *
      V1[4]) + (P3[1] * (T3[10] * V1[4] - T3[12] * V1[2]) + P3[2] * (T3[11] *
      V1[2] - T3[10] * V1[3])))))) + P2[3] * (P1[0] * (P3[1] * (T3[16] * V1[5]
      - T3[17] * V1[4]) + (P3[2] * (T3[17] * V1[3] - T3[15] * V1[5]) + P3[3] *
      (T3[15] * V1[4] - T3[16] * V1[3]))) + (P1[1] * (P3[0] * (T3[17] * V1[4] -
      T3[16] * V1[5]) + (P3[2] * (T3[14] * V1[5] - T3[17] * V1[2]) + P3[3] *
      (T3[16] * V1[2] - T3[14] * V1[4]))) + (P1[2] * (P3[0] * (T3[15] * V1[5] -
      T3[17] * V1[3]) + (P3[1] * (T3[17] * V1[2] - T3[14] * V1[5]) + P3[3] *
      (T3[14] * V1[3] - T3[15] * V1[2]))) + P1[3] * (P3[0] * (T3[16] * V1[3] -
      T3[15] * V1[4]) + (P3[1] * (T3[14] * V1[4] - T3[16] * V1[2]) + P3[2] *
      (T3[15] * V1[2] - T3[14] * V1[3])))))))));
  TMP27 = -1. * (P1[0] * (P2[0] * (P3[1] * (T3[14] * V2[4] - T3[10] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[14] * V2[3]) + P3[3] * (T3[10] * V2[3] -
      T3[6] * V2[4]))) + (P2[1] * (P3[0] * (T3[10] * V2[5] - T3[14] * V2[4]) +
      (P3[2] * (T3[14] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] -
      T3[10] * V2[2]))) + (P2[2] * (P3[0] * (T3[14] * V2[3] - T3[6] * V2[5]) +
      (P3[1] * (T3[2] * V2[5] - T3[14] * V2[2]) + P3[3] * (T3[6] * V2[2] -
      T3[2] * V2[3]))) + P2[3] * (P3[0] * (T3[6] * V2[4] - T3[10] * V2[3]) +
      (P3[1] * (T3[10] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] -
      T3[6] * V2[2])))))) + (P1[1] * (P2[0] * (P3[1] * (T3[11] * V2[5] - T3[15]
      * V2[4]) + (P3[2] * (T3[15] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] *
      V2[4] - T3[11] * V2[3]))) + (P2[1] * (P3[0] * (T3[15] * V2[4] - T3[11] *
      V2[5]) + (P3[2] * (T3[3] * V2[5] - T3[15] * V2[2]) + P3[3] * (T3[11] *
      V2[2] - T3[3] * V2[4]))) + (P2[2] * (P3[0] * (T3[7] * V2[5] - T3[15] *
      V2[3]) + (P3[1] * (T3[15] * V2[2] - T3[3] * V2[5]) + P3[3] * (T3[3] *
      V2[3] - T3[7] * V2[2]))) + P2[3] * (P3[0] * (T3[11] * V2[3] - T3[7] *
      V2[4]) + (P3[1] * (T3[3] * V2[4] - T3[11] * V2[2]) + P3[2] * (T3[7] *
      V2[2] - T3[3] * V2[3])))))) + (P1[2] * (P2[0] * (P3[1] * (T3[12] * V2[5]
      - T3[16] * V2[4]) + (P3[2] * (T3[16] * V2[3] - T3[8] * V2[5]) + P3[3] *
      (T3[8] * V2[4] - T3[12] * V2[3]))) + (P2[1] * (P3[0] * (T3[16] * V2[4] -
      T3[12] * V2[5]) + (P3[2] * (T3[4] * V2[5] - T3[16] * V2[2]) + P3[3] *
      (T3[12] * V2[2] - T3[4] * V2[4]))) + (P2[2] * (P3[0] * (T3[8] * V2[5] -
      T3[16] * V2[3]) + (P3[1] * (T3[16] * V2[2] - T3[4] * V2[5]) + P3[3] *
      (T3[4] * V2[3] - T3[8] * V2[2]))) + P2[3] * (P3[0] * (T3[12] * V2[3] -
      T3[8] * V2[4]) + (P3[1] * (T3[4] * V2[4] - T3[12] * V2[2]) + P3[2] *
      (T3[8] * V2[2] - T3[4] * V2[3])))))) + P1[3] * (P2[0] * (P3[1] * (T3[13]
      * V2[5] - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[9] * V2[5]) +
      P3[3] * (T3[9] * V2[4] - T3[13] * V2[3]))) + (P2[1] * (P3[0] * (T3[17] *
      V2[4] - T3[13] * V2[5]) + (P3[2] * (T3[5] * V2[5] - T3[17] * V2[2]) +
      P3[3] * (T3[13] * V2[2] - T3[5] * V2[4]))) + (P2[2] * (P3[0] * (T3[9] *
      V2[5] - T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[5] * V2[5]) +
      P3[3] * (T3[5] * V2[3] - T3[9] * V2[2]))) + P2[3] * (P3[0] * (T3[13] *
      V2[3] - T3[9] * V2[4]) + (P3[1] * (T3[5] * V2[4] - T3[13] * V2[2]) +
      P3[2] * (T3[9] * V2[2] - T3[5] * V2[3])))))))));
  TMP23 = -1. * (P1[0] * (P2[0] * (P3[1] * (T3[5] * V2[4] - T3[4] * V2[5]) +
      (P3[2] * (T3[3] * V2[5] - T3[5] * V2[3]) + P3[3] * (T3[4] * V2[3] - T3[3]
      * V2[4]))) + (P2[1] * (P3[0] * (T3[4] * V2[5] - T3[5] * V2[4]) + (P3[2] *
      (T3[5] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] - T3[4] *
      V2[2]))) + (P2[2] * (P3[0] * (T3[5] * V2[3] - T3[3] * V2[5]) + (P3[1] *
      (T3[2] * V2[5] - T3[5] * V2[2]) + P3[3] * (T3[3] * V2[2] - T3[2] *
      V2[3]))) + P2[3] * (P3[0] * (T3[3] * V2[4] - T3[4] * V2[3]) + (P3[1] *
      (T3[4] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] - T3[3] *
      V2[2])))))) + (P1[1] * (P2[0] * (P3[1] * (T3[8] * V2[5] - T3[9] * V2[4])
      + (P3[2] * (T3[9] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] * V2[4] -
      T3[8] * V2[3]))) + (P2[1] * (P3[0] * (T3[9] * V2[4] - T3[8] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[9] * V2[2]) + P3[3] * (T3[8] * V2[2] - T3[6]
      * V2[4]))) + (P2[2] * (P3[0] * (T3[7] * V2[5] - T3[9] * V2[3]) + (P3[1] *
      (T3[9] * V2[2] - T3[6] * V2[5]) + P3[3] * (T3[6] * V2[3] - T3[7] *
      V2[2]))) + P2[3] * (P3[0] * (T3[8] * V2[3] - T3[7] * V2[4]) + (P3[1] *
      (T3[6] * V2[4] - T3[8] * V2[2]) + P3[2] * (T3[7] * V2[2] - T3[6] *
      V2[3])))))) + (P1[2] * (P2[0] * (P3[1] * (T3[12] * V2[5] - T3[13] *
      V2[4]) + (P3[2] * (T3[13] * V2[3] - T3[11] * V2[5]) + P3[3] * (T3[11] *
      V2[4] - T3[12] * V2[3]))) + (P2[1] * (P3[0] * (T3[13] * V2[4] - T3[12] *
      V2[5]) + (P3[2] * (T3[10] * V2[5] - T3[13] * V2[2]) + P3[3] * (T3[12] *
      V2[2] - T3[10] * V2[4]))) + (P2[2] * (P3[0] * (T3[11] * V2[5] - T3[13] *
      V2[3]) + (P3[1] * (T3[13] * V2[2] - T3[10] * V2[5]) + P3[3] * (T3[10] *
      V2[3] - T3[11] * V2[2]))) + P2[3] * (P3[0] * (T3[12] * V2[3] - T3[11] *
      V2[4]) + (P3[1] * (T3[10] * V2[4] - T3[12] * V2[2]) + P3[2] * (T3[11] *
      V2[2] - T3[10] * V2[3])))))) + P1[3] * (P2[0] * (P3[1] * (T3[16] * V2[5]
      - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[15] * V2[5]) + P3[3] *
      (T3[15] * V2[4] - T3[16] * V2[3]))) + (P2[1] * (P3[0] * (T3[17] * V2[4] -
      T3[16] * V2[5]) + (P3[2] * (T3[14] * V2[5] - T3[17] * V2[2]) + P3[3] *
      (T3[16] * V2[2] - T3[14] * V2[4]))) + (P2[2] * (P3[0] * (T3[15] * V2[5] -
      T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[14] * V2[5]) + P3[3] *
      (T3[14] * V2[3] - T3[15] * V2[2]))) + P2[3] * (P3[0] * (T3[16] * V2[3] -
      T3[15] * V2[4]) + (P3[1] * (T3[14] * V2[4] - T3[16] * V2[2]) + P3[2] *
      (T3[15] * V2[2] - T3[14] * V2[3])))))))));
  TMP28 = -1. * (P1[0] * (P2[0] * (P3[1] * (T3[14] * V1[4] - T3[10] * V1[5]) +
      (P3[2] * (T3[6] * V1[5] - T3[14] * V1[3]) + P3[3] * (T3[10] * V1[3] -
      T3[6] * V1[4]))) + (P2[1] * (P3[0] * (T3[10] * V1[5] - T3[14] * V1[4]) +
      (P3[2] * (T3[14] * V1[2] - T3[2] * V1[5]) + P3[3] * (T3[2] * V1[4] -
      T3[10] * V1[2]))) + (P2[2] * (P3[0] * (T3[14] * V1[3] - T3[6] * V1[5]) +
      (P3[1] * (T3[2] * V1[5] - T3[14] * V1[2]) + P3[3] * (T3[6] * V1[2] -
      T3[2] * V1[3]))) + P2[3] * (P3[0] * (T3[6] * V1[4] - T3[10] * V1[3]) +
      (P3[1] * (T3[10] * V1[2] - T3[2] * V1[4]) + P3[2] * (T3[2] * V1[3] -
      T3[6] * V1[2])))))) + (P1[1] * (P2[0] * (P3[1] * (T3[11] * V1[5] - T3[15]
      * V1[4]) + (P3[2] * (T3[15] * V1[3] - T3[7] * V1[5]) + P3[3] * (T3[7] *
      V1[4] - T3[11] * V1[3]))) + (P2[1] * (P3[0] * (T3[15] * V1[4] - T3[11] *
      V1[5]) + (P3[2] * (T3[3] * V1[5] - T3[15] * V1[2]) + P3[3] * (T3[11] *
      V1[2] - T3[3] * V1[4]))) + (P2[2] * (P3[0] * (T3[7] * V1[5] - T3[15] *
      V1[3]) + (P3[1] * (T3[15] * V1[2] - T3[3] * V1[5]) + P3[3] * (T3[3] *
      V1[3] - T3[7] * V1[2]))) + P2[3] * (P3[0] * (T3[11] * V1[3] - T3[7] *
      V1[4]) + (P3[1] * (T3[3] * V1[4] - T3[11] * V1[2]) + P3[2] * (T3[7] *
      V1[2] - T3[3] * V1[3])))))) + (P1[2] * (P2[0] * (P3[1] * (T3[12] * V1[5]
      - T3[16] * V1[4]) + (P3[2] * (T3[16] * V1[3] - T3[8] * V1[5]) + P3[3] *
      (T3[8] * V1[4] - T3[12] * V1[3]))) + (P2[1] * (P3[0] * (T3[16] * V1[4] -
      T3[12] * V1[5]) + (P3[2] * (T3[4] * V1[5] - T3[16] * V1[2]) + P3[3] *
      (T3[12] * V1[2] - T3[4] * V1[4]))) + (P2[2] * (P3[0] * (T3[8] * V1[5] -
      T3[16] * V1[3]) + (P3[1] * (T3[16] * V1[2] - T3[4] * V1[5]) + P3[3] *
      (T3[4] * V1[3] - T3[8] * V1[2]))) + P2[3] * (P3[0] * (T3[12] * V1[3] -
      T3[8] * V1[4]) + (P3[1] * (T3[4] * V1[4] - T3[12] * V1[2]) + P3[2] *
      (T3[8] * V1[2] - T3[4] * V1[3])))))) + P1[3] * (P2[0] * (P3[1] * (T3[13]
      * V1[5] - T3[17] * V1[4]) + (P3[2] * (T3[17] * V1[3] - T3[9] * V1[5]) +
      P3[3] * (T3[9] * V1[4] - T3[13] * V1[3]))) + (P2[1] * (P3[0] * (T3[17] *
      V1[4] - T3[13] * V1[5]) + (P3[2] * (T3[5] * V1[5] - T3[17] * V1[2]) +
      P3[3] * (T3[13] * V1[2] - T3[5] * V1[4]))) + (P2[2] * (P3[0] * (T3[9] *
      V1[5] - T3[17] * V1[3]) + (P3[1] * (T3[17] * V1[2] - T3[5] * V1[5]) +
      P3[3] * (T3[5] * V1[3] - T3[9] * V1[2]))) + P2[3] * (P3[0] * (T3[13] *
      V1[3] - T3[9] * V1[4]) + (P3[1] * (T3[5] * V1[4] - T3[13] * V1[2]) +
      P3[2] * (T3[9] * V1[2] - T3[5] * V1[3])))))))));
  TMP29 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[14] * V2[4] - T3[10] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[14] * V2[3]) + P3[3] * (T3[10] * V2[3] -
      T3[6] * V2[4]))) + (P1[1] * (P3[0] * (T3[10] * V2[5] - T3[14] * V2[4]) +
      (P3[2] * (T3[14] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] -
      T3[10] * V2[2]))) + (P1[2] * (P3[0] * (T3[14] * V2[3] - T3[6] * V2[5]) +
      (P3[1] * (T3[2] * V2[5] - T3[14] * V2[2]) + P3[3] * (T3[6] * V2[2] -
      T3[2] * V2[3]))) + P1[3] * (P3[0] * (T3[6] * V2[4] - T3[10] * V2[3]) +
      (P3[1] * (T3[10] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] -
      T3[6] * V2[2])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[11] * V2[5] - T3[15]
      * V2[4]) + (P3[2] * (T3[15] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] *
      V2[4] - T3[11] * V2[3]))) + (P1[1] * (P3[0] * (T3[15] * V2[4] - T3[11] *
      V2[5]) + (P3[2] * (T3[3] * V2[5] - T3[15] * V2[2]) + P3[3] * (T3[11] *
      V2[2] - T3[3] * V2[4]))) + (P1[2] * (P3[0] * (T3[7] * V2[5] - T3[15] *
      V2[3]) + (P3[1] * (T3[15] * V2[2] - T3[3] * V2[5]) + P3[3] * (T3[3] *
      V2[3] - T3[7] * V2[2]))) + P1[3] * (P3[0] * (T3[11] * V2[3] - T3[7] *
      V2[4]) + (P3[1] * (T3[3] * V2[4] - T3[11] * V2[2]) + P3[2] * (T3[7] *
      V2[2] - T3[3] * V2[3])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[12] * V2[5]
      - T3[16] * V2[4]) + (P3[2] * (T3[16] * V2[3] - T3[8] * V2[5]) + P3[3] *
      (T3[8] * V2[4] - T3[12] * V2[3]))) + (P1[1] * (P3[0] * (T3[16] * V2[4] -
      T3[12] * V2[5]) + (P3[2] * (T3[4] * V2[5] - T3[16] * V2[2]) + P3[3] *
      (T3[12] * V2[2] - T3[4] * V2[4]))) + (P1[2] * (P3[0] * (T3[8] * V2[5] -
      T3[16] * V2[3]) + (P3[1] * (T3[16] * V2[2] - T3[4] * V2[5]) + P3[3] *
      (T3[4] * V2[3] - T3[8] * V2[2]))) + P1[3] * (P3[0] * (T3[12] * V2[3] -
      T3[8] * V2[4]) + (P3[1] * (T3[4] * V2[4] - T3[12] * V2[2]) + P3[2] *
      (T3[8] * V2[2] - T3[4] * V2[3])))))) + P2[3] * (P1[0] * (P3[1] * (T3[13]
      * V2[5] - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[9] * V2[5]) +
      P3[3] * (T3[9] * V2[4] - T3[13] * V2[3]))) + (P1[1] * (P3[0] * (T3[17] *
      V2[4] - T3[13] * V2[5]) + (P3[2] * (T3[5] * V2[5] - T3[17] * V2[2]) +
      P3[3] * (T3[13] * V2[2] - T3[5] * V2[4]))) + (P1[2] * (P3[0] * (T3[9] *
      V2[5] - T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[5] * V2[5]) +
      P3[3] * (T3[5] * V2[3] - T3[9] * V2[2]))) + P1[3] * (P3[0] * (T3[13] *
      V2[3] - T3[9] * V2[4]) + (P3[1] * (T3[5] * V2[4] - T3[13] * V2[2]) +
      P3[2] * (T3[9] * V2[2] - T3[5] * V2[3])))))))));
  TMP14 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP17 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP30 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[14] * V1[4] - T3[10] * V1[5]) +
      (P3[2] * (T3[6] * V1[5] - T3[14] * V1[3]) + P3[3] * (T3[10] * V1[3] -
      T3[6] * V1[4]))) + (P1[1] * (P3[0] * (T3[10] * V1[5] - T3[14] * V1[4]) +
      (P3[2] * (T3[14] * V1[2] - T3[2] * V1[5]) + P3[3] * (T3[2] * V1[4] -
      T3[10] * V1[2]))) + (P1[2] * (P3[0] * (T3[14] * V1[3] - T3[6] * V1[5]) +
      (P3[1] * (T3[2] * V1[5] - T3[14] * V1[2]) + P3[3] * (T3[6] * V1[2] -
      T3[2] * V1[3]))) + P1[3] * (P3[0] * (T3[6] * V1[4] - T3[10] * V1[3]) +
      (P3[1] * (T3[10] * V1[2] - T3[2] * V1[4]) + P3[2] * (T3[2] * V1[3] -
      T3[6] * V1[2])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[11] * V1[5] - T3[15]
      * V1[4]) + (P3[2] * (T3[15] * V1[3] - T3[7] * V1[5]) + P3[3] * (T3[7] *
      V1[4] - T3[11] * V1[3]))) + (P1[1] * (P3[0] * (T3[15] * V1[4] - T3[11] *
      V1[5]) + (P3[2] * (T3[3] * V1[5] - T3[15] * V1[2]) + P3[3] * (T3[11] *
      V1[2] - T3[3] * V1[4]))) + (P1[2] * (P3[0] * (T3[7] * V1[5] - T3[15] *
      V1[3]) + (P3[1] * (T3[15] * V1[2] - T3[3] * V1[5]) + P3[3] * (T3[3] *
      V1[3] - T3[7] * V1[2]))) + P1[3] * (P3[0] * (T3[11] * V1[3] - T3[7] *
      V1[4]) + (P3[1] * (T3[3] * V1[4] - T3[11] * V1[2]) + P3[2] * (T3[7] *
      V1[2] - T3[3] * V1[3])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[12] * V1[5]
      - T3[16] * V1[4]) + (P3[2] * (T3[16] * V1[3] - T3[8] * V1[5]) + P3[3] *
      (T3[8] * V1[4] - T3[12] * V1[3]))) + (P1[1] * (P3[0] * (T3[16] * V1[4] -
      T3[12] * V1[5]) + (P3[2] * (T3[4] * V1[5] - T3[16] * V1[2]) + P3[3] *
      (T3[12] * V1[2] - T3[4] * V1[4]))) + (P1[2] * (P3[0] * (T3[8] * V1[5] -
      T3[16] * V1[3]) + (P3[1] * (T3[16] * V1[2] - T3[4] * V1[5]) + P3[3] *
      (T3[4] * V1[3] - T3[8] * V1[2]))) + P1[3] * (P3[0] * (T3[12] * V1[3] -
      T3[8] * V1[4]) + (P3[1] * (T3[4] * V1[4] - T3[12] * V1[2]) + P3[2] *
      (T3[8] * V1[2] - T3[4] * V1[3])))))) + P2[3] * (P1[0] * (P3[1] * (T3[13]
      * V1[5] - T3[17] * V1[4]) + (P3[2] * (T3[17] * V1[3] - T3[9] * V1[5]) +
      P3[3] * (T3[9] * V1[4] - T3[13] * V1[3]))) + (P1[1] * (P3[0] * (T3[17] *
      V1[4] - T3[13] * V1[5]) + (P3[2] * (T3[5] * V1[5] - T3[17] * V1[2]) +
      P3[3] * (T3[13] * V1[2] - T3[5] * V1[4]))) + (P1[2] * (P3[0] * (T3[9] *
      V1[5] - T3[17] * V1[3]) + (P3[1] * (T3[17] * V1[2] - T3[5] * V1[5]) +
      P3[3] * (T3[5] * V1[3] - T3[9] * V1[2]))) + P1[3] * (P3[0] * (T3[13] *
      V1[3] - T3[9] * V1[4]) + (P3[1] * (T3[5] * V1[4] - T3[13] * V1[2]) +
      P3[2] * (T3[9] * V1[2] - T3[5] * V1[3])))))))));
  vertex = COUP * - 1. * (TMP14 * (+cI * (TMP24 + TMP26 + TMP28 + TMP30)) +
      TMP17 * (+cI * (TMP23 + TMP25 + TMP27 + TMP29)));
}


void FFV2P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] + F2[3]
      * F1[5]);
  V3[3] = denom * - cI * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] * F1[2] - F2[4]
      * F1[3]);
  V3[4] = denom * - cI * (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4] *
      F1[3] + F2[3] * F1[4]));
  V3[5] = denom * - cI * (F2[5] * F1[3] + F2[2] * F1[4] - F2[4] * F1[2] - F2[3]
      * F1[5]);
}


void FFT1_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP22; 
  complex<double> TMP10; 
  double P2[4]; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP22 = (F2[4] * F1[4] + F2[5] * F1[5] - F2[2] * F1[2] - F2[3] * F1[3]); 
  vertex = COUP * - TMP22 * (+cI * (TMP10 + TMP11)); 
}

void FFT1_2_3_5_0(complex<double> F1[], complex<double> F2[], complex<double>
    T3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> tmp; 
  FFT1_0(F1, F2, T3, COUP1, vertex); 
  FFT2_0(F1, F2, T3, COUP2, tmp); 
  vertex = vertex + tmp; 
  FFT3_0(F1, F2, T3, COUP3, tmp); 
  vertex = vertex + tmp; 
  FFT5_0(F1, F2, T3, COUP4, tmp); 
  vertex = vertex + tmp; 
}

void VVT2_0(complex<double> V1[], complex<double> V2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP36; 
  complex<double> TMP35; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP11 = (P1[0] * - 1. * (T3[3] * P2[1] + T3[4] * P2[2] + T3[5] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[8] * P2[2] + T3[9] * P2[3]
      - T3[6] * P2[0]) + (P1[2] * (T3[11] * P2[1] + T3[12] * P2[2] + T3[13] *
      P2[3] - T3[10] * P2[0]) + P1[3] * (T3[15] * P2[1] + T3[16] * P2[2] +
      T3[17] * P2[3] - T3[14] * P2[0]))));
  TMP10 = (P1[0] * - 1. * (T3[6] * P2[1] + T3[10] * P2[2] + T3[14] * P2[3] -
      T3[2] * P2[0]) + (P1[1] * (T3[7] * P2[1] + T3[11] * P2[2] + T3[15] *
      P2[3] - T3[3] * P2[0]) + (P1[2] * (T3[8] * P2[1] + T3[12] * P2[2] +
      T3[16] * P2[3] - T3[4] * P2[0]) + P1[3] * (T3[9] * P2[1] + T3[13] * P2[2]
      + T3[17] * P2[3] - T3[5] * P2[0]))));
  TMP35 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP36 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  vertex = COUP * (TMP10 * (-cI * (TMP36) + cI * (TMP35)) + TMP11 * (-cI *
      (TMP36) + cI * (TMP35)));
}


void VVT13_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP38; 
  complex<double> TMP36; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP15; 
  complex<double> TMP32; 
  complex<double> denom; 
  complex<double> TMP13; 
  complex<double> TMP34; 
  complex<double> TMP29; 
  double OM1; 
  complex<double> TMP9; 
  complex<double> TMP35; 
  double P1[4]; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP29 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP35 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP32 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP36 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP34 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP9 * (TMP12 * (OM1 * P1[0] * (-cI * (TMP29 + TMP32) + cI *
      (TMP35 + TMP36)) + (V2[3] * (+cI * (T3[6] + T3[3])) + (V2[4] * (+cI *
      (T3[10] + T3[4])) + (V2[5] * (+cI * (T3[14] + T3[5])) - 2. * cI * (T3[2]
      * V2[2]))))) + TMP38 * (OM1 * P1[0] * (-cI * (TMP34) + cI * (TMP13)) +
      (P2[1] * - 1. * (+cI * (T3[3] + T3[6])) + (P2[2] * - 1. * (+cI * (T3[4] +
      T3[10])) + (P2[3] * - 1. * (+cI * (T3[5] + T3[14])) + 2. * cI * (P2[0] *
      T3[2])))))) + P3[0] * (TMP12 * (+cI * (TMP29 + TMP32)) - TMP38 * (+cI *
      (TMP13 + TMP15))));
  V1[3] = denom * (TMP9 * (TMP12 * (OM1 * P1[1] * (-cI * (TMP29 + TMP32) + cI *
      (TMP35 + TMP36)) + (V2[2] * - 1. * (+cI * (T3[3] + T3[6])) + (V2[4] *
      (+cI * (T3[11] + T3[8])) + (V2[5] * (+cI * (T3[15] + T3[9])) + 2. * cI *
      (T3[7] * V2[3]))))) + TMP38 * (OM1 * P1[1] * (-cI * (TMP34) + cI *
      (TMP13)) + (P2[0] * (+cI * (T3[6] + T3[3])) + (P2[2] * - 1. * (+cI *
      (T3[8] + T3[11])) + (P2[3] * - 1. * (+cI * (T3[9] + T3[15])) - 2. * cI *
      (P2[1] * T3[7])))))) + P3[1] * (TMP12 * (+cI * (TMP29 + TMP32)) - TMP38 *
      (+cI * (TMP13 + TMP15))));
  V1[4] = denom * (TMP9 * (TMP12 * (OM1 * P1[2] * (-cI * (TMP29 + TMP32) + cI *
      (TMP35 + TMP36)) + (V2[2] * - 1. * (+cI * (T3[4] + T3[10])) + (V2[3] *
      (+cI * (T3[8] + T3[11])) + (V2[5] * (+cI * (T3[16] + T3[13])) + 2. * cI *
      (T3[12] * V2[4]))))) + TMP38 * (OM1 * P1[2] * (-cI * (TMP34) + cI *
      (TMP13)) + (P2[0] * (+cI * (T3[10] + T3[4])) + (P2[1] * - 1. * (+cI *
      (T3[11] + T3[8])) + (P2[3] * - 1. * (+cI * (T3[13] + T3[16])) - 2. * cI *
      (P2[2] * T3[12])))))) + P3[2] * (TMP12 * (+cI * (TMP29 + TMP32)) - TMP38
      * (+cI * (TMP13 + TMP15))));
  V1[5] = denom * (TMP9 * (TMP12 * (OM1 * P1[3] * (-cI * (TMP29 + TMP32) + cI *
      (TMP35 + TMP36)) + (V2[2] * - 1. * (+cI * (T3[5] + T3[14])) + (V2[3] *
      (+cI * (T3[9] + T3[15])) + (V2[4] * (+cI * (T3[13] + T3[16])) + 2. * cI *
      (T3[17] * V2[5]))))) + TMP38 * (OM1 * P1[3] * (-cI * (TMP34) + cI *
      (TMP13)) + (P2[0] * (+cI * (T3[14] + T3[5])) + (P2[1] * - 1. * (+cI *
      (T3[15] + T3[9])) + (P2[2] * - 1. * (+cI * (T3[16] + T3[13])) - 2. * cI *
      (P2[3] * T3[17])))))) + P3[3] * (TMP12 * (+cI * (TMP29 + TMP32)) - TMP38
      * (+cI * (TMP13 + TMP15))));
}


void FFT4_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP5; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP6; 
  complex<double> TMP4; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP5 = -1. * (P2[0] * (F1[2] * (F2[4] * (T3[2] + T3[14]) + F2[5] * (T3[6] +
      cI * (T3[10]))) + (F1[3] * (F2[4] * (T3[6] - cI * (T3[10])) + F2[5] *
      (T3[2] - T3[14])) + (F1[4] * (F2[2] * (T3[14] - T3[2]) + F2[3] * (T3[6] +
      cI * (T3[10]))) + F1[5] * (F2[2] * (T3[6] - cI * (T3[10])) - F2[3] *
      (T3[2] + T3[14]))))) + (P2[1] * (F1[2] * (F2[4] * - 1. * (T3[3] + T3[15])
      - F2[5] * (T3[7] + cI * (T3[11]))) + (F1[3] * (F2[4] * (+cI * (T3[11]) -
      T3[7]) + F2[5] * (T3[15] - T3[3])) + (F1[4] * (F2[2] * (T3[3] - T3[15]) -
      F2[3] * (T3[7] + cI * (T3[11]))) + F1[5] * (F2[2] * (+cI * (T3[11]) -
      T3[7]) + F2[3] * (T3[3] + T3[15]))))) + (P2[2] * (F1[2] * (F2[4] * - 1. *
      (T3[4] + T3[16]) - F2[5] * (T3[8] + cI * (T3[12]))) + (F1[3] * (F2[4] *
      (+cI * (T3[12]) - T3[8]) + F2[5] * (T3[16] - T3[4])) + (F1[4] * (F2[2] *
      (T3[4] - T3[16]) - F2[3] * (T3[8] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (+cI * (T3[12]) - T3[8]) + F2[3] * (T3[4] + T3[16]))))) + P2[3] * (F1[2]
      * (F2[4] * - 1. * (T3[5] + T3[17]) - F2[5] * (T3[9] + cI * (T3[13]))) +
      (F1[3] * (F2[4] * (+cI * (T3[13]) - T3[9]) + F2[5] * (T3[17] - T3[5])) +
      (F1[4] * (F2[2] * (T3[5] - T3[17]) - F2[3] * (T3[9] + cI * (T3[13]))) +
      F1[5] * (F2[2] * (+cI * (T3[13]) - T3[9]) + F2[3] * (T3[5] +
      T3[17]))))))));
  TMP4 = -1. * (P1[0] * (F1[2] * (F2[4] * (T3[2] + T3[14]) + F2[5] * (T3[6] +
      cI * (T3[10]))) + (F1[3] * (F2[4] * (T3[6] - cI * (T3[10])) + F2[5] *
      (T3[2] - T3[14])) + (F1[4] * (F2[2] * (T3[14] - T3[2]) + F2[3] * (T3[6] +
      cI * (T3[10]))) + F1[5] * (F2[2] * (T3[6] - cI * (T3[10])) - F2[3] *
      (T3[2] + T3[14]))))) + (P1[1] * (F1[2] * (F2[4] * - 1. * (T3[3] + T3[15])
      - F2[5] * (T3[7] + cI * (T3[11]))) + (F1[3] * (F2[4] * (+cI * (T3[11]) -
      T3[7]) + F2[5] * (T3[15] - T3[3])) + (F1[4] * (F2[2] * (T3[3] - T3[15]) -
      F2[3] * (T3[7] + cI * (T3[11]))) + F1[5] * (F2[2] * (+cI * (T3[11]) -
      T3[7]) + F2[3] * (T3[3] + T3[15]))))) + (P1[2] * (F1[2] * (F2[4] * - 1. *
      (T3[4] + T3[16]) - F2[5] * (T3[8] + cI * (T3[12]))) + (F1[3] * (F2[4] *
      (+cI * (T3[12]) - T3[8]) + F2[5] * (T3[16] - T3[4])) + (F1[4] * (F2[2] *
      (T3[4] - T3[16]) - F2[3] * (T3[8] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (+cI * (T3[12]) - T3[8]) + F2[3] * (T3[4] + T3[16]))))) + P1[3] * (F1[2]
      * (F2[4] * - 1. * (T3[5] + T3[17]) - F2[5] * (T3[9] + cI * (T3[13]))) +
      (F1[3] * (F2[4] * (+cI * (T3[13]) - T3[9]) + F2[5] * (T3[17] - T3[5])) +
      (F1[4] * (F2[2] * (T3[5] - T3[17]) - F2[3] * (T3[9] + cI * (T3[13]))) +
      F1[5] * (F2[2] * (+cI * (T3[13]) - T3[9]) + F2[3] * (T3[5] +
      T3[17]))))))));
  TMP7 = -1. * (P2[0] * (F1[2] * (F2[4] * (T3[2] + T3[5]) + F2[5] * (T3[3] + cI
      * (T3[4]))) + (F1[3] * (F2[4] * (T3[3] - cI * (T3[4])) + F2[5] * (T3[2] -
      T3[5])) + (F1[4] * (F2[2] * (T3[5] - T3[2]) + F2[3] * (T3[3] + cI *
      (T3[4]))) + F1[5] * (F2[2] * (T3[3] - cI * (T3[4])) - F2[3] * (T3[2] +
      T3[5]))))) + (P2[1] * (F1[2] * (F2[4] * - 1. * (T3[6] + T3[9]) - F2[5] *
      (T3[7] + cI * (T3[8]))) + (F1[3] * (F2[4] * (+cI * (T3[8]) - T3[7]) +
      F2[5] * (T3[9] - T3[6])) + (F1[4] * (F2[2] * (T3[6] - T3[9]) - F2[3] *
      (T3[7] + cI * (T3[8]))) + F1[5] * (F2[2] * (+cI * (T3[8]) - T3[7]) +
      F2[3] * (T3[6] + T3[9]))))) + (P2[2] * (F1[2] * (F2[4] * - 1. * (T3[10] +
      T3[13]) - F2[5] * (T3[11] + cI * (T3[12]))) + (F1[3] * (F2[4] * (+cI *
      (T3[12]) - T3[11]) + F2[5] * (T3[13] - T3[10])) + (F1[4] * (F2[2] *
      (T3[10] - T3[13]) - F2[3] * (T3[11] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (+cI * (T3[12]) - T3[11]) + F2[3] * (T3[10] + T3[13]))))) + P2[3] *
      (F1[2] * (F2[4] * - 1. * (T3[14] + T3[17]) - F2[5] * (T3[15] + cI *
      (T3[16]))) + (F1[3] * (F2[4] * (+cI * (T3[16]) - T3[15]) + F2[5] *
      (T3[17] - T3[14])) + (F1[4] * (F2[2] * (T3[14] - T3[17]) - F2[3] *
      (T3[15] + cI * (T3[16]))) + F1[5] * (F2[2] * (+cI * (T3[16]) - T3[15]) +
      F2[3] * (T3[14] + T3[17]))))))));
  TMP6 = -1. * (P1[0] * (F1[2] * (F2[4] * (T3[2] + T3[5]) + F2[5] * (T3[3] + cI
      * (T3[4]))) + (F1[3] * (F2[4] * (T3[3] - cI * (T3[4])) + F2[5] * (T3[2] -
      T3[5])) + (F1[4] * (F2[2] * (T3[5] - T3[2]) + F2[3] * (T3[3] + cI *
      (T3[4]))) + F1[5] * (F2[2] * (T3[3] - cI * (T3[4])) - F2[3] * (T3[2] +
      T3[5]))))) + (P1[1] * (F1[2] * (F2[4] * - 1. * (T3[6] + T3[9]) - F2[5] *
      (T3[7] + cI * (T3[8]))) + (F1[3] * (F2[4] * (+cI * (T3[8]) - T3[7]) +
      F2[5] * (T3[9] - T3[6])) + (F1[4] * (F2[2] * (T3[6] - T3[9]) - F2[3] *
      (T3[7] + cI * (T3[8]))) + F1[5] * (F2[2] * (+cI * (T3[8]) - T3[7]) +
      F2[3] * (T3[6] + T3[9]))))) + (P1[2] * (F1[2] * (F2[4] * - 1. * (T3[10] +
      T3[13]) - F2[5] * (T3[11] + cI * (T3[12]))) + (F1[3] * (F2[4] * (+cI *
      (T3[12]) - T3[11]) + F2[5] * (T3[13] - T3[10])) + (F1[4] * (F2[2] *
      (T3[10] - T3[13]) - F2[3] * (T3[11] + cI * (T3[12]))) + F1[5] * (F2[2] *
      (+cI * (T3[12]) - T3[11]) + F2[3] * (T3[10] + T3[13]))))) + P1[3] *
      (F1[2] * (F2[4] * - 1. * (T3[14] + T3[17]) - F2[5] * (T3[15] + cI *
      (T3[16]))) + (F1[3] * (F2[4] * (+cI * (T3[16]) - T3[15]) + F2[5] *
      (T3[17] - T3[14])) + (F1[4] * (F2[2] * (T3[14] - T3[17]) - F2[3] *
      (T3[15] + cI * (T3[16]))) + F1[5] * (F2[2] * (+cI * (T3[16]) - T3[15]) +
      F2[3] * (T3[14] + T3[17]))))))));
  vertex = COUP * (-cI * (TMP5 + TMP7) + cI * (TMP4 + TMP6)); 
}


void FFV8_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
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
  F1[2] = denom * 4. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3]
      - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] -
      V3[2])))) + (+1./4. * (M1 * (F2[5] * (V3[3] + cI * (V3[4])) + 4. * (F2[4]
      * 1./4. * (V3[2] + V3[5])))) + F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) +
      (P1[1] * - 1. * (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * 4. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (+1./4. * (M1 * (F2[5] * (V3[2] - V3[5]) + 4. *
      (F2[4] * 1./4. * (V3[3] - cI * (V3[4]))))) + F2[3] * (P1[0] * - 1. *
      (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI
      * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 4. * (V3[5] - V3[2]) + 4. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * 4. * (+cI * (V3[4]) - V3[3]) + 4. * (F2[3] *
      (V3[2] + V3[5])))));
}


void FFT3_1(complex<double> F2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + T3[0]; 
  F1[1] = +F2[1] + T3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (F2[3] * (P1[0] * (P1[3] * (T3[9] + T3[15] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[6] - T3[3]) + (P1[1] *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8])) + (P1[2]
      * (T3[8] + T3[11] + cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) +
      (P1[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[0] * (T3[6]
      + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] +
      T3[8])) - P2[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P1[1] *
      (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13] + cI * (T3[3] + T3[15]
      + T3[6] + T3[9])) + (P1[3] * (+2. * (T3[7]) + cI * (T3[11] + T3[8]) - 2.
      * (T3[17]) - T3[5] - T3[14]) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] +
      T3[16] + T3[10] + T3[13]) + (P2[0] * - 1. * (T3[14] + T3[5] + 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[17]) +
      2. * cI * (T3[12])) + (P1[2] * - 1. * (+cI * (T3[4] + T3[16] + T3[10] +
      T3[13])) + (P2[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P2[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[0] * - 1. * (+cI *
      (T3[14] + T3[5]) + 2. * cI * (T3[2])) + P2[3] * (+cI * (T3[5] + T3[14]) +
      2. * cI * (T3[17]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P2[1] * - 1.
      * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8] + T3[11] + 2.
      * cI * (T3[12]))))))))) + (F2[2] * (P1[0] * (P1[1] * (T3[15] + T3[9] + cI
      * (T3[10] + T3[4]) - 2. * (T3[3] + T3[6])) + (P1[2] * - 1. * (+2. *
      (T3[4] + T3[10]) + cI * (T3[6] + T3[3]) - T3[16] - T3[13]) + (P1[3] * 2.
      * (T3[17] + T3[2] - T3[5] - T3[14]) + (P2[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[0] * - 1. *
      (T3[14] + T3[5] - 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] - 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] - 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * - 1. * (T3[3] + T3[6] + cI * (T3[13] + T3[16]) - 2. * (T3[9] +
      T3[15])) + (P1[2] * (+2. * (T3[13] + T3[16]) + cI * (T3[9] + T3[15]) -
      T3[4] - T3[10]) + (P2[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] - 2.
      * (T3[17])) + (P2[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P2[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. * (T3[8] + T3[11] -
      cI * (T3[12]) + cI * (T3[7])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] *
      - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + (P2[1] * (+cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] + cI * (T3[6] + T3[3]))
      + (P2[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] + T3[15])) + (P1[2] *
      (+2. * (T3[12]) + cI * (T3[8] + T3[11])) + (P2[1] * - 1. * (T3[11] +
      T3[8] + 2. * cI * (T3[7])) - P2[2] * (+2. * (T3[12]) + cI * (T3[8] +
      T3[11]))))))))) + M1 * (F2[4] * (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] *
      (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (T3[14] * (P1[0] + P2[3] - P2[0] - P1[3]) + (T3[5] * (P2[3] +
      P1[0] - P1[3] - P2[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] *
      (P2[3] - P1[3]))))))))) + F2[5] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P1[1]) + cI
      * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) + cI * (P2[1]) -
      P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P2[1] - P1[1]))))))))))));
  F1[3] = denom * cI * (F2[2] * (P1[0] * (P1[3] * (+cI * (T3[13] + T3[10] +
      T3[16] + T3[4]) - T3[9] - T3[6] - T3[15] - T3[3]) + (P1[1] * (T3[14] +
      T3[5] + cI * (T3[11] + T3[8]) - 2. * (T3[7] + T3[2])) + (P1[2] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12] + T3[2]) + cI * (T3[14] + T3[5])) +
      (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[0] * (+cI *
      (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P1[1]
      * (P1[2] * (T3[4] + T3[10] - cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])
      - T3[16] - T3[13]) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[7]) - cI *
      (T3[11] + T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P2[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[12]) +
      2. * cI * (T3[17])) + (P1[2] * (-cI * (T3[4] + T3[10]) + cI * (T3[16] +
      T3[13])) + (P2[1] * (-cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6])) +
      (P2[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] + T3[10])) + (P2[0] *
      (-2. * cI * (T3[2]) + cI * (T3[14] + T3[5])) + P2[3] * (-2. * cI *
      (T3[17]) + cI * (T3[5] + T3[14]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15]
      - cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[1] *
      (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] + T3[11] - 2. *
      cI * (T3[12]))))))))) + (F2[3] * (P1[0] * (P1[1] * (T3[15] + T3[9] + 2. *
      (T3[3] + T3[6]) + cI * (T3[10] + T3[4])) + (P1[2] * (T3[16] + T3[13] + 2.
      * (T3[4] + T3[10]) - cI * (T3[6] + T3[3])) + (P1[3] * 2. * (T3[5] +
      T3[17] + T3[2] + T3[14]) + (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * -
      1. * (T3[14] + T3[5] + 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] + 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * - 1. * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) + cI * (T3[13] +
      T3[16])) + (P1[2] * - 1. * (T3[4] + T3[10] - cI * (T3[9] + T3[15]) + 2. *
      (T3[13] + T3[16])) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] *
      (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] + 2.
      * (T3[17])) + (P2[0] * - 1. * (T3[14] + T3[5] + 2. * (T3[2])) + P2[3] *
      (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[1] * (P1[2] * - 2. * (T3[8]
      + T3[11] - cI * (T3[7]) + cI * (T3[12])) + (P2[0] * - 1. * (T3[6] + T3[3]
      + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) +
      (P2[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + P2[2] * (T3[8] +
      T3[11] + 2. * cI * (T3[12]))))))) + P1[2] * (P2[0] * (+cI * (T3[6] +
      T3[3]) - T3[10] - T3[4]) + (P2[3] * (T3[13] + T3[16] - cI * (T3[9] +
      T3[15])) + (P1[2] * (+cI * (T3[8] + T3[11]) - 2. * (T3[12])) + (P2[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) - P2[2] * (+cI * (T3[8] + T3[11]) -
      2. * (T3[12]))))))))) + M1 * (F2[4] * (P1[0] * (+cI * (T3[10] + T3[4]) -
      T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) +
      (P2[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI * (P1[1]) +
      cI * (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) -
      P2[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P1[1] - P2[1]))))))))) + F2[5] * (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (T3[14] * (P1[0] + P1[3] - P2[0] - P2[3]) + (T3[5] * (P1[3] + P1[0] -
      P2[3] - P2[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))))));
  F1[4] = denom * cI * (F2[5] * (P1[0] * (P1[3] * - 1. * (T3[9] + T3[6] +
      T3[15] + T3[3] + cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P1[1] * - 1.
      * (+2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8]) - T3[14] - T3[5]) +
      (P1[2] * - 1. * (T3[8] + T3[11] - cI * (T3[14] + T3[5]) + 2. * cI *
      (T3[12] + T3[2])) + (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8])) + P2[2] * (T3[8] + T3[11] + 2. * cI * (T3[12])))))))))
      + (P1[1] * (P1[2] * (T3[4] + T3[10] - cI * (T3[15] + T3[9]) + cI * (T3[3]
      + T3[6]) - T3[16] - T3[13]) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[7]) +
      cI * (T3[11] + T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15]
      - T3[9]) + (P2[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P2[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - 2. * cI * (T3[17]) + cI * (T3[5] + T3[14]) +
      2. * cI * (T3[12])) + (P1[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] +
      T3[10])) + (P2[1] * (-cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])) +
      (P2[2] * (-cI * (T3[4] + T3[10]) + cI * (T3[16] + T3[13])) + (P2[0] * -
      1. * (-2. * cI * (T3[2]) + cI * (T3[14] + T3[5])) - P2[3] * (-2. * cI *
      (T3[17]) + cI * (T3[5] + T3[14]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15]
      + cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] +
      T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8]
      + T3[11] + 2. * cI * (T3[12]))))))))) + (F2[4] * (P1[0] * (P1[1] * - 1. *
      (T3[15] + T3[9] + 2. * (T3[3] + T3[6]) - cI * (T3[10] + T3[4])) + (P1[2]
      * - 1. * (T3[16] + T3[13] + 2. * (T3[4] + T3[10]) + cI * (T3[6] + T3[3]))
      + (P1[3] * - 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (P1[0] * (T3[14] + T3[5] + 2. * (T3[2])) + (P2[0] * - 1. * (T3[14] +
      T3[5] + 2. * (T3[2])) + P2[3] * (T3[5] + T3[14] + 2. * (T3[17]))))))))) +
      (P1[3] * (P1[1] * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) - cI * (T3[13] +
      T3[16])) + (P1[2] * (T3[4] + T3[10] + 2. * (T3[13] + T3[16]) + cI *
      (T3[9] + T3[15])) + (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[3] * (T3[5] +
      T3[14] + 2. * (T3[17])) + (P2[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. *
      (T3[8] + T3[11] - cI * (T3[12]) + cI * (T3[7])) + (P2[0] * (T3[6] + T3[3]
      - cI * (T3[10] + T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      (P2[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] +
      T3[11] - 2. * cI * (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] + cI
      * (T3[6] + T3[3])) + (P2[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] +
      T3[15])) + (P1[2] * (+2. * (T3[12]) + cI * (T3[8] + T3[11])) + (P2[1] * -
      1. * (T3[11] + T3[8] + 2. * cI * (T3[7])) - P2[2] * (+2. * (T3[12]) + cI
      * (T3[8] + T3[11]))))))))) + M1 * (F2[2] * (P1[1] * (T3[15] + T3[9] -
      T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[1] *
      (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] -
      T3[13]) + (T3[14] * (P2[0] + P2[3] - P1[0] - P1[3]) + (T3[5] * (P2[3] +
      P2[0] - P1[3] - P1[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] *
      (P1[3] - P2[3]))))))))) + F2[3] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P2[1]) + cI * (P1[1]) -
      P2[2]) + (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] *
      (P1[1] - P2[1]))))))))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (P1[3] * (T3[6] + T3[3] - cI * (T3[10]
      + T3[4]) + cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) - cI * (T3[11] + T3[8])) + (P1[2]
      * (+cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2]) - T3[8] - T3[11])
      + (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[0] * (+cI *
      (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P1[1]
      * (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13] - cI * (T3[3] + T3[15] +
      T3[6] + T3[9])) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[17]) + cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4]
      + T3[16] + T3[10] + T3[13]) + (P2[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[2] * (P1[3] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17])) + (P1[2] * - 1. * (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) +
      (P2[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P2[2] * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P2[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P1[3] * (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P2[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] *
      (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI *
      (T3[11] + T3[8]) - 2. * (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F2[5] * (P1[0] * (P1[1] * - 1. * (+2. * (T3[3] +
      T3[6]) + cI * (T3[10] + T3[4]) - T3[15] - T3[9]) + (P1[2] * (T3[16] +
      T3[13] + cI * (T3[6] + T3[3]) - 2. * (T3[4] + T3[10])) + (P1[3] * 2. *
      (T3[17] + T3[2] - T3[5] - T3[14]) + (P2[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[0] * - 1. *
      (T3[14] + T3[5] - 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] - 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] - 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * (+2. * (T3[9] + T3[15]) + cI * (T3[13] + T3[16]) - T3[3] -
      T3[6]) + (P1[2] * - 1. * (T3[4] + T3[10] + cI * (T3[9] + T3[15]) - 2. *
      (T3[13] + T3[16])) + (P2[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] - 2.
      * (T3[17])) + (P2[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P2[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. * (T3[8] + T3[11] -
      cI * (T3[7]) + cI * (T3[12])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] +
      T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P2[1] * - 1. * (+2. *
      (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] - cI * (T3[6] + T3[3]))
      + (P2[3] * (+cI * (T3[9] + T3[15]) - T3[13] - T3[16]) + (P1[2] * - 1. *
      (+cI * (T3[8] + T3[11]) - 2. * (T3[12])) + (P2[1] * - 1. * (T3[11] +
      T3[8] - 2. * cI * (T3[7])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M1 * (F2[2] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6]
      - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[0] *
      (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] * (+cI * (T3[13] +
      T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) - P2[2])
      + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] * (P1[1] -
      P2[1]))))))))) + F2[3] * (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (T3[14] * (P1[0] + P2[3] - P2[0] - P1[3]) + (T3[5] * (P2[3] + P1[0] -
      P1[3] - P2[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))))));
}


void VVT6_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP59; 
  complex<double> TMP61; 
  double P1[4]; 
  complex<double> TMP57; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP55; 
  complex<double> denom; 
  complex<double> TMP64; 
  double OM1; 
  complex<double> TMP9; 
  complex<double> TMP63; 
  complex<double> TMP38; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP55 = -1. * (P1[0] * (P2[0] * (P3[1] * (T3[5] * V2[4] - T3[4] * V2[5]) +
      (P3[2] * (T3[3] * V2[5] - T3[5] * V2[3]) + P3[3] * (T3[4] * V2[3] - T3[3]
      * V2[4]))) + (P2[1] * (P3[0] * (T3[4] * V2[5] - T3[5] * V2[4]) + (P3[2] *
      (T3[5] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] - T3[4] *
      V2[2]))) + (P2[2] * (P3[0] * (T3[5] * V2[3] - T3[3] * V2[5]) + (P3[1] *
      (T3[2] * V2[5] - T3[5] * V2[2]) + P3[3] * (T3[3] * V2[2] - T3[2] *
      V2[3]))) + P2[3] * (P3[0] * (T3[3] * V2[4] - T3[4] * V2[3]) + (P3[1] *
      (T3[4] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] - T3[3] *
      V2[2])))))) + (P1[1] * (P2[0] * (P3[1] * (T3[8] * V2[5] - T3[9] * V2[4])
      + (P3[2] * (T3[9] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] * V2[4] -
      T3[8] * V2[3]))) + (P2[1] * (P3[0] * (T3[9] * V2[4] - T3[8] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[9] * V2[2]) + P3[3] * (T3[8] * V2[2] - T3[6]
      * V2[4]))) + (P2[2] * (P3[0] * (T3[7] * V2[5] - T3[9] * V2[3]) + (P3[1] *
      (T3[9] * V2[2] - T3[6] * V2[5]) + P3[3] * (T3[6] * V2[3] - T3[7] *
      V2[2]))) + P2[3] * (P3[0] * (T3[8] * V2[3] - T3[7] * V2[4]) + (P3[1] *
      (T3[6] * V2[4] - T3[8] * V2[2]) + P3[2] * (T3[7] * V2[2] - T3[6] *
      V2[3])))))) + (P1[2] * (P2[0] * (P3[1] * (T3[12] * V2[5] - T3[13] *
      V2[4]) + (P3[2] * (T3[13] * V2[3] - T3[11] * V2[5]) + P3[3] * (T3[11] *
      V2[4] - T3[12] * V2[3]))) + (P2[1] * (P3[0] * (T3[13] * V2[4] - T3[12] *
      V2[5]) + (P3[2] * (T3[10] * V2[5] - T3[13] * V2[2]) + P3[3] * (T3[12] *
      V2[2] - T3[10] * V2[4]))) + (P2[2] * (P3[0] * (T3[11] * V2[5] - T3[13] *
      V2[3]) + (P3[1] * (T3[13] * V2[2] - T3[10] * V2[5]) + P3[3] * (T3[10] *
      V2[3] - T3[11] * V2[2]))) + P2[3] * (P3[0] * (T3[12] * V2[3] - T3[11] *
      V2[4]) + (P3[1] * (T3[10] * V2[4] - T3[12] * V2[2]) + P3[2] * (T3[11] *
      V2[2] - T3[10] * V2[3])))))) + P1[3] * (P2[0] * (P3[1] * (T3[16] * V2[5]
      - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[15] * V2[5]) + P3[3] *
      (T3[15] * V2[4] - T3[16] * V2[3]))) + (P2[1] * (P3[0] * (T3[17] * V2[4] -
      T3[16] * V2[5]) + (P3[2] * (T3[14] * V2[5] - T3[17] * V2[2]) + P3[3] *
      (T3[16] * V2[2] - T3[14] * V2[4]))) + (P2[2] * (P3[0] * (T3[15] * V2[5] -
      T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[14] * V2[5]) + P3[3] *
      (T3[14] * V2[3] - T3[15] * V2[2]))) + P2[3] * (P3[0] * (T3[16] * V2[3] -
      T3[15] * V2[4]) + (P3[1] * (T3[14] * V2[4] - T3[16] * V2[2]) + P3[2] *
      (T3[15] * V2[2] - T3[14] * V2[3])))))))));
  TMP57 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[5] * V2[4] - T3[4] * V2[5]) +
      (P3[2] * (T3[3] * V2[5] - T3[5] * V2[3]) + P3[3] * (T3[4] * V2[3] - T3[3]
      * V2[4]))) + (P1[1] * (P3[0] * (T3[4] * V2[5] - T3[5] * V2[4]) + (P3[2] *
      (T3[5] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] - T3[4] *
      V2[2]))) + (P1[2] * (P3[0] * (T3[5] * V2[3] - T3[3] * V2[5]) + (P3[1] *
      (T3[2] * V2[5] - T3[5] * V2[2]) + P3[3] * (T3[3] * V2[2] - T3[2] *
      V2[3]))) + P1[3] * (P3[0] * (T3[3] * V2[4] - T3[4] * V2[3]) + (P3[1] *
      (T3[4] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] - T3[3] *
      V2[2])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[8] * V2[5] - T3[9] * V2[4])
      + (P3[2] * (T3[9] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] * V2[4] -
      T3[8] * V2[3]))) + (P1[1] * (P3[0] * (T3[9] * V2[4] - T3[8] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[9] * V2[2]) + P3[3] * (T3[8] * V2[2] - T3[6]
      * V2[4]))) + (P1[2] * (P3[0] * (T3[7] * V2[5] - T3[9] * V2[3]) + (P3[1] *
      (T3[9] * V2[2] - T3[6] * V2[5]) + P3[3] * (T3[6] * V2[3] - T3[7] *
      V2[2]))) + P1[3] * (P3[0] * (T3[8] * V2[3] - T3[7] * V2[4]) + (P3[1] *
      (T3[6] * V2[4] - T3[8] * V2[2]) + P3[2] * (T3[7] * V2[2] - T3[6] *
      V2[3])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[12] * V2[5] - T3[13] *
      V2[4]) + (P3[2] * (T3[13] * V2[3] - T3[11] * V2[5]) + P3[3] * (T3[11] *
      V2[4] - T3[12] * V2[3]))) + (P1[1] * (P3[0] * (T3[13] * V2[4] - T3[12] *
      V2[5]) + (P3[2] * (T3[10] * V2[5] - T3[13] * V2[2]) + P3[3] * (T3[12] *
      V2[2] - T3[10] * V2[4]))) + (P1[2] * (P3[0] * (T3[11] * V2[5] - T3[13] *
      V2[3]) + (P3[1] * (T3[13] * V2[2] - T3[10] * V2[5]) + P3[3] * (T3[10] *
      V2[3] - T3[11] * V2[2]))) + P1[3] * (P3[0] * (T3[12] * V2[3] - T3[11] *
      V2[4]) + (P3[1] * (T3[10] * V2[4] - T3[12] * V2[2]) + P3[2] * (T3[11] *
      V2[2] - T3[10] * V2[3])))))) + P2[3] * (P1[0] * (P3[1] * (T3[16] * V2[5]
      - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[15] * V2[5]) + P3[3] *
      (T3[15] * V2[4] - T3[16] * V2[3]))) + (P1[1] * (P3[0] * (T3[17] * V2[4] -
      T3[16] * V2[5]) + (P3[2] * (T3[14] * V2[5] - T3[17] * V2[2]) + P3[3] *
      (T3[16] * V2[2] - T3[14] * V2[4]))) + (P1[2] * (P3[0] * (T3[15] * V2[5] -
      T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[14] * V2[5]) + P3[3] *
      (T3[14] * V2[3] - T3[15] * V2[2]))) + P1[3] * (P3[0] * (T3[16] * V2[3] -
      T3[15] * V2[4]) + (P3[1] * (T3[14] * V2[4] - T3[16] * V2[2]) + P3[2] *
      (T3[15] * V2[2] - T3[14] * V2[3])))))))));
  TMP59 = -1. * (P1[0] * (P2[0] * (P3[1] * (T3[14] * V2[4] - T3[10] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[14] * V2[3]) + P3[3] * (T3[10] * V2[3] -
      T3[6] * V2[4]))) + (P2[1] * (P3[0] * (T3[10] * V2[5] - T3[14] * V2[4]) +
      (P3[2] * (T3[14] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] -
      T3[10] * V2[2]))) + (P2[2] * (P3[0] * (T3[14] * V2[3] - T3[6] * V2[5]) +
      (P3[1] * (T3[2] * V2[5] - T3[14] * V2[2]) + P3[3] * (T3[6] * V2[2] -
      T3[2] * V2[3]))) + P2[3] * (P3[0] * (T3[6] * V2[4] - T3[10] * V2[3]) +
      (P3[1] * (T3[10] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] -
      T3[6] * V2[2])))))) + (P1[1] * (P2[0] * (P3[1] * (T3[11] * V2[5] - T3[15]
      * V2[4]) + (P3[2] * (T3[15] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] *
      V2[4] - T3[11] * V2[3]))) + (P2[1] * (P3[0] * (T3[15] * V2[4] - T3[11] *
      V2[5]) + (P3[2] * (T3[3] * V2[5] - T3[15] * V2[2]) + P3[3] * (T3[11] *
      V2[2] - T3[3] * V2[4]))) + (P2[2] * (P3[0] * (T3[7] * V2[5] - T3[15] *
      V2[3]) + (P3[1] * (T3[15] * V2[2] - T3[3] * V2[5]) + P3[3] * (T3[3] *
      V2[3] - T3[7] * V2[2]))) + P2[3] * (P3[0] * (T3[11] * V2[3] - T3[7] *
      V2[4]) + (P3[1] * (T3[3] * V2[4] - T3[11] * V2[2]) + P3[2] * (T3[7] *
      V2[2] - T3[3] * V2[3])))))) + (P1[2] * (P2[0] * (P3[1] * (T3[12] * V2[5]
      - T3[16] * V2[4]) + (P3[2] * (T3[16] * V2[3] - T3[8] * V2[5]) + P3[3] *
      (T3[8] * V2[4] - T3[12] * V2[3]))) + (P2[1] * (P3[0] * (T3[16] * V2[4] -
      T3[12] * V2[5]) + (P3[2] * (T3[4] * V2[5] - T3[16] * V2[2]) + P3[3] *
      (T3[12] * V2[2] - T3[4] * V2[4]))) + (P2[2] * (P3[0] * (T3[8] * V2[5] -
      T3[16] * V2[3]) + (P3[1] * (T3[16] * V2[2] - T3[4] * V2[5]) + P3[3] *
      (T3[4] * V2[3] - T3[8] * V2[2]))) + P2[3] * (P3[0] * (T3[12] * V2[3] -
      T3[8] * V2[4]) + (P3[1] * (T3[4] * V2[4] - T3[12] * V2[2]) + P3[2] *
      (T3[8] * V2[2] - T3[4] * V2[3])))))) + P1[3] * (P2[0] * (P3[1] * (T3[13]
      * V2[5] - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[9] * V2[5]) +
      P3[3] * (T3[9] * V2[4] - T3[13] * V2[3]))) + (P2[1] * (P3[0] * (T3[17] *
      V2[4] - T3[13] * V2[5]) + (P3[2] * (T3[5] * V2[5] - T3[17] * V2[2]) +
      P3[3] * (T3[13] * V2[2] - T3[5] * V2[4]))) + (P2[2] * (P3[0] * (T3[9] *
      V2[5] - T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[5] * V2[5]) +
      P3[3] * (T3[5] * V2[3] - T3[9] * V2[2]))) + P2[3] * (P3[0] * (T3[13] *
      V2[3] - T3[9] * V2[4]) + (P3[1] * (T3[5] * V2[4] - T3[13] * V2[2]) +
      P3[2] * (T3[9] * V2[2] - T3[5] * V2[3])))))))));
  TMP61 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[14] * V2[4] - T3[10] * V2[5]) +
      (P3[2] * (T3[6] * V2[5] - T3[14] * V2[3]) + P3[3] * (T3[10] * V2[3] -
      T3[6] * V2[4]))) + (P1[1] * (P3[0] * (T3[10] * V2[5] - T3[14] * V2[4]) +
      (P3[2] * (T3[14] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[4] -
      T3[10] * V2[2]))) + (P1[2] * (P3[0] * (T3[14] * V2[3] - T3[6] * V2[5]) +
      (P3[1] * (T3[2] * V2[5] - T3[14] * V2[2]) + P3[3] * (T3[6] * V2[2] -
      T3[2] * V2[3]))) + P1[3] * (P3[0] * (T3[6] * V2[4] - T3[10] * V2[3]) +
      (P3[1] * (T3[10] * V2[2] - T3[2] * V2[4]) + P3[2] * (T3[2] * V2[3] -
      T3[6] * V2[2])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[11] * V2[5] - T3[15]
      * V2[4]) + (P3[2] * (T3[15] * V2[3] - T3[7] * V2[5]) + P3[3] * (T3[7] *
      V2[4] - T3[11] * V2[3]))) + (P1[1] * (P3[0] * (T3[15] * V2[4] - T3[11] *
      V2[5]) + (P3[2] * (T3[3] * V2[5] - T3[15] * V2[2]) + P3[3] * (T3[11] *
      V2[2] - T3[3] * V2[4]))) + (P1[2] * (P3[0] * (T3[7] * V2[5] - T3[15] *
      V2[3]) + (P3[1] * (T3[15] * V2[2] - T3[3] * V2[5]) + P3[3] * (T3[3] *
      V2[3] - T3[7] * V2[2]))) + P1[3] * (P3[0] * (T3[11] * V2[3] - T3[7] *
      V2[4]) + (P3[1] * (T3[3] * V2[4] - T3[11] * V2[2]) + P3[2] * (T3[7] *
      V2[2] - T3[3] * V2[3])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[12] * V2[5]
      - T3[16] * V2[4]) + (P3[2] * (T3[16] * V2[3] - T3[8] * V2[5]) + P3[3] *
      (T3[8] * V2[4] - T3[12] * V2[3]))) + (P1[1] * (P3[0] * (T3[16] * V2[4] -
      T3[12] * V2[5]) + (P3[2] * (T3[4] * V2[5] - T3[16] * V2[2]) + P3[3] *
      (T3[12] * V2[2] - T3[4] * V2[4]))) + (P1[2] * (P3[0] * (T3[8] * V2[5] -
      T3[16] * V2[3]) + (P3[1] * (T3[16] * V2[2] - T3[4] * V2[5]) + P3[3] *
      (T3[4] * V2[3] - T3[8] * V2[2]))) + P1[3] * (P3[0] * (T3[12] * V2[3] -
      T3[8] * V2[4]) + (P3[1] * (T3[4] * V2[4] - T3[12] * V2[2]) + P3[2] *
      (T3[8] * V2[2] - T3[4] * V2[3])))))) + P2[3] * (P1[0] * (P3[1] * (T3[13]
      * V2[5] - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[3] - T3[9] * V2[5]) +
      P3[3] * (T3[9] * V2[4] - T3[13] * V2[3]))) + (P1[1] * (P3[0] * (T3[17] *
      V2[4] - T3[13] * V2[5]) + (P3[2] * (T3[5] * V2[5] - T3[17] * V2[2]) +
      P3[3] * (T3[13] * V2[2] - T3[5] * V2[4]))) + (P1[2] * (P3[0] * (T3[9] *
      V2[5] - T3[17] * V2[3]) + (P3[1] * (T3[17] * V2[2] - T3[5] * V2[5]) +
      P3[3] * (T3[5] * V2[3] - T3[9] * V2[2]))) + P1[3] * (P3[0] * (T3[13] *
      V2[3] - T3[9] * V2[4]) + (P3[1] * (T3[5] * V2[4] - T3[13] * V2[2]) +
      P3[2] * (T3[9] * V2[2] - T3[5] * V2[3])))))))));
  TMP63 = -1. * (P1[0] * (P1[1] * (P2[2] * (P3[3] * - 1. * (T3[7] + T3[2]) +
      (P3[0] * T3[5] + P3[1] * T3[9])) + (P2[3] * (P3[2] * (T3[7] + T3[2]) +
      (-P3[0] * T3[4] - P3[1] * T3[8])) + (P2[0] * (P3[3] * T3[4] - P3[2] *
      T3[5]) + P2[1] * (P3[3] * T3[8] - P3[2] * T3[9])))) + (P1[2] * (P2[1] *
      (P3[3] * (T3[12] + T3[2]) + (-P3[0] * T3[5] - P3[2] * T3[13])) + (P2[3] *
      (P3[1] * - 1. * (T3[12] + T3[2]) + (P3[0] * T3[3] + P3[2] * T3[11])) +
      (P2[0] * (P3[1] * T3[5] - P3[3] * T3[3]) + P2[2] * (P3[1] * T3[13] -
      P3[3] * T3[11])))) + (P1[3] * (P2[1] * (P3[2] * - 1. * (T3[17] + T3[2]) +
      (P3[0] * T3[4] + P3[3] * T3[16])) + (P2[2] * (P3[1] * (T3[17] + T3[2]) +
      (-P3[0] * T3[3] - P3[3] * T3[15])) + (P2[0] * (P3[2] * T3[3] - P3[1] *
      T3[4]) + P2[3] * (P3[2] * T3[15] - P3[1] * T3[16])))) + P1[0] * (P2[1] *
      (P3[2] * T3[5] - P3[3] * T3[4]) + (P2[2] * (P3[3] * T3[3] - P3[1] *
      T3[5]) + P2[3] * (P3[1] * T3[4] - P3[2] * T3[3])))))) + (P1[1] * (P1[2] *
      (P2[0] * (P3[3] * (T3[7] - T3[12]) + (P3[2] * T3[13] - P3[1] * T3[9])) +
      (P2[3] * (P3[0] * (T3[12] - T3[7]) + (P3[1] * T3[6] - P3[2] * T3[10])) +
      (P2[1] * (P3[0] * T3[9] - P3[3] * T3[6]) + P2[2] * (P3[3] * T3[10] -
      P3[0] * T3[13])))) + (P1[3] * (P2[0] * (P3[2] * (T3[17] - T3[7]) + (P3[1]
      * T3[8] - P3[3] * T3[16])) + (P2[2] * (P3[0] * (T3[7] - T3[17]) + (P3[3]
      * T3[14] - P3[1] * T3[6])) + (P2[1] * (P3[2] * T3[6] - P3[0] * T3[8]) +
      P2[3] * (P3[0] * T3[16] - P3[2] * T3[14])))) + P1[1] * (P2[0] * (P3[2] *
      T3[9] - P3[3] * T3[8]) + (P2[2] * (P3[3] * T3[6] - P3[0] * T3[9]) + P2[3]
      * (P3[0] * T3[8] - P3[2] * T3[6]))))) + (P1[2] * (P1[3] * (P2[0] * (P3[1]
      * (T3[12] - T3[17]) + (P3[3] * T3[15] - P3[2] * T3[11])) + (P2[1] *
      (P3[0] * (T3[17] - T3[12]) + (P3[2] * T3[10] - P3[3] * T3[14])) + (P2[2]
      * (P3[0] * T3[11] - P3[1] * T3[10]) + P2[3] * (P3[1] * T3[14] - P3[0] *
      T3[15])))) + P1[2] * (P2[0] * (P3[3] * T3[11] - P3[1] * T3[13]) + (P2[1]
      * (P3[0] * T3[13] - P3[3] * T3[10]) + P2[3] * (P3[1] * T3[10] - P3[0] *
      T3[11])))) + P1[3] * P1[3] * (P2[0] * (P3[1] * T3[16] - P3[2] * T3[15]) +
      (P2[1] * (P3[2] * T3[14] - P3[0] * T3[16]) + P2[2] * (P3[0] * T3[15] -
      P3[1] * T3[14]))))));
  TMP64 = -1. * (P1[0] * (P1[1] * (P2[2] * (P3[3] * - 1. * (T3[7] + T3[2]) +
      (P3[0] * T3[14] + P3[1] * T3[15])) + (P2[3] * (P3[2] * (T3[7] + T3[2]) +
      (-P3[0] * T3[10] - P3[1] * T3[11])) + (P2[0] * (P3[3] * T3[10] - P3[2] *
      T3[14]) + P2[1] * (P3[3] * T3[11] - P3[2] * T3[15])))) + (P1[2] * (P2[1]
      * (P3[3] * (T3[12] + T3[2]) + (-P3[0] * T3[14] - P3[2] * T3[16])) +
      (P2[3] * (P3[1] * - 1. * (T3[12] + T3[2]) + (P3[0] * T3[6] + P3[2] *
      T3[8])) + (P2[0] * (P3[1] * T3[14] - P3[3] * T3[6]) + P2[2] * (P3[1] *
      T3[16] - P3[3] * T3[8])))) + (P1[3] * (P2[1] * (P3[2] * - 1. * (T3[17] +
      T3[2]) + (P3[0] * T3[10] + P3[3] * T3[13])) + (P2[2] * (P3[1] * (T3[17] +
      T3[2]) + (-P3[0] * T3[6] - P3[3] * T3[9])) + (P2[0] * (P3[2] * T3[6] -
      P3[1] * T3[10]) + P2[3] * (P3[2] * T3[9] - P3[1] * T3[13])))) + P1[0] *
      (P2[1] * (P3[2] * T3[14] - P3[3] * T3[10]) + (P2[2] * (P3[3] * T3[6] -
      P3[1] * T3[14]) + P2[3] * (P3[1] * T3[10] - P3[2] * T3[6])))))) + (P1[1]
      * (P1[2] * (P2[0] * (P3[3] * (T3[7] - T3[12]) + (P3[2] * T3[16] - P3[1] *
      T3[15])) + (P2[3] * (P3[0] * (T3[12] - T3[7]) + (P3[1] * T3[3] - P3[2] *
      T3[4])) + (P2[1] * (P3[0] * T3[15] - P3[3] * T3[3]) + P2[2] * (P3[3] *
      T3[4] - P3[0] * T3[16])))) + (P1[3] * (P2[0] * (P3[2] * (T3[17] - T3[7])
      + (P3[1] * T3[11] - P3[3] * T3[13])) + (P2[2] * (P3[0] * (T3[7] - T3[17])
      + (P3[3] * T3[5] - P3[1] * T3[3])) + (P2[1] * (P3[2] * T3[3] - P3[0] *
      T3[11]) + P2[3] * (P3[0] * T3[13] - P3[2] * T3[5])))) + P1[1] * (P2[0] *
      (P3[2] * T3[15] - P3[3] * T3[11]) + (P2[2] * (P3[3] * T3[3] - P3[0] *
      T3[15]) + P2[3] * (P3[0] * T3[11] - P3[2] * T3[3]))))) + (P1[2] * (P1[3]
      * (P2[0] * (P3[1] * (T3[12] - T3[17]) + (P3[3] * T3[9] - P3[2] * T3[8]))
      + (P2[1] * (P3[0] * (T3[17] - T3[12]) + (P3[2] * T3[4] - P3[3] * T3[5]))
      + (P2[2] * (P3[0] * T3[8] - P3[1] * T3[4]) + P2[3] * (P3[1] * T3[5] -
      P3[0] * T3[9])))) + P1[2] * (P2[0] * (P3[3] * T3[8] - P3[1] * T3[16]) +
      (P2[1] * (P3[0] * T3[16] - P3[3] * T3[4]) + P2[3] * (P3[1] * T3[4] -
      P3[0] * T3[8])))) + P1[3] * P1[3] * (P2[0] * (P3[1] * T3[13] - P3[2] *
      T3[9]) + (P2[1] * (P3[2] * T3[5] - P3[0] * T3[13]) + P2[2] * (P3[0] *
      T3[9] - P3[1] * T3[5]))))));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP38 * (P3[1] * (P1[2] * (P2[0] * - 1. * (+cI * (T3[5] +
      T3[14])) + (P2[1] * (+cI * (T3[9] + T3[15])) + (P2[2] * 2. * (+cI *
      (T3[13] + T3[16])) + 2. * (P2[3] * (-cI * (T3[12]) + cI * (T3[17])))))) +
      (P1[3] * (P2[0] * (+cI * (T3[4] + T3[10])) + (P2[1] * - 1. * (+cI *
      (T3[8] + T3[11])) + (P2[2] * 2. * (-cI * (T3[12]) + cI * (T3[17])) - 2. *
      (P2[3] * (+cI * (T3[16] + T3[13])))))) + (P1[0] * (P2[2] * - 1. * (+cI *
      (T3[5] + T3[14])) + P2[3] * (+cI * (T3[4] + T3[10]))) + P1[1] * (P2[2] *
      (+cI * (T3[9] + T3[15])) - P2[3] * (+cI * (T3[8] + T3[11])))))) + (P3[2]
      * (P1[1] * (P2[0] * (+cI * (T3[5] + T3[14])) + (P2[1] * - 2. * (+cI *
      (T3[9] + T3[15])) + (P2[2] * - 1. * (+cI * (T3[13] + T3[16])) + 2. *
      (P2[3] * (-cI * (T3[17]) + cI * (T3[7])))))) + (P1[3] * (P2[0] * - 1. *
      (+cI * (T3[3] + T3[6])) + (P2[1] * 2. * (-cI * (T3[17]) + cI * (T3[7])) +
      (P2[2] * (+cI * (T3[11] + T3[8])) + 2. * (P2[3] * (+cI * (T3[15] +
      T3[9])))))) + (P1[0] * (P2[1] * (+cI * (T3[5] + T3[14])) - P2[3] * (+cI *
      (T3[3] + T3[6]))) + P1[2] * (P2[1] * - 1. * (+cI * (T3[13] + T3[16])) +
      P2[3] * (+cI * (T3[11] + T3[8])))))) + (P3[3] * (P1[1] * (P2[0] * - 1. *
      (+cI * (T3[4] + T3[10])) + (P2[1] * 2. * (+cI * (T3[8] + T3[11])) +
      (P2[2] * 2. * (-cI * (T3[7]) + cI * (T3[12])) + P2[3] * (+cI * (T3[16] +
      T3[13]))))) + (P1[2] * (P2[0] * (+cI * (T3[3] + T3[6])) + (P2[1] * 2. *
      (-cI * (T3[7]) + cI * (T3[12])) + (P2[2] * - 2. * (+cI * (T3[11] +
      T3[8])) - P2[3] * (+cI * (T3[15] + T3[9]))))) + (P1[0] * (P2[1] * - 1. *
      (+cI * (T3[4] + T3[10])) + P2[2] * (+cI * (T3[3] + T3[6]))) + P1[3] *
      (P2[1] * (+cI * (T3[16] + T3[13])) - P2[2] * (+cI * (T3[15] + T3[9]))))))
      + OM1 * P1[0] * (+cI * (TMP63 + TMP64))))) + (OM1 * P1[0] * TMP9 * (+cI *
      (TMP55 + TMP57 + TMP59 + TMP61)) - P3[0] * (+cI * (TMP55 + TMP57 + TMP59
      + TMP61))));
  V1[3] = denom * (TMP38 * (P3[0] * (P1[2] * (P2[0] * - 1. * (+cI * (T3[5] +
      T3[14])) + (P2[1] * (+cI * (T3[9] + T3[15])) + (P2[2] * 2. * (+cI *
      (T3[13] + T3[16])) + 2. * (P2[3] * (-cI * (T3[12]) + cI * (T3[17])))))) +
      (P1[3] * (P2[0] * (+cI * (T3[4] + T3[10])) + (P2[1] * - 1. * (+cI *
      (T3[8] + T3[11])) + (P2[2] * 2. * (-cI * (T3[12]) + cI * (T3[17])) - 2. *
      (P2[3] * (+cI * (T3[16] + T3[13])))))) + (P1[0] * (P2[2] * - 1. * (+cI *
      (T3[5] + T3[14])) + P2[3] * (+cI * (T3[4] + T3[10]))) + P1[1] * (P2[2] *
      (+cI * (T3[9] + T3[15])) - P2[3] * (+cI * (T3[8] + T3[11])))))) + (P3[2]
      * (P1[0] * (P2[0] * 2. * (+cI * (T3[5] + T3[14])) + (P2[1] * - 1. * (+cI
      * (T3[9] + T3[15])) + (P2[2] * - 1. * (+cI * (T3[13] + T3[16])) - 2. *
      (P2[3] * (+cI * (T3[2] + T3[17])))))) + (P1[3] * (P2[0] * - 2. * (+cI *
      (T3[17] + T3[2])) + (P2[1] * (+cI * (T3[6] + T3[3])) + (P2[2] * (+cI *
      (T3[10] + T3[4])) + 2. * (P2[3] * (+cI * (T3[14] + T3[5])))))) + (P1[1] *
      (P2[0] * - 1. * (+cI * (T3[9] + T3[15])) + P2[3] * (+cI * (T3[6] +
      T3[3]))) + P1[2] * (P2[0] * - 1. * (+cI * (T3[13] + T3[16])) + P2[3] *
      (+cI * (T3[10] + T3[4])))))) + (P3[3] * (P1[0] * (P2[0] * - 2. * (+cI *
      (T3[4] + T3[10])) + (P2[1] * (+cI * (T3[8] + T3[11])) + (P2[2] * 2. *
      (+cI * (T3[2] + T3[12])) + P2[3] * (+cI * (T3[16] + T3[13]))))) + (P1[2]
      * (P2[0] * 2. * (+cI * (T3[12] + T3[2])) + (P2[1] * - 1. * (+cI * (T3[6]
      + T3[3])) + (P2[2] * - 2. * (+cI * (T3[10] + T3[4])) - P2[3] * (+cI *
      (T3[14] + T3[5]))))) + (P1[1] * (P2[0] * (+cI * (T3[8] + T3[11])) - P2[2]
      * (+cI * (T3[6] + T3[3]))) + P1[3] * (P2[0] * (+cI * (T3[16] + T3[13])) -
      P2[2] * (+cI * (T3[14] + T3[5])))))) + OM1 * P1[1] * (+cI * (TMP63 +
      TMP64))))) + (OM1 * P1[1] * TMP9 * (+cI * (TMP55 + TMP57 + TMP59 +
      TMP61)) - P3[1] * (+cI * (TMP55 + TMP57 + TMP59 + TMP61))));
  V1[4] = denom * (TMP38 * (P3[0] * (P1[1] * (P2[0] * (+cI * (T3[5] + T3[14]))
      + (P2[1] * - 2. * (+cI * (T3[9] + T3[15])) + (P2[2] * - 1. * (+cI *
      (T3[13] + T3[16])) + 2. * (P2[3] * (-cI * (T3[17]) + cI * (T3[7])))))) +
      (P1[3] * (P2[0] * - 1. * (+cI * (T3[3] + T3[6])) + (P2[1] * 2. * (-cI *
      (T3[17]) + cI * (T3[7])) + (P2[2] * (+cI * (T3[11] + T3[8])) + 2. *
      (P2[3] * (+cI * (T3[15] + T3[9])))))) + (P1[0] * (P2[1] * (+cI * (T3[5] +
      T3[14])) - P2[3] * (+cI * (T3[3] + T3[6]))) + P1[2] * (P2[1] * - 1. *
      (+cI * (T3[13] + T3[16])) + P2[3] * (+cI * (T3[11] + T3[8])))))) + (P3[1]
      * (P1[0] * (P2[0] * - 2. * (+cI * (T3[5] + T3[14])) + (P2[1] * (+cI *
      (T3[9] + T3[15])) + (P2[2] * (+cI * (T3[13] + T3[16])) + 2. * (P2[3] *
      (+cI * (T3[2] + T3[17])))))) + (P1[3] * (P2[0] * 2. * (+cI * (T3[17] +
      T3[2])) + (P2[1] * - 1. * (+cI * (T3[6] + T3[3])) + (P2[2] * - 1. * (+cI
      * (T3[10] + T3[4])) - 2. * (P2[3] * (+cI * (T3[14] + T3[5])))))) + (P1[1]
      * (P2[0] * (+cI * (T3[9] + T3[15])) - P2[3] * (+cI * (T3[6] + T3[3]))) +
      P1[2] * (P2[0] * (+cI * (T3[13] + T3[16])) - P2[3] * (+cI * (T3[10] +
      T3[4])))))) + (P3[3] * (P1[0] * (P2[0] * 2. * (+cI * (T3[3] + T3[6])) +
      (P2[1] * - 2. * (+cI * (T3[2] + T3[7])) + (P2[2] * - 1. * (+cI * (T3[11]
      + T3[8])) - P2[3] * (+cI * (T3[15] + T3[9]))))) + (P1[1] * (P2[0] * - 2.
      * (+cI * (T3[7] + T3[2])) + (P2[1] * 2. * (+cI * (T3[6] + T3[3])) +
      (P2[2] * (+cI * (T3[10] + T3[4])) + P2[3] * (+cI * (T3[14] + T3[5]))))) +
      (P1[2] * (P2[0] * - 1. * (+cI * (T3[11] + T3[8])) + P2[1] * (+cI *
      (T3[10] + T3[4]))) + P1[3] * (P2[0] * - 1. * (+cI * (T3[15] + T3[9])) +
      P2[1] * (+cI * (T3[14] + T3[5])))))) + OM1 * P1[2] * (+cI * (TMP63 +
      TMP64))))) + (OM1 * P1[2] * TMP9 * (+cI * (TMP55 + TMP57 + TMP59 +
      TMP61)) - P3[2] * (+cI * (TMP55 + TMP57 + TMP59 + TMP61))));
  V1[5] = denom * (TMP38 * (P3[0] * (P1[1] * (P2[0] * - 1. * (+cI * (T3[4] +
      T3[10])) + (P2[1] * 2. * (+cI * (T3[8] + T3[11])) + (P2[2] * 2. * (-cI *
      (T3[7]) + cI * (T3[12])) + P2[3] * (+cI * (T3[16] + T3[13]))))) + (P1[2]
      * (P2[0] * (+cI * (T3[3] + T3[6])) + (P2[1] * 2. * (-cI * (T3[7]) + cI *
      (T3[12])) + (P2[2] * - 2. * (+cI * (T3[11] + T3[8])) - P2[3] * (+cI *
      (T3[15] + T3[9]))))) + (P1[0] * (P2[1] * - 1. * (+cI * (T3[4] + T3[10]))
      + P2[2] * (+cI * (T3[3] + T3[6]))) + P1[3] * (P2[1] * (+cI * (T3[16] +
      T3[13])) - P2[2] * (+cI * (T3[15] + T3[9])))))) + (P3[1] * (P1[0] *
      (P2[0] * 2. * (+cI * (T3[4] + T3[10])) + (P2[1] * - 1. * (+cI * (T3[8] +
      T3[11])) + (P2[2] * - 2. * (+cI * (T3[2] + T3[12])) - P2[3] * (+cI *
      (T3[16] + T3[13]))))) + (P1[2] * (P2[0] * - 2. * (+cI * (T3[12] + T3[2]))
      + (P2[1] * (+cI * (T3[6] + T3[3])) + (P2[2] * 2. * (+cI * (T3[10] +
      T3[4])) + P2[3] * (+cI * (T3[14] + T3[5]))))) + (P1[1] * (P2[0] * - 1. *
      (+cI * (T3[8] + T3[11])) + P2[2] * (+cI * (T3[6] + T3[3]))) + P1[3] *
      (P2[0] * - 1. * (+cI * (T3[16] + T3[13])) + P2[2] * (+cI * (T3[14] +
      T3[5])))))) + (P3[2] * (P1[0] * (P2[0] * - 2. * (+cI * (T3[3] + T3[6])) +
      (P2[1] * 2. * (+cI * (T3[2] + T3[7])) + (P2[2] * (+cI * (T3[11] + T3[8]))
      + P2[3] * (+cI * (T3[15] + T3[9]))))) + (P1[1] * (P2[0] * 2. * (+cI *
      (T3[7] + T3[2])) + (P2[1] * - 2. * (+cI * (T3[6] + T3[3])) + (P2[2] * -
      1. * (+cI * (T3[10] + T3[4])) - P2[3] * (+cI * (T3[14] + T3[5]))))) +
      (P1[2] * (P2[0] * (+cI * (T3[11] + T3[8])) - P2[1] * (+cI * (T3[10] +
      T3[4]))) + P1[3] * (P2[0] * (+cI * (T3[15] + T3[9])) - P2[1] * (+cI *
      (T3[14] + T3[5])))))) + OM1 * P1[3] * (+cI * (TMP63 + TMP64))))) + (OM1 *
      P1[3] * TMP9 * (+cI * (TMP55 + TMP57 + TMP59 + TMP61)) - P3[3] * (+cI *
      (TMP55 + TMP57 + TMP59 + TMP61))));
}


void VVT10_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> TMP47; 
  complex<double> TMP9; 
  complex<double> TMP25; 
  complex<double> TMP37; 
  double P2[4]; 
  complex<double> TMP46; 
  double OM3; 
  complex<double> TMP49; 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP30; 
  complex<double> TMP16; 
  complex<double> denom; 
  complex<double> TMP48; 
  complex<double> TMP38; 
  complex<double> TMP26; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP25 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP30 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP46 = (P1[0] * P1[0] - P1[1] * P1[1] - P1[2] * P1[2] - P1[3] * P1[3]); 
  TMP47 = (P2[0] * P2[0] - P2[1] * P2[1] - P2[2] * P2[2] - P2[3] * P2[3]); 
  TMP48 = (P1[0] * V1[2] - P1[1] * V1[3] - P1[2] * V1[4] - P1[3] * V1[5]); 
  TMP49 = (P2[0] * V2[2] - P2[1] * V2[3] - P2[2] * V2[4] - P2[3] * V2[5]); 
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 1./3. * (OM3 * (P3[0] * (P3[0] * (OM3 * (TMP12 * (TMP9 *
      (TMP25 * - 2. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 2.
      * (-2. * cI * (TMP30) + cI * (TMP37)) + 2. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2. * (+cI * (TMP30 + TMP37)) - 2. * cI * (TMP16 *
      TMP25)) - 2. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 *
      2. * (+cI * (TMP26 + TMP38)) - 2. * cI * (TMP16 * TMP25)) - 2. * cI *
      (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1. * (-2. * cI * (TMP16)
      + cI * (TMP46 + TMP47)) + (-2. * cI * (TMP26 * TMP30) - cI * (TMP38 *
      TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * (+cI * (TMP46 + TMP47)) +
      (+cI * (TMP37 * TMP47 + TMP12 * TMP48))) + (TMP25 * - 1. * (+cI * (TMP12
      * TMP46 + TMP9 * TMP47)) + TMP30 * (+cI * (TMP38 * TMP46 + TMP9 *
      TMP49)))))) + (TMP12 * (TMP26 * (P1[0] * 3. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (P2[0] * - 6. * (+cI * (TMP30 + TMP37)) - 3. * cI * (TMP9 *
      V1[2]))) + (TMP25 * (P1[0] * 6. * (-cI * (TMP16) + cI * (TMP9)) + 6. *
      (P2[0] * (+cI * (TMP16 + TMP9)))) + 3. * (V2[2] * (-cI * (TMP9 * TMP30) +
      cI * (TMP16 * TMP37))))) + (TMP9 * (TMP30 * (P1[0] * - 6. * (+cI * (TMP26
      + TMP38)) + 3. * (P2[0] * (-cI * (TMP38) + 2. * cI * (TMP26)))) + TMP16 *
      (TMP25 * 6. * (-cI * (P2[0]) + cI * (P1[0])) + 3. * cI * (V1[2] *
      TMP38))) + 3. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] + P2[0])))))) +
      (TMP12 * (TMP9 * (TMP25 * - 1. * (-2. * cI * (TMP16) + cI * (TMP9 +
      TMP12)) + (TMP26 * (-2. * cI * (TMP30) + cI * (TMP37)) + cI * (TMP30 *
      TMP38))) + (TMP12 * (TMP26 * (+cI * (TMP30 + TMP37)) - cI * (TMP16 *
      TMP25)) - cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 * (+cI *
      (TMP26 + TMP38)) - cI * (TMP16 * TMP25)) - cI * (TMP16 * TMP37 *
      TMP38)))) + (TMP16 * (TMP25 * (P1[0] * 3. * (-cI * (P1[0]) + 2. * cI *
      (P2[0])) + (+1./1. * cI * (TMP46 + TMP47) - 2./1. * cI * (TMP16) - 3. *
      cI * (P2[0] * P2[0]))) + (TMP37 * 1./1. * (+cI * (TMP49) - 3./1. * cI *
      (P2[0] * V2[2])) + (+2./1. * cI * (TMP26 * TMP30) + TMP38 * 1./1. * (+cI
      * (TMP48) - 3./1. * cI * (P1[0] * V1[2]))))) + (TMP26 * (TMP30 * (P1[0] *
      3. * (-2. * cI * (P2[0]) + cI * (P1[0])) + (-1./1. * cI * (TMP46 + TMP47)
      + 3. * cI * (P2[0] * P2[0]))) + (+1./1. * (TMP12 * (-cI * (TMP48) + 3./1.
      * cI * (P1[0] * V1[2]))) + TMP37 * 1./1. * (-cI * (TMP47) + 3./1. * cI *
      (P2[0] * P2[0])))) + (TMP9 * 1./1. * (TMP30 * (-cI * (TMP49) + 3./1. * cI
      * (P2[0] * V2[2])) + 1./1. * (TMP25 * 1./1. * (+cI * (TMP47) - 3./1. * cI
      * (P2[0] * P2[0])))) + (+1./1. * (TMP30 * TMP38 * (-cI * (TMP46) + 3./1.
      * cI * (P1[0] * P1[0]))) + TMP12 * 1./1. * TMP25 * (+cI * (TMP46) - 3./1.
      * cI * (P1[0] * P1[0])))))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[1] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[3]) + P2[1] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[3] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[1] * (-cI * (TMP16) + cI * (TMP9)) + P2[1] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] + P2[1]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[1]) + cI * (P1[1]))) +
      cI * (V1[3] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[1] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[1] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[1] * (TMP12 * (TMP26 * (P1[0] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[2]) + P2[0] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[2] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[0] * (-cI * (TMP16) + cI * (TMP9)) + P2[0] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] +
      P2[0]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[0]) + cI *
      (P1[0]))) + cI * (V1[2] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[0] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[0] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[0] * (P1[1] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[3] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[1] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[0] * (P2[1] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[3] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[1] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[1] * V2[2] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[1] * 1./2. * V1[2] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[2] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[4]) + P2[2] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[4] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[2] * (-cI * (TMP16) + cI * (TMP9)) + P2[2] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] + P2[2]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[2]) + cI * (P1[2]))) +
      cI * (V1[4] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[2] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[2] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[2] * (TMP12 * (TMP26 * (P1[0] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[2]) + P2[0] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[2] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[0] * (-cI * (TMP16) + cI * (TMP9)) + P2[0] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] +
      P2[0]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[0]) + cI *
      (P1[0]))) + cI * (V1[2] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[0] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[0] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[0] * (P1[2] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[4] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[2] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[0] * (P2[2] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[4] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[2] * V2[2] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * 1./2. * V1[2] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[3] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[5]) + P2[3] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[5] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[3] * (-cI * (TMP16) + cI * (TMP9)) + P2[3] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[3]) + cI * (P1[3]))) +
      cI * (V1[5] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[3] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[3] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[3] * (TMP12 * (TMP26 * (P1[0] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[2]) + P2[0] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[2] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[0] * (-cI * (TMP16) + cI * (TMP9)) + P2[0] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] +
      P2[0]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[0]) + cI *
      (P1[0]))) + cI * (V1[2] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[0] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[0] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[0] * (P1[3] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[5] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[3] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[0] * (P2[3] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[5] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[3] * V2[2] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * 1./2. * V1[2] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[1] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[3]) + P2[1] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[3] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[1] * (-cI * (TMP16) + cI * (TMP9)) + P2[1] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] + P2[1]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[1]) + cI * (P1[1]))) +
      cI * (V1[3] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[1] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[1] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[1] * (TMP12 * (TMP26 * (P1[0] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[2]) + P2[0] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[2] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[0] * (-cI * (TMP16) + cI * (TMP9)) + P2[0] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] +
      P2[0]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[0]) + cI *
      (P1[0]))) + cI * (V1[2] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[0] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[0] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[0] * (P1[1] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[3] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[1] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[0] * (P2[1] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[3] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[1] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[1] * V2[2] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[1] * 1./2. * V1[2] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[7] = denom * 1./3. * (OM3 * (P3[1] * (P3[1] * (OM3 * (TMP12 * (TMP9 *
      (TMP25 * - 2. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 2.
      * (-2. * cI * (TMP30) + cI * (TMP37)) + 2. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2. * (+cI * (TMP30 + TMP37)) - 2. * cI * (TMP16 *
      TMP25)) - 2. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 *
      2. * (+cI * (TMP26 + TMP38)) - 2. * cI * (TMP16 * TMP25)) - 2. * cI *
      (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1. * (-2. * cI * (TMP16)
      + cI * (TMP46 + TMP47)) + (-2. * cI * (TMP26 * TMP30) - cI * (TMP38 *
      TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * (+cI * (TMP46 + TMP47)) +
      (+cI * (TMP37 * TMP47 + TMP12 * TMP48))) + (TMP25 * - 1. * (+cI * (TMP12
      * TMP46 + TMP9 * TMP47)) + TMP30 * (+cI * (TMP38 * TMP46 + TMP9 *
      TMP49)))))) + (TMP12 * (TMP26 * (P1[1] * 3. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (P2[1] * - 6. * (+cI * (TMP30 + TMP37)) - 3. * cI * (TMP9 *
      V1[3]))) + (TMP25 * (P1[1] * 6. * (-cI * (TMP16) + cI * (TMP9)) + 6. *
      (P2[1] * (+cI * (TMP16 + TMP9)))) + 3. * (V2[3] * (-cI * (TMP9 * TMP30) +
      cI * (TMP16 * TMP37))))) + (TMP9 * (TMP30 * (P1[1] * - 6. * (+cI * (TMP26
      + TMP38)) + 3. * (P2[1] * (-cI * (TMP38) + 2. * cI * (TMP26)))) + TMP16 *
      (TMP25 * 6. * (-cI * (P2[1]) + cI * (P1[1])) + 3. * cI * (V1[3] *
      TMP38))) + 3. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] + P2[1])))))) +
      (TMP12 * (TMP9 * (TMP25 * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) +
      (TMP26 * (-cI * (TMP37) + 2. * cI * (TMP30)) - cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * - 1. * (+cI * (TMP30 + TMP37)) + cI * (TMP16 * TMP25))
      + cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 * - 1. * (+cI *
      (TMP26 + TMP38)) + cI * (TMP16 * TMP25)) + cI * (TMP16 * TMP37 *
      TMP38)))) + (TMP16 * (TMP25 * (P1[1] * 3. * (-cI * (P1[1]) + 2. * cI *
      (P2[1])) + (-1./1. * cI * (TMP46 + TMP47) + 2./1. * cI * (TMP16) - 3. *
      cI * (P2[1] * P2[1]))) + (TMP37 * - 1./1. * (+cI * (TMP49) + 3./1. * cI *
      (P2[1] * V2[3])) + (-2./1. * cI * (TMP26 * TMP30) + TMP38 * - 1./1. *
      (+cI * (TMP48) + 3./1. * cI * (P1[1] * V1[3]))))) + (TMP26 * (TMP30 *
      (P1[1] * 3. * (-2. * cI * (P2[1]) + cI * (P1[1])) + (+1./1. * cI * (TMP46
      + TMP47) + 3. * cI * (P2[1] * P2[1]))) + (+1./1. * (TMP12 * (+cI *
      (TMP48) + 3./1. * cI * (P1[1] * V1[3]))) + TMP37 * 1./1. * (+cI * (TMP47)
      + 3./1. * cI * (P2[1] * P2[1])))) + (TMP9 * 1./1. * (TMP30 * (+cI *
      (TMP49) + 3./1. * cI * (P2[1] * V2[3])) + 1./1. * (TMP25 * - 1./1. * (+cI
      * (TMP47) + 3./1. * cI * (P2[1] * P2[1])))) + (+1./1. * (TMP30 * TMP38 *
      (+cI * (TMP46) + 3./1. * cI * (P1[1] * P1[1]))) + TMP12 * - 1./1. * TMP25
      * (+cI * (TMP46) + 3./1. * cI * (P1[1] * P1[1])))))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[2] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[4]) + P2[2] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[4] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[2] * (-cI * (TMP16) + cI * (TMP9)) + P2[2] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] + P2[2]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[2]) + cI * (P1[2]))) +
      cI * (V1[4] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[2] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[2] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[2] * (TMP12 * (TMP26 * (P1[1] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[3]) + P2[1] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[3] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[1] * (-cI * (TMP16) + cI * (TMP9)) + P2[1] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] +
      P2[1]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[1]) + cI *
      (P1[1]))) + cI * (V1[3] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[1] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[1] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[1] * (P1[2] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[4] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[2] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[1] * (P2[2] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[4] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[2] * V2[3] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * 1./2. * V1[3] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[3] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[5]) + P2[3] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[5] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[3] * (-cI * (TMP16) + cI * (TMP9)) + P2[3] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[3]) + cI * (P1[3]))) +
      cI * (V1[5] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[3] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[3] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[3] * (TMP12 * (TMP26 * (P1[1] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[3]) + P2[1] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[3] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[1] * (-cI * (TMP16) + cI * (TMP9)) + P2[1] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] +
      P2[1]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[1]) + cI *
      (P1[1]))) + cI * (V1[3] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[1] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[1] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[1] * (P1[3] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[5] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[3] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[1] * (P2[3] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[5] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[3] * V2[3] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * 1./2. * V1[3] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[2] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[4]) + P2[2] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[4] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[2] * (-cI * (TMP16) + cI * (TMP9)) + P2[2] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] + P2[2]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[2]) + cI * (P1[2]))) +
      cI * (V1[4] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[2] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[2] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[2] * (TMP12 * (TMP26 * (P1[0] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[2]) + P2[0] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[2] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[0] * (-cI * (TMP16) + cI * (TMP9)) + P2[0] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] +
      P2[0]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[0]) + cI *
      (P1[0]))) + cI * (V1[2] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[0] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[0] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[0] * (P1[2] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[4] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[2] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[0] * (P2[2] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[4] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[2] * V2[2] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * 1./2. * V1[2] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[2] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[4]) + P2[2] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[4] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[2] * (-cI * (TMP16) + cI * (TMP9)) + P2[2] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] + P2[2]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[2]) + cI * (P1[2]))) +
      cI * (V1[4] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[2] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[2] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[2] * (TMP12 * (TMP26 * (P1[1] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[3]) + P2[1] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[3] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[1] * (-cI * (TMP16) + cI * (TMP9)) + P2[1] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] +
      P2[1]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[1]) + cI *
      (P1[1]))) + cI * (V1[3] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[1] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[1] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[1] * (P1[2] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[4] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[2] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[1] * (P2[2] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[4] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[2] * V2[3] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[2] * 1./2. * V1[3] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[12] = denom * 1./3. * (OM3 * (P3[2] * (P3[2] * (OM3 * (TMP12 * (TMP9 *
      (TMP25 * - 2. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 2.
      * (-2. * cI * (TMP30) + cI * (TMP37)) + 2. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2. * (+cI * (TMP30 + TMP37)) - 2. * cI * (TMP16 *
      TMP25)) - 2. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 *
      2. * (+cI * (TMP26 + TMP38)) - 2. * cI * (TMP16 * TMP25)) - 2. * cI *
      (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1. * (-2. * cI * (TMP16)
      + cI * (TMP46 + TMP47)) + (-2. * cI * (TMP26 * TMP30) - cI * (TMP38 *
      TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * (+cI * (TMP46 + TMP47)) +
      (+cI * (TMP37 * TMP47 + TMP12 * TMP48))) + (TMP25 * - 1. * (+cI * (TMP12
      * TMP46 + TMP9 * TMP47)) + TMP30 * (+cI * (TMP38 * TMP46 + TMP9 *
      TMP49)))))) + (TMP12 * (TMP26 * (P1[2] * 3. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (P2[2] * - 6. * (+cI * (TMP30 + TMP37)) - 3. * cI * (TMP9 *
      V1[4]))) + (TMP25 * (P1[2] * 6. * (-cI * (TMP16) + cI * (TMP9)) + 6. *
      (P2[2] * (+cI * (TMP16 + TMP9)))) + 3. * (V2[4] * (-cI * (TMP9 * TMP30) +
      cI * (TMP16 * TMP37))))) + (TMP9 * (TMP30 * (P1[2] * - 6. * (+cI * (TMP26
      + TMP38)) + 3. * (P2[2] * (-cI * (TMP38) + 2. * cI * (TMP26)))) + TMP16 *
      (TMP25 * 6. * (-cI * (P2[2]) + cI * (P1[2])) + 3. * cI * (V1[4] *
      TMP38))) + 3. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] + P2[2])))))) +
      (TMP12 * (TMP9 * (TMP25 * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) +
      (TMP26 * (-cI * (TMP37) + 2. * cI * (TMP30)) - cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * - 1. * (+cI * (TMP30 + TMP37)) + cI * (TMP16 * TMP25))
      + cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 * - 1. * (+cI *
      (TMP26 + TMP38)) + cI * (TMP16 * TMP25)) + cI * (TMP16 * TMP37 *
      TMP38)))) + (TMP16 * (TMP25 * (P1[2] * 3. * (-cI * (P1[2]) + 2. * cI *
      (P2[2])) + (-1./1. * cI * (TMP46 + TMP47) + 2./1. * cI * (TMP16) - 3. *
      cI * (P2[2] * P2[2]))) + (TMP37 * - 1./1. * (+cI * (TMP49) + 3./1. * cI *
      (P2[2] * V2[4])) + (-2./1. * cI * (TMP26 * TMP30) + TMP38 * - 1./1. *
      (+cI * (TMP48) + 3./1. * cI * (P1[2] * V1[4]))))) + (TMP26 * (TMP30 *
      (P1[2] * 3. * (-2. * cI * (P2[2]) + cI * (P1[2])) + (+1./1. * cI * (TMP46
      + TMP47) + 3. * cI * (P2[2] * P2[2]))) + (+1./1. * (TMP12 * (+cI *
      (TMP48) + 3./1. * cI * (P1[2] * V1[4]))) + TMP37 * 1./1. * (+cI * (TMP47)
      + 3./1. * cI * (P2[2] * P2[2])))) + (TMP9 * 1./1. * (TMP30 * (+cI *
      (TMP49) + 3./1. * cI * (P2[2] * V2[4])) + 1./1. * (TMP25 * - 1./1. * (+cI
      * (TMP47) + 3./1. * cI * (P2[2] * P2[2])))) + (+1./1. * (TMP30 * TMP38 *
      (+cI * (TMP46) + 3./1. * cI * (P1[2] * P1[2]))) + TMP12 * - 1./1. * TMP25
      * (+cI * (TMP46) + 3./1. * cI * (P1[2] * P1[2])))))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[3] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[5]) + P2[3] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[5] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[3] * (-cI * (TMP16) + cI * (TMP9)) + P2[3] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[3]) + cI * (P1[3]))) +
      cI * (V1[5] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[3] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[3] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[3] * (TMP12 * (TMP26 * (P1[2] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[4]) + P2[2] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[4] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[2] * (-cI * (TMP16) + cI * (TMP9)) + P2[2] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] +
      P2[2]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[2]) + cI *
      (P1[2]))) + cI * (V1[4] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[2] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[2] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[2] * (P1[3] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[5] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[3] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[2] * (P2[3] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[5] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[3] * V2[4] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * 1./2. * V1[4] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[3] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[5]) + P2[3] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[5] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[3] * (-cI * (TMP16) + cI * (TMP9)) + P2[3] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[3]) + cI * (P1[3]))) +
      cI * (V1[5] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[3] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[3] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[3] * (TMP12 * (TMP26 * (P1[0] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[2]) + P2[0] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[2] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[0] * (-cI * (TMP16) + cI * (TMP9)) + P2[0] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[0] +
      P2[0]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[0]) + cI *
      (P1[0]))) + cI * (V1[2] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[0] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[0] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[0] * (P1[3] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[5] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[3] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[0] * (P2[3] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[5] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[3] * V2[2] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * 1./2. * V1[2] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[3] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[5]) + P2[3] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[5] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[3] * (-cI * (TMP16) + cI * (TMP9)) + P2[3] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[3]) + cI * (P1[3]))) +
      cI * (V1[5] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[3] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[3] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[3] * (TMP12 * (TMP26 * (P1[1] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[3]) + P2[1] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[3] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[1] * (-cI * (TMP16) + cI * (TMP9)) + P2[1] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[1] +
      P2[1]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[1]) + cI *
      (P1[1]))) + cI * (V1[3] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[1] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[1] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[1] * (P1[3] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[5] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[3] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[1] * (P2[3] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[5] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[3] * V2[3] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * 1./2. * V1[3] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP12 * (TMP9 * (TMP25 * -
      2./3. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 4./3. *
      (+1./2. * cI * (TMP37) - cI * (TMP30)) + 2./3. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2./3. * (+cI * (TMP30 + TMP37)) - 2./3. * cI * (TMP16 *
      TMP25)) - 2./3. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30
      * 2./3. * (+cI * (TMP26 + TMP38)) - 2./3. * cI * (TMP16 * TMP25)) - 2./3.
      * cI * (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1./3. * (-2. * cI
      * (TMP16) + cI * (TMP46 + TMP47)) + (-2./3. * cI * (TMP26 * TMP30) -
      1./3. * cI * (TMP38 * TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * 1./3.
      * (+cI * (TMP46 + TMP47)) + (+1./3. * cI * (TMP37 * TMP47 + TMP12 *
      TMP48))) + (TMP25 * - 1./3. * (+cI * (TMP12 * TMP46 + TMP9 * TMP47)) +
      1./3. * (TMP30 * (+cI * (TMP38 * TMP46 + TMP9 * TMP49))))))) + (TMP12 *
      (TMP26 * (P1[3] * 1./2. * (-cI * (TMP37) + 2. * cI * (TMP30)) + (-1./2. *
      cI * (TMP9 * V1[5]) + P2[3] * - 1. * (+cI * (TMP30 + TMP37)))) + (+1./2.
      * (V2[5] * (-cI * (TMP9 * TMP30) + cI * (TMP16 * TMP37))) + TMP25 *
      (P1[3] * (-cI * (TMP16) + cI * (TMP9)) + P2[3] * (+cI * (TMP16 +
      TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3]))) +
      TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[3]) + cI * (P1[3]))) +
      cI * (V1[5] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[3] * (-cI * (TMP38) +
      2. * cI * (TMP26)) + 2. * (P1[3] * - 1. * (+cI * (TMP26 + TMP38)))))))))
      + P3[3] * (TMP12 * (TMP26 * (P1[2] * 1./2. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (-1./2. * cI * (TMP9 * V1[4]) + P2[2] * - 1. * (+cI * (TMP30 +
      TMP37)))) + (+1./2. * (V2[4] * (-cI * (TMP9 * TMP30) + cI * (TMP16 *
      TMP37))) + TMP25 * (P1[2] * (-cI * (TMP16) + cI * (TMP9)) + P2[2] * (+cI
      * (TMP16 + TMP9))))) + (+1./2. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[2] +
      P2[2]))) + TMP9 * 1./2. * (TMP16 * (+2. * (TMP25 * (-cI * (P2[2]) + cI *
      (P1[2]))) + cI * (V1[4] * TMP38)) + 2. * (TMP30 * 1./2. * (P2[2] * (-cI *
      (TMP38) + 2. * cI * (TMP26)) + 2. * (P1[2] * - 1. * (+cI * (TMP26 +
      TMP38))))))))) + (P1[2] * (P1[3] * (TMP25 * - 1. * (+cI * (TMP16 +
      TMP12)) + TMP30 * (+cI * (TMP26 + TMP38))) + (+1./2. * (V1[5] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))) + P2[3] * (-cI * (TMP26 * TMP30)
      + cI * (TMP16 * TMP25)))) + (P2[2] * (P2[3] * (TMP25 * - 1. * (+cI *
      (TMP16 + TMP9)) + TMP26 * (+cI * (TMP30 + TMP37))) + (+1./2. * (V2[5] *
      (-cI * (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * (-cI * (TMP26 *
      TMP30) + cI * (TMP16 * TMP25)))) + (+1./2. * (P2[3] * V2[4] * (-cI *
      (TMP16 * TMP37) + cI * (TMP9 * TMP30))) + P1[3] * 1./2. * V1[4] * (-cI *
      (TMP16 * TMP38) + cI * (TMP12 * TMP26))))));
  T3[17] = denom * 1./3. * (OM3 * (P3[3] * (P3[3] * (OM3 * (TMP12 * (TMP9 *
      (TMP25 * - 2. * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) + (TMP26 * 2.
      * (-2. * cI * (TMP30) + cI * (TMP37)) + 2. * cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * 2. * (+cI * (TMP30 + TMP37)) - 2. * cI * (TMP16 *
      TMP25)) - 2. * cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 *
      2. * (+cI * (TMP26 + TMP38)) - 2. * cI * (TMP16 * TMP25)) - 2. * cI *
      (TMP16 * TMP37 * TMP38))) + (TMP16 * (TMP25 * - 1. * (-2. * cI * (TMP16)
      + cI * (TMP46 + TMP47)) + (-2. * cI * (TMP26 * TMP30) - cI * (TMP38 *
      TMP48 + TMP37 * TMP49))) + (TMP26 * (TMP30 * (+cI * (TMP46 + TMP47)) +
      (+cI * (TMP37 * TMP47 + TMP12 * TMP48))) + (TMP25 * - 1. * (+cI * (TMP12
      * TMP46 + TMP9 * TMP47)) + TMP30 * (+cI * (TMP38 * TMP46 + TMP9 *
      TMP49)))))) + (TMP12 * (TMP26 * (P1[3] * 3. * (-cI * (TMP37) + 2. * cI *
      (TMP30)) + (P2[3] * - 6. * (+cI * (TMP30 + TMP37)) - 3. * cI * (TMP9 *
      V1[5]))) + (TMP25 * (P1[3] * 6. * (-cI * (TMP16) + cI * (TMP9)) + 6. *
      (P2[3] * (+cI * (TMP16 + TMP9)))) + 3. * (V2[5] * (-cI * (TMP9 * TMP30) +
      cI * (TMP16 * TMP37))))) + (TMP9 * (TMP30 * (P1[3] * - 6. * (+cI * (TMP26
      + TMP38)) + 3. * (P2[3] * (-cI * (TMP38) + 2. * cI * (TMP26)))) + TMP16 *
      (TMP25 * 6. * (-cI * (P2[3]) + cI * (P1[3])) + 3. * cI * (V1[5] *
      TMP38))) + 3. * (TMP16 * TMP37 * TMP38 * (+cI * (P1[3] + P2[3])))))) +
      (TMP12 * (TMP9 * (TMP25 * (-2. * cI * (TMP16) + cI * (TMP9 + TMP12)) +
      (TMP26 * (-cI * (TMP37) + 2. * cI * (TMP30)) - cI * (TMP30 * TMP38))) +
      (TMP12 * (TMP26 * - 1. * (+cI * (TMP30 + TMP37)) + cI * (TMP16 * TMP25))
      + cI * (TMP16 * TMP37 * TMP38))) + TMP9 * (TMP9 * (TMP30 * - 1. * (+cI *
      (TMP26 + TMP38)) + cI * (TMP16 * TMP25)) + cI * (TMP16 * TMP37 *
      TMP38)))) + (TMP16 * (TMP25 * (P1[3] * 3. * (-cI * (P1[3]) + 2. * cI *
      (P2[3])) + (-1./1. * cI * (TMP46 + TMP47) + 2./1. * cI * (TMP16) - 3. *
      cI * (P2[3] * P2[3]))) + (TMP37 * - 1./1. * (+cI * (TMP49) + 3./1. * cI *
      (P2[3] * V2[5])) + (-2./1. * cI * (TMP26 * TMP30) + TMP38 * - 1./1. *
      (+cI * (TMP48) + 3./1. * cI * (P1[3] * V1[5]))))) + (TMP26 * (TMP30 *
      (P1[3] * 3. * (-2. * cI * (P2[3]) + cI * (P1[3])) + (+1./1. * cI * (TMP46
      + TMP47) + 3. * cI * (P2[3] * P2[3]))) + (+1./1. * (TMP12 * (+cI *
      (TMP48) + 3./1. * cI * (P1[3] * V1[5]))) + TMP37 * 1./1. * (+cI * (TMP47)
      + 3./1. * cI * (P2[3] * P2[3])))) + (TMP9 * 1./1. * (TMP30 * (+cI *
      (TMP49) + 3./1. * cI * (P2[3] * V2[5])) + 1./1. * (TMP25 * - 1./1. * (+cI
      * (TMP47) + 3./1. * cI * (P2[3] * P2[3])))) + (+1./1. * (TMP30 * TMP38 *
      (+cI * (TMP46) + 3./1. * cI * (P1[3] * P1[3]))) + TMP12 * - 1./1. * TMP25
      * (+cI * (TMP46) + 3./1. * cI * (P1[3] * P1[3])))))));
}

void VVT10_11_12_13_2_3_6_7_8_9_3(complex<double> V1[], complex<double> V2[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
//   double P2[4]; 
//   double OM3; 
//   double P1[4]; 
  complex<double> Ttmp[18]; 
  complex<double> denom; 
  int i; 
  VVT10_3(V1, V2, COUP1, M3, W3, T3); 
  VVT11_3(V1, V2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT12_3(V1, V2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT13_3(V1, V2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT2_3(V1, V2, COUP5, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT3_3(V1, V2, COUP6, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT6_3(V1, V2, COUP7, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT7_3(V1, V2, COUP8, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT8_3(V1, V2, COUP9, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT9_3(V1, V2, COUP10, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}

void FFV6_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
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
  F2[4] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3]
      + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (+1./2. * (M2 * (F1[3] * (V3[3] - cI * (V3[4])) + 2. * (F1[2]
      * 1./2. * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) +
      (P2[1] * - 1. * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./2. * (M2 * (F1[3] * (V3[2] - V3[5]) + 2.
      * (F1[2] * 1./2. * (V3[3] + cI * (V3[4]))))) + F1[5] * (P2[0] * - 1. *
      (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI
      * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void VVT12_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP34; 
  double P1[4]; 
  complex<double> TMP36; 
  double P2[4]; 
  complex<double> TMP16; 
  complex<double> TMP15; 
  complex<double> TMP26; 
  complex<double> TMP32; 
  complex<double> denom; 
  complex<double> TMP13; 
  complex<double> TMP29; 
  double OM1; 
  complex<double> TMP35; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP29 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP35 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP32 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP36 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP34 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP16 * (OM1 * P1[0] * (-cI * (TMP29 + TMP32) + cI * (TMP35
      + TMP36)) + (V2[3] * (+cI * (T3[6] + T3[3])) + (V2[4] * (+cI * (T3[10] +
      T3[4])) + (V2[5] * (+cI * (T3[14] + T3[5])) - 2. * cI * (T3[2] *
      V2[2]))))) + (TMP26 * (OM1 * P1[0] * (-cI * (TMP34) + cI * (TMP13)) +
      (P2[1] * - 1. * (+cI * (T3[3] + T3[6])) + (P2[2] * - 1. * (+cI * (T3[4] +
      T3[10])) + (P2[3] * - 1. * (+cI * (T3[5] + T3[14])) + 2. * cI * (P2[0] *
      T3[2]))))) + (P2[0] * (+cI * (TMP29 + TMP32)) - V2[2] * (+cI * (TMP13 +
      TMP15)))));
  V1[3] = denom * (TMP16 * (OM1 * P1[1] * (-cI * (TMP29 + TMP32) + cI * (TMP35
      + TMP36)) + (V2[2] * - 1. * (+cI * (T3[3] + T3[6])) + (V2[4] * (+cI *
      (T3[11] + T3[8])) + (V2[5] * (+cI * (T3[15] + T3[9])) + 2. * cI * (T3[7]
      * V2[3]))))) + (TMP26 * (OM1 * P1[1] * (-cI * (TMP34) + cI * (TMP13)) +
      (P2[0] * (+cI * (T3[6] + T3[3])) + (P2[2] * - 1. * (+cI * (T3[8] +
      T3[11])) + (P2[3] * - 1. * (+cI * (T3[9] + T3[15])) - 2. * cI * (P2[1] *
      T3[7]))))) + (P2[1] * (+cI * (TMP29 + TMP32)) - V2[3] * (+cI * (TMP13 +
      TMP15)))));
  V1[4] = denom * (TMP16 * (OM1 * P1[2] * (-cI * (TMP29 + TMP32) + cI * (TMP35
      + TMP36)) + (V2[2] * - 1. * (+cI * (T3[4] + T3[10])) + (V2[3] * (+cI *
      (T3[8] + T3[11])) + (V2[5] * (+cI * (T3[16] + T3[13])) + 2. * cI *
      (T3[12] * V2[4]))))) + (TMP26 * (OM1 * P1[2] * (-cI * (TMP34) + cI *
      (TMP13)) + (P2[0] * (+cI * (T3[10] + T3[4])) + (P2[1] * - 1. * (+cI *
      (T3[11] + T3[8])) + (P2[3] * - 1. * (+cI * (T3[13] + T3[16])) - 2. * cI *
      (P2[2] * T3[12]))))) + (P2[2] * (+cI * (TMP29 + TMP32)) - V2[4] * (+cI *
      (TMP13 + TMP15)))));
  V1[5] = denom * (TMP16 * (OM1 * P1[3] * (-cI * (TMP29 + TMP32) + cI * (TMP35
      + TMP36)) + (V2[2] * - 1. * (+cI * (T3[5] + T3[14])) + (V2[3] * (+cI *
      (T3[9] + T3[15])) + (V2[4] * (+cI * (T3[13] + T3[16])) + 2. * cI *
      (T3[17] * V2[5]))))) + (TMP26 * (OM1 * P1[3] * (-cI * (TMP34) + cI *
      (TMP13)) + (P2[0] * (+cI * (T3[14] + T3[5])) + (P2[1] * - 1. * (+cI *
      (T3[15] + T3[9])) + (P2[2] * - 1. * (+cI * (T3[16] + T3[13])) - 2. * cI *
      (P2[3] * T3[17]))))) + (P2[3] * (+cI * (TMP29 + TMP32)) - V2[5] * (+cI *
      (TMP13 + TMP15)))));
}


void FFV5_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP0; 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - cI * TMP0; 
}

void FFV5_7_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV5_0(F1, F2, V3, COUP1, vertex); 
  FFV7_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV5_6_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV5_0(F1, F2, V3, COUP1, vertex); 
  FFV6_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV5_8_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFV5_0(F1, F2, V3, COUP1, vertex); 
  FFV8_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVT7_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP16; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP25 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * TMP25 * (OM3 * (P3[0] * (P3[0] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[0] * TMP9 + P1[0] *
      TMP12))) + 1./3. * cI * (TMP9 * TMP12)) + (-1./3. * cI * (TMP16) + cI *
      (P1[0] * P2[0])));
  T3[6] = denom * TMP25 * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[1] * TMP9 + P1[1] * TMP12))) -
      P3[1] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[1] +
      P1[1] * P2[0])));
  T3[10] = denom * TMP25 * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] * TMP12))) -
      P3[2] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[2] +
      P1[2] * P2[0])));
  T3[14] = denom * TMP25 * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[3] +
      P1[3] * P2[0])));
  T3[3] = denom * TMP25 * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[1] * TMP12 + P2[1] * TMP9))) -
      P3[1] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[1] * P2[0] +
      P1[0] * P2[1])));
  T3[7] = denom * 2. * TMP25 * (OM3 * (P3[1] * (P3[1] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[1] * TMP9 + P1[1] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[1] * P2[1]) + 1./3.
      * cI * (TMP16)));
  T3[11] = denom * TMP25 * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] * TMP12))) -
      P3[2] * (+cI * (P1[1] * TMP12 + P2[1] * TMP9))) + (+cI * (P1[1] * P2[2] +
      P1[2] * P2[1])));
  T3[15] = denom * TMP25 * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[1] * TMP12 + P2[1] * TMP9))) + (+cI * (P1[1] * P2[3] +
      P1[3] * P2[1])));
  T3[4] = denom * TMP25 * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[2] * TMP12 + P2[2] * TMP9))) -
      P3[2] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[2] * P2[0] +
      P1[0] * P2[2])));
  T3[8] = denom * TMP25 * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[2] * TMP12 + P2[2] * TMP9))) -
      P3[2] * (+cI * (P2[1] * TMP9 + P1[1] * TMP12))) + (+cI * (P1[2] * P2[1] +
      P1[1] * P2[2])));
  T3[12] = denom * 2. * TMP25 * (OM3 * (P3[2] * (P3[2] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[2] * P2[2]) + 1./3.
      * cI * (TMP16)));
  T3[16] = denom * TMP25 * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[2] * TMP12 + P2[2] * TMP9))) + (+cI * (P1[2] * P2[3] +
      P1[3] * P2[2])));
  T3[5] = denom * TMP25 * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[3] * P2[0] +
      P1[0] * P2[3])));
  T3[9] = denom * TMP25 * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[1] * TMP9 + P1[1] * TMP12))) + (+cI * (P1[3] * P2[1] +
      P1[1] * P2[3])));
  T3[13] = denom * TMP25 * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[2] * TMP9 + P1[2] * TMP12))) + (+cI * (P1[3] * P2[2] +
      P1[2] * P2[3])));
  T3[17] = denom * 2. * TMP25 * (OM3 * (P3[3] * (P3[3] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[3] * P2[3]) + 1./3.
      * cI * (TMP16)));
}


void VVT11_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP36; 
  complex<double> denom; 
  double OM1; 
  complex<double> TMP35; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP36 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP35 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * - cI * (OM1 * - P1[0] * (TMP35 + TMP36) + (V2[3] * - 1. *
      (T3[6] + T3[3]) + (V2[4] * - 1. * (T3[10] + T3[4]) + (V2[5] * - 1. *
      (T3[14] + T3[5]) + 2. * (T3[2] * V2[2])))));
  V1[3] = denom * cI * (OM1 * P1[1] * (TMP35 + TMP36) + (V2[2] * - 1. * (T3[3]
      + T3[6]) + (V2[4] * (T3[11] + T3[8]) + (V2[5] * (T3[15] + T3[9]) + 2. *
      (T3[7] * V2[3])))));
  V1[4] = denom * cI * (OM1 * P1[2] * (TMP35 + TMP36) + (V2[2] * - 1. * (T3[4]
      + T3[10]) + (V2[3] * (T3[8] + T3[11]) + (V2[5] * (T3[16] + T3[13]) + 2. *
      (T3[12] * V2[4])))));
  V1[5] = denom * cI * (OM1 * P1[3] * (TMP35 + TMP36) + (V2[2] * - 1. * (T3[5]
      + T3[14]) + (V2[3] * (T3[9] + T3[15]) + (V2[4] * (T3[13] + T3[16]) + 2. *
      (T3[17] * V2[5])))));
}


void FFT2_1(complex<double> F2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + T3[0]; 
  F1[1] = +F2[1] + T3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (F2[3] * (P1[0] * (P1[3] * (T3[9] + T3[15] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[6] - T3[3]) + (P1[1] *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8])) + (P1[2]
      * (T3[8] + T3[11] + cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) +
      (P1[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[0] * (T3[6]
      + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] +
      T3[8])) - P2[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P1[1] *
      (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13] + cI * (T3[3] + T3[15]
      + T3[6] + T3[9])) + (P1[3] * (+2. * (T3[7]) + cI * (T3[11] + T3[8]) - 2.
      * (T3[17]) - T3[5] - T3[14]) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] +
      T3[16] + T3[10] + T3[13]) + (P2[0] * - 1. * (T3[14] + T3[5] + 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[17]) +
      2. * cI * (T3[12])) + (P1[2] * - 1. * (+cI * (T3[4] + T3[16] + T3[10] +
      T3[13])) + (P2[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P2[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[0] * - 1. * (+cI *
      (T3[14] + T3[5]) + 2. * cI * (T3[2])) + P2[3] * (+cI * (T3[5] + T3[14]) +
      2. * cI * (T3[17]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P2[1] * - 1.
      * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8] + T3[11] + 2.
      * cI * (T3[12]))))))))) + (F2[2] * (P1[0] * (P1[1] * (T3[15] + T3[9] + cI
      * (T3[10] + T3[4]) - 2. * (T3[3] + T3[6])) + (P1[2] * - 1. * (+2. *
      (T3[4] + T3[10]) + cI * (T3[6] + T3[3]) - T3[16] - T3[13]) + (P1[3] * 2.
      * (T3[17] + T3[2] - T3[5] - T3[14]) + (P2[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[0] * - 1. *
      (T3[14] + T3[5] - 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] - 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] - 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * - 1. * (T3[3] + T3[6] + cI * (T3[13] + T3[16]) - 2. * (T3[9] +
      T3[15])) + (P1[2] * (+2. * (T3[13] + T3[16]) + cI * (T3[9] + T3[15]) -
      T3[4] - T3[10]) + (P2[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] - 2.
      * (T3[17])) + (P2[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P2[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. * (T3[8] + T3[11] -
      cI * (T3[12]) + cI * (T3[7])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] *
      - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + (P2[1] * (+cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] + cI * (T3[6] + T3[3]))
      + (P2[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] + T3[15])) + (P1[2] *
      (+2. * (T3[12]) + cI * (T3[8] + T3[11])) + (P2[1] * - 1. * (T3[11] +
      T3[8] + 2. * cI * (T3[7])) - P2[2] * (+2. * (T3[12]) + cI * (T3[8] +
      T3[11]))))))))) + M1 * (F2[4] * (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (T3[14] * (P2[0] + P1[3] - P1[0] - P2[3]) + (T3[5] * (P1[3] +
      P2[0] - P2[3] - P1[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] *
      (P1[3] - P2[3]))))))))) + F2[5] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P2[1]) + cI * (P1[1]) -
      P2[2]) + (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] *
      (P1[1] - P2[1]))))))))))));
  F1[3] = denom * cI * (F2[2] * (P1[0] * (P1[3] * (+cI * (T3[13] + T3[10] +
      T3[16] + T3[4]) - T3[9] - T3[6] - T3[15] - T3[3]) + (P1[1] * (T3[14] +
      T3[5] + cI * (T3[11] + T3[8]) - 2. * (T3[7] + T3[2])) + (P1[2] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12] + T3[2]) + cI * (T3[14] + T3[5])) +
      (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[0] * (+cI *
      (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P1[1]
      * (P1[2] * (T3[4] + T3[10] - cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])
      - T3[16] - T3[13]) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[7]) - cI *
      (T3[11] + T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P2[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[12]) +
      2. * cI * (T3[17])) + (P1[2] * (-cI * (T3[4] + T3[10]) + cI * (T3[16] +
      T3[13])) + (P2[1] * (-cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6])) +
      (P2[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] + T3[10])) + (P2[0] *
      (-2. * cI * (T3[2]) + cI * (T3[14] + T3[5])) + P2[3] * (-2. * cI *
      (T3[17]) + cI * (T3[5] + T3[14]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15]
      - cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[1] *
      (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] + T3[11] - 2. *
      cI * (T3[12]))))))))) + (F2[3] * (P1[0] * (P1[1] * (T3[15] + T3[9] + 2. *
      (T3[3] + T3[6]) + cI * (T3[10] + T3[4])) + (P1[2] * (T3[16] + T3[13] + 2.
      * (T3[4] + T3[10]) - cI * (T3[6] + T3[3])) + (P1[3] * 2. * (T3[5] +
      T3[17] + T3[2] + T3[14]) + (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * -
      1. * (T3[14] + T3[5] + 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] + 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * - 1. * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) + cI * (T3[13] +
      T3[16])) + (P1[2] * - 1. * (T3[4] + T3[10] - cI * (T3[9] + T3[15]) + 2. *
      (T3[13] + T3[16])) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] *
      (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] + 2.
      * (T3[17])) + (P2[0] * - 1. * (T3[14] + T3[5] + 2. * (T3[2])) + P2[3] *
      (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[1] * (P1[2] * - 2. * (T3[8]
      + T3[11] - cI * (T3[7]) + cI * (T3[12])) + (P2[0] * - 1. * (T3[6] + T3[3]
      + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) +
      (P2[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + P2[2] * (T3[8] +
      T3[11] + 2. * cI * (T3[12]))))))) + P1[2] * (P2[0] * (+cI * (T3[6] +
      T3[3]) - T3[10] - T3[4]) + (P2[3] * (T3[13] + T3[16] - cI * (T3[9] +
      T3[15])) + (P1[2] * (+cI * (T3[8] + T3[11]) - 2. * (T3[12])) + (P2[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) - P2[2] * (+cI * (T3[8] + T3[11]) -
      2. * (T3[12]))))))))) + M1 * (F2[4] * (P1[0] * (T3[6] + T3[3] - cI *
      (T3[10] + T3[4])) + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) +
      (P2[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] +
      T3[15] - cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P2[1]) + cI *
      (P1[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) - P1[2])
      + (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] * (P2[1] -
      P1[1]))))))))) + F2[5] * (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) +
      (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[1] * (T3[3] + T3[6] -
      T3[15] - T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (T3[14] *
      (P2[0] + P2[3] - P1[0] - P1[3]) + (T3[5] * (P2[3] + P2[0] - P1[3] -
      P1[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] * (P1[3] -
      P2[3]))))))))))));
  F1[4] = denom * - cI * (F2[5] * (P1[0] * (P1[3] * - 1. * (T3[9] + T3[6] +
      T3[15] + T3[3] + cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P1[1] * - 1.
      * (+2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8]) - T3[14] - T3[5]) +
      (P1[2] * - 1. * (T3[8] + T3[11] - cI * (T3[14] + T3[5]) + 2. * cI *
      (T3[12] + T3[2])) + (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8])) + P2[2] * (T3[8] + T3[11] + 2. * cI * (T3[12])))))))))
      + (P1[1] * (P1[2] * (T3[4] + T3[10] - cI * (T3[15] + T3[9]) + cI * (T3[3]
      + T3[6]) - T3[16] - T3[13]) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[7]) +
      cI * (T3[11] + T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15]
      - T3[9]) + (P2[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P2[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - 2. * cI * (T3[17]) + cI * (T3[5] + T3[14]) +
      2. * cI * (T3[12])) + (P1[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] +
      T3[10])) + (P2[1] * (-cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])) +
      (P2[2] * (-cI * (T3[4] + T3[10]) + cI * (T3[16] + T3[13])) + (P2[0] * -
      1. * (-2. * cI * (T3[2]) + cI * (T3[14] + T3[5])) - P2[3] * (-2. * cI *
      (T3[17]) + cI * (T3[5] + T3[14]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15]
      + cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] +
      T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8]
      + T3[11] + 2. * cI * (T3[12]))))))))) + (F2[4] * (P1[0] * (P1[1] * - 1. *
      (T3[15] + T3[9] + 2. * (T3[3] + T3[6]) - cI * (T3[10] + T3[4])) + (P1[2]
      * - 1. * (T3[16] + T3[13] + 2. * (T3[4] + T3[10]) + cI * (T3[6] + T3[3]))
      + (P1[3] * - 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (P1[0] * (T3[14] + T3[5] + 2. * (T3[2])) + (P2[0] * - 1. * (T3[14] +
      T3[5] + 2. * (T3[2])) + P2[3] * (T3[5] + T3[14] + 2. * (T3[17]))))))))) +
      (P1[3] * (P1[1] * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) - cI * (T3[13] +
      T3[16])) + (P1[2] * (T3[4] + T3[10] + 2. * (T3[13] + T3[16]) + cI *
      (T3[9] + T3[15])) + (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[3] * (T3[5] +
      T3[14] + 2. * (T3[17])) + (P2[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. *
      (T3[8] + T3[11] - cI * (T3[12]) + cI * (T3[7])) + (P2[0] * (T3[6] + T3[3]
      - cI * (T3[10] + T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      (P2[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] +
      T3[11] - 2. * cI * (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] + cI
      * (T3[6] + T3[3])) + (P2[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] +
      T3[15])) + (P1[2] * (+2. * (T3[12]) + cI * (T3[8] + T3[11])) + (P2[1] * -
      1. * (T3[11] + T3[8] + 2. * cI * (T3[7])) - P2[2] * (+2. * (T3[12]) + cI
      * (T3[8] + T3[11]))))))))) + M1 * (F2[2] * (P1[1] * (T3[3] + T3[6] -
      T3[15] - T3[9]) + (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] *
      (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] + T3[13] - T3[4] -
      T3[10]) + (T3[14] * (P1[0] + P1[3] - P2[0] - P2[3]) + (T3[5] * (P1[3] +
      P1[0] - P2[3] - P2[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] *
      (P2[3] - P1[3]))))))))) + F2[3] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P1[1]) + cI
      * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) + cI * (P2[1]) -
      P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P2[1] - P1[1]))))))))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (P1[3] * (T3[6] + T3[3] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * -
      1. * (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) - cI * (T3[11] + T3[8])) +
      (P1[2] * (+cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2]) - T3[8] -
      T3[11]) + (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[0] *
      (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI
      * (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P1[1]
      * (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13] - cI * (T3[3] + T3[15] +
      T3[6] + T3[9])) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[17]) + cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4]
      + T3[16] + T3[10] + T3[13]) + (P2[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[2] * (P1[3] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17])) + (P1[2] * - 1. * (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) +
      (P2[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P2[2] * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P2[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P1[3] * (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P2[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] *
      (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI *
      (T3[11] + T3[8]) - 2. * (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F2[5] * (P1[0] * (P1[1] * - 1. * (+2. * (T3[3] +
      T3[6]) + cI * (T3[10] + T3[4]) - T3[15] - T3[9]) + (P1[2] * (T3[16] +
      T3[13] + cI * (T3[6] + T3[3]) - 2. * (T3[4] + T3[10])) + (P1[3] * 2. *
      (T3[17] + T3[2] - T3[5] - T3[14]) + (P2[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[0] * - 1. *
      (T3[14] + T3[5] - 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] - 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] - 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * (+2. * (T3[9] + T3[15]) + cI * (T3[13] + T3[16]) - T3[3] -
      T3[6]) + (P1[2] * - 1. * (T3[4] + T3[10] + cI * (T3[9] + T3[15]) - 2. *
      (T3[13] + T3[16])) + (P2[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] - 2.
      * (T3[17])) + (P2[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P2[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. * (T3[8] + T3[11] -
      cI * (T3[7]) + cI * (T3[12])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] +
      T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P2[1] * - 1. * (+2. *
      (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] - cI * (T3[6] + T3[3]))
      + (P2[3] * (+cI * (T3[9] + T3[15]) - T3[13] - T3[16]) + (P1[2] * - 1. *
      (+cI * (T3[8] + T3[11]) - 2. * (T3[12])) + (P2[1] * - 1. * (T3[11] +
      T3[8] - 2. * cI * (T3[7])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M1 * (F2[2] * (P1[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0] *
      (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI
      * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) -
      P1[2]) + (T3[8] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) - P1[2]) + (T3[12]
      * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] * (P2[1] -
      P1[1]))))))))) + F2[3] * (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (T3[14] * (P2[0] + P1[3] - P1[0] - P2[3]) + (T3[5] * (P1[3] +
      P2[0] - P2[3] - P1[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] *
      (P1[3] - P2[3]))))))))))));
}


void VVT8_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP16; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  double OM3; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP25 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP30 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[0] * - 1./3. *
      (+cI * (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[0] * TMP9 +
      P1[0] * TMP12))) + 1./3. * cI * (P3[0] * TMP26 * TMP30)) + TMP26 * TMP30
      * (TMP12 * (-cI * (P1[0]) + 2./3. * cI * (P3[0] * OM3 * TMP9)) - cI *
      (P2[0] * TMP9))) + 1./3. * (TMP12 * TMP9 * (-cI * (TMP16 * TMP25) + cI *
      (TMP26 * TMP30)))) + (TMP16 * (TMP25 * (-cI * (P1[0] * P2[0]) + 1./3. *
      cI * (TMP16)) - 1./3. * cI * (TMP26 * TMP30)) + cI * (P1[0] * P2[0] *
      TMP26 * TMP30)));
  T3[6] = denom * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[1] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[1] * TMP9 + P1[1]
      * TMP12))) + 2./3. * cI * (P3[1] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[1]) + 4./3. * cI * (P3[1] * OM3 * TMP9)) - cI *
      (P2[1] * TMP9))) + P3[1] * (P1[0] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[0] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[0] * P2[1] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[1] * P2[0] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[10] = denom * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[2] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[2] * TMP9 + P1[2]
      * TMP12))) + 2./3. * cI * (P3[2] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[2]) + 4./3. * cI * (P3[2] * OM3 * TMP9)) - cI *
      (P2[2] * TMP9))) + P3[2] * (P1[0] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[0] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[0] * P2[2] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[2] * P2[0] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[14] = denom * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[3] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[3] * TMP9 + P1[3]
      * TMP12))) + 2./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[3]) + 4./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + P3[3] * (P1[0] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[0] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[0] * P2[3] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[3] * P2[0] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[3] = denom * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[1] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P1[1] * TMP12 + P2[1]
      * TMP9))) + 2./3. * cI * (P3[1] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[1]) + 4./3. * cI * (P3[1] * OM3 * TMP9)) - cI *
      (P2[1] * TMP9))) + P3[1] * (P1[0] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[0] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[0] * P2[1] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[1] * P2[0] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (TMP16 * (TMP25 * (P3[1] * - 1./3. *
      (+cI * (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[1] * TMP9 +
      P1[1] * TMP12))) + 1./3. * cI * (P3[1] * TMP26 * TMP30)) + TMP26 * TMP30
      * (TMP12 * (-cI * (P1[1]) + 2./3. * cI * (P3[1] * OM3 * TMP9)) - cI *
      (P2[1] * TMP9))) + 1./3. * (TMP12 * TMP9 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)))) + (TMP16 * (TMP25 * - 1. * (+cI * (P1[1] * P2[1]) +
      1./3. * cI * (TMP16)) + 1./3. * cI * (TMP26 * TMP30)) + cI * (P1[1] *
      P2[1] * TMP26 * TMP30)));
  T3[11] = denom * (OM3 * (P3[1] * (TMP16 * (TMP25 * (P3[2] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[2] * TMP9 + P1[2]
      * TMP12))) + 2./3. * cI * (P3[2] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[2]) + 4./3. * cI * (P3[2] * OM3 * TMP9)) - cI *
      (P2[2] * TMP9))) + P3[2] * (P1[1] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[1] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[1] * P2[2] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[2] * P2[1] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[15] = denom * (OM3 * (P3[1] * (TMP16 * (TMP25 * (P3[3] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[3] * TMP9 + P1[3]
      * TMP12))) + 2./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[3]) + 4./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + P3[3] * (P1[1] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[1] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[1] * P2[3] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[3] * P2[1] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[4] = denom * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[2] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P1[2] * TMP12 + P2[2]
      * TMP9))) + 2./3. * cI * (P3[2] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[2]) + 4./3. * cI * (P3[2] * OM3 * TMP9)) - cI *
      (P2[2] * TMP9))) + P3[2] * (P1[0] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[0] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[0] * P2[2] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[2] * P2[0] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[8] = denom * (OM3 * (P3[1] * (TMP16 * (TMP25 * (P3[2] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P1[2] * TMP12 + P2[2]
      * TMP9))) + 2./3. * cI * (P3[2] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[2]) + 4./3. * cI * (P3[2] * OM3 * TMP9)) - cI *
      (P2[2] * TMP9))) + P3[2] * (P1[1] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[1] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[1] * P2[2] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[2] * P2[1] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (TMP16 * (TMP25 * (P3[2] * - 1./3. *
      (+cI * (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[2] * TMP9 +
      P1[2] * TMP12))) + 1./3. * cI * (P3[2] * TMP26 * TMP30)) + TMP26 * TMP30
      * (TMP12 * (-cI * (P1[2]) + 2./3. * cI * (P3[2] * OM3 * TMP9)) - cI *
      (P2[2] * TMP9))) + 1./3. * (TMP12 * TMP9 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)))) + (TMP16 * (TMP25 * - 1. * (+cI * (P1[2] * P2[2]) +
      1./3. * cI * (TMP16)) + 1./3. * cI * (TMP26 * TMP30)) + cI * (P1[2] *
      P2[2] * TMP26 * TMP30)));
  T3[16] = denom * (OM3 * (P3[2] * (TMP16 * (TMP25 * (P3[3] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[3] * TMP9 + P1[3]
      * TMP12))) + 2./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[3]) + 4./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + P3[3] * (P1[2] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[2] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[2] * P2[3] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[3] * P2[2] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[5] = denom * (OM3 * (P3[0] * (TMP16 * (TMP25 * (P3[3] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P1[3] * TMP12 + P2[3]
      * TMP9))) + 2./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[3]) + 4./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + P3[3] * (P1[0] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[0] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[0] * P2[3] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[3] * P2[0] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[9] = denom * (OM3 * (P3[1] * (TMP16 * (TMP25 * (P3[3] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P1[3] * TMP12 + P2[3]
      * TMP9))) + 2./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[3]) + 4./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + P3[3] * (P1[1] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[1] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[1] * P2[3] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[3] * P2[1] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[13] = denom * (OM3 * (P3[2] * (TMP16 * (TMP25 * (P3[3] * - 2./3. * (+cI *
      (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P1[3] * TMP12 + P2[3]
      * TMP9))) + 2./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30 *
      (TMP12 * (-cI * (P1[3]) + 4./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + P3[3] * (P1[2] * TMP12 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)) + P2[2] * TMP9 * (-cI * (TMP26 * TMP30) + cI * (TMP16 *
      TMP25)))) + (P1[2] * P2[3] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30)) + P1[3] * P2[2] * (-cI * (TMP16 * TMP25) + cI * (TMP26 *
      TMP30))));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (TMP16 * (TMP25 * (P3[3] * - 1./3. *
      (+cI * (TMP16) + 2. * cI * (OM3 * TMP9 * TMP12)) + (+cI * (P2[3] * TMP9 +
      P1[3] * TMP12))) + 1./3. * cI * (P3[3] * TMP26 * TMP30)) + TMP26 * TMP30
      * (TMP12 * (-cI * (P1[3]) + 2./3. * cI * (P3[3] * OM3 * TMP9)) - cI *
      (P2[3] * TMP9))) + 1./3. * (TMP12 * TMP9 * (-cI * (TMP26 * TMP30) + cI *
      (TMP16 * TMP25)))) + (TMP16 * (TMP25 * - 1. * (+cI * (P1[3] * P2[3]) +
      1./3. * cI * (TMP16)) + 1./3. * cI * (TMP26 * TMP30)) + cI * (P1[3] *
      P2[3] * TMP26 * TMP30)));
}


void FFV8_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
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
  F2[4] = denom * - 4. * cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + (+1./4. * (M2 * (F1[3] * (+cI * (V3[4]) - V3[3]) + 4. *
      (F1[2] * - 1./4. * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI *
      (V3[4])) + (P2[1] * - 1. * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] +
      V3[5])) + P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * - 4. * cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./4. * (M2 * (F1[3] * (V3[5] - V3[2]) + 4.
      * (F1[2] * - 1./4. * (V3[3] + cI * (V3[4]))))) + F1[5] * (P2[0] * - 1. *
      (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI
      * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void VVT2_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP16; 
  complex<double> TMP66; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP65; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP65 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP66 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (P3[0] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[0] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[0] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)))) + (P1[0] * P2[0] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP66) + cI * (TMP65)))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[1] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[1] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[1] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[1] * (-cI * (TMP65) + cI * (TMP66)) + P1[1] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI
      * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI
      * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[1] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[1] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[1] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[1] * (-cI * (TMP65) + cI * (TMP66)) + P1[1] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (P3[1] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[1] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[1] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP66) + cI * (TMP65)))) + (P1[1] * P2[1] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI * (TMP66)))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI
      * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI
      * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (P3[2] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP66) + cI * (TMP65)))) + (P1[2] * P2[2] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI * (TMP66)))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI
      * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[2] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[2] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[2] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[2] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 * (-cI
      * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[2] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[2] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[2] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[2] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (P3[3] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP66) + cI * (TMP65)))) + (P1[3] * P2[3] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI * (TMP66)))));
}


void VVT5_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP54; 
  complex<double> denom; 
  complex<double> TMP53; 
  complex<double> TMP9; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP53 = -1. * (P1[0] * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] *
      V1[4]))) + (P1[1] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] *
      V1[2]))) + (P1[2] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + P1[3] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] *
      (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] *
      V1[2]))))));
  TMP54 = -1. * (P1[0] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))) + (P1[1] * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + (P1[2] * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (TMP37 * (OM3 * P3[0] * (TMP12 * (P1[1] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) +
      P1[3] * (P3[1] * V2[4] - P3[2] * V2[3]))) - 1./3. * (P3[0] * TMP54)) +
      (P2[0] * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      1./3. * (TMP54))) + TMP38 * (OM3 * P3[0] * (TMP9 * (P2[1] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3]
      * (P3[1] * V1[4] - P3[2] * V1[3]))) - 1./3. * (P3[0] * TMP53)) + (P1[0] *
      (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + 1./3. *
      (TMP53))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + 2./3. * (P3[1] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) +
      2./3. * (P3[1] * TMP53))) + P3[1] * (TMP12 * TMP37 * (P1[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 * TMP9 * (P2[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP37 * (P2[0] * (P1[0] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] - P3[0] * V2[5]) +
      P1[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + P2[1] * (P1[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P1[0] * (P2[0] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[4] - P3[2] * V1[2]))) + P1[1] * (P2[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP53))) + P3[2] * (TMP12 * TMP37 * (P1[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 * TMP9 * (P2[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP37 * (P2[0] * (P1[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P1[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P2[2] * (P1[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P1[0] * (P2[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P1[2] * (P2[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP53))) + P3[3] * (TMP12 * TMP37 * (P1[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 * TMP9 * (P2[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP37 * (P2[0] * (P1[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P2[3] * (P1[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P1[0] * (P2[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P1[3] * (P2[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + 2./3. * (P3[1] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) +
      2./3. * (P3[1] * TMP53))) + P3[1] * (TMP12 * TMP37 * (P1[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 * TMP9 * (P2[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP37 * (P2[0] * (P1[0] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] - P3[0] * V2[5]) +
      P1[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + P2[1] * (P1[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P1[0] * (P2[0] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[4] - P3[2] * V1[2]))) + P1[1] * (P2[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[7] = denom * 2. * cI * (TMP37 * (OM3 * P3[1] * (TMP12 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + 1./3. * (P3[1] * TMP54)) + (P2[1] *
      (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + 1./3. *
      (TMP54))) + TMP38 * (OM3 * P3[1] * (TMP9 * (P2[0] * (P3[3] * V1[4] -
      P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] *
      (P3[2] * V1[2] - P3[0] * V1[4]))) + 1./3. * (P3[1] * TMP53)) + (P1[1] *
      (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) + 1./3. *
      (TMP53))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP53))) + P3[2] * (TMP12 * TMP37 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 * TMP9 * (P2[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP37 * (P2[1] * (P1[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P1[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P2[2] * (P1[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P1[1] * (P2[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P1[2] * (P2[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP53))) + P3[3] * (TMP12 * TMP37 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 * TMP9 * (P2[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP37 * (P2[1] * (P1[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P2[3] * (P1[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P1[1] * (P2[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P1[3] * (P2[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP53))) + P3[2] * (TMP12 * TMP37 * (P1[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 * TMP9 * (P2[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP37 * (P2[0] * (P1[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P1[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P2[2] * (P1[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P1[0] * (P2[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P1[2] * (P2[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP53))) + P3[2] * (TMP12 * TMP37 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 * TMP9 * (P2[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP37 * (P2[1] * (P1[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P1[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P2[2] * (P1[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P1[1] * (P2[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P1[2] * (P2[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[12] = denom * 2. * cI * (TMP37 * (OM3 * P3[2] * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 1./3. * (P3[2] * TMP54)) + (P2[2] *
      (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + 1./3. *
      (TMP54))) + TMP38 * (OM3 * P3[2] * (TMP9 * (P2[0] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[3] - P3[1] * V1[2]))) + 1./3. * (P3[2] * TMP53)) + (P1[2] *
      (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) + 1./3. *
      (TMP53))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP53))) + P3[3] * (TMP12 * TMP37 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + TMP38 * TMP9 * (P2[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))))) + (TMP37 * (P2[2] * (P1[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P2[3] * (P1[0] * (P3[3] *
      V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[1] * V2[2] - P3[0] * V2[3])))) + TMP38 * (P1[2] * (P2[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P1[3] * (P2[0] * (P3[3] * V1[3] -
      P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] *
      (P3[1] * V1[2] - P3[0] * V1[3]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP53))) + P3[3] * (TMP12 * TMP37 * (P1[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P1[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 * TMP9 * (P2[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))))) + (TMP37 * (P2[0] * (P1[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P2[3] * (P1[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P1[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P1[0] * (P2[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P1[3] * (P2[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP53))) + P3[3] * (TMP12 * TMP37 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 * TMP9 * (P2[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))))) + (TMP37 * (P2[1] * (P1[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P2[3] * (P1[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P1[1] * (P2[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P1[3] * (P2[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP54)) + TMP38 *
      (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP53))) + P3[3] * (TMP12 * TMP37 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + TMP38 * TMP9 * (P2[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))))) + (TMP37 * (P2[2] * (P1[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P2[3] * (P1[0] * (P3[3] *
      V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[1] * V2[2] - P3[0] * V2[3])))) + TMP38 * (P1[2] * (P2[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P1[3] * (P2[0] * (P3[3] * V1[3] -
      P3[1] * V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] *
      (P3[1] * V1[2] - P3[0] * V1[3]))))));
  T3[17] = denom * 2. * cI * (TMP37 * (OM3 * P3[3] * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 1./3. * (P3[3] * TMP54)) + (P2[3] *
      (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] -
      P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + 1./3. *
      (TMP54))) + TMP38 * (OM3 * P3[3] * (TMP9 * (P2[0] * (P3[2] * V1[3] -
      P3[1] * V1[4]) + (P2[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] *
      (P3[1] * V1[2] - P3[0] * V1[3]))) + 1./3. * (P3[3] * TMP53)) + (P1[3] *
      (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] * V1[2] -
      P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) + 1./3. *
      (TMP53))));
}


void FFT1_2(complex<double> F1[], complex<double> T3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP13; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + T3[0]; 
  F2[1] = +F1[1] + T3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (TMP13 * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) - F1[2] * M2)) + TMP15 * (F1[4] * (P2[0] - P2[3]) +
      (F1[5] * (+cI * (P2[2]) - P2[1]) - F1[2] * M2)));
  F2[3] = denom * - cI * (TMP13 * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * -
      1. * (P2[0] + P2[3]) + F1[3] * M2)) + TMP15 * (F1[4] * (P2[1] + cI *
      (P2[2])) + (F1[5] * - 1. * (P2[0] + P2[3]) + F1[3] * M2)));
  F2[4] = denom * cI * (TMP13 * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] * (+cI
      * (P2[2]) - P2[1]) + F1[4] * M2)) + TMP15 * (F1[2] * - 1. * (P2[0] +
      P2[3]) + (F1[3] * (+cI * (P2[2]) - P2[1]) + F1[4] * M2)));
  F2[5] = denom * - cI * (TMP13 * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) - F1[5] * M2)) + TMP15 * (F1[2] * (P2[1] + cI * (P2[2]))
      + (F1[3] * (P2[0] - P2[3]) - F1[5] * M2)));
}

void FFT1_2_3_5_2(complex<double> F1[], complex<double> T3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P1[4]; 
//   double P2[4]; 
  int i; 
  complex<double> Ftmp[6]; 
  FFT1_2(F1, T3, COUP1, M2, W2, F2); 
  FFT2_2(F1, T3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFT3_2(F1, T3, COUP3, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFT5_2(F1, T3, COUP4, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFT1_2_4_5_2(complex<double> F1[], complex<double> T3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P1[4]; 
//   double P2[4]; 
  int i; 
  complex<double> Ftmp[6]; 
  FFT1_2(F1, T3, COUP1, M2, W2, F2); 
  FFT2_2(F1, T3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFT4_2(F1, T3, COUP3, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFT5_2(F1, T3, COUP4, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFT4_2(complex<double> F1[], complex<double> T3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + T3[0]; 
  F2[1] = +F1[1] + T3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[3] * (P2[0] * (P2[3] * (T3[9] + T3[6] + T3[15] +
      T3[3] - cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P2[1] * - 1. *
      (T3[14] + T3[5] + cI * (T3[11] + T3[8]) - 2. * (T3[7] + T3[2])) + (P2[2]
      * (T3[8] + T3[11] - 2. * cI * (T3[12] + T3[2]) + cI * (T3[14] + T3[5])) +
      (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P1[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0] * (+cI * (T3[10] + T3[4]) -
      T3[6] - T3[3]) + (P1[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P1[2]
      * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[16] + T3[13] - cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6]) - T3[4] -
      T3[10]) + (P2[3] * (+2. * (T3[17]) + cI * (T3[11] + T3[8]) - 2. * (T3[7])
      - T3[5] - T3[14]) + (P1[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P1[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] + T3[9] - T3[3] -
      T3[6]) + (P1[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P1[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P2[2] * (P2[3] * (-2. * cI * (T3[17]) +
      cI * (T3[5] + T3[14]) + 2. * cI * (T3[12]) - T3[8] - T3[11]) + (P1[1] *
      (-cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])) + (P1[2] * (-cI * (T3[4]
      + T3[10]) + cI * (T3[16] + T3[13])) + (P2[2] * (-cI * (T3[16] + T3[13]) +
      cI * (T3[4] + T3[10])) + (P1[0] * - 1. * (-2. * cI * (T3[2]) + cI *
      (T3[14] + T3[5])) - P1[3] * (-2. * cI * (T3[17]) + cI * (T3[5] +
      T3[14]))))))) + P2[3] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3])
      + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. * (+cI * (T3[11] +
      T3[8]) - 2. * (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F1[2] * (P2[0] * (P2[1] * (T3[15] + T3[9] + 2. *
      (T3[3] + T3[6]) + cI * (T3[10] + T3[4])) + (P2[2] * (T3[16] + T3[13] + 2.
      * (T3[4] + T3[10]) - cI * (T3[6] + T3[3])) + (P1[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (P2[3] * 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P1[0] *
      (T3[14] + T3[5] + 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] + 2. *
      (T3[17])) - P2[0] * (T3[14] + T3[5] + 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * - 1. * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) + cI * (T3[13] +
      T3[16])) + (P2[2] * - 1. * (T3[4] + T3[10] - cI * (T3[9] + T3[15]) + 2. *
      (T3[13] + T3[16])) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[2] *
      (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * - 1. * (T3[14] + T3[5] + 2.
      * (T3[2])) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[17])) - P2[3] * (T3[5] +
      T3[14] + 2. * (T3[17])))))))) + (P2[1] * (P1[0] * - 1. * (T3[6] + T3[3] +
      cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[7]) + cI *
      (T3[12])) + (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P1[2] *
      (T3[8] + T3[11] + 2. * cI * (T3[12])) - P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8]))))))) + P2[2] * (P1[0] * (+cI * (T3[6] + T3[3]) - T3[10]
      - T3[4]) + (P1[3] * (T3[13] + T3[16] - cI * (T3[9] + T3[15])) + (P1[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) + (P1[2] * - 1. * (+cI * (T3[8] +
      T3[11]) - 2. * (T3[12])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M2 * (F1[4] * (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (T3[14] * (P1[0] + P1[3] - P2[0] - P2[3]) + (T3[5] * (P1[3] + P1[0] -
      P2[3] - P2[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))) + F1[5] * (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4]))
      + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0] * (+cI *
      (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI *
      (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) -
      P1[2]) + (T3[8] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) - P1[2]) + (T3[12]
      * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] * (P2[1] -
      P1[1]))))))))))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (P2[3] * (T3[9] + T3[15] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[6] - T3[3]) + (P2[1] *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8])) + (P2[2]
      * (T3[8] + T3[11] + cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) +
      (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI
      * (T3[10] + T3[4])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] +
      T3[8])) - P1[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P2[1] *
      (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13] + cI * (T3[3] + T3[15]
      + T3[6] + T3[9])) + (P2[3] * (+2. * (T3[7]) + cI * (T3[11] + T3[8]) - 2.
      * (T3[17]) - T3[5] - T3[14]) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[0] * - 1. * (T3[14] + T3[5] + 2. * (T3[2]))
      + P1[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[2] * (P2[3] *
      (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[17]) + 2. * cI *
      (T3[12])) + (P1[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P1[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[2] * - 1. * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P1[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P1[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P2[3] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4]))
      + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P2[3] *
      (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P1[1] * - 1. * (+2. *
      (T3[7]) + cI * (T3[11] + T3[8])) - P1[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))))) + (F1[3] * (P2[0] * (P2[1] * - 1. * (T3[15] + T3[9] + cI
      * (T3[10] + T3[4]) - 2. * (T3[3] + T3[6])) + (P2[2] * (+2. * (T3[4] +
      T3[10]) + cI * (T3[6] + T3[3]) - T3[16] - T3[13]) + (P1[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (P2[3] * 2. * (T3[5] + T3[14] - T3[17] - T3[2]) + (P1[0] * - 1. * (T3[14]
      + T3[5] - 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] - 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * (T3[3] + T3[6] + cI * (T3[13] + T3[16]) - 2. * (T3[9] + T3[15]))
      + (P2[2] * - 1. * (+2. * (T3[13] + T3[16]) + cI * (T3[9] + T3[15]) -
      T3[4] - T3[10]) + (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P1[2] *
      (T3[16] + T3[13] - T3[4] - T3[10]) + (P1[0] * - 1. * (T3[14] + T3[5] - 2.
      * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[1] * (P1[0] * (+cI * (T3[10]
      + T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] +
      T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[12]) + cI *
      (T3[7])) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      (P1[2] * (T3[8] + T3[11] - 2. * cI * (T3[12])) + P2[1] * (+cI * (T3[11] +
      T3[8]) - 2. * (T3[7]))))))) + P2[2] * (P1[0] * - 1. * (T3[10] + T3[4] +
      cI * (T3[6] + T3[3])) + (P1[3] * (T3[13] + T3[16] + cI * (T3[9] +
      T3[15])) + (P1[1] * (T3[11] + T3[8] + 2. * cI * (T3[7])) + (P1[2] * (+2.
      * (T3[12]) + cI * (T3[8] + T3[11])) - P2[2] * (+2. * (T3[12]) + cI *
      (T3[8] + T3[11]))))))))) + M2 * (F1[4] * (P1[0] * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] -
      cI * (P1[1]) + cI * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) +
      cI * (P2[1]) - P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) +
      2. * (T3[7] * (P2[1] - P1[1]))))))))) + F1[5] * (P1[1] * (T3[3] + T3[15]
      + T3[6] + T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] *
      - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16]
      + T3[10] + T3[13]) + (T3[14] * (P2[0] + P1[3] - P1[0] - P2[3]) + (T3[5] *
      (P1[3] + P2[0] - P2[3] - P1[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. *
      (T3[17] * (P1[3] - P2[3]))))))))))));
  F2[4] = denom * - cI * (F1[5] * (P2[0] * (P2[3] * (T3[6] + T3[3] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[1] * -
      1. * (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) - cI * (T3[11] + T3[8])) +
      (P2[2] * (+cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2]) - T3[8] -
      T3[11]) + (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P1[3] *
      (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] - cI
      * (T3[10] + T3[4])) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P2[1]
      * (P2[2] * (T3[4] + T3[16] + T3[10] + T3[13] - cI * (T3[3] + T3[15] +
      T3[6] + T3[9])) + (P2[3] * (T3[5] + T3[14] + 2. * (T3[17]) + cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] *
      (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) - P1[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[2] *
      (P2[3] * - 1. * (T3[8] + T3[11] - 2. * cI * (T3[12]) + cI * (T3[5] +
      T3[14]) + 2. * cI * (T3[17])) + (P1[1] * (+cI * (T3[3] + T3[15] + T3[6] +
      T3[9])) + (P1[2] * (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[2] *
      - 1. * (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P1[0] * - 1. * (+cI
      * (T3[14] + T3[5]) + 2. * cI * (T3[2])) + P1[3] * (+cI * (T3[5] + T3[14])
      + 2. * cI * (T3[17]))))))) + P2[3] * (P1[0] * (+cI * (T3[10] + T3[4]) -
      T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) +
      (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. *
      (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. *
      cI * (T3[12]))))))))) + (F1[4] * (P2[0] * (P2[1] * (+2. * (T3[3] + T3[6])
      + cI * (T3[10] + T3[4]) - T3[15] - T3[9]) + (P2[2] * - 1. * (T3[16] +
      T3[13] + cI * (T3[6] + T3[3]) - 2. * (T3[4] + T3[10])) + (P1[1] * (T3[15]
      + T3[9] - T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (P2[3] * 2. * (T3[5] + T3[14] - T3[17] - T3[2]) + (P1[0] * - 1. * (T3[14]
      + T3[5] - 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] - 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * - 1. * (+2. * (T3[9] + T3[15]) + cI * (T3[13] + T3[16]) - T3[3]
      - T3[6]) + (P2[2] * (T3[4] + T3[10] + cI * (T3[9] + T3[15]) - 2. *
      (T3[13] + T3[16])) + (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P1[2] *
      (T3[16] + T3[13] - T3[4] - T3[10]) + (P1[0] * - 1. * (T3[14] + T3[5] - 2.
      * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[1] * (P1[0] * - 1. * (T3[6]
      + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[7]) + cI
      * (T3[12])) + (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P1[2] *
      (T3[8] + T3[11] + 2. * cI * (T3[12])) - P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8]))))))) + P2[2] * (P1[0] * (+cI * (T3[6] + T3[3]) - T3[10]
      - T3[4]) + (P1[3] * (T3[13] + T3[16] - cI * (T3[9] + T3[15])) + (P1[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) + (P1[2] * - 1. * (+cI * (T3[8] +
      T3[11]) - 2. * (T3[12])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M2 * (F1[2] * (P1[1] * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. *
      (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] +
      T3[10] + T3[13]) + (T3[14] * (P2[0] + P1[3] - P1[0] - P2[3]) + (T3[5] *
      (P1[3] + P2[0] - P2[3] - P1[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. *
      (T3[17] * (P1[3] - P2[3]))))))))) + F1[3] * (P1[0] * (+cI * (T3[10] +
      T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] +
      T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] *
      (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI *
      (P1[1]) + cI * (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. *
      (T3[7] * (P1[1] - P2[1]))))))))))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * (P2[3] * (T3[9] + T3[6] + T3[15] +
      T3[3] + cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P2[1] * (+2. * (T3[7]
      + T3[2]) + cI * (T3[11] + T3[8]) - T3[14] - T3[5]) + (P2[2] * (T3[8] +
      T3[11] - cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) + (P1[0] *
      (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15]
      + cI * (T3[13] + T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) -
      P1[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[16] + T3[13] - cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9]) - T3[4] -
      T3[10]) + (P2[3] * - 1. * (T3[5] + T3[14] + 2. * (T3[7]) + cI * (T3[11] +
      T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15] - T3[9]) +
      (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] + T3[9] -
      T3[3] - T3[6]) + (P1[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P1[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[2] * (P2[3] * - 1. * (T3[8]
      + T3[11] - 2. * cI * (T3[17]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[12])) + (P1[1] * (-cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6])) +
      (P1[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] + T3[10])) + (P2[2] *
      (-cI * (T3[4] + T3[10]) + cI * (T3[16] + T3[13])) + (P1[0] * (-2. * cI *
      (T3[2]) + cI * (T3[14] + T3[5])) + P1[3] * (-2. * cI * (T3[17]) + cI *
      (T3[5] + T3[14]))))))) + P2[3] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P1[1] * (+2.
      * (T3[7]) + cI * (T3[11] + T3[8])) + P1[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))))) + (F1[5] * (P2[0] * (P2[1] * - 1. * (T3[15] + T3[9] + 2.
      * (T3[3] + T3[6]) - cI * (T3[10] + T3[4])) + (P2[2] * - 1. * (T3[16] +
      T3[13] + 2. * (T3[4] + T3[10]) + cI * (T3[6] + T3[3])) + (P1[1] * (T3[3]
      + T3[15] + T3[6] + T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (P2[3] * - 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P1[0] * - 1. *
      (T3[14] + T3[5] + 2. * (T3[2])) + (P1[3] * (T3[5] + T3[14] + 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] + 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) - cI * (T3[13] + T3[16]))
      + (P2[2] * (T3[4] + T3[10] + 2. * (T3[13] + T3[16]) + cI * (T3[9] +
      T3[15])) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[2] * -
      1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] + 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[1] * (P1[0] * (T3[6] + T3[3]
      - cI * (T3[10] + T3[4])) + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P2[2] * 2. * (T3[8] + T3[11] - cI * (T3[12]) + cI * (T3[7])) +
      (P1[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + (P1[2] * - 1. * (T3[8]
      + T3[11] - 2. * cI * (T3[12])) - P2[1] * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7]))))))) + P2[2] * (P1[0] * (T3[10] + T3[4] + cI * (T3[6] + T3[3]))
      + (P1[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] + T3[15])) + (P1[1] * -
      1. * (T3[11] + T3[8] + 2. * cI * (T3[7])) + (P1[2] * - 1. * (+2. *
      (T3[12]) + cI * (T3[8] + T3[11])) + P2[2] * (+2. * (T3[12]) + cI * (T3[8]
      + T3[11]))))))))) + M2 * (F1[2] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P1[1]) + cI
      * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) + cI * (P2[1]) -
      P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P2[1] - P1[1]))))))))) + F1[3] * (P1[1] * (T3[15] + T3[9] - T3[3] -
      T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[1] * (T3[3] +
      T3[6] - T3[15] - T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) +
      (T3[14] * (P2[0] + P2[3] - P1[0] - P1[3]) + (T3[5] * (P2[3] + P2[0] -
      P1[3] - P1[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] * (P1[3] -
      P2[3]))))))))))));
}


void FFT3_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP9; 
  complex<double> TMP8; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = -1. * (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI *
      (P3[2]))) + (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] -
      P3[3])) + (F1[4] * (F2[2] * (P3[3] - P3[0]) + F2[3] * (P3[1] + cI *
      (P3[2]))) + F1[5] * (F2[2] * (P3[1] - cI * (P3[2])) - F2[3] * (P3[0] +
      P3[3])))));
  TMP11 = -1. * (F1[2] * (F2[4] * (P2[0] + P2[3]) + F2[5] * (P2[1] + cI *
      (P2[2]))) + (F1[3] * (F2[4] * (P2[1] - cI * (P2[2])) + F2[5] * (P2[0] -
      P2[3])) + (F1[4] * (F2[2] * (P2[3] - P2[0]) + F2[3] * (P2[1] + cI *
      (P2[2]))) + F1[5] * (F2[2] * (P2[1] - cI * (P2[2])) - F2[3] * (P2[0] +
      P2[3])))));
  TMP8 = -1. * (F1[2] * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI *
      (P1[2]))) + (F1[3] * (F2[4] * (P1[1] - cI * (P1[2])) + F2[5] * (P1[0] -
      P1[3])) + (F1[4] * (F2[2] * (P1[3] - P1[0]) + F2[3] * (P1[1] + cI *
      (P1[2]))) + F1[5] * (F2[2] * (P1[1] - cI * (P1[2])) - F2[3] * (P1[0] +
      P1[3])))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[2] * F2[4] +
      F1[3] * F2[5] + 2./3. * (P3[0] * OM3 * TMP10) - F1[4] * F2[2] - F1[5] *
      F2[3]) + (TMP9 * (F1[2] * F2[4] + F1[3] * F2[5] + 2./3. * (P3[0] * OM3 *
      TMP10) - F1[4] * F2[2] - F1[5] * F2[3]) + (P3[0] * 1./3. * (TMP8 - TMP11)
      + TMP10 * (P2[0] - P1[0])))) + 1./3. * (TMP10 * (TMP9 - TMP12))) + (P1[0]
      * (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) +
      (P2[0] * (F1[2] * F2[4] + F1[3] * F2[5] - F1[4] * F2[2] - F1[5] * F2[3])
      + (-1./3. * (TMP8) + 1./3. * (TMP11)))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP12 * (F1[2] * F2[5] + F1[3] * F2[4]
      + F1[4] * F2[3] + F1[5] * F2[2] - 4./3. * (P3[1] * OM3 * TMP10)) + (TMP9
      * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2] -
      4./3. * (P3[1] * OM3 * TMP10)) + (P3[1] * 2./3. * (TMP8 - TMP11) + TMP10
      * (P2[1] - P1[1])))) + P3[1] * (TMP12 * (F1[4] * F2[2] + F1[5] * F2[3] -
      F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] * F2[5] -
      F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) + (F1[2] *
      (F2[4] * (P2[1] - P1[1]) + F2[5] * (P1[0] - P2[0])) + (F1[3] * (F2[4] *
      (P1[0] - P2[0]) + F2[5] * (P2[1] - P1[1])) + (F1[4] * (F2[2] * (P1[1] -
      P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] * (P1[0] - P2[0]) +
      F2[3] * (P1[1] - P2[1]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5] +
      F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 *
      (F1[2] * F2[4] + F1[3] * F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 *
      (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P2[2] - P1[2]) + F2[5] * (-cI *
      (P2[0]) + cI * (P1[0]))) + (F1[3] * (F2[4] * (-cI * (P1[0]) + cI *
      (P2[0])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[0]) + cI * (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0])
      + cI * (P2[0])) + F2[3] * (P1[2] - P2[2]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5] *
      F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) +
      (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) -
      F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) + TMP10
      * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[4] * F2[2] + F1[5] * F2[3] -
      F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] * F2[5] -
      F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) + (F1[2] *
      F2[4] * (P1[0] + P2[3] - P1[3] - P2[0]) + (F1[3] * F2[5] * (P2[3] + P2[0]
      - P1[3] - P1[0]) + (F1[4] * F2[2] * (P1[3] + P1[0] - P2[3] - P2[0]) +
      F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP12 * (F1[2] * F2[5] + F1[3] * F2[4]
      + F1[4] * F2[3] + F1[5] * F2[2] - 4./3. * (P3[1] * OM3 * TMP10)) + (TMP9
      * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2] -
      4./3. * (P3[1] * OM3 * TMP10)) + (P3[1] * 2./3. * (TMP8 - TMP11) + TMP10
      * (P2[1] - P1[1])))) + P3[1] * (TMP12 * (F1[4] * F2[2] + F1[5] * F2[3] -
      F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] * F2[5] -
      F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) + (F1[2] *
      (F2[4] * (P2[1] - P1[1]) + F2[5] * (P1[0] - P2[0])) + (F1[3] * (F2[4] *
      (P1[0] - P2[0]) + F2[5] * (P2[1] - P1[1])) + (F1[4] * (F2[2] * (P1[1] -
      P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] * (P1[0] - P2[0]) +
      F2[3] * (P1[1] - P2[1]))))));
  T3[7] = denom * 2. * cI * (OM3 * (P3[1] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2] - 2./3. * (P3[1] * OM3 * TMP10)) +
      (TMP9 * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] *
      F2[2] - 2./3. * (P3[1] * OM3 * TMP10)) + (P3[1] * 1./3. * (TMP8 - TMP11)
      + TMP10 * (P2[1] - P1[1])))) + 1./3. * (TMP10 * (TMP12 - TMP9))) + (P1[1]
      * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) +
      (P2[1] * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] *
      F2[2]) + (-1./3. * (TMP11) + 1./3. * (TMP8)))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5] +
      F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 *
      - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) +
      TMP10 * (P2[1] - P1[1])))) + (F1[2] * F2[5] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (F1[3] * F2[4] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (F1[4] * F2[3] * (P1[2] - cI * (P2[1]) + cI * (P1[1])
      - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) -
      P2[2])))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5] *
      F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) +
      (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) -
      F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) + TMP10
      * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[2] * F2[5] + F1[3] * F2[4] +
      F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 * - 1. * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + TMP10 * (P2[1] - P1[1])))) +
      (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P1[3] - P2[3])) + (F1[3] *
      (F2[4] * (P1[3] - P2[3]) + F2[5] * (P2[1] - P1[1])) + (F1[4] * (F2[2] *
      (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] * (P1[3] -
      P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5] +
      F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 *
      (F1[2] * F2[4] + F1[3] * F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 *
      (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P2[2] - P1[2]) + F2[5] * (-cI *
      (P2[0]) + cI * (P1[0]))) + (F1[3] * (F2[4] * (-cI * (P1[0]) + cI *
      (P2[0])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[0]) + cI * (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0])
      + cI * (P2[0])) + F2[3] * (P1[2] - P2[2]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5] +
      F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 *
      - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) +
      TMP10 * (P2[1] - P1[1])))) + (F1[2] * F2[5] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (F1[3] * F2[4] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (F1[4] * F2[3] * (P1[2] - cI * (P2[1]) + cI * (P1[1])
      - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) -
      P2[2])))));
  T3[12] = denom * 2. * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (-cI * (F1[2] *
      F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 2./3. *
      (P3[2] * OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) +
      cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 2./3. * (P3[2] * OM3 * TMP10)) +
      (P3[2] * 1./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + 1./3. *
      (TMP10 * (TMP12 - TMP9))) + (P1[2] * (-cI * (F1[3] * F2[4] + F1[5] *
      F2[2]) + cI * (F1[2] * F2[5] + F1[4] * F2[3])) + (P2[2] * (-cI * (F1[2] *
      F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2])) + (-1./3.
      * (TMP11) + 1./3. * (TMP8)))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2])
      + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10)
      - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) +
      TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (-cI * (F1[3] * F2[4] +
      F1[5] * F2[2]) + cI * (F1[2] * F2[5] + F1[4] * F2[3])) + (TMP9 * (-cI *
      (F1[2] * F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2])) +
      TMP10 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P1[2] - P2[2]) + F2[5] *
      (-cI * (P2[3]) + cI * (P1[3]))) + (F1[3] * (F2[4] * (-cI * (P1[3]) + cI *
      (P2[3])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] * (F2[2] * (-cI * (P1[3])
      + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2])
      + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10)
      - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) +
      TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[4] * F2[2] + F1[5] *
      F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] *
      F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) +
      (F1[2] * F2[4] * (P1[0] + P2[3] - P1[3] - P2[0]) + (F1[3] * F2[5] *
      (P2[0] + P2[3] - P1[0] - P1[3]) + (F1[4] * F2[2] * (P1[0] + P1[3] - P2[0]
      - P2[3]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2])
      + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10)
      - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) +
      TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 * - 1. * (F1[2] * F2[5] +
      F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + TMP10 * (P2[1] -
      P1[1])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P1[3] - P2[3]))
      + (F1[3] * (F2[4] * (P1[3] - P2[3]) + F2[5] * (P2[1] - P1[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] *
      (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2])
      + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10)
      - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) +
      TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (-cI * (F1[3] * F2[4] +
      F1[5] * F2[2]) + cI * (F1[2] * F2[5] + F1[4] * F2[3])) + (TMP9 * (-cI *
      (F1[2] * F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2])) +
      TMP10 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P1[2] - P2[2]) + F2[5] *
      (-cI * (P2[3]) + cI * (P1[3]))) + (F1[3] * (F2[4] * (-cI * (P1[3]) + cI *
      (P2[3])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] * (F2[2] * (-cI * (P1[3])
      + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[17] = denom * 2. * cI * (OM3 * (P3[3] * (TMP12 * - 1. * (F1[3] * F2[5] +
      F1[5] * F2[3] + 2./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] *
      F2[2]) + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 2./3. * (P3[3] * OM3 *
      TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 1./3. * (TMP8 - TMP11)
      + TMP10 * (P2[3] - P1[3])))) + 1./3. * (TMP10 * (TMP12 - TMP9))) + (P1[3]
      * (F1[2] * F2[4] + F1[4] * F2[2] - F1[3] * F2[5] - F1[5] * F2[3]) +
      (P2[3] * (F1[3] * F2[5] + F1[5] * F2[3] - F1[2] * F2[4] - F1[4] * F2[2])
      + (-1./3. * (TMP11) + 1./3. * (TMP8)))));
}


void FFT5_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP16; 
  double OM3; 
  complex<double> TMP14; 
  complex<double> denom; 
  complex<double> TMP9; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP14 = (F1[2] * F2[2] + F1[3] * F2[3] + F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * TMP14 * (OM3 * (P3[0] * (P3[0] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[0] * TMP9 + P1[0] *
      TMP12))) + 1./3. * cI * (TMP9 * TMP12)) + (-1./3. * cI * (TMP16) + cI *
      (P1[0] * P2[0])));
  T3[6] = denom * TMP14 * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[1] * TMP9 + P1[1] * TMP12))) -
      P3[1] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[1] +
      P1[1] * P2[0])));
  T3[10] = denom * TMP14 * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] * TMP12))) -
      P3[2] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[2] +
      P1[2] * P2[0])));
  T3[14] = denom * TMP14 * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[3] +
      P1[3] * P2[0])));
  T3[3] = denom * TMP14 * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[1] * TMP12 + P2[1] * TMP9))) -
      P3[1] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[1] * P2[0] +
      P1[0] * P2[1])));
  T3[7] = denom * 2. * TMP14 * (OM3 * (P3[1] * (P3[1] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[1] * TMP9 + P1[1] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[1] * P2[1]) + 1./3.
      * cI * (TMP16)));
  T3[11] = denom * TMP14 * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] * TMP12))) -
      P3[2] * (+cI * (P1[1] * TMP12 + P2[1] * TMP9))) + (+cI * (P1[1] * P2[2] +
      P1[2] * P2[1])));
  T3[15] = denom * TMP14 * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[1] * TMP12 + P2[1] * TMP9))) + (+cI * (P1[1] * P2[3] +
      P1[3] * P2[1])));
  T3[4] = denom * TMP14 * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[2] * TMP12 + P2[2] * TMP9))) -
      P3[2] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[2] * P2[0] +
      P1[0] * P2[2])));
  T3[8] = denom * TMP14 * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[2] * TMP12 + P2[2] * TMP9))) -
      P3[2] * (+cI * (P2[1] * TMP9 + P1[1] * TMP12))) + (+cI * (P1[2] * P2[1] +
      P1[1] * P2[2])));
  T3[12] = denom * 2. * TMP14 * (OM3 * (P3[2] * (P3[2] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[2] * P2[2]) + 1./3.
      * cI * (TMP16)));
  T3[16] = denom * TMP14 * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[2] * TMP12 + P2[2] * TMP9))) + (+cI * (P1[2] * P2[3] +
      P1[3] * P2[2])));
  T3[5] = denom * TMP14 * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[3] * P2[0] +
      P1[0] * P2[3])));
  T3[9] = denom * TMP14 * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[1] * TMP9 + P1[1] * TMP12))) + (+cI * (P1[3] * P2[1] +
      P1[1] * P2[3])));
  T3[13] = denom * TMP14 * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[2] * TMP9 + P1[2] * TMP12))) + (+cI * (P1[3] * P2[2] +
      P1[2] * P2[3])));
  T3[17] = denom * 2. * TMP14 * (OM3 * (P3[3] * (P3[3] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[3] * P2[3]) + 1./3.
      * cI * (TMP16)));
}


void VVT3_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP71; 
  complex<double> denom; 
  complex<double> TMP73; 
  complex<double> TMP72; 
  double OM1; 
  complex<double> TMP74; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP74 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[10] * V2[5] - T3[14] * V2[4]) +
      (P3[2] * (T3[14] * V2[3] - T3[6] * V2[5]) + P3[3] * (T3[6] * V2[4] -
      T3[10] * V2[3]))) + (P1[1] * (P3[0] * (T3[14] * V2[4] - T3[10] * V2[5]) +
      (P3[2] * (T3[2] * V2[5] - T3[14] * V2[2]) + P3[3] * (T3[10] * V2[2] -
      T3[2] * V2[4]))) + (P1[2] * (P3[0] * (T3[6] * V2[5] - T3[14] * V2[3]) +
      (P3[1] * (T3[14] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[3] -
      T3[6] * V2[2]))) + P1[3] * (P3[0] * (T3[10] * V2[3] - T3[6] * V2[4]) +
      (P3[1] * (T3[2] * V2[4] - T3[10] * V2[2]) + P3[2] * (T3[6] * V2[2] -
      T3[2] * V2[3])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[15] * V2[4] - T3[11]
      * V2[5]) + (P3[2] * (T3[7] * V2[5] - T3[15] * V2[3]) + P3[3] * (T3[11] *
      V2[3] - T3[7] * V2[4]))) + (P1[1] * (P3[0] * (T3[11] * V2[5] - T3[15] *
      V2[4]) + (P3[2] * (T3[15] * V2[2] - T3[3] * V2[5]) + P3[3] * (T3[3] *
      V2[4] - T3[11] * V2[2]))) + (P1[2] * (P3[0] * (T3[15] * V2[3] - T3[7] *
      V2[5]) + (P3[1] * (T3[3] * V2[5] - T3[15] * V2[2]) + P3[3] * (T3[7] *
      V2[2] - T3[3] * V2[3]))) + P1[3] * (P3[0] * (T3[7] * V2[4] - T3[11] *
      V2[3]) + (P3[1] * (T3[11] * V2[2] - T3[3] * V2[4]) + P3[2] * (T3[3] *
      V2[3] - T3[7] * V2[2])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[16] * V2[4]
      - T3[12] * V2[5]) + (P3[2] * (T3[8] * V2[5] - T3[16] * V2[3]) + P3[3] *
      (T3[12] * V2[3] - T3[8] * V2[4]))) + (P1[1] * (P3[0] * (T3[12] * V2[5] -
      T3[16] * V2[4]) + (P3[2] * (T3[16] * V2[2] - T3[4] * V2[5]) + P3[3] *
      (T3[4] * V2[4] - T3[12] * V2[2]))) + (P1[2] * (P3[0] * (T3[16] * V2[3] -
      T3[8] * V2[5]) + (P3[1] * (T3[4] * V2[5] - T3[16] * V2[2]) + P3[3] *
      (T3[8] * V2[2] - T3[4] * V2[3]))) + P1[3] * (P3[0] * (T3[8] * V2[4] -
      T3[12] * V2[3]) + (P3[1] * (T3[12] * V2[2] - T3[4] * V2[4]) + P3[2] *
      (T3[4] * V2[3] - T3[8] * V2[2])))))) + P2[3] * (P1[0] * (P3[1] * (T3[17]
      * V2[4] - T3[13] * V2[5]) + (P3[2] * (T3[9] * V2[5] - T3[17] * V2[3]) +
      P3[3] * (T3[13] * V2[3] - T3[9] * V2[4]))) + (P1[1] * (P3[0] * (T3[13] *
      V2[5] - T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[2] - T3[5] * V2[5]) +
      P3[3] * (T3[5] * V2[4] - T3[13] * V2[2]))) + (P1[2] * (P3[0] * (T3[17] *
      V2[3] - T3[9] * V2[5]) + (P3[1] * (T3[5] * V2[5] - T3[17] * V2[2]) +
      P3[3] * (T3[9] * V2[2] - T3[5] * V2[3]))) + P1[3] * (P3[0] * (T3[9] *
      V2[4] - T3[13] * V2[3]) + (P3[1] * (T3[13] * V2[2] - T3[5] * V2[4]) +
      P3[2] * (T3[5] * V2[3] - T3[9] * V2[2])))))))));
  TMP73 = -1. * (P1[0] * (P1[1] * (P3[2] * (V2[5] * (T3[7] + T3[2]) + (-T3[15]
      * V2[3] - T3[14] * V2[2])) + (P3[3] * (V2[4] * - 1. * (T3[7] + T3[2]) +
      (T3[11] * V2[3] + T3[10] * V2[2])) + (P3[0] * (T3[14] * V2[4] - T3[10] *
      V2[5]) + P3[1] * (T3[15] * V2[4] - T3[11] * V2[5])))) + (P1[2] * (P3[1] *
      (V2[5] * - 1. * (T3[12] + T3[2]) + (T3[16] * V2[4] + T3[14] * V2[2])) +
      (P3[3] * (V2[3] * (T3[12] + T3[2]) + (-T3[8] * V2[4] - T3[6] * V2[2])) +
      (P3[0] * (T3[6] * V2[5] - T3[14] * V2[3]) + P3[2] * (T3[8] * V2[5] -
      T3[16] * V2[3])))) + (P1[3] * (P3[1] * (V2[4] * (T3[17] + T3[2]) +
      (-T3[13] * V2[5] - T3[10] * V2[2])) + (P3[2] * (V2[3] * - 1. * (T3[17] +
      T3[2]) + (T3[9] * V2[5] + T3[6] * V2[2])) + (P3[0] * (T3[10] * V2[3] -
      T3[6] * V2[4]) + P3[3] * (T3[13] * V2[3] - T3[9] * V2[4])))) + P1[0] *
      (P3[1] * (T3[10] * V2[5] - T3[14] * V2[4]) + (P3[2] * (T3[14] * V2[3] -
      T3[6] * V2[5]) + P3[3] * (T3[6] * V2[4] - T3[10] * V2[3])))))) + (P1[1] *
      (P1[2] * (P3[0] * (V2[5] * (T3[12] - T3[7]) + (T3[15] * V2[3] - T3[16] *
      V2[4])) + (P3[3] * (V2[2] * (T3[7] - T3[12]) + (T3[4] * V2[4] - T3[3] *
      V2[3])) + (P3[1] * (T3[3] * V2[5] - T3[15] * V2[2]) + P3[2] * (T3[16] *
      V2[2] - T3[4] * V2[5])))) + (P1[3] * (P3[0] * (V2[4] * (T3[7] - T3[17]) +
      (T3[13] * V2[5] - T3[11] * V2[3])) + (P3[2] * (V2[2] * (T3[17] - T3[7]) +
      (T3[3] * V2[3] - T3[5] * V2[5])) + (P3[1] * (T3[11] * V2[2] - T3[3] *
      V2[4]) + P3[3] * (T3[5] * V2[4] - T3[13] * V2[2])))) + P1[1] * (P3[0] *
      (T3[11] * V2[5] - T3[15] * V2[4]) + (P3[2] * (T3[15] * V2[2] - T3[3] *
      V2[5]) + P3[3] * (T3[3] * V2[4] - T3[11] * V2[2]))))) + (P1[2] * (P1[3] *
      (P3[0] * (V2[3] * (T3[17] - T3[12]) + (T3[8] * V2[4] - T3[9] * V2[5])) +
      (P3[1] * (V2[2] * (T3[12] - T3[17]) + (T3[5] * V2[5] - T3[4] * V2[4])) +
      (P3[2] * (T3[4] * V2[3] - T3[8] * V2[2]) + P3[3] * (T3[9] * V2[2] - T3[5]
      * V2[3])))) + P1[2] * (P3[0] * (T3[16] * V2[3] - T3[8] * V2[5]) + (P3[1]
      * (T3[4] * V2[5] - T3[16] * V2[2]) + P3[3] * (T3[8] * V2[2] - T3[4] *
      V2[3])))) + P1[3] * P1[3] * (P3[0] * (T3[9] * V2[4] - T3[13] * V2[3]) +
      (P3[1] * (T3[13] * V2[2] - T3[5] * V2[4]) + P3[2] * (T3[5] * V2[3] -
      T3[9] * V2[2]))))));
  TMP72 = -1. * (P2[0] * (P1[0] * (P3[1] * (T3[4] * V2[5] - T3[5] * V2[4]) +
      (P3[2] * (T3[5] * V2[3] - T3[3] * V2[5]) + P3[3] * (T3[3] * V2[4] - T3[4]
      * V2[3]))) + (P1[1] * (P3[0] * (T3[5] * V2[4] - T3[4] * V2[5]) + (P3[2] *
      (T3[2] * V2[5] - T3[5] * V2[2]) + P3[3] * (T3[4] * V2[2] - T3[2] *
      V2[4]))) + (P1[2] * (P3[0] * (T3[3] * V2[5] - T3[5] * V2[3]) + (P3[1] *
      (T3[5] * V2[2] - T3[2] * V2[5]) + P3[3] * (T3[2] * V2[3] - T3[3] *
      V2[2]))) + P1[3] * (P3[0] * (T3[4] * V2[3] - T3[3] * V2[4]) + (P3[1] *
      (T3[2] * V2[4] - T3[4] * V2[2]) + P3[2] * (T3[3] * V2[2] - T3[2] *
      V2[3])))))) + (P2[1] * (P1[0] * (P3[1] * (T3[9] * V2[4] - T3[8] * V2[5])
      + (P3[2] * (T3[7] * V2[5] - T3[9] * V2[3]) + P3[3] * (T3[8] * V2[3] -
      T3[7] * V2[4]))) + (P1[1] * (P3[0] * (T3[8] * V2[5] - T3[9] * V2[4]) +
      (P3[2] * (T3[9] * V2[2] - T3[6] * V2[5]) + P3[3] * (T3[6] * V2[4] - T3[8]
      * V2[2]))) + (P1[2] * (P3[0] * (T3[9] * V2[3] - T3[7] * V2[5]) + (P3[1] *
      (T3[6] * V2[5] - T3[9] * V2[2]) + P3[3] * (T3[7] * V2[2] - T3[6] *
      V2[3]))) + P1[3] * (P3[0] * (T3[7] * V2[4] - T3[8] * V2[3]) + (P3[1] *
      (T3[8] * V2[2] - T3[6] * V2[4]) + P3[2] * (T3[6] * V2[3] - T3[7] *
      V2[2])))))) + (P2[2] * (P1[0] * (P3[1] * (T3[13] * V2[4] - T3[12] *
      V2[5]) + (P3[2] * (T3[11] * V2[5] - T3[13] * V2[3]) + P3[3] * (T3[12] *
      V2[3] - T3[11] * V2[4]))) + (P1[1] * (P3[0] * (T3[12] * V2[5] - T3[13] *
      V2[4]) + (P3[2] * (T3[13] * V2[2] - T3[10] * V2[5]) + P3[3] * (T3[10] *
      V2[4] - T3[12] * V2[2]))) + (P1[2] * (P3[0] * (T3[13] * V2[3] - T3[11] *
      V2[5]) + (P3[1] * (T3[10] * V2[5] - T3[13] * V2[2]) + P3[3] * (T3[11] *
      V2[2] - T3[10] * V2[3]))) + P1[3] * (P3[0] * (T3[11] * V2[4] - T3[12] *
      V2[3]) + (P3[1] * (T3[12] * V2[2] - T3[10] * V2[4]) + P3[2] * (T3[10] *
      V2[3] - T3[11] * V2[2])))))) + P2[3] * (P1[0] * (P3[1] * (T3[17] * V2[4]
      - T3[16] * V2[5]) + (P3[2] * (T3[15] * V2[5] - T3[17] * V2[3]) + P3[3] *
      (T3[16] * V2[3] - T3[15] * V2[4]))) + (P1[1] * (P3[0] * (T3[16] * V2[5] -
      T3[17] * V2[4]) + (P3[2] * (T3[17] * V2[2] - T3[14] * V2[5]) + P3[3] *
      (T3[14] * V2[4] - T3[16] * V2[2]))) + (P1[2] * (P3[0] * (T3[17] * V2[3] -
      T3[15] * V2[5]) + (P3[1] * (T3[14] * V2[5] - T3[17] * V2[2]) + P3[3] *
      (T3[15] * V2[2] - T3[14] * V2[3]))) + P1[3] * (P3[0] * (T3[15] * V2[4] -
      T3[16] * V2[3]) + (P3[1] * (T3[16] * V2[2] - T3[14] * V2[4]) + P3[2] *
      (T3[14] * V2[3] - T3[15] * V2[2])))))))));
  TMP71 = -1. * (P1[0] * (P1[1] * (P3[2] * (V2[5] * (T3[7] + T3[2]) + (-T3[9] *
      V2[3] - T3[5] * V2[2])) + (P3[3] * (V2[4] * - 1. * (T3[7] + T3[2]) +
      (T3[8] * V2[3] + T3[4] * V2[2])) + (P3[0] * (T3[5] * V2[4] - T3[4] *
      V2[5]) + P3[1] * (T3[9] * V2[4] - T3[8] * V2[5])))) + (P1[2] * (P3[1] *
      (V2[5] * - 1. * (T3[12] + T3[2]) + (T3[13] * V2[4] + T3[5] * V2[2])) +
      (P3[3] * (V2[3] * (T3[12] + T3[2]) + (-T3[11] * V2[4] - T3[3] * V2[2])) +
      (P3[0] * (T3[3] * V2[5] - T3[5] * V2[3]) + P3[2] * (T3[11] * V2[5] -
      T3[13] * V2[3])))) + (P1[3] * (P3[1] * (V2[4] * (T3[17] + T3[2]) +
      (-T3[16] * V2[5] - T3[4] * V2[2])) + (P3[2] * (V2[3] * - 1. * (T3[17] +
      T3[2]) + (T3[15] * V2[5] + T3[3] * V2[2])) + (P3[0] * (T3[4] * V2[3] -
      T3[3] * V2[4]) + P3[3] * (T3[16] * V2[3] - T3[15] * V2[4])))) + P1[0] *
      (P3[1] * (T3[4] * V2[5] - T3[5] * V2[4]) + (P3[2] * (T3[5] * V2[3] -
      T3[3] * V2[5]) + P3[3] * (T3[3] * V2[4] - T3[4] * V2[3])))))) + (P1[1] *
      (P1[2] * (P3[0] * (V2[5] * (T3[12] - T3[7]) + (T3[9] * V2[3] - T3[13] *
      V2[4])) + (P3[3] * (V2[2] * (T3[7] - T3[12]) + (T3[10] * V2[4] - T3[6] *
      V2[3])) + (P3[1] * (T3[6] * V2[5] - T3[9] * V2[2]) + P3[2] * (T3[13] *
      V2[2] - T3[10] * V2[5])))) + (P1[3] * (P3[0] * (V2[4] * (T3[7] - T3[17])
      + (T3[16] * V2[5] - T3[8] * V2[3])) + (P3[2] * (V2[2] * (T3[17] - T3[7])
      + (T3[6] * V2[3] - T3[14] * V2[5])) + (P3[1] * (T3[8] * V2[2] - T3[6] *
      V2[4]) + P3[3] * (T3[14] * V2[4] - T3[16] * V2[2])))) + P1[1] * (P3[0] *
      (T3[8] * V2[5] - T3[9] * V2[4]) + (P3[2] * (T3[9] * V2[2] - T3[6] *
      V2[5]) + P3[3] * (T3[6] * V2[4] - T3[8] * V2[2]))))) + (P1[2] * (P1[3] *
      (P3[0] * (V2[3] * (T3[17] - T3[12]) + (T3[11] * V2[4] - T3[15] * V2[5]))
      + (P3[1] * (V2[2] * (T3[12] - T3[17]) + (T3[14] * V2[5] - T3[10] *
      V2[4])) + (P3[2] * (T3[10] * V2[3] - T3[11] * V2[2]) + P3[3] * (T3[15] *
      V2[2] - T3[14] * V2[3])))) + P1[2] * (P3[0] * (T3[13] * V2[3] - T3[11] *
      V2[5]) + (P3[1] * (T3[10] * V2[5] - T3[13] * V2[2]) + P3[3] * (T3[11] *
      V2[2] - T3[10] * V2[3])))) + P1[3] * P1[3] * (P3[0] * (T3[15] * V2[4] -
      T3[16] * V2[3]) + (P3[1] * (T3[16] * V2[2] - T3[14] * V2[4]) + P3[2] *
      (T3[14] * V2[3] - T3[15] * V2[2]))))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * 2. * cI * (P3[1] * (V2[4] * (P1[0] * - 1./2. * (T3[5] +
      T3[14]) + (P1[1] * 1./2. * (T3[9] + T3[15]) + (P1[2] * 1./2. * (T3[13] +
      T3[16]) + (P2[0] * 1./2. * (T3[5] + T3[14]) + (P2[1] * - 1./2. * (T3[9] +
      T3[15]) + (P2[2] * - 1./2. * (T3[13] + T3[16]) + T3[17] * (P1[3] -
      P2[3]))))))) + V2[5] * (P1[0] * 1./2. * (T3[4] + T3[10]) + (P1[1] * -
      1./2. * (T3[8] + T3[11]) + (P1[3] * - 1./2. * (T3[16] + T3[13]) + (P2[0]
      * - 1./2. * (T3[4] + T3[10]) + (P2[1] * 1./2. * (T3[8] + T3[11]) + (P2[3]
      * 1./2. * (T3[16] + T3[13]) + T3[12] * (P2[2] - P1[2])))))))) + (P3[2] *
      (V2[3] * (P1[0] * 1./2. * (T3[5] + T3[14]) + (P1[1] * - 1./2. * (T3[9] +
      T3[15]) + (P1[2] * - 1./2. * (T3[13] + T3[16]) + (P2[0] * - 1./2. *
      (T3[5] + T3[14]) + (P2[1] * 1./2. * (T3[9] + T3[15]) + (P2[2] * 1./2. *
      (T3[13] + T3[16]) + T3[17] * (P2[3] - P1[3]))))))) + V2[5] * (P1[0] * -
      1./2. * (T3[3] + T3[6]) + (P1[2] * 1./2. * (T3[11] + T3[8]) + (P1[3] *
      1./2. * (T3[15] + T3[9]) + (P2[0] * 1./2. * (T3[3] + T3[6]) + (P2[2] * -
      1./2. * (T3[11] + T3[8]) + (P2[3] * - 1./2. * (T3[15] + T3[9]) + T3[7] *
      (P1[1] - P2[1])))))))) + (+1./2. * (OM1 * P1[0] * (TMP71 + TMP73 - TMP72
      - TMP74)) + P3[3] * (V2[3] * (P1[0] * - 1./2. * (T3[4] + T3[10]) + (P1[1]
      * 1./2. * (T3[8] + T3[11]) + (P1[3] * 1./2. * (T3[16] + T3[13]) + (P2[0]
      * 1./2. * (T3[4] + T3[10]) + (P2[1] * - 1./2. * (T3[8] + T3[11]) + (P2[3]
      * - 1./2. * (T3[16] + T3[13]) + T3[12] * (P1[2] - P2[2]))))))) + V2[4] *
      (P1[0] * 1./2. * (T3[3] + T3[6]) + (P1[2] * - 1./2. * (T3[11] + T3[8]) +
      (P1[3] * - 1./2. * (T3[15] + T3[9]) + (P2[0] * - 1./2. * (T3[3] + T3[6])
      + (P2[2] * 1./2. * (T3[11] + T3[8]) + (P2[3] * 1./2. * (T3[15] + T3[9]) +
      T3[7] * (P2[1] - P1[1])))))))))));
  V1[3] = denom * - 2. * cI * (P3[0] * (V2[4] * (P1[0] * 1./2. * (T3[5] +
      T3[14]) + (P1[1] * - 1./2. * (T3[9] + T3[15]) + (P1[2] * - 1./2. *
      (T3[13] + T3[16]) + (P2[0] * - 1./2. * (T3[5] + T3[14]) + (P2[1] * 1./2.
      * (T3[9] + T3[15]) + (P2[2] * 1./2. * (T3[13] + T3[16]) + T3[17] * (P2[3]
      - P1[3]))))))) + V2[5] * (P1[0] * - 1./2. * (T3[4] + T3[10]) + (P1[1] *
      1./2. * (T3[8] + T3[11]) + (P1[3] * 1./2. * (T3[16] + T3[13]) + (P2[0] *
      1./2. * (T3[4] + T3[10]) + (P2[1] * - 1./2. * (T3[8] + T3[11]) + (P2[3] *
      - 1./2. * (T3[16] + T3[13]) + T3[12] * (P1[2] - P2[2])))))))) + (P3[2] *
      (V2[2] * (P1[0] * - 1./2. * (T3[5] + T3[14]) + (P1[1] * 1./2. * (T3[9] +
      T3[15]) + (P1[2] * 1./2. * (T3[13] + T3[16]) + (P2[0] * 1./2. * (T3[5] +
      T3[14]) + (P2[1] * - 1./2. * (T3[9] + T3[15]) + (P2[2] * - 1./2. *
      (T3[13] + T3[16]) + T3[17] * (P1[3] - P2[3]))))))) + V2[5] * (P1[1] * -
      1./2. * (T3[6] + T3[3]) + (P1[2] * - 1./2. * (T3[10] + T3[4]) + (P1[3] *
      - 1./2. * (T3[14] + T3[5]) + (P2[1] * 1./2. * (T3[6] + T3[3]) + (P2[2] *
      1./2. * (T3[10] + T3[4]) + (P2[3] * 1./2. * (T3[14] + T3[5]) + T3[2] *
      (P1[0] - P2[0])))))))) + (+1./2. * (OM1 * P1[1] * (TMP72 + TMP74 - TMP71
      - TMP73)) + P3[3] * (V2[2] * (P1[0] * 1./2. * (T3[4] + T3[10]) + (P1[1] *
      - 1./2. * (T3[8] + T3[11]) + (P1[3] * - 1./2. * (T3[16] + T3[13]) +
      (P2[0] * - 1./2. * (T3[4] + T3[10]) + (P2[1] * 1./2. * (T3[8] + T3[11]) +
      (P2[3] * 1./2. * (T3[16] + T3[13]) + T3[12] * (P2[2] - P1[2]))))))) +
      V2[4] * (P1[1] * 1./2. * (T3[6] + T3[3]) + (P1[2] * 1./2. * (T3[10] +
      T3[4]) + (P1[3] * 1./2. * (T3[14] + T3[5]) + (P2[1] * - 1./2. * (T3[6] +
      T3[3]) + (P2[2] * - 1./2. * (T3[10] + T3[4]) + (P2[3] * - 1./2. * (T3[14]
      + T3[5]) + T3[2] * (P2[0] - P1[0])))))))))));
  V1[4] = denom * - 2. * cI * (P3[0] * (V2[3] * (P1[0] * - 1./2. * (T3[5] +
      T3[14]) + (P1[1] * 1./2. * (T3[9] + T3[15]) + (P1[2] * 1./2. * (T3[13] +
      T3[16]) + (P2[0] * 1./2. * (T3[5] + T3[14]) + (P2[1] * - 1./2. * (T3[9] +
      T3[15]) + (P2[2] * - 1./2. * (T3[13] + T3[16]) + T3[17] * (P1[3] -
      P2[3]))))))) + V2[5] * (P1[0] * 1./2. * (T3[3] + T3[6]) + (P1[2] * -
      1./2. * (T3[11] + T3[8]) + (P1[3] * - 1./2. * (T3[15] + T3[9]) + (P2[0] *
      - 1./2. * (T3[3] + T3[6]) + (P2[2] * 1./2. * (T3[11] + T3[8]) + (P2[3] *
      1./2. * (T3[15] + T3[9]) + T3[7] * (P2[1] - P1[1])))))))) + (P3[1] *
      (V2[2] * (P1[0] * 1./2. * (T3[5] + T3[14]) + (P1[1] * - 1./2. * (T3[9] +
      T3[15]) + (P1[2] * - 1./2. * (T3[13] + T3[16]) + (P2[0] * - 1./2. *
      (T3[5] + T3[14]) + (P2[1] * 1./2. * (T3[9] + T3[15]) + (P2[2] * 1./2. *
      (T3[13] + T3[16]) + T3[17] * (P2[3] - P1[3]))))))) + V2[5] * (P1[1] *
      1./2. * (T3[6] + T3[3]) + (P1[2] * 1./2. * (T3[10] + T3[4]) + (P1[3] *
      1./2. * (T3[14] + T3[5]) + (P2[1] * - 1./2. * (T3[6] + T3[3]) + (P2[2] *
      - 1./2. * (T3[10] + T3[4]) + (P2[3] * - 1./2. * (T3[14] + T3[5]) + T3[2]
      * (P2[0] - P1[0])))))))) + (+1./2. * (OM1 * P1[2] * (TMP72 + TMP74 -
      TMP71 - TMP73)) + P3[3] * (V2[2] * (P1[0] * - 1./2. * (T3[3] + T3[6]) +
      (P1[2] * 1./2. * (T3[11] + T3[8]) + (P1[3] * 1./2. * (T3[15] + T3[9]) +
      (P2[0] * 1./2. * (T3[3] + T3[6]) + (P2[2] * - 1./2. * (T3[11] + T3[8]) +
      (P2[3] * - 1./2. * (T3[15] + T3[9]) + T3[7] * (P1[1] - P2[1]))))))) +
      V2[3] * (P1[1] * - 1./2. * (T3[6] + T3[3]) + (P1[2] * - 1./2. * (T3[10] +
      T3[4]) + (P1[3] * - 1./2. * (T3[14] + T3[5]) + (P2[1] * 1./2. * (T3[6] +
      T3[3]) + (P2[2] * 1./2. * (T3[10] + T3[4]) + (P2[3] * 1./2. * (T3[14] +
      T3[5]) + T3[2] * (P1[0] - P2[0])))))))))));
  V1[5] = denom * - 2. * cI * (P3[0] * (V2[3] * (P1[0] * 1./2. * (T3[4] +
      T3[10]) + (P1[1] * - 1./2. * (T3[8] + T3[11]) + (P1[3] * - 1./2. *
      (T3[16] + T3[13]) + (P2[0] * - 1./2. * (T3[4] + T3[10]) + (P2[1] * 1./2.
      * (T3[8] + T3[11]) + (P2[3] * 1./2. * (T3[16] + T3[13]) + T3[12] * (P2[2]
      - P1[2]))))))) + V2[4] * (P1[0] * - 1./2. * (T3[3] + T3[6]) + (P1[2] *
      1./2. * (T3[11] + T3[8]) + (P1[3] * 1./2. * (T3[15] + T3[9]) + (P2[0] *
      1./2. * (T3[3] + T3[6]) + (P2[2] * - 1./2. * (T3[11] + T3[8]) + (P2[3] *
      - 1./2. * (T3[15] + T3[9]) + T3[7] * (P1[1] - P2[1])))))))) + (P3[1] *
      (V2[2] * (P1[0] * - 1./2. * (T3[4] + T3[10]) + (P1[1] * 1./2. * (T3[8] +
      T3[11]) + (P1[3] * 1./2. * (T3[16] + T3[13]) + (P2[0] * 1./2. * (T3[4] +
      T3[10]) + (P2[1] * - 1./2. * (T3[8] + T3[11]) + (P2[3] * - 1./2. *
      (T3[16] + T3[13]) + T3[12] * (P1[2] - P2[2]))))))) + V2[4] * (P1[1] * -
      1./2. * (T3[6] + T3[3]) + (P1[2] * - 1./2. * (T3[10] + T3[4]) + (P1[3] *
      - 1./2. * (T3[14] + T3[5]) + (P2[1] * 1./2. * (T3[6] + T3[3]) + (P2[2] *
      1./2. * (T3[10] + T3[4]) + (P2[3] * 1./2. * (T3[14] + T3[5]) + T3[2] *
      (P1[0] - P2[0])))))))) + (+1./2. * (OM1 * P1[3] * (TMP72 + TMP74 - TMP71
      - TMP73)) + P3[2] * (V2[2] * (P1[0] * 1./2. * (T3[3] + T3[6]) + (P1[2] *
      - 1./2. * (T3[11] + T3[8]) + (P1[3] * - 1./2. * (T3[15] + T3[9]) + (P2[0]
      * - 1./2. * (T3[3] + T3[6]) + (P2[2] * 1./2. * (T3[11] + T3[8]) + (P2[3]
      * 1./2. * (T3[15] + T3[9]) + T3[7] * (P2[1] - P1[1]))))))) + V2[3] *
      (P1[1] * 1./2. * (T3[6] + T3[3]) + (P1[2] * 1./2. * (T3[10] + T3[4]) +
      (P1[3] * 1./2. * (T3[14] + T3[5]) + (P2[1] * - 1./2. * (T3[6] + T3[3]) +
      (P2[2] * - 1./2. * (T3[10] + T3[4]) + (P2[3] * - 1./2. * (T3[14] + T3[5])
      + T3[2] * (P2[0] - P1[0])))))))))));
}


void VVT10_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP39; 
  double P3[4]; 
  complex<double> TMP9; 
  complex<double> TMP43; 
  double P2[4]; 
  complex<double> TMP26; 
  complex<double> TMP15; 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP45; 
  complex<double> TMP16; 
  complex<double> denom; 
  complex<double> TMP13; 
  complex<double> TMP38; 
  complex<double> TMP44; 
  complex<double> TMP40; 
  double OM1; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP43 = (P2[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P2[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P2[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P2[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP44 = (P2[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P2[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P2[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P2[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP45 = (P1[0] * (P1[1] * - 1. * (T3[3] + T3[6]) + (P1[2] * - 1. * (T3[4] +
      T3[10]) + (P1[3] * - 1. * (T3[5] + T3[14]) + P1[0] * T3[2]))) + (P1[1] *
      (P1[2] * (T3[8] + T3[11]) + (P1[3] * (T3[9] + T3[15]) + P1[1] * T3[7])) +
      (P1[2] * (P1[3] * (T3[13] + T3[16]) + P1[2] * T3[12]) + P1[3] * P1[3] *
      T3[17])));
  TMP40 = (P2[0] * (P2[1] * - 1. * (T3[6] + T3[3]) + (P2[2] * - 1. * (T3[10] +
      T3[4]) + (P2[3] * - 1. * (T3[14] + T3[5]) + P2[0] * T3[2]))) + (P2[1] *
      (P2[2] * (T3[11] + T3[8]) + (P2[3] * (T3[15] + T3[9]) + P2[1] * T3[7])) +
      (P2[2] * (P2[3] * (T3[16] + T3[13]) + P2[2] * T3[12]) + P2[3] * P2[3] *
      T3[17])));
  TMP39 = (P1[0] * (P1[1] * - 1. * (T3[6] + T3[3]) + (P1[2] * - 1. * (T3[10] +
      T3[4]) + (P1[3] * - 1. * (T3[14] + T3[5]) + P1[0] * T3[2]))) + (P1[1] *
      (P1[2] * (T3[11] + T3[8]) + (P1[3] * (T3[15] + T3[9]) + P1[1] * T3[7])) +
      (P1[2] * (P1[3] * (T3[16] + T3[13]) + P1[2] * T3[12]) + P1[3] * P1[3] *
      T3[17])));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * 1./2. * (TMP16 * (TMP38 * (P1[0] * (OM1 * (-cI * (TMP45) + cI
      * (TMP39)) + 2. * cI * (T3[2])) + (P1[1] * - 1. * (+cI * (T3[3] + T3[6]))
      + (P1[2] * - 1. * (+cI * (T3[4] + T3[10])) - P1[3] * (+cI * (T3[5] +
      T3[14]))))) + (V2[2] * 2. * (-cI * (TMP13 + TMP15) + cI * (TMP39 +
      TMP40)) + P3[0] * (+cI * (TMP43 + TMP44)))) + (TMP26 * (TMP12 * (P1[0] *
      (OM1 * (-cI * (TMP39) + cI * (TMP45)) - 2. * cI * (T3[2])) + (P1[1] *
      (+cI * (T3[3] + T3[6])) + (P1[2] * (+cI * (T3[4] + T3[10])) + P1[3] *
      (+cI * (T3[5] + T3[14]))))) + (P2[0] * 2. * (-cI * (TMP39 + TMP40) + cI *
      (TMP13 + TMP15)) - 2. * cI * (P3[0] * TMP40))) + (P2[0] * (TMP9 * - 1. *
      (+cI * (TMP43 + TMP44)) - 2. * cI * (TMP38 * TMP39)) + 2. * (V2[2] * (+cI
      * (TMP12 * TMP39 + TMP9 * TMP40))))));
  V1[3] = denom * 1./2. * (TMP16 * (TMP38 * (P1[1] * (OM1 * (-cI * (TMP45) + cI
      * (TMP39)) - 2. * cI * (T3[7])) + (P1[0] * (+cI * (T3[6] + T3[3])) +
      (P1[2] * - 1. * (+cI * (T3[8] + T3[11])) - P1[3] * (+cI * (T3[9] +
      T3[15]))))) + (V2[3] * 2. * (-cI * (TMP13 + TMP15) + cI * (TMP39 +
      TMP40)) + P3[1] * (+cI * (TMP43 + TMP44)))) + (TMP26 * (TMP12 * (P1[1] *
      (OM1 * (-cI * (TMP39) + cI * (TMP45)) + 2. * cI * (T3[7])) + (P1[0] * -
      1. * (+cI * (T3[6] + T3[3])) + (P1[2] * (+cI * (T3[8] + T3[11])) + P1[3]
      * (+cI * (T3[9] + T3[15]))))) + (P2[1] * 2. * (-cI * (TMP39 + TMP40) + cI
      * (TMP13 + TMP15)) - 2. * cI * (P3[1] * TMP40))) + (P2[1] * (TMP9 * - 1.
      * (+cI * (TMP43 + TMP44)) - 2. * cI * (TMP38 * TMP39)) + 2. * (V2[3] *
      (+cI * (TMP12 * TMP39 + TMP9 * TMP40))))));
  V1[4] = denom * 1./2. * (TMP16 * (TMP38 * (P1[2] * (OM1 * (-cI * (TMP45) + cI
      * (TMP39)) - 2. * cI * (T3[12])) + (P1[0] * (+cI * (T3[10] + T3[4])) +
      (P1[1] * - 1. * (+cI * (T3[11] + T3[8])) - P1[3] * (+cI * (T3[13] +
      T3[16]))))) + (V2[4] * 2. * (-cI * (TMP13 + TMP15) + cI * (TMP39 +
      TMP40)) + P3[2] * (+cI * (TMP43 + TMP44)))) + (TMP26 * (TMP12 * (P1[2] *
      (OM1 * (-cI * (TMP39) + cI * (TMP45)) + 2. * cI * (T3[12])) + (P1[0] * -
      1. * (+cI * (T3[10] + T3[4])) + (P1[1] * (+cI * (T3[11] + T3[8])) + P1[3]
      * (+cI * (T3[13] + T3[16]))))) + (P2[2] * 2. * (-cI * (TMP39 + TMP40) +
      cI * (TMP13 + TMP15)) - 2. * cI * (P3[2] * TMP40))) + (P2[2] * (TMP9 * -
      1. * (+cI * (TMP43 + TMP44)) - 2. * cI * (TMP38 * TMP39)) + 2. * (V2[4] *
      (+cI * (TMP12 * TMP39 + TMP9 * TMP40))))));
  V1[5] = denom * 1./2. * (TMP16 * (TMP38 * (P1[3] * (OM1 * (-cI * (TMP45) + cI
      * (TMP39)) - 2. * cI * (T3[17])) + (P1[0] * (+cI * (T3[14] + T3[5])) +
      (P1[1] * - 1. * (+cI * (T3[15] + T3[9])) - P1[2] * (+cI * (T3[16] +
      T3[13]))))) + (V2[5] * 2. * (-cI * (TMP13 + TMP15) + cI * (TMP39 +
      TMP40)) + P3[3] * (+cI * (TMP43 + TMP44)))) + (TMP26 * (TMP12 * (P1[3] *
      (OM1 * (-cI * (TMP39) + cI * (TMP45)) + 2. * cI * (T3[17])) + (P1[0] * -
      1. * (+cI * (T3[14] + T3[5])) + (P1[1] * (+cI * (T3[15] + T3[9])) + P1[2]
      * (+cI * (T3[16] + T3[13]))))) + (P2[3] * 2. * (-cI * (TMP39 + TMP40) +
      cI * (TMP13 + TMP15)) - 2. * cI * (P3[3] * TMP40))) + (P2[3] * (TMP9 * -
      1. * (+cI * (TMP43 + TMP44)) - 2. * cI * (TMP38 * TMP39)) + 2. * (V2[5] *
      (+cI * (TMP12 * TMP39 + TMP9 * TMP40))))));
}

void VVT10_11_12_13_2_3_6_7_8_9_1(complex<double> V2[], complex<double> T3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M1, double W1, complex<double> V1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
  complex<double> Vtmp[6]; 
//   double P2[4]; 
//   double P1[4]; 
  complex<double> denom; 
  int i; 
//   double OM1; 
  VVT10_1(V2, T3, COUP1, M1, W1, V1); 
  VVT11_1(V2, T3, COUP2, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT12_1(V2, T3, COUP3, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT13_1(V2, T3, COUP4, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT2_1(V2, T3, COUP5, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT3_1(V2, T3, COUP6, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT6_1(V2, T3, COUP7, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT7_1(V2, T3, COUP8, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT8_1(V2, T3, COUP9, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVT9_1(V2, T3, COUP10, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
}

void FFT2_2(complex<double> F1[], complex<double> T3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + T3[0]; 
  F2[1] = +F1[1] + T3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[3] * (P2[0] * (P2[3] * (T3[9] + T3[6] + T3[15] +
      T3[3] - cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P2[1] * - 1. *
      (T3[14] + T3[5] + cI * (T3[11] + T3[8]) - 2. * (T3[7] + T3[2])) + (P2[2]
      * (T3[8] + T3[11] - 2. * cI * (T3[12] + T3[2]) + cI * (T3[14] + T3[5])) +
      (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P1[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0] * (+cI * (T3[10] + T3[4]) -
      T3[6] - T3[3]) + (P1[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P1[2]
      * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[16] + T3[13] - cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6]) - T3[4] -
      T3[10]) + (P2[3] * (+2. * (T3[17]) + cI * (T3[11] + T3[8]) - 2. * (T3[7])
      - T3[5] - T3[14]) + (P1[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P1[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] + T3[9] - T3[3] -
      T3[6]) + (P1[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P1[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P2[2] * (P2[3] * (-2. * cI * (T3[17]) +
      cI * (T3[5] + T3[14]) + 2. * cI * (T3[12]) - T3[8] - T3[11]) + (P1[1] *
      (-cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])) + (P1[2] * (-cI * (T3[4]
      + T3[10]) + cI * (T3[16] + T3[13])) + (P2[2] * (-cI * (T3[16] + T3[13]) +
      cI * (T3[4] + T3[10])) + (P1[0] * - 1. * (-2. * cI * (T3[2]) + cI *
      (T3[14] + T3[5])) - P1[3] * (-2. * cI * (T3[17]) + cI * (T3[5] +
      T3[14]))))))) + P2[3] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3])
      + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. * (+cI * (T3[11] +
      T3[8]) - 2. * (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F1[2] * (P2[0] * (P2[1] * (T3[15] + T3[9] + 2. *
      (T3[3] + T3[6]) + cI * (T3[10] + T3[4])) + (P2[2] * (T3[16] + T3[13] + 2.
      * (T3[4] + T3[10]) - cI * (T3[6] + T3[3])) + (P1[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (P2[3] * 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P1[0] *
      (T3[14] + T3[5] + 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] + 2. *
      (T3[17])) - P2[0] * (T3[14] + T3[5] + 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * - 1. * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) + cI * (T3[13] +
      T3[16])) + (P2[2] * - 1. * (T3[4] + T3[10] - cI * (T3[9] + T3[15]) + 2. *
      (T3[13] + T3[16])) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[2] *
      (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * - 1. * (T3[14] + T3[5] + 2.
      * (T3[2])) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[17])) - P2[3] * (T3[5] +
      T3[14] + 2. * (T3[17])))))))) + (P2[1] * (P1[0] * - 1. * (T3[6] + T3[3] +
      cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[7]) + cI *
      (T3[12])) + (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P1[2] *
      (T3[8] + T3[11] + 2. * cI * (T3[12])) - P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8]))))))) + P2[2] * (P1[0] * (+cI * (T3[6] + T3[3]) - T3[10]
      - T3[4]) + (P1[3] * (T3[13] + T3[16] - cI * (T3[9] + T3[15])) + (P1[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) + (P1[2] * - 1. * (+cI * (T3[8] +
      T3[11]) - 2. * (T3[12])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M2 * (F1[4] * (P1[1] * (T3[15] + T3[9] - T3[3] -
      T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[1] * (T3[3] +
      T3[6] - T3[15] - T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) +
      (T3[14] * (P2[0] + P2[3] - P1[0] - P1[3]) + (T3[5] * (P2[3] + P2[0] -
      P1[3] - P1[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] * (P1[3] -
      P2[3]))))))))) + F1[5] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] -
      T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[0] *
      (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] * (+cI * (T3[13] +
      T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) - P2[2])
      + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] * (P1[1] -
      P2[1]))))))))))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (P2[3] * (T3[9] + T3[15] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[6] - T3[3]) + (P2[1] *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8])) + (P2[2]
      * (T3[8] + T3[11] + cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) +
      (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI
      * (T3[10] + T3[4])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] +
      T3[8])) - P1[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P2[1] *
      (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13] + cI * (T3[3] + T3[15]
      + T3[6] + T3[9])) + (P2[3] * (+2. * (T3[7]) + cI * (T3[11] + T3[8]) - 2.
      * (T3[17]) - T3[5] - T3[14]) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[0] * - 1. * (T3[14] + T3[5] + 2. * (T3[2]))
      + P1[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[2] * (P2[3] *
      (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[17]) + 2. * cI *
      (T3[12])) + (P1[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P1[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[2] * - 1. * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P1[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P1[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P2[3] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4]))
      + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P2[3] *
      (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P1[1] * - 1. * (+2. *
      (T3[7]) + cI * (T3[11] + T3[8])) - P1[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))))) + (F1[3] * (P2[0] * (P2[1] * - 1. * (T3[15] + T3[9] + cI
      * (T3[10] + T3[4]) - 2. * (T3[3] + T3[6])) + (P2[2] * (+2. * (T3[4] +
      T3[10]) + cI * (T3[6] + T3[3]) - T3[16] - T3[13]) + (P1[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (P2[3] * 2. * (T3[5] + T3[14] - T3[17] - T3[2]) + (P1[0] * - 1. * (T3[14]
      + T3[5] - 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] - 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * (T3[3] + T3[6] + cI * (T3[13] + T3[16]) - 2. * (T3[9] + T3[15]))
      + (P2[2] * - 1. * (+2. * (T3[13] + T3[16]) + cI * (T3[9] + T3[15]) -
      T3[4] - T3[10]) + (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P1[2] *
      (T3[16] + T3[13] - T3[4] - T3[10]) + (P1[0] * - 1. * (T3[14] + T3[5] - 2.
      * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[1] * (P1[0] * (+cI * (T3[10]
      + T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] +
      T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[12]) + cI *
      (T3[7])) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      (P1[2] * (T3[8] + T3[11] - 2. * cI * (T3[12])) + P2[1] * (+cI * (T3[11] +
      T3[8]) - 2. * (T3[7]))))))) + P2[2] * (P1[0] * - 1. * (T3[10] + T3[4] +
      cI * (T3[6] + T3[3])) + (P1[3] * (T3[13] + T3[16] + cI * (T3[9] +
      T3[15])) + (P1[1] * (T3[11] + T3[8] + 2. * cI * (T3[7])) + (P1[2] * (+2.
      * (T3[12]) + cI * (T3[8] + T3[11])) - P2[2] * (+2. * (T3[12]) + cI *
      (T3[8] + T3[11]))))))))) + M2 * (F1[4] * (P1[0] * - 1. * (T3[6] + T3[3] +
      cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * -
      1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P1[2] - cI *
      (P2[1]) + cI * (P1[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P2[1]) + cI *
      (P1[1]) - P2[2]) + (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. *
      (T3[7] * (P1[1] - P2[1]))))))))) + F1[5] * (P1[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] +
      T3[16] + T3[10] + T3[13]) + (T3[14] * (P1[0] + P2[3] - P2[0] - P1[3]) +
      (T3[5] * (P2[3] + P1[0] - P1[3] - P2[0]) + (T3[2] * 2. * (P1[0] - P2[0])
      + 2. * (T3[17] * (P2[3] - P1[3]))))))))))));
  F2[4] = denom * cI * (F1[5] * (P2[0] * (P2[3] * (T3[6] + T3[3] - cI * (T3[10]
      + T3[4]) + cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[1] * - 1. *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) - cI * (T3[11] + T3[8])) + (P2[2]
      * (+cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2]) - T3[8] - T3[11])
      + (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] +
      T3[15] - cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10]
      + T3[4])) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      P1[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[4] + T3[16] + T3[10] + T3[13] - cI * (T3[3] + T3[15] + T3[6] +
      T3[9])) + (P2[3] * (T3[5] + T3[14] + 2. * (T3[17]) + cI * (T3[11] +
      T3[8]) - 2. * (T3[7])) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P1[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[2] * (P2[3] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17])) + (P1[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P1[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[2] * - 1. * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P1[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P1[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P2[3] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] -
      T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[3] *
      (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. * (+cI *
      (T3[11] + T3[8]) - 2. * (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F1[4] * (P2[0] * (P2[1] * (+2. * (T3[3] + T3[6]) + cI
      * (T3[10] + T3[4]) - T3[15] - T3[9]) + (P2[2] * - 1. * (T3[16] + T3[13] +
      cI * (T3[6] + T3[3]) - 2. * (T3[4] + T3[10])) + (P1[1] * (T3[15] + T3[9]
      - T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[3] *
      2. * (T3[5] + T3[14] - T3[17] - T3[2]) + (P1[0] * - 1. * (T3[14] + T3[5]
      - 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) +
      P2[0] * (T3[14] + T3[5] - 2. * (T3[2]))))))))) + (P2[3] * (P2[1] * - 1. *
      (+2. * (T3[9] + T3[15]) + cI * (T3[13] + T3[16]) - T3[3] - T3[6]) +
      (P2[2] * (T3[4] + T3[10] + cI * (T3[9] + T3[15]) - 2. * (T3[13] +
      T3[16])) + (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P1[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P1[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[1] * (P1[0] * - 1. * (T3[6]
      + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[7]) + cI
      * (T3[12])) + (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P1[2] *
      (T3[8] + T3[11] + 2. * cI * (T3[12])) - P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8]))))))) + P2[2] * (P1[0] * (+cI * (T3[6] + T3[3]) - T3[10]
      - T3[4]) + (P1[3] * (T3[13] + T3[16] - cI * (T3[9] + T3[15])) + (P1[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) + (P1[2] * - 1. * (+cI * (T3[8] +
      T3[11]) - 2. * (T3[12])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M2 * (F1[2] * (P1[1] * - 1. * (T3[3] + T3[15] + T3[6]
      + T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] *
      (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (T3[14] * (P1[0] + P2[3] - P2[0] - P1[3]) + (T3[5] * (P2[3] +
      P1[0] - P1[3] - P2[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] *
      (P2[3] - P1[3]))))))))) + F1[3] * (P1[0] * (T3[6] + T3[3] - cI * (T3[10]
      + T3[4])) + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0]
      * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] -
      cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P2[1]) + cI * (P1[1])
      - P1[2]) + (T3[8] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) - P1[2]) +
      (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] * (P2[1] -
      P1[1]))))))))))));
  F2[5] = denom * - cI * (F1[4] * (P2[0] * (P2[3] * (T3[9] + T3[6] + T3[15] +
      T3[3] + cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P2[1] * (+2. * (T3[7]
      + T3[2]) + cI * (T3[11] + T3[8]) - T3[14] - T3[5]) + (P2[2] * (T3[8] +
      T3[11] - cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) + (P1[0] *
      (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15]
      + cI * (T3[13] + T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) -
      P1[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[16] + T3[13] - cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9]) - T3[4] -
      T3[10]) + (P2[3] * - 1. * (T3[5] + T3[14] + 2. * (T3[7]) + cI * (T3[11] +
      T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15] - T3[9]) +
      (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] + T3[9] -
      T3[3] - T3[6]) + (P1[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P1[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[2] * (P2[3] * - 1. * (T3[8]
      + T3[11] - 2. * cI * (T3[17]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[12])) + (P1[1] * (-cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6])) +
      (P1[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] + T3[10])) + (P2[2] *
      (-cI * (T3[4] + T3[10]) + cI * (T3[16] + T3[13])) + (P1[0] * (-2. * cI *
      (T3[2]) + cI * (T3[14] + T3[5])) + P1[3] * (-2. * cI * (T3[17]) + cI *
      (T3[5] + T3[14]))))))) + P2[3] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P1[1] * (+2.
      * (T3[7]) + cI * (T3[11] + T3[8])) + P1[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))))) + (F1[5] * (P2[0] * (P2[1] * - 1. * (T3[15] + T3[9] + 2.
      * (T3[3] + T3[6]) - cI * (T3[10] + T3[4])) + (P2[2] * - 1. * (T3[16] +
      T3[13] + 2. * (T3[4] + T3[10]) + cI * (T3[6] + T3[3])) + (P1[1] * (T3[3]
      + T3[15] + T3[6] + T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (P2[3] * - 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P1[0] * - 1. *
      (T3[14] + T3[5] + 2. * (T3[2])) + (P1[3] * (T3[5] + T3[14] + 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] + 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) - cI * (T3[13] + T3[16]))
      + (P2[2] * (T3[4] + T3[10] + 2. * (T3[13] + T3[16]) + cI * (T3[9] +
      T3[15])) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[2] * -
      1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] + 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[1] * (P1[0] * (T3[6] + T3[3]
      - cI * (T3[10] + T3[4])) + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P2[2] * 2. * (T3[8] + T3[11] - cI * (T3[12]) + cI * (T3[7])) +
      (P1[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + (P1[2] * - 1. * (T3[8]
      + T3[11] - 2. * cI * (T3[12])) - P2[1] * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7]))))))) + P2[2] * (P1[0] * (T3[10] + T3[4] + cI * (T3[6] + T3[3]))
      + (P1[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] + T3[15])) + (P1[1] * -
      1. * (T3[11] + T3[8] + 2. * cI * (T3[7])) + (P1[2] * - 1. * (+2. *
      (T3[12]) + cI * (T3[8] + T3[11])) + P2[2] * (+2. * (T3[12]) + cI * (T3[8]
      + T3[11]))))))))) + M2 * (F1[2] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P2[1]) + cI * (P1[1]) -
      P2[2]) + (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] *
      (P1[1] - P2[1]))))))))) + F1[3] * (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (T3[14] * (P1[0] + P1[3] - P2[0] - P2[3]) + (T3[5] * (P1[3] + P1[0] -
      P2[3] - P2[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))))));
}


void VVT12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  double OM3; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP16; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP25 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP30 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (P3[0] * (OM3 * (TMP12 * 2./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 2./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-2./3. * cI * (TMP26 * TMP30) + 2./3.
      * cI * (TMP16 * TMP25))) + (P1[0] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[2] * TMP37 + V1[2] * TMP38)) + (+cI * (TMP12 *
      V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) + (TMP12 * 1./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 1./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37))))) + (TMP16 * (-2./3. * cI * (TMP25) + cI
      * (V2[2] * V1[2])) + (TMP26 * (-cI * (P2[0] * V1[2]) + 2./3. * cI *
      (TMP30)) + P1[0] * (-cI * (V2[2] * TMP30) + cI * (P2[0] * TMP25)))));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[1] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[3] * TMP37 + V1[3] * TMP38)) + (+cI * (TMP12 *
      V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) + P3[1] * (P1[0] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V1[2] * TMP38 + V2[2] *
      TMP37)) + (+cI * (TMP12 * V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) +
      (P1[0] * (-cI * (V2[3] * TMP30) + cI * (P2[1] * TMP25)) + (P1[1] * (-cI *
      (V2[2] * TMP30) + cI * (P2[0] * TMP25)) + (TMP16 * (+cI * (V2[3] * V1[2]
      + V2[2] * V1[3])) - TMP26 * (+cI * (P2[0] * V1[3] + P2[1] * V1[2]))))));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[2] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[4] * TMP37 + V1[4] * TMP38)) + (+cI * (TMP12 *
      V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) + P3[2] * (P1[0] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V1[2] * TMP38 + V2[2] *
      TMP37)) + (+cI * (TMP12 * V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) +
      (P1[0] * (-cI * (V2[4] * TMP30) + cI * (P2[2] * TMP25)) + (P1[2] * (-cI *
      (V2[2] * TMP30) + cI * (P2[0] * TMP25)) + (TMP16 * (+cI * (V2[4] * V1[2]
      + V2[2] * V1[4])) - TMP26 * (+cI * (P2[0] * V1[4] + P2[2] * V1[2]))))));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[5] * TMP37 + V1[5] * TMP38)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + P3[3] * (P1[0] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V1[2] * TMP38 + V2[2] *
      TMP37)) + (+cI * (TMP12 * V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) +
      (P1[0] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)) + (P1[3] * (-cI *
      (V2[2] * TMP30) + cI * (P2[0] * TMP25)) + (TMP16 * (+cI * (V2[5] * V1[2]
      + V2[2] * V1[5])) - TMP26 * (+cI * (P2[0] * V1[5] + P2[3] * V1[2]))))));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[1] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V1[3] * TMP38 + V2[3] * TMP37)) + (+cI * (TMP12 *
      V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) + P3[1] * (P1[0] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V2[2] * TMP37 + V1[2] *
      TMP38)) + (+cI * (TMP12 * V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) +
      (P1[0] * (-cI * (V2[3] * TMP30) + cI * (P2[1] * TMP25)) + (P1[1] * (-cI *
      (V2[2] * TMP30) + cI * (P2[0] * TMP25)) + (TMP16 * (+cI * (V2[2] * V1[3]
      + V2[3] * V1[2])) - TMP26 * (+cI * (P2[1] * V1[2] + P2[0] * V1[3]))))));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (P3[1] * (OM3 * (TMP12 * 2./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 2./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-2./3. * cI * (TMP26 * TMP30) + 2./3.
      * cI * (TMP16 * TMP25))) + (P1[1] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[3] * TMP37 + V1[3] * TMP38)) + (+cI * (TMP12 *
      V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) + (TMP12 * 1./3. * (-cI *
      (TMP9 * TMP25) + cI * (TMP26 * TMP37)) + 1./3. * (TMP38 * (-cI * (TMP16 *
      TMP37) + cI * (TMP9 * TMP30))))) + (TMP16 * (+cI * (V2[3] * V1[3]) +
      2./3. * cI * (TMP25)) + (TMP26 * - 1. * (+cI * (P2[1] * V1[3]) + 2./3. *
      cI * (TMP30)) + P1[1] * (-cI * (V2[3] * TMP30) + cI * (P2[1] * TMP25)))));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[2] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[4] * TMP37 + V1[4] * TMP38)) + (+cI * (TMP12 *
      V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) + P3[2] * (P1[1] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V1[3] * TMP38 + V2[3] *
      TMP37)) + (+cI * (TMP12 * V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) +
      (P1[1] * (-cI * (V2[4] * TMP30) + cI * (P2[2] * TMP25)) + (P1[2] * (-cI *
      (V2[3] * TMP30) + cI * (P2[1] * TMP25)) + (TMP16 * (+cI * (V2[4] * V1[3]
      + V2[3] * V1[4])) - TMP26 * (+cI * (P2[1] * V1[4] + P2[2] * V1[3]))))));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[5] * TMP37 + V1[5] * TMP38)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + P3[3] * (P1[1] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V1[3] * TMP38 + V2[3] *
      TMP37)) + (+cI * (TMP12 * V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) +
      (P1[1] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)) + (P1[3] * (-cI *
      (V2[3] * TMP30) + cI * (P2[1] * TMP25)) + (TMP16 * (+cI * (V2[5] * V1[3]
      + V2[3] * V1[5])) - TMP26 * (+cI * (P2[1] * V1[5] + P2[3] * V1[3]))))));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[2] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V1[4] * TMP38 + V2[4] * TMP37)) + (+cI * (TMP12 *
      V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) + P3[2] * (P1[0] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V2[2] * TMP37 + V1[2] *
      TMP38)) + (+cI * (TMP12 * V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) +
      (P1[0] * (-cI * (V2[4] * TMP30) + cI * (P2[2] * TMP25)) + (P1[2] * (-cI *
      (V2[2] * TMP30) + cI * (P2[0] * TMP25)) + (TMP16 * (+cI * (V2[2] * V1[4]
      + V2[4] * V1[2])) - TMP26 * (+cI * (P2[2] * V1[2] + P2[0] * V1[4]))))));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[2] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V1[4] * TMP38 + V2[4] * TMP37)) + (+cI * (TMP12 *
      V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) + P3[2] * (P1[1] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V2[3] * TMP37 + V1[3] *
      TMP38)) + (+cI * (TMP12 * V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) +
      (P1[1] * (-cI * (V2[4] * TMP30) + cI * (P2[2] * TMP25)) + (P1[2] * (-cI *
      (V2[3] * TMP30) + cI * (P2[1] * TMP25)) + (TMP16 * (+cI * (V2[3] * V1[4]
      + V2[4] * V1[3])) - TMP26 * (+cI * (P2[2] * V1[3] + P2[1] * V1[4]))))));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (P3[2] * (OM3 * (TMP12 * 2./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 2./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-2./3. * cI * (TMP26 * TMP30) + 2./3.
      * cI * (TMP16 * TMP25))) + (P1[2] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[4] * TMP37 + V1[4] * TMP38)) + (+cI * (TMP12 *
      V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) + (TMP12 * 1./3. * (-cI *
      (TMP9 * TMP25) + cI * (TMP26 * TMP37)) + 1./3. * (TMP38 * (-cI * (TMP16 *
      TMP37) + cI * (TMP9 * TMP30))))) + (TMP16 * (+cI * (V2[4] * V1[4]) +
      2./3. * cI * (TMP25)) + (TMP26 * - 1. * (+cI * (P2[2] * V1[4]) + 2./3. *
      cI * (TMP30)) + P1[2] * (-cI * (V2[4] * TMP30) + cI * (P2[2] * TMP25)))));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[5] * TMP37 + V1[5] * TMP38)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + P3[3] * (P1[2] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V1[4] * TMP38 + V2[4] *
      TMP37)) + (+cI * (TMP12 * V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) +
      (P1[2] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)) + (P1[3] * (-cI *
      (V2[4] * TMP30) + cI * (P2[2] * TMP25)) + (TMP16 * (+cI * (V2[5] * V1[4]
      + V2[4] * V1[5])) - TMP26 * (+cI * (P2[2] * V1[5] + P2[3] * V1[4]))))));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V1[5] * TMP38 + V2[5] * TMP37)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + P3[3] * (P1[0] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[0] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V2[2] * TMP37 + V1[2] *
      TMP38)) + (+cI * (TMP12 * V1[2] * TMP26 + TMP9 * V2[2] * TMP30)))))) +
      (P1[0] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)) + (P1[3] * (-cI *
      (V2[2] * TMP30) + cI * (P2[0] * TMP25)) + (TMP16 * (+cI * (V2[2] * V1[5]
      + V2[5] * V1[2])) - TMP26 * (+cI * (P2[3] * V1[2] + P2[0] * V1[5]))))));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V1[5] * TMP38 + V2[5] * TMP37)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + P3[3] * (P1[1] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[1] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V2[3] * TMP37 + V1[3] *
      TMP38)) + (+cI * (TMP12 * V1[3] * TMP26 + TMP9 * V2[3] * TMP30)))))) +
      (P1[1] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)) + (P1[3] * (-cI *
      (V2[3] * TMP30) + cI * (P2[1] * TMP25)) + (TMP16 * (+cI * (V2[3] * V1[5]
      + V2[5] * V1[3])) - TMP26 * (+cI * (P2[3] * V1[3] + P2[1] * V1[5]))))));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * (OM3 * (TMP12 * 4./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 4./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-4./3. * cI * (TMP26 * TMP30) + 4./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V1[5] * TMP38 + V2[5] * TMP37)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + P3[3] * (P1[2] * (-cI *
      (TMP12 * TMP25) + cI * (TMP30 * TMP38)) + (P2[2] * (-cI * (TMP9 * TMP25)
      + cI * (TMP26 * TMP37)) + (TMP16 * - 1. * (+cI * (V2[4] * TMP37 + V1[4] *
      TMP38)) + (+cI * (TMP12 * V1[4] * TMP26 + TMP9 * V2[4] * TMP30)))))) +
      (P1[2] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)) + (P1[3] * (-cI *
      (V2[4] * TMP30) + cI * (P2[2] * TMP25)) + (TMP16 * (+cI * (V2[4] * V1[5]
      + V2[5] * V1[4])) - TMP26 * (+cI * (P2[3] * V1[4] + P2[2] * V1[5]))))));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (P3[3] * (OM3 * (TMP12 * 2./3. * (-cI *
      (TMP26 * TMP37) + cI * (TMP9 * TMP25)) + 2./3. * (TMP38 * (-cI * (TMP9 *
      TMP30) + cI * (TMP16 * TMP37)))) + (-2./3. * cI * (TMP26 * TMP30) + 2./3.
      * cI * (TMP16 * TMP25))) + (P1[3] * (-cI * (TMP12 * TMP25) + cI * (TMP30
      * TMP38)) + (P2[3] * (-cI * (TMP9 * TMP25) + cI * (TMP26 * TMP37)) +
      (TMP16 * - 1. * (+cI * (V2[5] * TMP37 + V1[5] * TMP38)) + (+cI * (TMP12 *
      V1[5] * TMP26 + TMP9 * V2[5] * TMP30)))))) + (TMP12 * 1./3. * (-cI *
      (TMP9 * TMP25) + cI * (TMP26 * TMP37)) + 1./3. * (TMP38 * (-cI * (TMP16 *
      TMP37) + cI * (TMP9 * TMP30))))) + (TMP16 * (+cI * (V2[5] * V1[5]) +
      2./3. * cI * (TMP25)) + (TMP26 * - 1. * (+cI * (P2[3] * V1[5]) + 2./3. *
      cI * (TMP30)) + P1[3] * (-cI * (V2[5] * TMP30) + cI * (P2[3] * TMP25)))));
}

void FFV5_6_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV5_2(F1, V3, COUP1, M2, W2, F2); 
  FFV6_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV5_8_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFV5_2(F1, V3, COUP1, M2, W2, F2); 
  FFV8_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void VVT11_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP37; 
  complex<double> TMP38; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP25; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP25 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (P3[0] * (P3[0] * 1./3. * (+cI * (TMP25) + 2. *
      cI * (OM3 * TMP37 * TMP38)) + (-cI * (V2[2] * TMP37 + V1[2] * TMP38))) +
      1./3. * cI * (TMP37 * TMP38)) + (-1./3. * cI * (TMP25) + cI * (V2[2] *
      V1[2])));
  T3[6] = denom * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V2[3] * TMP37 + V1[3] * TMP38))) - P3[1]
      * (+cI * (V1[2] * TMP38 + V2[2] * TMP37))) + (+cI * (V2[3] * V1[2] +
      V2[2] * V1[3])));
  T3[10] = denom * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V2[4] * TMP37 + V1[4] * TMP38))) - P3[2]
      * (+cI * (V1[2] * TMP38 + V2[2] * TMP37))) + (+cI * (V2[4] * V1[2] +
      V2[2] * V1[4])));
  T3[14] = denom * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V2[5] * TMP37 + V1[5] * TMP38))) - P3[3]
      * (+cI * (V1[2] * TMP38 + V2[2] * TMP37))) + (+cI * (V2[5] * V1[2] +
      V2[2] * V1[5])));
  T3[3] = denom * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V1[3] * TMP38 + V2[3] * TMP37))) - P3[1]
      * (+cI * (V2[2] * TMP37 + V1[2] * TMP38))) + (+cI * (V2[2] * V1[3] +
      V2[3] * V1[2])));
  T3[7] = denom * 2. * (OM3 * (P3[1] * (P3[1] * 1./3. * (+cI * (TMP25) + 2. *
      cI * (OM3 * TMP37 * TMP38)) + (-cI * (V2[3] * TMP37 + V1[3] * TMP38))) -
      1./3. * cI * (TMP37 * TMP38)) + (+cI * (V2[3] * V1[3]) + 1./3. * cI *
      (TMP25)));
  T3[11] = denom * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V2[4] * TMP37 + V1[4] * TMP38))) - P3[2]
      * (+cI * (V1[3] * TMP38 + V2[3] * TMP37))) + (+cI * (V2[4] * V1[3] +
      V2[3] * V1[4])));
  T3[15] = denom * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V2[5] * TMP37 + V1[5] * TMP38))) - P3[3]
      * (+cI * (V1[3] * TMP38 + V2[3] * TMP37))) + (+cI * (V2[5] * V1[3] +
      V2[3] * V1[5])));
  T3[4] = denom * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V1[4] * TMP38 + V2[4] * TMP37))) - P3[2]
      * (+cI * (V2[2] * TMP37 + V1[2] * TMP38))) + (+cI * (V2[2] * V1[4] +
      V2[4] * V1[2])));
  T3[8] = denom * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V1[4] * TMP38 + V2[4] * TMP37))) - P3[2]
      * (+cI * (V2[3] * TMP37 + V1[3] * TMP38))) + (+cI * (V2[3] * V1[4] +
      V2[4] * V1[3])));
  T3[12] = denom * 2. * (OM3 * (P3[2] * (P3[2] * 1./3. * (+cI * (TMP25) + 2. *
      cI * (OM3 * TMP37 * TMP38)) + (-cI * (V2[4] * TMP37 + V1[4] * TMP38))) -
      1./3. * cI * (TMP37 * TMP38)) + (+cI * (V2[4] * V1[4]) + 1./3. * cI *
      (TMP25)));
  T3[16] = denom * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V2[5] * TMP37 + V1[5] * TMP38))) - P3[3]
      * (+cI * (V1[4] * TMP38 + V2[4] * TMP37))) + (+cI * (V2[5] * V1[4] +
      V2[4] * V1[5])));
  T3[5] = denom * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V1[5] * TMP38 + V2[5] * TMP37))) - P3[3]
      * (+cI * (V2[2] * TMP37 + V1[2] * TMP38))) + (+cI * (V2[2] * V1[5] +
      V2[5] * V1[2])));
  T3[9] = denom * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V1[5] * TMP38 + V2[5] * TMP37))) - P3[3]
      * (+cI * (V2[3] * TMP37 + V1[3] * TMP38))) + (+cI * (V2[3] * V1[5] +
      V2[5] * V1[3])));
  T3[13] = denom * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP25) + 2. * cI *
      (OM3 * TMP37 * TMP38)) + (-cI * (V1[5] * TMP38 + V2[5] * TMP37))) - P3[3]
      * (+cI * (V2[4] * TMP37 + V1[4] * TMP38))) + (+cI * (V2[4] * V1[5] +
      V2[5] * V1[4])));
  T3[17] = denom * 2. * (OM3 * (P3[3] * (P3[3] * 1./3. * (+cI * (TMP25) + 2. *
      cI * (OM3 * TMP37 * TMP38)) + (-cI * (V2[5] * TMP37 + V1[5] * TMP38))) -
      1./3. * cI * (TMP37 * TMP38)) + (+cI * (V2[5] * V1[5]) + 1./3. * cI *
      (TMP25)));
}


void VVT1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP16; 
  complex<double> TMP66; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP65; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP65 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP66 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 8. * (OM3 * (P3[0] * (P3[0] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[0] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[0] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP65) + cI * (TMP66)))) + (P1[0] * P2[0] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP66) + cI * (TMP65)))));
  T3[6] = denom * 4. * (OM3 * (P3[0] * (P3[1] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[1] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[1] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[1] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[1] * (-cI * (TMP65) + cI * (TMP66)) + P1[1] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[10] = denom * 4. * (OM3 * (P3[0] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[14] = denom * 4. * (OM3 * (P3[0] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[3] = denom * 4. * (OM3 * (P3[0] * (P3[1] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[1] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[1] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[1] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[1] * (-cI * (TMP65) + cI * (TMP66)) + P1[1] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[7] = denom * 8. * (OM3 * (P3[1] * (P3[1] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[1] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[1] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP66) + cI * (TMP65)))) + (P1[1] * P2[1] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI * (TMP66)))));
  T3[11] = denom * 4. * (OM3 * (P3[1] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[15] = denom * 4. * (OM3 * (P3[1] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[4] = denom * 4. * (OM3 * (P3[0] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[8] = denom * 4. * (OM3 * (P3[1] * (P3[2] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[2] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[2] * (-cI * (TMP65) + cI * (TMP66)) + P1[2] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[12] = denom * 8. * (OM3 * (P3[2] * (P3[2] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[2] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[2] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP66) + cI * (TMP65)))) + (P1[2] * P2[2] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI * (TMP66)))));
  T3[16] = denom * 4. * (OM3 * (P3[2] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[2] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[2] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[2] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[2] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[5] = denom * 4. * (OM3 * (P3[0] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[0] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[0] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[0] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[0] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[9] = denom * 4. * (OM3 * (P3[1] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[1] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[1] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[1] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[1] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[13] = denom * 4. * (OM3 * (P3[2] * (P3[3] * (OM3 * 4./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 2./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + P3[3] * (P1[2] * TMP12 * (-cI *
      (TMP66) + cI * (TMP65)) + P2[2] * TMP9 * (-cI * (TMP66) + cI * (TMP65))))
      + (P1[2] * P2[3] * (-cI * (TMP65) + cI * (TMP66)) + P1[3] * P2[2] * (-cI
      * (TMP65) + cI * (TMP66))));
  T3[17] = denom * 8. * (OM3 * (P3[3] * (P3[3] * (OM3 * 2./3. * TMP12 * TMP9 *
      (-cI * (TMP65) + cI * (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI *
      (TMP66)))) + (P1[3] * TMP12 * (-cI * (TMP66) + cI * (TMP65)) + P2[3] *
      TMP9 * (-cI * (TMP66) + cI * (TMP65)))) + 1./3. * (TMP12 * TMP9 * (-cI *
      (TMP66) + cI * (TMP65)))) + (P1[3] * P2[3] * (-cI * (TMP65) + cI *
      (TMP66)) + 1./3. * (TMP16 * (-cI * (TMP65) + cI * (TMP66)))));
}

void VVT1_10_11_12_13_3_5_7_8_9_3(complex<double> V1[], complex<double> V2[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
//   double P2[4]; 
//   double OM3; 
//   double P1[4]; 
  complex<double> Ttmp[18]; 
  complex<double> denom; 
  int i; 
  VVT1_3(V1, V2, COUP1, M3, W3, T3); 
  VVT10_3(V1, V2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT11_3(V1, V2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT12_3(V1, V2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT13_3(V1, V2, COUP5, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT3_3(V1, V2, COUP6, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT5_3(V1, V2, COUP7, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT7_3(V1, V2, COUP8, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT8_3(V1, V2, COUP9, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  VVT9_3(V1, V2, COUP10, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}

void FFV8_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  complex<double> TMP0; 
  TMP1 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - 1. * (+cI * (TMP0) + 4. * cI * (TMP1)); 
}


void VVT2_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP13; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * - 2. * cI * (TMP13 * (P1[1] * (P2[3] * V2[4] - P2[2] * V2[5])
      + (P1[2] * (P2[1] * V2[5] - P2[3] * V2[3]) + P1[3] * (P2[2] * V2[3] -
      P2[1] * V2[4]))) + TMP15 * (P1[1] * (P2[3] * V2[4] - P2[2] * V2[5]) +
      (P1[2] * (P2[1] * V2[5] - P2[3] * V2[3]) + P1[3] * (P2[2] * V2[3] - P2[1]
      * V2[4]))));
  V1[3] = denom * 2. * cI * (TMP13 * (P1[0] * (P2[2] * V2[5] - P2[3] * V2[4]) +
      (P1[2] * (P2[3] * V2[2] - P2[0] * V2[5]) + P1[3] * (P2[0] * V2[4] - P2[2]
      * V2[2]))) + TMP15 * (P1[0] * (P2[2] * V2[5] - P2[3] * V2[4]) + (P1[2] *
      (P2[3] * V2[2] - P2[0] * V2[5]) + P1[3] * (P2[0] * V2[4] - P2[2] *
      V2[2]))));
  V1[4] = denom * 2. * cI * (TMP13 * (P1[0] * (P2[3] * V2[3] - P2[1] * V2[5]) +
      (P1[1] * (P2[0] * V2[5] - P2[3] * V2[2]) + P1[3] * (P2[1] * V2[2] - P2[0]
      * V2[3]))) + TMP15 * (P1[0] * (P2[3] * V2[3] - P2[1] * V2[5]) + (P1[1] *
      (P2[0] * V2[5] - P2[3] * V2[2]) + P1[3] * (P2[1] * V2[2] - P2[0] *
      V2[3]))));
  V1[5] = denom * 2. * cI * (TMP13 * (P1[0] * (P2[1] * V2[4] - P2[2] * V2[3]) +
      (P1[1] * (P2[2] * V2[2] - P2[0] * V2[4]) + P1[2] * (P2[0] * V2[3] - P2[1]
      * V2[2]))) + TMP15 * (P1[0] * (P2[1] * V2[4] - P2[2] * V2[3]) + (P1[1] *
      (P2[2] * V2[2] - P2[0] * V2[4]) + P1[2] * (P2[0] * V2[3] - P2[1] *
      V2[2]))));
}

void FFT1_2_4_5_0(complex<double> F1[], complex<double> F2[], complex<double>
    T3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> tmp; 
  FFT1_0(F1, F2, T3, COUP1, vertex); 
  FFT2_0(F1, F2, T3, COUP2, tmp); 
  vertex = vertex + tmp; 
  FFT4_0(F1, F2, T3, COUP3, tmp); 
  vertex = vertex + tmp; 
  FFT5_0(F1, F2, T3, COUP4, tmp); 
  vertex = vertex + tmp; 
}

void FFV6_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
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
  F1[2] = denom * - 2. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + (+1./2. * (M1 * (+2. * (F2[4] * - 1./2. * (V3[2] + V3[5]))
      - F2[5] * (V3[3] + cI * (V3[4])))) + F2[3] * (P1[0] * (V3[3] + cI *
      (V3[4])) + (P1[1] * - 1. * (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI *
      (V3[2] + V3[5])) + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * - 2. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1]
      * (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] *
      (+cI * (V3[4]) - V3[3])))) + (+1./2. * (M1 * (F2[5] * (V3[5] - V3[2]) +
      2. * (F2[4] * 1./2. * (+cI * (V3[4]) - V3[3])))) + F2[3] * (P1[0] * - 1.
      * (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] -
      cI * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * cI * (F2[4] * (P1[0] * - 1. * (V3[2] + V3[5]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[2]
      + V3[5])))) + (F2[5] * (P1[0] * - 1. * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 2. * (V3[5] - V3[2]) + 2. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + M1 * (F2[2] * 2. * (+cI * (V3[4]) - V3[3]) + 2. *
      (F2[3] * (V3[2] + V3[5])))));
}


void FFT5_1(complex<double> F2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP13; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + T3[0]; 
  F1[1] = +F2[1] + T3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (TMP13 * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) - F2[2] * M1)) + TMP15 * (F2[4] * (P1[0] + P1[3]) + (F2[5]
      * (P1[1] + cI * (P1[2])) - F2[2] * M1)));
  F1[3] = denom * cI * (TMP13 * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) + F2[3] * M1)) + TMP15 * (F2[4] * (+cI * (P1[2]) - P1[1])
      + (F2[5] * (P1[3] - P1[0]) + F2[3] * M1)));
  F1[4] = denom * cI * (TMP13 * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) + F2[4] * M1)) + TMP15 * (F2[2] * (P1[3] - P1[0]) + (F2[3] *
      (P1[1] + cI * (P1[2])) + F2[4] * M1)));
  F1[5] = denom * - cI * (TMP13 * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) - F2[5] * M1)) + TMP15 * (F2[2] * (+cI * (P1[2]) - P1[1])
      + (F2[3] * (P1[0] + P1[3]) - F2[5] * M1)));
}


void VVT3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP76; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP75; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP76 = -1. * (P2[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P2[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P2[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P2[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP75 = -1. * (P1[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (OM3 * P3[0] * (TMP12 * (P3[1] * (V2[4] * V1[5] -
      V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] *
      (V2[3] * V1[4] - V2[4] * V1[3]))) + (TMP9 * (P3[1] * (V2[5] * V1[4] -
      V2[4] * V1[5]) + (P3[2] * (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] *
      (V2[4] * V1[3] - V2[3] * V1[4]))) + 1./3. * (P3[0] * (TMP76 - TMP75)))) +
      (P1[0] * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] *
      V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] * V1[3]))) +
      (P2[0] * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[3] *
      V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] * V1[4]))) +
      (-1./3. * (TMP76) + 1./3. * (TMP75)))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP12 * (P3[0] * (V2[5] * V1[4] - V2[4]
      * V1[5]) + (P3[2] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] *
      V1[2] - V2[2] * V1[4]))) + (TMP9 * (P3[0] * (V2[4] * V1[5] - V2[5] *
      V1[4]) + (P3[2] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[4] - V2[4] * V1[2]))) + 2./3. * (P3[1] * (TMP75 - TMP76)))) + P3[1] *
      (TMP12 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[3] *
      V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] * V1[4]))) + TMP9
      * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] * V1[3]))))) + (P3[2] *
      (V1[5] * (V2[2] * (P1[0] - P2[0]) + V2[3] * (P1[1] - P2[1])) + V2[5] *
      (V1[2] * (P2[0] - P1[0]) + V1[3] * (P2[1] - P1[1]))) + (P3[3] * (V1[4] *
      (V2[2] * (P2[0] - P1[0]) + V2[3] * (P2[1] - P1[1])) + V2[4] * (V1[2] *
      (P1[0] - P2[0]) + V1[3] * (P1[1] - P2[1]))) + (P3[0] * (V1[4] * V2[5] *
      (P1[0] - P2[0]) + V1[5] * V2[4] * (P2[0] - P1[0])) + P3[1] * (V1[4] *
      V2[5] * (P1[1] - P2[1]) + V1[5] * V2[4] * (P2[1] - P1[1]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP12 * (P3[0] * (V2[3] * V1[5] - V2[5]
      * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP9 * (P3[0] * (V2[5] * V1[3] - V2[3] *
      V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 2./3. * (P3[2] * (TMP75 - TMP76)))) + P3[2] *
      (TMP12 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[3] *
      V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] * V1[4]))) + TMP9
      * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] * V1[3]))))) + (P3[1] *
      (V1[5] * (V2[2] * (P2[0] - P1[0]) + V2[4] * (P2[2] - P1[2])) + V2[5] *
      (V1[2] * (P1[0] - P2[0]) + V1[4] * (P1[2] - P2[2]))) + (P3[3] * (V1[3] *
      (V2[2] * (P1[0] - P2[0]) + V2[4] * (P1[2] - P2[2])) + V2[3] * (V1[2] *
      (P2[0] - P1[0]) + V1[4] * (P2[2] - P1[2]))) + (P3[0] * (V1[3] * V2[5] *
      (P2[0] - P1[0]) + V1[5] * V2[3] * (P1[0] - P2[0])) + P3[2] * (V1[3] *
      V2[5] * (P2[2] - P1[2]) + V1[5] * V2[3] * (P1[2] - P2[2]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP12 * (P3[0] * (V2[4] * V1[3] - V2[3]
      * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[4] - V2[4] *
      V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + 2./3. * (P3[3] * (TMP75 - TMP76)))) + P3[3] *
      (TMP12 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[3] *
      V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] * V1[4]))) + TMP9
      * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] * V1[3]))))) + (P3[1] *
      (V1[4] * (V2[2] * (P1[0] - P2[0]) + V2[5] * (P1[3] - P2[3])) + V2[4] *
      (V1[2] * (P2[0] - P1[0]) + V1[5] * (P2[3] - P1[3]))) + (P3[2] * (V1[3] *
      (V2[2] * (P2[0] - P1[0]) + V2[5] * (P2[3] - P1[3])) + V2[3] * (V1[2] *
      (P1[0] - P2[0]) + V1[5] * (P1[3] - P2[3]))) + (P3[0] * (V1[3] * V2[4] *
      (P1[0] - P2[0]) + V1[4] * V2[3] * (P2[0] - P1[0])) + P3[3] * (V1[3] *
      V2[4] * (P1[3] - P2[3]) + V1[4] * V2[3] * (P2[3] - P1[3]))))));
  T3[6] = denom * - cI * (OM3 * (P3[0] * (TMP12 * (P3[0] * (V2[4] * V1[5] -
      V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[4] - V2[4] * V1[2]))) + (TMP9 * (P3[0] * (V2[5] * V1[4] -
      V2[4] * V1[5]) + (P3[2] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]))) + 2./3. * (P3[1] * (TMP76 - TMP75)))) +
      P3[1] * (TMP12 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + TMP9 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))))) + (P3[2] * (V1[5] * (V2[2] * (P2[0] - P1[0]) + V2[3] * (P2[1]
      - P1[1])) + V2[5] * (V1[2] * (P1[0] - P2[0]) + V1[3] * (P1[1] - P2[1])))
      + (P3[3] * (V1[4] * (V2[2] * (P1[0] - P2[0]) + V2[3] * (P1[1] - P2[1])) +
      V2[4] * (V1[2] * (P2[0] - P1[0]) + V1[3] * (P2[1] - P1[1]))) + (P3[0] *
      (V1[4] * V2[5] * (P2[0] - P1[0]) + V1[5] * V2[4] * (P1[0] - P2[0])) +
      P3[1] * (V1[4] * V2[5] * (P2[1] - P1[1]) + V1[5] * V2[4] * (P1[1] -
      P2[1]))))));
  T3[7] = denom * 2. * cI * (OM3 * P3[1] * (TMP12 * (P3[0] * (V2[5] * V1[4] -
      V2[4] * V1[5]) + (P3[2] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[4] * V1[2] - V2[2] * V1[4]))) + (TMP9 * (P3[0] * (V2[4] * V1[5] -
      V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[4] - V2[4] * V1[2]))) + 1./3. * (P3[1] * (TMP75 - TMP76)))) +
      (P1[1] * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] * V1[4]))) +
      (P2[1] * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] *
      V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] * V1[2]))) +
      (-1./3. * (TMP76) + 1./3. * (TMP75)))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP12 * (P3[0] * (V2[3] * V1[5] - V2[5]
      * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + (TMP9 * (P3[0] * (V2[5] * V1[3] - V2[3] *
      V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + 2./3. * (P3[2] * (TMP75 - TMP76)))) + P3[2] *
      (TMP12 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] * V1[4]))) + TMP9
      * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[2] -
      V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] * V1[2]))))) + (P3[0] *
      (V1[5] * (V2[3] * (P1[1] - P2[1]) + V2[4] * (P2[2] - P1[2])) + V2[5] *
      (V1[3] * (P2[1] - P1[1]) + V1[4] * (P1[2] - P2[2]))) + (P3[3] * (V1[2] *
      (V2[3] * (P2[1] - P1[1]) + V2[4] * (P1[2] - P2[2])) + V2[2] * (V1[3] *
      (P1[1] - P2[1]) + V1[4] * (P2[2] - P1[2]))) + (P3[1] * (V1[2] * V2[5] *
      (P1[1] - P2[1]) + V1[5] * V2[2] * (P2[1] - P1[1])) + P3[2] * (V1[2] *
      V2[5] * (P2[2] - P1[2]) + V1[5] * V2[2] * (P1[2] - P2[2]))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP12 * (P3[0] * (V2[4] * V1[3] - V2[3]
      * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] *
      V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[4] - V2[4] *
      V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] *
      V1[3] - V2[3] * V1[2]))) + 2./3. * (P3[3] * (TMP75 - TMP76)))) + P3[3] *
      (TMP12 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] * V1[4]))) + TMP9
      * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] * (V2[5] * V1[2] -
      V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] * V1[2]))))) + (P3[0] *
      (V1[4] * (V2[3] * (P2[1] - P1[1]) + V2[5] * (P1[3] - P2[3])) + V2[4] *
      (V1[3] * (P1[1] - P2[1]) + V1[5] * (P2[3] - P1[3]))) + (P3[2] * (V1[2] *
      (V2[3] * (P1[1] - P2[1]) + V2[5] * (P2[3] - P1[3])) + V2[2] * (V1[3] *
      (P2[1] - P1[1]) + V1[5] * (P1[3] - P2[3]))) + (P3[1] * (V1[2] * V2[4] *
      (P2[1] - P1[1]) + V1[4] * V2[2] * (P1[1] - P2[1])) + P3[3] * (V1[2] *
      V2[4] * (P1[3] - P2[3]) + V1[4] * V2[2] * (P2[3] - P1[3]))))));
  T3[10] = denom * - cI * (OM3 * (P3[0] * (TMP12 * (P3[0] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[5] -
      V2[5] * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 2./3. * (P3[2] * (TMP76 - TMP75)))) +
      P3[2] * (TMP12 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + TMP9 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))))) + (P3[1] * (V1[5] * (V2[2] * (P1[0] - P2[0]) + V2[4] * (P1[2]
      - P2[2])) + V2[5] * (V1[2] * (P2[0] - P1[0]) + V1[4] * (P2[2] - P1[2])))
      + (P3[3] * (V1[3] * (V2[2] * (P2[0] - P1[0]) + V2[4] * (P2[2] - P1[2])) +
      V2[3] * (V1[2] * (P1[0] - P2[0]) + V1[4] * (P1[2] - P2[2]))) + (P3[0] *
      (V1[3] * V2[5] * (P1[0] - P2[0]) + V1[5] * V2[3] * (P2[0] - P1[0])) +
      P3[2] * (V1[3] * V2[5] * (P1[2] - P2[2]) + V1[5] * V2[3] * (P2[2] -
      P1[2]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP12 * (P3[0] * (V2[3] * V1[5] -
      V2[5] * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + (TMP9 * (P3[0] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + 2./3. * (P3[2] * (TMP75 - TMP76)))) +
      P3[2] * (TMP12 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + TMP9 * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))))) + (P3[0] * (V1[5] * (V2[3] * (P1[1] - P2[1]) + V2[4] * (P2[2]
      - P1[2])) + V2[5] * (V1[3] * (P2[1] - P1[1]) + V1[4] * (P1[2] - P2[2])))
      + (P3[3] * (V1[2] * (V2[3] * (P2[1] - P1[1]) + V2[4] * (P1[2] - P2[2])) +
      V2[2] * (V1[3] * (P1[1] - P2[1]) + V1[4] * (P2[2] - P1[2]))) + (P3[1] *
      (V1[2] * V2[5] * (P1[1] - P2[1]) + V1[5] * V2[2] * (P2[1] - P1[1])) +
      P3[2] * (V1[2] * V2[5] * (P2[2] - P1[2]) + V1[5] * V2[2] * (P1[2] -
      P2[2]))))));
  T3[12] = denom * 2. * cI * (OM3 * P3[2] * (TMP12 * (P3[0] * (V2[3] * V1[5] -
      V2[5] * V1[3]) + (P3[1] * (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + (TMP9 * (P3[0] * (V2[5] * V1[3] -
      V2[3] * V1[5]) + (P3[1] * (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + 1./3. * (P3[2] * (TMP75 - TMP76)))) +
      (P1[2] * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] * (V2[5] *
      V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] * V1[2]))) +
      (P2[2] * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] * (V2[2] *
      V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] * V1[3]))) +
      (-1./3. * (TMP76) + 1./3. * (TMP75)))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP12 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 2./3. * (P3[3] * (TMP75 - TMP76)))) +
      P3[3] * (TMP12 * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + TMP9 * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))) + (P3[0] * (V1[3] * (V2[4] * (P1[2] - P2[2]) + V2[5] * (P2[3]
      - P1[3])) + V2[3] * (V1[4] * (P2[2] - P1[2]) + V1[5] * (P1[3] - P2[3])))
      + (P3[1] * (V1[2] * (V2[4] * (P2[2] - P1[2]) + V2[5] * (P1[3] - P2[3])) +
      V2[2] * (V1[4] * (P1[2] - P2[2]) + V1[5] * (P2[3] - P1[3]))) + (P3[2] *
      (V1[2] * V2[3] * (P1[2] - P2[2]) + V1[3] * V2[2] * (P2[2] - P1[2])) +
      P3[3] * (V1[2] * V2[3] * (P2[3] - P1[3]) + V1[3] * V2[2] * (P1[3] -
      P2[3]))))));
  T3[14] = denom * - cI * (OM3 * (P3[0] * (TMP12 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + (TMP9 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + 2./3. * (P3[3] * (TMP76 - TMP75)))) +
      P3[3] * (TMP12 * (P3[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P3[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + TMP9 * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))))) + (P3[1] * (V1[4] * (V2[2] * (P2[0] - P1[0]) + V2[5] * (P2[3]
      - P1[3])) + V2[4] * (V1[2] * (P1[0] - P2[0]) + V1[5] * (P1[3] - P2[3])))
      + (P3[2] * (V1[3] * (V2[2] * (P1[0] - P2[0]) + V2[5] * (P1[3] - P2[3])) +
      V2[3] * (V1[2] * (P2[0] - P1[0]) + V1[5] * (P2[3] - P1[3]))) + (P3[0] *
      (V1[3] * V2[4] * (P2[0] - P1[0]) + V1[4] * V2[3] * (P1[0] - P2[0])) +
      P3[3] * (V1[3] * V2[4] * (P2[3] - P1[3]) + V1[4] * V2[3] * (P1[3] -
      P2[3]))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP12 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 2./3. * (P3[3] * (TMP75 - TMP76)))) +
      P3[3] * (TMP12 * (P3[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + TMP9 * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))))) + (P3[0] * (V1[4] * (V2[3] * (P2[1] - P1[1]) + V2[5] * (P1[3]
      - P2[3])) + V2[4] * (V1[3] * (P1[1] - P2[1]) + V1[5] * (P2[3] - P1[3])))
      + (P3[2] * (V1[2] * (V2[3] * (P1[1] - P2[1]) + V2[5] * (P2[3] - P1[3])) +
      V2[2] * (V1[3] * (P2[1] - P1[1]) + V1[5] * (P1[3] - P2[3]))) + (P3[1] *
      (V1[2] * V2[4] * (P2[1] - P1[1]) + V1[4] * V2[2] * (P1[1] - P2[1])) +
      P3[3] * (V1[2] * V2[4] * (P1[3] - P2[3]) + V1[4] * V2[2] * (P2[3] -
      P1[3]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP12 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 2./3. * (P3[3] * (TMP75 - TMP76)))) +
      P3[3] * (TMP12 * (P3[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P3[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + TMP9 * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))) + (P3[0] * (V1[3] * (V2[4] * (P1[2] - P2[2]) + V2[5] * (P2[3]
      - P1[3])) + V2[3] * (V1[4] * (P2[2] - P1[2]) + V1[5] * (P1[3] - P2[3])))
      + (P3[1] * (V1[2] * (V2[4] * (P2[2] - P1[2]) + V2[5] * (P1[3] - P2[3])) +
      V2[2] * (V1[4] * (P1[2] - P2[2]) + V1[5] * (P2[3] - P1[3]))) + (P3[2] *
      (V1[2] * V2[3] * (P1[2] - P2[2]) + V1[3] * V2[2] * (P2[2] - P1[2])) +
      P3[3] * (V1[2] * V2[3] * (P2[3] - P1[3]) + V1[3] * V2[2] * (P1[3] -
      P2[3]))))));
  T3[17] = denom * 2. * cI * (OM3 * P3[3] * (TMP12 * (P3[0] * (V2[4] * V1[3] -
      V2[3] * V1[4]) + (P3[1] * (V2[2] * V1[4] - V2[4] * V1[2]) + P3[2] *
      (V2[3] * V1[2] - V2[2] * V1[3]))) + (TMP9 * (P3[0] * (V2[3] * V1[4] -
      V2[4] * V1[3]) + (P3[1] * (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] *
      (V2[2] * V1[3] - V2[3] * V1[2]))) + 1./3. * (P3[3] * (TMP75 - TMP76)))) +
      (P1[3] * (P3[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P3[1] * (V2[2] *
      V1[4] - V2[4] * V1[2]) + P3[2] * (V2[3] * V1[2] - V2[2] * V1[3]))) +
      (P2[3] * (P3[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] * (V2[4] *
      V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] * V1[3] - V2[3] * V1[2]))) +
      (-1./3. * (TMP76) + 1./3. * (TMP75)))));
}


void FFV2_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP50; 
  TMP50 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * - cI * TMP50; 
}

void FFV5_6_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV5_1(F2, V3, COUP1, M1, W1, F1); 
  FFV6_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV5_8_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV5_1(F2, V3, COUP1, M1, W1, F1); 
  FFV8_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}


void VVT9_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  double OM3; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP9; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP30 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 2./3. * P3[0] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[0] + P1[0]))) + (+1./3. * cI * (TMP12
      + TMP9))) + P3[0] * (-cI * (TMP9 * V2[2]) + 1./3. * cI * (P3[0] *
      TMP26))) + P3[0] * TMP38 * (-cI * (TMP12 * V1[2]) + 1./3. * cI * (P3[0] *
      TMP30))) + (TMP37 * (-1./3. * cI * (TMP26) + cI * (P1[0] * V2[2])) +
      TMP38 * (-1./3. * cI * (TMP30) + cI * (P2[0] * V1[2]))));
  T3[6] = denom * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 4./3. * P3[1] * (+cI
      * (TMP12 + TMP9)) + (-cI * (P2[1] + P1[1]))) - P3[1] * (+cI * (P2[0] +
      P1[0]))) + (P3[0] * (-cI * (TMP9 * V2[3]) + 2./3. * cI * (P3[1] * TMP26))
      - cI * (P3[1] * TMP9 * V2[2]))) + TMP38 * (P3[0] * (-cI * (TMP12 * V1[3])
      + 2./3. * cI * (P3[1] * TMP30)) - cI * (P3[1] * TMP12 * V1[2]))) + (TMP37
      * (+cI * (P1[0] * V2[3] + P1[1] * V2[2])) + TMP38 * (+cI * (P2[0] * V1[3]
      + P2[1] * V1[2]))));
  T3[10] = denom * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 4./3. * P3[2] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI * (P2[0]
      + P1[0]))) + (P3[0] * (-cI * (TMP9 * V2[4]) + 2./3. * cI * (P3[2] *
      TMP26)) - cI * (P3[2] * TMP9 * V2[2]))) + TMP38 * (P3[0] * (-cI * (TMP12
      * V1[4]) + 2./3. * cI * (P3[2] * TMP30)) - cI * (P3[2] * TMP12 * V1[2])))
      + (TMP37 * (+cI * (P1[0] * V2[4] + P1[2] * V2[2])) + TMP38 * (+cI *
      (P2[0] * V1[4] + P2[2] * V1[2]))));
  T3[14] = denom * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 4./3. * P3[3] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI * (P2[0]
      + P1[0]))) + (P3[0] * (-cI * (TMP9 * V2[5]) + 2./3. * cI * (P3[3] *
      TMP26)) - cI * (P3[3] * TMP9 * V2[2]))) + TMP38 * (P3[0] * (-cI * (TMP12
      * V1[5]) + 2./3. * cI * (P3[3] * TMP30)) - cI * (P3[3] * TMP12 * V1[2])))
      + (TMP37 * (+cI * (P1[0] * V2[5] + P1[3] * V2[2])) + TMP38 * (+cI *
      (P2[0] * V1[5] + P2[3] * V1[2]))));
  T3[3] = denom * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 4./3. * P3[1] * (+cI
      * (TMP12 + TMP9)) + (-cI * (P2[1] + P1[1]))) - P3[1] * (+cI * (P2[0] +
      P1[0]))) + (P3[0] * (-cI * (TMP9 * V2[3]) + 2./3. * cI * (P3[1] * TMP26))
      - cI * (P3[1] * TMP9 * V2[2]))) + TMP38 * (P3[0] * (-cI * (TMP12 * V1[3])
      + 2./3. * cI * (P3[1] * TMP30)) - cI * (P3[1] * TMP12 * V1[2]))) + (TMP37
      * (+cI * (P1[1] * V2[2] + P1[0] * V2[3])) + TMP38 * (+cI * (P2[1] * V1[2]
      + P2[0] * V1[3]))));
  T3[7] = denom * 2. * (OM3 * (TMP37 * (TMP38 * (P3[1] * (OM3 * 2./3. * P3[1] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[1] + P1[1]))) + (-1./3. * cI * (TMP12
      + TMP9))) + P3[1] * (-cI * (TMP9 * V2[3]) + 1./3. * cI * (P3[1] *
      TMP26))) + P3[1] * TMP38 * (-cI * (TMP12 * V1[3]) + 1./3. * cI * (P3[1] *
      TMP30))) + (TMP37 * (+cI * (P1[1] * V2[3]) + 1./3. * cI * (TMP26)) +
      TMP38 * (+cI * (P2[1] * V1[3]) + 1./3. * cI * (TMP30))));
  T3[11] = denom * (OM3 * (TMP37 * (TMP38 * (P3[1] * (OM3 * 4./3. * P3[2] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI * (P2[1]
      + P1[1]))) + (P3[1] * (-cI * (TMP9 * V2[4]) + 2./3. * cI * (P3[2] *
      TMP26)) - cI * (P3[2] * TMP9 * V2[3]))) + TMP38 * (P3[1] * (-cI * (TMP12
      * V1[4]) + 2./3. * cI * (P3[2] * TMP30)) - cI * (P3[2] * TMP12 * V1[3])))
      + (TMP37 * (+cI * (P1[1] * V2[4] + P1[2] * V2[3])) + TMP38 * (+cI *
      (P2[1] * V1[4] + P2[2] * V1[3]))));
  T3[15] = denom * (OM3 * (TMP37 * (TMP38 * (P3[1] * (OM3 * 4./3. * P3[3] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI * (P2[1]
      + P1[1]))) + (P3[1] * (-cI * (TMP9 * V2[5]) + 2./3. * cI * (P3[3] *
      TMP26)) - cI * (P3[3] * TMP9 * V2[3]))) + TMP38 * (P3[1] * (-cI * (TMP12
      * V1[5]) + 2./3. * cI * (P3[3] * TMP30)) - cI * (P3[3] * TMP12 * V1[3])))
      + (TMP37 * (+cI * (P1[1] * V2[5] + P1[3] * V2[3])) + TMP38 * (+cI *
      (P2[1] * V1[5] + P2[3] * V1[3]))));
  T3[4] = denom * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 4./3. * P3[2] * (+cI
      * (TMP12 + TMP9)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI * (P2[0] +
      P1[0]))) + (P3[0] * (-cI * (TMP9 * V2[4]) + 2./3. * cI * (P3[2] * TMP26))
      - cI * (P3[2] * TMP9 * V2[2]))) + TMP38 * (P3[0] * (-cI * (TMP12 * V1[4])
      + 2./3. * cI * (P3[2] * TMP30)) - cI * (P3[2] * TMP12 * V1[2]))) + (TMP37
      * (+cI * (P1[2] * V2[2] + P1[0] * V2[4])) + TMP38 * (+cI * (P2[2] * V1[2]
      + P2[0] * V1[4]))));
  T3[8] = denom * (OM3 * (TMP37 * (TMP38 * (P3[1] * (OM3 * 4./3. * P3[2] * (+cI
      * (TMP12 + TMP9)) + (-cI * (P2[2] + P1[2]))) - P3[2] * (+cI * (P2[1] +
      P1[1]))) + (P3[1] * (-cI * (TMP9 * V2[4]) + 2./3. * cI * (P3[2] * TMP26))
      - cI * (P3[2] * TMP9 * V2[3]))) + TMP38 * (P3[1] * (-cI * (TMP12 * V1[4])
      + 2./3. * cI * (P3[2] * TMP30)) - cI * (P3[2] * TMP12 * V1[3]))) + (TMP37
      * (+cI * (P1[2] * V2[3] + P1[1] * V2[4])) + TMP38 * (+cI * (P2[2] * V1[3]
      + P2[1] * V1[4]))));
  T3[12] = denom * 2. * (OM3 * (TMP37 * (TMP38 * (P3[2] * (OM3 * 2./3. * P3[2]
      * (+cI * (TMP12 + TMP9)) + (-cI * (P2[2] + P1[2]))) + (-1./3. * cI *
      (TMP12 + TMP9))) + P3[2] * (-cI * (TMP9 * V2[4]) + 1./3. * cI * (P3[2] *
      TMP26))) + P3[2] * TMP38 * (-cI * (TMP12 * V1[4]) + 1./3. * cI * (P3[2] *
      TMP30))) + (TMP37 * (+cI * (P1[2] * V2[4]) + 1./3. * cI * (TMP26)) +
      TMP38 * (+cI * (P2[2] * V1[4]) + 1./3. * cI * (TMP30))));
  T3[16] = denom * (OM3 * (TMP37 * (TMP38 * (P3[2] * (OM3 * 4./3. * P3[3] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI * (P2[2]
      + P1[2]))) + (P3[2] * (-cI * (TMP9 * V2[5]) + 2./3. * cI * (P3[3] *
      TMP26)) - cI * (P3[3] * TMP9 * V2[4]))) + TMP38 * (P3[2] * (-cI * (TMP12
      * V1[5]) + 2./3. * cI * (P3[3] * TMP30)) - cI * (P3[3] * TMP12 * V1[4])))
      + (TMP37 * (+cI * (P1[2] * V2[5] + P1[3] * V2[4])) + TMP38 * (+cI *
      (P2[2] * V1[5] + P2[3] * V1[4]))));
  T3[5] = denom * (OM3 * (TMP37 * (TMP38 * (P3[0] * (OM3 * 4./3. * P3[3] * (+cI
      * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI * (P2[0] +
      P1[0]))) + (P3[0] * (-cI * (TMP9 * V2[5]) + 2./3. * cI * (P3[3] * TMP26))
      - cI * (P3[3] * TMP9 * V2[2]))) + TMP38 * (P3[0] * (-cI * (TMP12 * V1[5])
      + 2./3. * cI * (P3[3] * TMP30)) - cI * (P3[3] * TMP12 * V1[2]))) + (TMP37
      * (+cI * (P1[3] * V2[2] + P1[0] * V2[5])) + TMP38 * (+cI * (P2[3] * V1[2]
      + P2[0] * V1[5]))));
  T3[9] = denom * (OM3 * (TMP37 * (TMP38 * (P3[1] * (OM3 * 4./3. * P3[3] * (+cI
      * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI * (P2[1] +
      P1[1]))) + (P3[1] * (-cI * (TMP9 * V2[5]) + 2./3. * cI * (P3[3] * TMP26))
      - cI * (P3[3] * TMP9 * V2[3]))) + TMP38 * (P3[1] * (-cI * (TMP12 * V1[5])
      + 2./3. * cI * (P3[3] * TMP30)) - cI * (P3[3] * TMP12 * V1[3]))) + (TMP37
      * (+cI * (P1[3] * V2[3] + P1[1] * V2[5])) + TMP38 * (+cI * (P2[3] * V1[3]
      + P2[1] * V1[5]))));
  T3[13] = denom * (OM3 * (TMP37 * (TMP38 * (P3[2] * (OM3 * 4./3. * P3[3] *
      (+cI * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) - P3[3] * (+cI * (P2[2]
      + P1[2]))) + (P3[2] * (-cI * (TMP9 * V2[5]) + 2./3. * cI * (P3[3] *
      TMP26)) - cI * (P3[3] * TMP9 * V2[4]))) + TMP38 * (P3[2] * (-cI * (TMP12
      * V1[5]) + 2./3. * cI * (P3[3] * TMP30)) - cI * (P3[3] * TMP12 * V1[4])))
      + (TMP37 * (+cI * (P1[3] * V2[4] + P1[2] * V2[5])) + TMP38 * (+cI *
      (P2[3] * V1[4] + P2[2] * V1[5]))));
  T3[17] = denom * 2. * (OM3 * (TMP37 * (TMP38 * (P3[3] * (OM3 * 2./3. * P3[3]
      * (+cI * (TMP12 + TMP9)) + (-cI * (P2[3] + P1[3]))) + (-1./3. * cI *
      (TMP12 + TMP9))) + P3[3] * (-cI * (TMP9 * V2[5]) + 1./3. * cI * (P3[3] *
      TMP26))) + P3[3] * TMP38 * (-cI * (TMP12 * V1[5]) + 1./3. * cI * (P3[3] *
      TMP30))) + (TMP37 * (+cI * (P1[3] * V2[5]) + 1./3. * cI * (TMP26)) +
      TMP38 * (+cI * (P2[3] * V1[5]) + 1./3. * cI * (TMP30))));
}


void FFT1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP17; 
  double P3[4]; 
  complex<double> TMP16; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP9; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP17 = (F1[4] * F2[4] + F1[5] * F2[5] - F1[2] * F2[2] - F1[3] * F2[3]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * TMP17 * (OM3 * (P3[0] * (P3[0] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[0] * TMP9 + P1[0] *
      TMP12))) + 1./3. * cI * (TMP9 * TMP12)) + (-1./3. * cI * (TMP16) + cI *
      (P1[0] * P2[0])));
  T3[6] = denom * TMP17 * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[1] * TMP9 + P1[1] * TMP12))) -
      P3[1] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[1] +
      P1[1] * P2[0])));
  T3[10] = denom * TMP17 * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] * TMP12))) -
      P3[2] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[2] +
      P1[2] * P2[0])));
  T3[14] = denom * TMP17 * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[0] * TMP12 + P2[0] * TMP9))) + (+cI * (P1[0] * P2[3] +
      P1[3] * P2[0])));
  T3[3] = denom * TMP17 * (OM3 * (P3[0] * (P3[1] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[1] * TMP12 + P2[1] * TMP9))) -
      P3[1] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[1] * P2[0] +
      P1[0] * P2[1])));
  T3[7] = denom * 2. * TMP17 * (OM3 * (P3[1] * (P3[1] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[1] * TMP9 + P1[1] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[1] * P2[1]) + 1./3.
      * cI * (TMP16)));
  T3[11] = denom * TMP17 * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] * TMP12))) -
      P3[2] * (+cI * (P1[1] * TMP12 + P2[1] * TMP9))) + (+cI * (P1[1] * P2[2] +
      P1[2] * P2[1])));
  T3[15] = denom * TMP17 * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[1] * TMP12 + P2[1] * TMP9))) + (+cI * (P1[1] * P2[3] +
      P1[3] * P2[1])));
  T3[4] = denom * TMP17 * (OM3 * (P3[0] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[2] * TMP12 + P2[2] * TMP9))) -
      P3[2] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[2] * P2[0] +
      P1[0] * P2[2])));
  T3[8] = denom * TMP17 * (OM3 * (P3[1] * (P3[2] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[2] * TMP12 + P2[2] * TMP9))) -
      P3[2] * (+cI * (P2[1] * TMP9 + P1[1] * TMP12))) + (+cI * (P1[2] * P2[1] +
      P1[1] * P2[2])));
  T3[12] = denom * 2. * TMP17 * (OM3 * (P3[2] * (P3[2] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[2] * TMP9 + P1[2] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[2] * P2[2]) + 1./3.
      * cI * (TMP16)));
  T3[16] = denom * TMP17 * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] * TMP12))) -
      P3[3] * (+cI * (P1[2] * TMP12 + P2[2] * TMP9))) + (+cI * (P1[2] * P2[3] +
      P1[3] * P2[2])));
  T3[5] = denom * TMP17 * (OM3 * (P3[0] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[0] * TMP9 + P1[0] * TMP12))) + (+cI * (P1[3] * P2[0] +
      P1[0] * P2[3])));
  T3[9] = denom * TMP17 * (OM3 * (P3[1] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[1] * TMP9 + P1[1] * TMP12))) + (+cI * (P1[3] * P2[1] +
      P1[1] * P2[3])));
  T3[13] = denom * TMP17 * (OM3 * (P3[2] * (P3[3] * 2./3. * (+cI * (TMP16) + 2.
      * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P1[3] * TMP12 + P2[3] * TMP9))) -
      P3[3] * (+cI * (P2[2] * TMP9 + P1[2] * TMP12))) + (+cI * (P1[3] * P2[2] +
      P1[2] * P2[3])));
  T3[17] = denom * 2. * TMP17 * (OM3 * (P3[3] * (P3[3] * 1./3. * (+cI * (TMP16)
      + 2. * cI * (OM3 * TMP9 * TMP12)) + (-cI * (P2[3] * TMP9 + P1[3] *
      TMP12))) - 1./3. * cI * (TMP9 * TMP12)) + (+cI * (P1[3] * P2[3]) + 1./3.
      * cI * (TMP16)));
}

void FFT1_2_4_5_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M3, double W3, complex<double> T3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ttmp[18]; 
//   double P1[4]; 
//   double P2[4]; 
//   double P3[4]; 
//   double OM3; 
  complex<double> denom; 
  int i; 
  FFT1_3(F1, F2, COUP1, M3, W3, T3); 
  FFT2_3(F1, F2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  FFT4_3(F1, F2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  FFT5_3(F1, F2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}
void FFT1_2_3_5_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M3, double W3, complex<double> T3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ttmp[18]; 
//   double P1[4]; 
//   double P2[4]; 
//   double P3[4]; 
//   double OM3; 
  complex<double> denom; 
  int i; 
  FFT1_3(F1, F2, COUP1, M3, W3, T3); 
  FFT2_3(F1, F2, COUP2, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  FFT3_3(F1, F2, COUP3, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
  FFT5_3(F1, F2, COUP4, M3, W3, Ttmp); 
  i = 2; 
  while (i < 18)
  {
    T3[i] = T3[i] + Ttmp[i]; 
    i++; 
  }
}

void FFT4_1(complex<double> F2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + T3[0]; 
  F1[1] = +F2[1] + T3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[3] * (P1[0] * (P1[3] * (T3[9] + T3[15] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[6] - T3[3]) + (P1[1] *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8])) + (P1[2]
      * (T3[8] + T3[11] + cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) +
      (P1[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[0] * (T3[6]
      + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] +
      T3[8])) - P2[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P1[1] *
      (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13] + cI * (T3[3] + T3[15]
      + T3[6] + T3[9])) + (P1[3] * (+2. * (T3[7]) + cI * (T3[11] + T3[8]) - 2.
      * (T3[17]) - T3[5] - T3[14]) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] +
      T3[16] + T3[10] + T3[13]) + (P2[0] * - 1. * (T3[14] + T3[5] + 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[17]) +
      2. * cI * (T3[12])) + (P1[2] * - 1. * (+cI * (T3[4] + T3[16] + T3[10] +
      T3[13])) + (P2[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P2[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[0] * - 1. * (+cI *
      (T3[14] + T3[5]) + 2. * cI * (T3[2])) + P2[3] * (+cI * (T3[5] + T3[14]) +
      2. * cI * (T3[17]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P2[1] * - 1.
      * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8] + T3[11] + 2.
      * cI * (T3[12]))))))))) + (F2[2] * (P1[0] * (P1[1] * (T3[15] + T3[9] + cI
      * (T3[10] + T3[4]) - 2. * (T3[3] + T3[6])) + (P1[2] * - 1. * (+2. *
      (T3[4] + T3[10]) + cI * (T3[6] + T3[3]) - T3[16] - T3[13]) + (P1[3] * 2.
      * (T3[17] + T3[2] - T3[5] - T3[14]) + (P2[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[0] * - 1. *
      (T3[14] + T3[5] - 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] - 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] - 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * - 1. * (T3[3] + T3[6] + cI * (T3[13] + T3[16]) - 2. * (T3[9] +
      T3[15])) + (P1[2] * (+2. * (T3[13] + T3[16]) + cI * (T3[9] + T3[15]) -
      T3[4] - T3[10]) + (P2[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] - 2.
      * (T3[17])) + (P2[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P2[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. * (T3[8] + T3[11] -
      cI * (T3[12]) + cI * (T3[7])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] *
      - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + (P2[1] * (+cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] + cI * (T3[6] + T3[3]))
      + (P2[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] + T3[15])) + (P1[2] *
      (+2. * (T3[12]) + cI * (T3[8] + T3[11])) + (P2[1] * - 1. * (T3[11] +
      T3[8] + 2. * cI * (T3[7])) - P2[2] * (+2. * (T3[12]) + cI * (T3[8] +
      T3[11]))))))))) + M1 * (F2[4] * (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] *
      (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (T3[14] * (P1[0] + P2[3] - P2[0] - P1[3]) + (T3[5] * (P2[3] +
      P1[0] - P1[3] - P2[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] *
      (P2[3] - P1[3]))))))))) + F2[5] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P1[1]) + cI
      * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) + cI * (P2[1]) -
      P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P2[1] - P1[1]))))))))))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (P1[3] * (+cI * (T3[13] + T3[10] +
      T3[16] + T3[4]) - T3[9] - T3[6] - T3[15] - T3[3]) + (P1[1] * (T3[14] +
      T3[5] + cI * (T3[11] + T3[8]) - 2. * (T3[7] + T3[2])) + (P1[2] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12] + T3[2]) + cI * (T3[14] + T3[5])) +
      (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[0] * (+cI *
      (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI *
      (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P1[1]
      * (P1[2] * (T3[4] + T3[10] - cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])
      - T3[16] - T3[13]) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[7]) - cI *
      (T3[11] + T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P2[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[12]) +
      2. * cI * (T3[17])) + (P1[2] * (-cI * (T3[4] + T3[10]) + cI * (T3[16] +
      T3[13])) + (P2[1] * (-cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6])) +
      (P2[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] + T3[10])) + (P2[0] *
      (-2. * cI * (T3[2]) + cI * (T3[14] + T3[5])) + P2[3] * (-2. * cI *
      (T3[17]) + cI * (T3[5] + T3[14]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15]
      - cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] +
      T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[1] *
      (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] + T3[11] - 2. *
      cI * (T3[12]))))))))) + (F2[3] * (P1[0] * (P1[1] * (T3[15] + T3[9] + 2. *
      (T3[3] + T3[6]) + cI * (T3[10] + T3[4])) + (P1[2] * (T3[16] + T3[13] + 2.
      * (T3[4] + T3[10]) - cI * (T3[6] + T3[3])) + (P1[3] * 2. * (T3[5] +
      T3[17] + T3[2] + T3[14]) + (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * -
      1. * (T3[14] + T3[5] + 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] + 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * - 1. * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) + cI * (T3[13] +
      T3[16])) + (P1[2] * - 1. * (T3[4] + T3[10] - cI * (T3[9] + T3[15]) + 2. *
      (T3[13] + T3[16])) + (P2[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] *
      (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] + 2.
      * (T3[17])) + (P2[0] * - 1. * (T3[14] + T3[5] + 2. * (T3[2])) + P2[3] *
      (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[1] * (P1[2] * - 2. * (T3[8]
      + T3[11] - cI * (T3[7]) + cI * (T3[12])) + (P2[0] * - 1. * (T3[6] + T3[3]
      + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) +
      (P2[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + P2[2] * (T3[8] +
      T3[11] + 2. * cI * (T3[12]))))))) + P1[2] * (P2[0] * (+cI * (T3[6] +
      T3[3]) - T3[10] - T3[4]) + (P2[3] * (T3[13] + T3[16] - cI * (T3[9] +
      T3[15])) + (P1[2] * (+cI * (T3[8] + T3[11]) - 2. * (T3[12])) + (P2[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) - P2[2] * (+cI * (T3[8] + T3[11]) -
      2. * (T3[12]))))))))) + M1 * (F2[4] * (P1[0] * (+cI * (T3[10] + T3[4]) -
      T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) +
      (P2[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI * (P1[1]) +
      cI * (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) -
      P2[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P1[1] - P2[1]))))))))) + F2[5] * (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (T3[14] * (P1[0] + P1[3] - P2[0] - P2[3]) + (T3[5] * (P1[3] + P1[0] -
      P2[3] - P2[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))))));
  F1[4] = denom * - cI * (F2[5] * (P1[0] * (P1[3] * - 1. * (T3[9] + T3[6] +
      T3[15] + T3[3] + cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P1[1] * - 1.
      * (+2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8]) - T3[14] - T3[5]) +
      (P1[2] * - 1. * (T3[8] + T3[11] - cI * (T3[14] + T3[5]) + 2. * cI *
      (T3[12] + T3[2])) + (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8])) + P2[2] * (T3[8] + T3[11] + 2. * cI * (T3[12])))))))))
      + (P1[1] * (P1[2] * (T3[4] + T3[10] - cI * (T3[15] + T3[9]) + cI * (T3[3]
      + T3[6]) - T3[16] - T3[13]) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[7]) +
      cI * (T3[11] + T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15]
      - T3[9]) + (P2[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P2[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) - P2[3] * (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P1[2] *
      (P1[3] * (T3[8] + T3[11] - 2. * cI * (T3[17]) + cI * (T3[5] + T3[14]) +
      2. * cI * (T3[12])) + (P1[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] +
      T3[10])) + (P2[1] * (-cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])) +
      (P2[2] * (-cI * (T3[4] + T3[10]) + cI * (T3[16] + T3[13])) + (P2[0] * -
      1. * (-2. * cI * (T3[2]) + cI * (T3[14] + T3[5])) - P2[3] * (-2. * cI *
      (T3[17]) + cI * (T3[5] + T3[14]))))))) + P1[3] * (P1[3] * (T3[9] + T3[15]
      + cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] +
      T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8]
      + T3[11] + 2. * cI * (T3[12]))))))))) + (F2[4] * (P1[0] * (P1[1] * - 1. *
      (T3[15] + T3[9] + 2. * (T3[3] + T3[6]) - cI * (T3[10] + T3[4])) + (P1[2]
      * - 1. * (T3[16] + T3[13] + 2. * (T3[4] + T3[10]) + cI * (T3[6] + T3[3]))
      + (P1[3] * - 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (P1[0] * (T3[14] + T3[5] + 2. * (T3[2])) + (P2[0] * - 1. * (T3[14] +
      T3[5] + 2. * (T3[2])) + P2[3] * (T3[5] + T3[14] + 2. * (T3[17]))))))))) +
      (P1[3] * (P1[1] * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) - cI * (T3[13] +
      T3[16])) + (P1[2] * (T3[4] + T3[10] + 2. * (T3[13] + T3[16]) + cI *
      (T3[9] + T3[15])) + (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[3] * (T3[5] +
      T3[14] + 2. * (T3[17])) + (P2[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. *
      (T3[8] + T3[11] - cI * (T3[12]) + cI * (T3[7])) + (P2[0] * (T3[6] + T3[3]
      - cI * (T3[10] + T3[4])) + (P2[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      (P2[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P2[2] * (T3[8] +
      T3[11] - 2. * cI * (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] + cI
      * (T3[6] + T3[3])) + (P2[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] +
      T3[15])) + (P1[2] * (+2. * (T3[12]) + cI * (T3[8] + T3[11])) + (P2[1] * -
      1. * (T3[11] + T3[8] + 2. * cI * (T3[7])) - P2[2] * (+2. * (T3[12]) + cI
      * (T3[8] + T3[11]))))))))) + M1 * (F2[2] * (P1[1] * (T3[15] + T3[9] -
      T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[1] *
      (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] -
      T3[13]) + (T3[14] * (P2[0] + P2[3] - P1[0] - P1[3]) + (T3[5] * (P2[3] +
      P2[0] - P1[3] - P1[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] *
      (P1[3] - P2[3]))))))))) + F2[3] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P2[1]) + cI * (P1[1]) -
      P2[2]) + (T3[12] * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] *
      (P1[1] - P2[1]))))))))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (P1[3] * (T3[6] + T3[3] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * -
      1. * (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) - cI * (T3[11] + T3[8])) +
      (P1[2] * (+cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2]) - T3[8] -
      T3[11]) + (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[0] *
      (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI
      * (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P1[1]
      * (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13] - cI * (T3[3] + T3[15] +
      T3[6] + T3[9])) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[17]) + cI * (T3[11]
      + T3[8]) - 2. * (T3[7])) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) +
      (P2[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4]
      + T3[16] + T3[10] + T3[13]) + (P2[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P2[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P1[2] * (P1[3] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17])) + (P1[2] * - 1. * (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) +
      (P2[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P2[2] * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P2[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P1[3] * (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P2[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] *
      (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[1] * - 1. * (+cI *
      (T3[11] + T3[8]) - 2. * (T3[7])) + P2[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F2[5] * (P1[0] * (P1[1] * - 1. * (+2. * (T3[3] +
      T3[6]) + cI * (T3[10] + T3[4]) - T3[15] - T3[9]) + (P1[2] * (T3[16] +
      T3[13] + cI * (T3[6] + T3[3]) - 2. * (T3[4] + T3[10])) + (P1[3] * 2. *
      (T3[17] + T3[2] - T3[5] - T3[14]) + (P2[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[0] * - 1. *
      (T3[14] + T3[5] - 2. * (T3[2])) + (P2[0] * (T3[14] + T3[5] - 2. *
      (T3[2])) + P2[3] * (T3[5] + T3[14] - 2. * (T3[17]))))))))) + (P1[3] *
      (P1[1] * (+2. * (T3[9] + T3[15]) + cI * (T3[13] + T3[16]) - T3[3] -
      T3[6]) + (P1[2] * - 1. * (T3[4] + T3[10] + cI * (T3[9] + T3[15]) - 2. *
      (T3[13] + T3[16])) + (P2[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P2[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P1[3] * - 1. * (T3[5] + T3[14] - 2.
      * (T3[17])) + (P2[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P2[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P1[1] * (P1[2] * 2. * (T3[8] + T3[11] -
      cI * (T3[7]) + cI * (T3[12])) + (P2[0] * (T3[6] + T3[3] + cI * (T3[10] +
      T3[4])) + (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P2[1] * - 1. * (+2. *
      (T3[7]) + cI * (T3[11] + T3[8])) - P2[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))) + P1[2] * (P2[0] * (T3[10] + T3[4] - cI * (T3[6] + T3[3]))
      + (P2[3] * (+cI * (T3[9] + T3[15]) - T3[13] - T3[16]) + (P1[2] * - 1. *
      (+cI * (T3[8] + T3[11]) - 2. * (T3[12])) + (P2[1] * - 1. * (T3[11] +
      T3[8] - 2. * cI * (T3[7])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M1 * (F2[2] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6]
      - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[0] *
      (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] * (+cI * (T3[13] +
      T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) - P2[2])
      + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] * (P1[1] -
      P2[1]))))))))) + F2[3] * (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P2[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (T3[14] * (P1[0] + P2[3] - P2[0] - P1[3]) + (T3[5] * (P2[3] + P1[0] -
      P1[3] - P2[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))))));
}

void FFV6_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP3 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (OM3 * 1./2. * P3[0] * (TMP2 - 2. * (TMP3)) +
      (-1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
      F2[3]));
  V3[3] = denom * 2. * cI * (OM3 * 1./2. * P3[1] * (TMP2 - 2. * (TMP3)) +
      (+1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
      F2[2]));
  V3[4] = denom * - 2. * cI * (OM3 * 1./2. * P3[2] * (+2. * (TMP3) - TMP2) +
      (-1./2. * cI * (F1[2] * F2[5]) + 1./2. * cI * (F1[3] * F2[4]) - cI *
      (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * - 2. * cI * (OM3 * 1./2. * P3[3] * (+2. * (TMP3) - TMP2) +
      (-1./2. * (F1[2] * F2[4]) + 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] +
      F1[5] * F2[3]));
}


void VVT13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  double OM3; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP30; 
  complex<double> TMP16; 
  complex<double> denom; 
  complex<double> TMP26; 
  complex<double> TMP25; 
  complex<double> TMP9; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP25 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP30 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * (TMP12 * (TMP37 * (TMP26 * - 1./3. * (+cI * (P3[0] *
      P3[0] * OM3) + - 1./1. * cI) - cI * (P1[0] * V2[2])) + TMP9 * (TMP25 *
      1./3. * (+cI * (P3[0] * P3[0] * OM3) + - 1./1. * cI) + cI * (V2[2] *
      V1[2]))) + TMP38 * (TMP37 * (TMP16 * 1./3. * (+cI * (P3[0] * P3[0] * OM3)
      + - 1./1. * cI) + cI * (P1[0] * P2[0])) + TMP9 * (TMP30 * - 1./3. * (+cI
      * (P3[0] * P3[0] * OM3) + - 1./1. * cI) - cI * (P2[0] * V1[2]))));
  T3[6] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[0] * V2[3] + P1[1] *
      V2[2]) + 2./3. * cI * (P3[0] * P3[1] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[3] * V1[2] + V2[2] * V1[3]) + 2./3. * cI * (P3[0] * P3[1] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[0] * P2[1] + P1[1] * P2[0]) +
      2./3. * cI * (P3[0] * P3[1] * OM3 * TMP16)) - TMP9 * (+cI * (P2[0] *
      V1[3] + P2[1] * V1[2]) + 2./3. * cI * (P3[0] * P3[1] * OM3 * TMP30))));
  T3[10] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[0] * V2[4] + P1[2] *
      V2[2]) + 2./3. * cI * (P3[0] * P3[2] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[4] * V1[2] + V2[2] * V1[4]) + 2./3. * cI * (P3[0] * P3[2] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[0] * P2[2] + P1[2] * P2[0]) +
      2./3. * cI * (P3[0] * P3[2] * OM3 * TMP16)) - TMP9 * (+cI * (P2[0] *
      V1[4] + P2[2] * V1[2]) + 2./3. * cI * (P3[0] * P3[2] * OM3 * TMP30))));
  T3[14] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[0] * V2[5] + P1[3] *
      V2[2]) + 2./3. * cI * (P3[0] * P3[3] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[5] * V1[2] + V2[2] * V1[5]) + 2./3. * cI * (P3[0] * P3[3] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[0] * P2[3] + P1[3] * P2[0]) +
      2./3. * cI * (P3[0] * P3[3] * OM3 * TMP16)) - TMP9 * (+cI * (P2[0] *
      V1[5] + P2[3] * V1[2]) + 2./3. * cI * (P3[0] * P3[3] * OM3 * TMP30))));
  T3[3] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[1] * V2[2] + P1[0] *
      V2[3]) + 2./3. * cI * (P3[0] * P3[1] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[2] * V1[3] + V2[3] * V1[2]) + 2./3. * cI * (P3[0] * P3[1] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[1] * P2[0] + P1[0] * P2[1]) +
      2./3. * cI * (P3[0] * P3[1] * OM3 * TMP16)) - TMP9 * (+cI * (P2[1] *
      V1[2] + P2[0] * V1[3]) + 2./3. * cI * (P3[0] * P3[1] * OM3 * TMP30))));
  T3[7] = denom * 2. * (TMP12 * (TMP37 * (TMP26 * - 1./3. * (+cI * (P3[1] *
      P3[1] * OM3) + 1./1. * cI) - cI * (P1[1] * V2[3])) + TMP9 * (TMP25 *
      1./3. * (+cI * (P3[1] * P3[1] * OM3) + 1./1. * cI) + cI * (V2[3] *
      V1[3]))) + TMP38 * (TMP37 * (TMP16 * 1./3. * (+cI * (P3[1] * P3[1] * OM3)
      + 1./1. * cI) + cI * (P1[1] * P2[1])) + TMP9 * (TMP30 * - 1./3. * (+cI *
      (P3[1] * P3[1] * OM3) + 1./1. * cI) - cI * (P2[1] * V1[3]))));
  T3[11] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[1] * V2[4] + P1[2] *
      V2[3]) + 2./3. * cI * (P3[1] * P3[2] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[4] * V1[3] + V2[3] * V1[4]) + 2./3. * cI * (P3[1] * P3[2] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[1] * P2[2] + P1[2] * P2[1]) +
      2./3. * cI * (P3[1] * P3[2] * OM3 * TMP16)) - TMP9 * (+cI * (P2[1] *
      V1[4] + P2[2] * V1[3]) + 2./3. * cI * (P3[1] * P3[2] * OM3 * TMP30))));
  T3[15] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[1] * V2[5] + P1[3] *
      V2[3]) + 2./3. * cI * (P3[1] * P3[3] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[5] * V1[3] + V2[3] * V1[5]) + 2./3. * cI * (P3[1] * P3[3] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[1] * P2[3] + P1[3] * P2[1]) +
      2./3. * cI * (P3[1] * P3[3] * OM3 * TMP16)) - TMP9 * (+cI * (P2[1] *
      V1[5] + P2[3] * V1[3]) + 2./3. * cI * (P3[1] * P3[3] * OM3 * TMP30))));
  T3[4] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[2] * V2[2] + P1[0] *
      V2[4]) + 2./3. * cI * (P3[0] * P3[2] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[2] * V1[4] + V2[4] * V1[2]) + 2./3. * cI * (P3[0] * P3[2] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[2] * P2[0] + P1[0] * P2[2]) +
      2./3. * cI * (P3[0] * P3[2] * OM3 * TMP16)) - TMP9 * (+cI * (P2[2] *
      V1[2] + P2[0] * V1[4]) + 2./3. * cI * (P3[0] * P3[2] * OM3 * TMP30))));
  T3[8] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[2] * V2[3] + P1[1] *
      V2[4]) + 2./3. * cI * (P3[1] * P3[2] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[3] * V1[4] + V2[4] * V1[3]) + 2./3. * cI * (P3[1] * P3[2] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[2] * P2[1] + P1[1] * P2[2]) +
      2./3. * cI * (P3[1] * P3[2] * OM3 * TMP16)) - TMP9 * (+cI * (P2[2] *
      V1[3] + P2[1] * V1[4]) + 2./3. * cI * (P3[1] * P3[2] * OM3 * TMP30))));
  T3[12] = denom * 2. * (TMP12 * (TMP37 * (TMP26 * - 1./3. * (+cI * (P3[2] *
      P3[2] * OM3) + 1./1. * cI) - cI * (P1[2] * V2[4])) + TMP9 * (TMP25 *
      1./3. * (+cI * (P3[2] * P3[2] * OM3) + 1./1. * cI) + cI * (V2[4] *
      V1[4]))) + TMP38 * (TMP37 * (TMP16 * 1./3. * (+cI * (P3[2] * P3[2] * OM3)
      + 1./1. * cI) + cI * (P1[2] * P2[2])) + TMP9 * (TMP30 * - 1./3. * (+cI *
      (P3[2] * P3[2] * OM3) + 1./1. * cI) - cI * (P2[2] * V1[4]))));
  T3[16] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[2] * V2[5] + P1[3] *
      V2[4]) + 2./3. * cI * (P3[2] * P3[3] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[5] * V1[4] + V2[4] * V1[5]) + 2./3. * cI * (P3[2] * P3[3] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[2] * P2[3] + P1[3] * P2[2]) +
      2./3. * cI * (P3[2] * P3[3] * OM3 * TMP16)) - TMP9 * (+cI * (P2[2] *
      V1[5] + P2[3] * V1[4]) + 2./3. * cI * (P3[2] * P3[3] * OM3 * TMP30))));
  T3[5] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[3] * V2[2] + P1[0] *
      V2[5]) + 2./3. * cI * (P3[0] * P3[3] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[2] * V1[5] + V2[5] * V1[2]) + 2./3. * cI * (P3[0] * P3[3] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[3] * P2[0] + P1[0] * P2[3]) +
      2./3. * cI * (P3[0] * P3[3] * OM3 * TMP16)) - TMP9 * (+cI * (P2[3] *
      V1[2] + P2[0] * V1[5]) + 2./3. * cI * (P3[0] * P3[3] * OM3 * TMP30))));
  T3[9] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[3] * V2[3] + P1[1] *
      V2[5]) + 2./3. * cI * (P3[1] * P3[3] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[3] * V1[5] + V2[5] * V1[3]) + 2./3. * cI * (P3[1] * P3[3] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[3] * P2[1] + P1[1] * P2[3]) +
      2./3. * cI * (P3[1] * P3[3] * OM3 * TMP16)) - TMP9 * (+cI * (P2[3] *
      V1[3] + P2[1] * V1[5]) + 2./3. * cI * (P3[1] * P3[3] * OM3 * TMP30))));
  T3[13] = denom * (TMP12 * (TMP37 * - 1. * (+cI * (P1[3] * V2[4] + P1[2] *
      V2[5]) + 2./3. * cI * (P3[2] * P3[3] * OM3 * TMP26)) + TMP9 * (+cI *
      (V2[4] * V1[5] + V2[5] * V1[4]) + 2./3. * cI * (P3[2] * P3[3] * OM3 *
      TMP25))) + TMP38 * (TMP37 * (+cI * (P1[3] * P2[2] + P1[2] * P2[3]) +
      2./3. * cI * (P3[2] * P3[3] * OM3 * TMP16)) - TMP9 * (+cI * (P2[3] *
      V1[4] + P2[2] * V1[5]) + 2./3. * cI * (P3[2] * P3[3] * OM3 * TMP30))));
  T3[17] = denom * 2. * (TMP12 * (TMP37 * (TMP26 * - 1./3. * (+cI * (P3[3] *
      P3[3] * OM3) + 1./1. * cI) - cI * (P1[3] * V2[5])) + TMP9 * (TMP25 *
      1./3. * (+cI * (P3[3] * P3[3] * OM3) + 1./1. * cI) + cI * (V2[5] *
      V1[5]))) + TMP38 * (TMP37 * (TMP16 * 1./3. * (+cI * (P3[3] * P3[3] * OM3)
      + 1./1. * cI) + cI * (P1[3] * P2[3])) + TMP9 * (TMP30 * - 1./3. * (+cI *
      (P3[3] * P3[3] * OM3) + 1./1. * cI) - cI * (P2[3] * V1[5]))));
}

void FFV5_6_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
//   double OM3; 
  int i; 
  complex<double> denom; 
  complex<double> Vtmp[6]; 
  FFV5_3(F1, F2, COUP1, M3, W3, V3); 
  FFV6_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV5_8_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P3[4]; 
//   double OM3; 
  int i; 
  complex<double> denom; 
  complex<double> Vtmp[6]; 
  FFV5_3(F1, F2, COUP1, M3, W3, V3); 
  FFV8_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void VVT4_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP52; 
  complex<double> TMP51; 
  complex<double> TMP9; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP51 = -1. * (P1[0] * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] *
      V2[4]))) + (P1[1] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[2] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[3] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))))));
  TMP52 = -1. * (P1[0] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] *
      V1[3]))) + (P1[1] * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + (P1[2] * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + P1[3] * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))))));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (TMP37 * (OM3 * P3[0] * (TMP9 * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3]))) - 1./3. * (P3[0] * TMP51)) + (P1[0] *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + 1./3. *
      (TMP51))) + TMP38 * (OM3 * P3[0] * (TMP12 * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))) - 1./3. * (P3[0] * TMP52)) + (P2[0] *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4]))) + 1./3. *
      (TMP52))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP9 * (P2[0] * (P3[3] * V2[4]
      - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] *
      (P3[2] * V2[2] - P3[0] * V2[4]))) + 2./3. * (P3[1] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) +
      2./3. * (P3[1] * TMP52))) + P3[1] * (TMP12 * TMP38 * (P1[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP37 * TMP9 * (P2[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP37 * (P1[0] * (P2[0] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] - P3[0] * V2[5]) +
      P2[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + P1[1] * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P2[0] * (P1[0] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[4] - P3[2] * V1[2]))) + P2[1] * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP9 * (P2[0] * (P3[1] * V2[5]
      - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP52))) + P3[2] * (TMP12 * TMP38 * (P1[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP37 * TMP9 * (P2[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP37 * (P1[0] * (P2[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P2[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P1[2] * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P2[0] * (P1[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P2[2] * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP9 * (P2[0] * (P3[2] * V2[3]
      - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP52))) + P3[3] * (TMP12 * TMP38 * (P1[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP37 * TMP9 * (P2[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP37 * (P1[0] * (P2[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P1[3] * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P2[0] * (P1[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P2[3] * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP9 * (P2[0] * (P3[3] * V2[4]
      - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] *
      (P3[2] * V2[2] - P3[0] * V2[4]))) + 2./3. * (P3[1] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4]))) +
      2./3. * (P3[1] * TMP52))) + P3[1] * (TMP12 * TMP38 * (P1[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP37 * TMP9 * (P2[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP37 * (P1[0] * (P2[0] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] - P3[0] * V2[5]) +
      P2[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + P1[1] * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P2[0] * (P1[0] * (P3[2] *
      V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[4] - P3[2] * V1[2]))) + P2[1] * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[7] = denom * 2. * cI * (TMP37 * (OM3 * P3[1] * (TMP9 * (P2[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + 1./3. * (P3[1] * TMP51)) + (P1[1] *
      (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + 1./3. *
      (TMP51))) + TMP38 * (OM3 * P3[1] * (TMP12 * (P1[0] * (P3[3] * V1[4] -
      P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] *
      (P3[2] * V1[2] - P3[0] * V1[4]))) + 1./3. * (P3[1] * TMP52)) + (P2[1] *
      (P1[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P1[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) + 1./3. *
      (TMP52))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP9 * (P2[0] * (P3[1] * V2[5]
      - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP52))) + P3[2] * (TMP12 * TMP38 * (P1[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP37 * TMP9 * (P2[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP37 * (P1[1] * (P2[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P2[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P1[2] * (P2[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P2[1] * (P1[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P2[2] * (P1[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP9 * (P2[0] * (P3[2] * V2[3]
      - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP52))) + P3[3] * (TMP12 * TMP38 * (P1[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP37 * TMP9 * (P2[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP37 * (P1[1] * (P2[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P1[3] * (P2[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P2[1] * (P1[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P2[3] * (P1[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP9 * (P2[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP52))) + P3[2] * (TMP12 * TMP38 * (P1[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP37 * TMP9 * (P2[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP37 * (P1[0] * (P2[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P2[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P1[2] * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P2[0] * (P1[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P2[2] * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP9 * (P2[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      2./3. * (P3[2] * TMP52))) + P3[2] * (TMP12 * TMP38 * (P1[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP37 * TMP9 * (P2[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP37 * (P1[1] * (P2[0] * (P3[3]
      * V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] - P3[3] * V2[2]) +
      P2[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + P1[2] * (P2[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P2[1] * (P1[0] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[1] * V1[2] - P3[0] * V1[3]))) + P2[2] * (P1[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[12] = denom * 2. * cI * (TMP37 * (OM3 * P3[2] * (TMP9 * (P2[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + 1./3. * (P3[2] * TMP51)) + (P1[2] *
      (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + 1./3. *
      (TMP51))) + TMP38 * (OM3 * P3[2] * (TMP12 * (P1[0] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[3] - P3[1] * V1[2]))) + 1./3. * (P3[2] * TMP52)) + (P2[2] *
      (P1[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) + 1./3. *
      (TMP52))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP37 * (TMP9 * (P2[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP52))) + P3[3] * (TMP12 * TMP38 * (P1[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + TMP37 * TMP9 * (P2[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))))) + (TMP37 * (P1[2] * (P2[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P1[3] * (P2[0] * (P3[3] *
      V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[1] * V2[2] - P3[0] * V2[3])))) + TMP38 * (P2[2] * (P1[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P2[3] * (P1[0] * (P3[3] * V1[3] -
      P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] *
      (P3[1] * V1[2] - P3[0] * V1[3]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP9 * (P2[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP52))) + P3[3] * (TMP12 * TMP38 * (P1[1] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3]
      * (P3[2] * V1[3] - P3[1] * V1[4]))) + TMP37 * TMP9 * (P2[1] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3]
      * (P3[2] * V2[3] - P3[1] * V2[4]))))) + (TMP37 * (P1[0] * (P2[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P1[3] * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3])))) + TMP38 * (P2[0] * (P1[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P2[3] * (P1[1] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] *
      (P3[1] * V1[4] - P3[2] * V1[3]))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP9 * (P2[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP52))) + P3[3] * (TMP12 * TMP38 * (P1[0] * (P3[3] *
      V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3]
      * (P3[2] * V1[2] - P3[0] * V1[4]))) + TMP37 * TMP9 * (P2[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))))) + (TMP37 * (P1[1] * (P2[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P1[3] * (P2[0] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[4] - P3[2] * V2[2])))) + TMP38 * (P2[1] * (P1[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P2[3] * (P1[0] * (P3[2] * V1[5] -
      P3[3] * V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] *
      (P3[0] * V1[4] - P3[2] * V1[2]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP37 * (TMP9 * (P2[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * TMP51)) + TMP38 *
      (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] * (P3[0] *
      V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      2./3. * (P3[3] * TMP52))) + P3[3] * (TMP12 * TMP38 * (P1[0] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + TMP37 * TMP9 * (P2[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))))) + (TMP37 * (P1[2] * (P2[0] * (P3[1]
      * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] - P3[0] * V2[4]) +
      P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + P1[3] * (P2[0] * (P3[3] *
      V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3]
      * (P3[1] * V2[2] - P3[0] * V2[3])))) + TMP38 * (P2[2] * (P1[0] * (P3[1] *
      V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2]
      * (P3[0] * V1[3] - P3[1] * V1[2]))) + P2[3] * (P1[0] * (P3[3] * V1[3] -
      P3[1] * V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] *
      (P3[1] * V1[2] - P3[0] * V1[3]))))));
  T3[17] = denom * 2. * cI * (TMP37 * (OM3 * P3[3] * (TMP9 * (P2[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + 1./3. * (P3[3] * TMP51)) + (P1[3] *
      (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] -
      P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + 1./3. *
      (TMP51))) + TMP38 * (OM3 * P3[3] * (TMP12 * (P1[0] * (P3[2] * V1[3] -
      P3[1] * V1[4]) + (P1[1] * (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] *
      (P3[1] * V1[2] - P3[0] * V1[3]))) + 1./3. * (P3[3] * TMP52)) + (P2[3] *
      (P1[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P1[1] * (P3[2] * V1[2] -
      P3[0] * V1[4]) + P1[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) + 1./3. *
      (TMP52))));
}


void FFV8_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP3 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - 4. * cI * (OM3 * - 1./4. * P3[0] * (TMP2 + 4. * (TMP3)) +
      (+1./4. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
      F2[3]));
  V3[3] = denom * - 4. * cI * (OM3 * - 1./4. * P3[1] * (TMP2 + 4. * (TMP3)) +
      (-1./4. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
      F2[2]));
  V3[4] = denom * 4. * cI * (OM3 * 1./4. * P3[2] * (TMP2 + 4. * (TMP3)) +
      (+1./4. * cI * (F1[2] * F2[5]) - 1./4. * cI * (F1[3] * F2[4]) - cI *
      (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * 4. * cI * (OM3 * 1./4. * P3[3] * (TMP2 + 4. * (TMP3)) +
      (+1./4. * (F1[2] * F2[4]) - 1./4. * (F1[3] * F2[5]) - F1[4] * F2[2] +
      F1[5] * F2[3]));
}


void VVT9_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP34; 
  complex<double> TMP38; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP15; 
  complex<double> TMP32; 
  complex<double> denom; 
  complex<double> TMP29; 
  double OM1; 
  complex<double> TMP9; 
  double P1[4]; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = T3[0].real(); 
  P3[1] = T3[1].real(); 
  P3[2] = T3[1].imag(); 
  P3[3] = T3[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP29 = (P1[0] * - 1. * (T3[6] * V2[3] + T3[10] * V2[4] + T3[14] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[11] * V2[4] + T3[15] *
      V2[5] - T3[3] * V2[2]) + (P1[2] * (T3[8] * V2[3] + T3[12] * V2[4] +
      T3[16] * V2[5] - T3[4] * V2[2]) + P1[3] * (T3[9] * V2[3] + T3[13] * V2[4]
      + T3[17] * V2[5] - T3[5] * V2[2]))));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP32 = (P1[0] * - 1. * (T3[3] * V2[3] + T3[4] * V2[4] + T3[5] * V2[5] -
      T3[2] * V2[2]) + (P1[1] * (T3[7] * V2[3] + T3[8] * V2[4] + T3[9] * V2[5]
      - T3[6] * V2[2]) + (P1[2] * (T3[11] * V2[3] + T3[12] * V2[4] + T3[13] *
      V2[5] - T3[10] * V2[2]) + P1[3] * (T3[15] * V2[3] + T3[16] * V2[4] +
      T3[17] * V2[5] - T3[14] * V2[2]))));
  TMP34 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * - cI * (TMP38 * (OM1 * - P1[0] * (TMP15 + TMP34) + (P2[1] * -
      1. * (T3[3] + T3[6]) + (P2[2] * - 1. * (T3[4] + T3[10]) + (P2[3] * - 1. *
      (T3[5] + T3[14]) + 2. * (P2[0] * T3[2]))))) + (OM1 * - P1[0] * TMP9 *
      (TMP29 + TMP32) + P3[0] * (TMP29 + TMP32)));
  V1[3] = denom * cI * (TMP38 * (OM1 * P1[1] * (TMP15 + TMP34) + (P2[0] * - 1.
      * (T3[6] + T3[3]) + (P2[2] * (T3[8] + T3[11]) + (P2[3] * (T3[9] + T3[15])
      + 2. * (P2[1] * T3[7]))))) + (OM1 * P1[1] * TMP9 * (TMP29 + TMP32) -
      P3[1] * (TMP29 + TMP32)));
  V1[4] = denom * cI * (TMP38 * (OM1 * P1[2] * (TMP15 + TMP34) + (P2[0] * - 1.
      * (T3[10] + T3[4]) + (P2[1] * (T3[11] + T3[8]) + (P2[3] * (T3[13] +
      T3[16]) + 2. * (P2[2] * T3[12]))))) + (OM1 * P1[2] * TMP9 * (TMP29 +
      TMP32) - P3[2] * (TMP29 + TMP32)));
  V1[5] = denom * cI * (TMP38 * (OM1 * P1[3] * (TMP15 + TMP34) + (P2[0] * - 1.
      * (T3[14] + T3[5]) + (P2[1] * (T3[15] + T3[9]) + (P2[2] * (T3[16] +
      T3[13]) + 2. * (P2[3] * T3[17]))))) + (OM1 * P1[3] * TMP9 * (TMP29 +
      TMP32) - P3[3] * (TMP29 + TMP32)));
}


void VVT6_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP37; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP54; 
  complex<double> denom; 
  complex<double> TMP53; 
  complex<double> TMP52; 
  complex<double> TMP51; 
  complex<double> TMP9; 
  complex<double> TMP38; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +V1[0] + V2[0]; 
  T3[1] = +V1[1] + V2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP51 = -1. * (P1[0] * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[1] * V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] *
      V2[4]))) + (P1[1] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[2] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + P1[3] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))))));
  TMP53 = -1. * (P1[0] * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] *
      V1[4]))) + (P1[1] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] *
      V1[2]))) + (P1[2] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + P1[3] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] *
      (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] *
      V1[2]))))));
  TMP52 = -1. * (P1[0] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] *
      V1[3]))) + (P1[1] * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + (P1[2] * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + P1[3] * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))))));
  TMP54 = -1. * (P1[0] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))) + (P1[1] * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] *
      V2[4]))) + (P1[2] * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + P1[3] * (P2[0] * (P3[2] * V2[3] - P3[1] * V2[4]) + (P2[1] *
      (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] * (P3[1] * V2[2] - P3[0] *
      V2[3]))))));
  TMP38 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP37 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (OM3 * P3[0] * (TMP37 * (TMP12 * (P1[1] * (P3[2]
      * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[3] - P3[1] * V2[5]) +
      P1[3] * (P3[1] * V2[4] - P3[2] * V2[3]))) + (TMP9 * (P2[1] * (P3[2] *
      V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3]
      * (P3[1] * V2[4] - P3[2] * V2[3]))) - 1./3. * (P3[0] * (TMP51 + TMP54))))
      + TMP38 * (TMP12 * (P1[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P1[2] *
      (P3[3] * V1[3] - P3[1] * V1[5]) + P1[3] * (P3[1] * V1[4] - P3[2] *
      V1[3]))) + (TMP9 * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] *
      (P3[3] * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] *
      V1[3]))) - 1./3. * (P3[0] * (TMP53 + TMP52))))) + (TMP37 * (P1[0] *
      (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + (P2[0] *
      (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + (+1./3. *
      (TMP51 + TMP54)))) + TMP38 * (P1[0] * (P2[1] * (P3[3] * V1[4] - P3[2] *
      V1[5]) + (P2[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] *
      V1[3] - P3[1] * V1[4]))) + (P2[0] * (P1[1] * (P3[3] * V1[4] - P3[2] *
      V1[5]) + (P1[2] * (P3[1] * V1[5] - P3[3] * V1[3]) + P1[3] * (P3[2] *
      V1[3] - P3[1] * V1[4]))) + (+1./3. * (TMP53 + TMP52))))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + (TMP9 * (P2[0] * (P3[3] * V2[4] -
      P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] *
      (P3[2] * V2[2] - P3[0] * V2[4]))) + 2./3. * (P3[1] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + (TMP9 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + 2./3. * (P3[1] * (TMP53 + TMP52))))) + P3[1] * (TMP12 * (TMP37
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4])))) + TMP9 *
      (TMP37 * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      TMP38 * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4])))))) +
      (TMP37 * (P1[0] * (P2[0] * 2. * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2]
      * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[1] * (P2[1] * 2. * (P3[2] * V2[5] - P3[3] * V2[4]) +
      (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2]
      * V2[3]))) + (P1[2] * (P2[0] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[1] *
      (P3[3] * V2[3] - P3[1] * V2[5])) + P1[3] * (P2[0] * (P3[0] * V2[4] -
      P3[2] * V2[2]) + P2[1] * (P3[1] * V2[4] - P3[2] * V2[3]))))) + TMP38 *
      (P1[0] * (P2[0] * 2. * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3]
      * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      (P1[1] * (P2[1] * 2. * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3]
      * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3]))) +
      (P1[2] * (P2[0] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[1] * (P3[3] *
      V1[3] - P3[1] * V1[5])) + P1[3] * (P2[0] * (P3[0] * V1[4] - P3[2] *
      V1[2]) + P2[1] * (P3[1] * V1[4] - P3[2] * V1[3])))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + (TMP9 * (P2[0] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + 2./3. * (P3[2] * (TMP53 + TMP52))))) + P3[2] * (TMP12 * (TMP37
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4])))) + TMP9 *
      (TMP37 * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      TMP38 * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4])))))) +
      (TMP37 * (P1[0] * (P2[0] * 2. * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1]
      * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + (P1[2] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      2. * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))) + (P1[1] * (P2[0] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[2] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[3] * (P2[0] * (P3[1] * V2[2] -
      P3[0] * V2[3]) + P2[2] * (P3[1] * V2[4] - P3[2] * V2[3]))))) + TMP38 *
      (P1[0] * (P2[0] * 2. * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0]
      * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      (P1[2] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * 2. * (P3[3]
      * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3]))) +
      (P1[1] * (P2[0] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[2] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[3] * (P2[0] * (P3[1] * V1[2] - P3[0] *
      V1[3]) + P2[2] * (P3[1] * V1[4] - P3[2] * V1[3])))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 2./3. * (P3[3] * (TMP53 + TMP52))))) + P3[3] * (TMP12 * (TMP37
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4])))) + TMP9 *
      (TMP37 * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      TMP38 * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4])))))) +
      (TMP37 * (P1[0] * (P2[0] * 2. * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1]
      * (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + (P1[3] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + 2. * (P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3])))) + (P1[1] * (P2[0] * (P3[2] * V2[2] - P3[0] * V2[4]) + P2[3] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[2] * (P2[0] * (P3[0] * V2[3] -
      P3[1] * V2[2]) + P2[3] * (P3[3] * V2[3] - P3[1] * V2[5]))))) + TMP38 *
      (P1[0] * (P2[0] * 2. * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2]
      * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      (P1[3] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + 2. * (P2[3] * (P3[1] * V1[4] - P3[2] * V1[3]))))
      + (P1[1] * (P2[0] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[3] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[2] * (P2[0] * (P3[0] * V1[3] - P3[1] *
      V1[2]) + P2[3] * (P3[3] * V1[3] - P3[1] * V1[5])))))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + (TMP9 * (P2[0] * (P3[3] * V2[4] -
      P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] *
      (P3[2] * V2[2] - P3[0] * V2[4]))) + 2./3. * (P3[1] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + (TMP9 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + 2./3. * (P3[1] * (TMP53 + TMP52))))) + P3[1] * (TMP12 * (TMP37
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4])))) + TMP9 *
      (TMP37 * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      TMP38 * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4])))))) +
      (TMP37 * (P1[0] * (P2[0] * 2. * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2]
      * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[1] * (P2[1] * 2. * (P3[2] * V2[5] - P3[3] * V2[4]) +
      (P2[2] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2]
      * V2[3]))) + (P1[2] * (P2[0] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[1] *
      (P3[3] * V2[3] - P3[1] * V2[5])) + P1[3] * (P2[0] * (P3[0] * V2[4] -
      P3[2] * V2[2]) + P2[1] * (P3[1] * V2[4] - P3[2] * V2[3]))))) + TMP38 *
      (P1[0] * (P2[0] * 2. * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3]
      * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      (P1[1] * (P2[1] * 2. * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3]
      * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3]))) +
      (P1[2] * (P2[0] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[1] * (P3[3] *
      V1[3] - P3[1] * V1[5])) + P1[3] * (P2[0] * (P3[0] * V1[4] - P3[2] *
      V1[2]) + P2[1] * (P3[1] * V1[4] - P3[2] * V1[3])))))));
  T3[7] = denom * 2. * cI * (OM3 * P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[3] *
      V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P1[3]
      * (P3[2] * V2[2] - P3[0] * V2[4]))) + (TMP9 * (P2[0] * (P3[3] * V2[4] -
      P3[2] * V2[5]) + (P2[2] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] *
      (P3[2] * V2[2] - P3[0] * V2[4]))) + 1./3. * (P3[1] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + (TMP9 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] *
      (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] *
      V1[4]))) + 1./3. * (P3[1] * (TMP53 + TMP52))))) + (TMP37 * (P1[1] *
      (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + (P2[1] *
      (P1[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P1[2] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P1[3] * (P3[0] * V2[4] - P3[2] * V2[2]))) + (+1./3. *
      (TMP51 + TMP54)))) + TMP38 * (P1[1] * (P2[0] * (P3[2] * V1[5] - P3[3] *
      V1[4]) + (P2[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] *
      V1[4] - P3[2] * V1[2]))) + (P2[1] * (P1[0] * (P3[2] * V1[5] - P3[3] *
      V1[4]) + (P1[2] * (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] *
      V1[4] - P3[2] * V1[2]))) + (+1./3. * (TMP53 + TMP52))))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + (TMP9 * (P2[0] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + 2./3. * (P3[2] * (TMP53 + TMP52))))) + P3[2] * (TMP12 * (TMP37
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 *
      (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4])))) + TMP9 *
      (TMP37 * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] *
      V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) +
      TMP38 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4])))))) +
      (TMP37 * (P1[1] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] * 2.
      * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + (P1[2] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      2. * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[0] * (P2[1] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[2] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[3] * (P2[1] * (P3[1] * V2[2] -
      P3[0] * V2[3]) + P2[2] * (P3[0] * V2[4] - P3[2] * V2[2]))))) + TMP38 *
      (P1[1] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * 2. * (P3[0]
      * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      (P1[2] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * 2. * (P3[3]
      * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      (P1[0] * (P2[1] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[2] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[3] * (P2[1] * (P3[1] * V1[2] - P3[0] *
      V1[3]) + P2[2] * (P3[0] * V1[4] - P3[2] * V1[2])))))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 2./3. * (P3[3] * (TMP53 + TMP52))))) + P3[3] * (TMP12 * (TMP37
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 *
      (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4])))) + TMP9 *
      (TMP37 * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] *
      V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) +
      TMP38 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4])))))) +
      (TMP37 * (P1[1] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] * 2.
      * (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + (P1[3] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + 2. * (P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2])))) + (P1[0] * (P2[1] * (P3[1] * V2[4] - P3[2] * V2[3]) + P2[3] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[2] * (P2[1] * (P3[0] * V2[3] -
      P3[1] * V2[2]) + P2[3] * (P3[3] * V2[2] - P3[0] * V2[5]))))) + TMP38 *
      (P1[1] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * 2. * (P3[2]
      * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      (P1[3] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + 2. * (P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))))
      + (P1[0] * (P2[1] * (P3[1] * V1[4] - P3[2] * V1[3]) + P2[3] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[2] * (P2[1] * (P3[0] * V1[3] - P3[1] *
      V1[2]) + P2[3] * (P3[3] * V1[2] - P3[0] * V1[5])))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + (TMP9 * (P2[0] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + 2./3. * (P3[2] * (TMP53 + TMP52))))) + P3[2] * (TMP12 * (TMP37
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4])))) + TMP9 *
      (TMP37 * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      TMP38 * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4])))))) +
      (TMP37 * (P1[0] * (P2[0] * 2. * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1]
      * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + (P1[2] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      2. * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3]))) + (P1[1] * (P2[0] * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[2] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[3] * (P2[0] * (P3[1] * V2[2] -
      P3[0] * V2[3]) + P2[2] * (P3[1] * V2[4] - P3[2] * V2[3]))))) + TMP38 *
      (P1[0] * (P2[0] * 2. * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0]
      * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      (P1[2] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * 2. * (P3[3]
      * V1[3] - P3[1] * V1[5]) + P2[3] * (P3[1] * V1[4] - P3[2] * V1[3]))) +
      (P1[1] * (P2[0] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[2] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[3] * (P2[0] * (P3[1] * V1[2] - P3[0] *
      V1[3]) + P2[2] * (P3[1] * V1[4] - P3[2] * V1[3])))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + (TMP9 * (P2[0] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 2./3. * (P3[2] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + 2./3. * (P3[2] * (TMP53 + TMP52))))) + P3[2] * (TMP12 * (TMP37
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 *
      (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4])))) + TMP9 *
      (TMP37 * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] *
      V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) +
      TMP38 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4])))))) +
      (TMP37 * (P1[1] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] * 2.
      * (P3[0] * V2[5] - P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3]))) + (P1[2] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      2. * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2]))) + (P1[0] * (P2[1] * (P3[3] * V2[3] - P3[1] * V2[5]) + P2[2] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[3] * (P2[1] * (P3[1] * V2[2] -
      P3[0] * V2[3]) + P2[2] * (P3[0] * V2[4] - P3[2] * V2[2]))))) + TMP38 *
      (P1[1] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * 2. * (P3[0]
      * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))) +
      (P1[2] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * 2. * (P3[3]
      * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))) +
      (P1[0] * (P2[1] * (P3[3] * V1[3] - P3[1] * V1[5]) + P2[2] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[3] * (P2[1] * (P3[1] * V1[2] - P3[0] *
      V1[3]) + P2[2] * (P3[0] * V1[4] - P3[2] * V1[2])))))));
  T3[12] = denom * 2. * cI * (OM3 * P3[2] * (TMP37 * (TMP12 * (P1[0] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P1[3]
      * (P3[0] * V2[3] - P3[1] * V2[2]))) + (TMP9 * (P2[0] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + (P2[1] * (P3[3] * V2[2] - P3[0] * V2[5]) + P2[3] *
      (P3[0] * V2[3] - P3[1] * V2[2]))) + 1./3. * (P3[2] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + (TMP9 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] *
      (P3[3] * V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] *
      V1[2]))) + 1./3. * (P3[2] * (TMP53 + TMP52))))) + (TMP37 * (P1[2] *
      (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P2[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + (P2[2] *
      (P1[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P1[1] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[1] * V2[2] - P3[0] * V2[3]))) + (+1./3. *
      (TMP51 + TMP54)))) + TMP38 * (P1[2] * (P2[0] * (P3[3] * V1[3] - P3[1] *
      V1[5]) + (P2[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P2[3] * (P3[1] *
      V1[2] - P3[0] * V1[3]))) + (P2[2] * (P1[0] * (P3[3] * V1[3] - P3[1] *
      V1[5]) + (P1[1] * (P3[0] * V1[5] - P3[3] * V1[2]) + P1[3] * (P3[1] *
      V1[2] - P3[0] * V1[3]))) + (+1./3. * (TMP53 + TMP52))))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 2./3. * (P3[3] * (TMP53 + TMP52))))) + P3[3] * (TMP12 * (TMP37
      * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) + TMP38 *
      (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2])))) + TMP9 *
      (TMP37 * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] *
      V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) +
      TMP38 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2])))))) +
      (TMP37 * (P1[2] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + 2. * (P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2])))) + (P1[3] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + 2. * (P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3])))) + (P1[0] * (P2[2] * (P3[1] * V2[4] - P3[2] * V2[3]) + P2[3] *
      (P3[3] * V2[3] - P3[1] * V2[5])) + P1[1] * (P2[2] * (P3[2] * V2[2] -
      P3[0] * V2[4]) + P2[3] * (P3[0] * V2[5] - P3[3] * V2[2]))))) + TMP38 *
      (P1[2] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + 2. * (P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))))
      + (P1[3] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + 2. * (P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))))
      + (P1[0] * (P2[2] * (P3[1] * V1[4] - P3[2] * V1[3]) + P2[3] * (P3[3] *
      V1[3] - P3[1] * V1[5])) + P1[1] * (P2[2] * (P3[2] * V1[2] - P3[0] *
      V1[4]) + P2[3] * (P3[0] * V1[5] - P3[3] * V1[2])))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 2./3. * (P3[3] * (TMP53 + TMP52))))) + P3[3] * (TMP12 * (TMP37
      * (P1[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[1] * V2[5] -
      P3[3] * V2[3]) + P1[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) + TMP38 *
      (P1[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[1] * V1[5] -
      P3[3] * V1[3]) + P1[3] * (P3[2] * V1[3] - P3[1] * V1[4])))) + TMP9 *
      (TMP37 * (P2[1] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[1] *
      V2[5] - P3[3] * V2[3]) + P2[3] * (P3[2] * V2[3] - P3[1] * V2[4]))) +
      TMP38 * (P2[1] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[1] *
      V1[5] - P3[3] * V1[3]) + P2[3] * (P3[2] * V1[3] - P3[1] * V1[4])))))) +
      (TMP37 * (P1[0] * (P2[0] * 2. * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1]
      * (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + (P1[3] * (P2[1] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[3] - P3[1] * V2[5]) + 2. * (P2[3] * (P3[1] * V2[4] - P3[2] *
      V2[3])))) + (P1[1] * (P2[0] * (P3[2] * V2[2] - P3[0] * V2[4]) + P2[3] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[2] * (P2[0] * (P3[0] * V2[3] -
      P3[1] * V2[2]) + P2[3] * (P3[3] * V2[3] - P3[1] * V2[5]))))) + TMP38 *
      (P1[0] * (P2[0] * 2. * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2]
      * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      (P1[3] * (P2[1] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[3] - P3[1] * V1[5]) + 2. * (P2[3] * (P3[1] * V1[4] - P3[2] * V1[3]))))
      + (P1[1] * (P2[0] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[3] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[2] * (P2[0] * (P3[0] * V1[3] - P3[1] *
      V1[2]) + P2[3] * (P3[3] * V1[3] - P3[1] * V1[5])))))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 2./3. * (P3[3] * (TMP53 + TMP52))))) + P3[3] * (TMP12 * (TMP37
      * (P1[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P1[2] * (P3[0] * V2[5] -
      P3[3] * V2[2]) + P1[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) + TMP38 *
      (P1[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P1[2] * (P3[0] * V1[5] -
      P3[3] * V1[2]) + P1[3] * (P3[2] * V1[2] - P3[0] * V1[4])))) + TMP9 *
      (TMP37 * (P2[0] * (P3[3] * V2[4] - P3[2] * V2[5]) + (P2[2] * (P3[0] *
      V2[5] - P3[3] * V2[2]) + P2[3] * (P3[2] * V2[2] - P3[0] * V2[4]))) +
      TMP38 * (P2[0] * (P3[3] * V1[4] - P3[2] * V1[5]) + (P2[2] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + P2[3] * (P3[2] * V1[2] - P3[0] * V1[4])))))) +
      (TMP37 * (P1[1] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] * 2.
      * (P3[2] * V2[2] - P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2]))) + (P1[3] * (P2[0] * (P3[2] * V2[5] - P3[3] * V2[4]) + (P2[2] *
      (P3[3] * V2[2] - P3[0] * V2[5]) + 2. * (P2[3] * (P3[0] * V2[4] - P3[2] *
      V2[2])))) + (P1[0] * (P2[1] * (P3[1] * V2[4] - P3[2] * V2[3]) + P2[3] *
      (P3[2] * V2[5] - P3[3] * V2[4])) + P1[2] * (P2[1] * (P3[0] * V2[3] -
      P3[1] * V2[2]) + P2[3] * (P3[3] * V2[2] - P3[0] * V2[5]))))) + TMP38 *
      (P1[1] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * 2. * (P3[2]
      * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))) +
      (P1[3] * (P2[0] * (P3[2] * V1[5] - P3[3] * V1[4]) + (P2[2] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + 2. * (P2[3] * (P3[0] * V1[4] - P3[2] * V1[2]))))
      + (P1[0] * (P2[1] * (P3[1] * V1[4] - P3[2] * V1[3]) + P2[3] * (P3[2] *
      V1[5] - P3[3] * V1[4])) + P1[2] * (P2[1] * (P3[0] * V1[3] - P3[1] *
      V1[2]) + P2[3] * (P3[3] * V1[2] - P3[0] * V1[5])))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 2./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 2./3. * (P3[3] * (TMP53 + TMP52))))) + P3[3] * (TMP12 * (TMP37
      * (P1[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P1[1] * (P3[3] * V2[2] -
      P3[0] * V2[5]) + P1[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) + TMP38 *
      (P1[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P1[1] * (P3[3] * V1[2] -
      P3[0] * V1[5]) + P1[3] * (P3[0] * V1[3] - P3[1] * V1[2])))) + TMP9 *
      (TMP37 * (P2[0] * (P3[1] * V2[5] - P3[3] * V2[3]) + (P2[1] * (P3[3] *
      V2[2] - P3[0] * V2[5]) + P2[3] * (P3[0] * V2[3] - P3[1] * V2[2]))) +
      TMP38 * (P2[0] * (P3[1] * V1[5] - P3[3] * V1[3]) + (P2[1] * (P3[3] *
      V1[2] - P3[0] * V1[5]) + P2[3] * (P3[0] * V1[3] - P3[1] * V1[2])))))) +
      (TMP37 * (P1[2] * (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] *
      (P3[2] * V2[2] - P3[0] * V2[4]) + 2. * (P2[2] * (P3[0] * V2[3] - P3[1] *
      V2[2])))) + (P1[3] * (P2[0] * (P3[3] * V2[3] - P3[1] * V2[5]) + (P2[1] *
      (P3[0] * V2[5] - P3[3] * V2[2]) + 2. * (P2[3] * (P3[1] * V2[2] - P3[0] *
      V2[3])))) + (P1[0] * (P2[2] * (P3[1] * V2[4] - P3[2] * V2[3]) + P2[3] *
      (P3[3] * V2[3] - P3[1] * V2[5])) + P1[1] * (P2[2] * (P3[2] * V2[2] -
      P3[0] * V2[4]) + P2[3] * (P3[0] * V2[5] - P3[3] * V2[2]))))) + TMP38 *
      (P1[2] * (P2[0] * (P3[1] * V1[4] - P3[2] * V1[3]) + (P2[1] * (P3[2] *
      V1[2] - P3[0] * V1[4]) + 2. * (P2[2] * (P3[0] * V1[3] - P3[1] * V1[2]))))
      + (P1[3] * (P2[0] * (P3[3] * V1[3] - P3[1] * V1[5]) + (P2[1] * (P3[0] *
      V1[5] - P3[3] * V1[2]) + 2. * (P2[3] * (P3[1] * V1[2] - P3[0] * V1[3]))))
      + (P1[0] * (P2[2] * (P3[1] * V1[4] - P3[2] * V1[3]) + P2[3] * (P3[3] *
      V1[3] - P3[1] * V1[5])) + P1[1] * (P2[2] * (P3[2] * V1[2] - P3[0] *
      V1[4]) + P2[3] * (P3[0] * V1[5] - P3[3] * V1[2])))))));
  T3[17] = denom * 2. * cI * (OM3 * P3[3] * (TMP37 * (TMP12 * (P1[0] * (P3[2] *
      V2[3] - P3[1] * V2[4]) + (P1[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P1[2]
      * (P3[1] * V2[2] - P3[0] * V2[3]))) + (TMP9 * (P2[0] * (P3[2] * V2[3] -
      P3[1] * V2[4]) + (P2[1] * (P3[0] * V2[4] - P3[2] * V2[2]) + P2[2] *
      (P3[1] * V2[2] - P3[0] * V2[3]))) + 1./3. * (P3[3] * (TMP51 + TMP54)))) +
      TMP38 * (TMP12 * (P1[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P1[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P1[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + (TMP9 * (P2[0] * (P3[2] * V1[3] - P3[1] * V1[4]) + (P2[1] *
      (P3[0] * V1[4] - P3[2] * V1[2]) + P2[2] * (P3[1] * V1[2] - P3[0] *
      V1[3]))) + 1./3. * (P3[3] * (TMP53 + TMP52))))) + (TMP37 * (P1[3] *
      (P2[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P2[1] * (P3[2] * V2[2] -
      P3[0] * V2[4]) + P2[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + (P2[3] *
      (P1[0] * (P3[1] * V2[4] - P3[2] * V2[3]) + (P1[1] * (P3[2] * V2[2] -
      P3[0] * V2[4]) + P1[2] * (P3[0] * V2[3] - P3[1] * V2[2]))) + (+1./3. *
      (TMP51 + TMP54)))) + TMP38 * (P1[3] * (P2[0] * (P3[1] * V1[4] - P3[2] *
      V1[3]) + (P2[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P2[2] * (P3[0] *
      V1[3] - P3[1] * V1[2]))) + (P2[3] * (P1[0] * (P3[1] * V1[4] - P3[2] *
      V1[3]) + (P1[1] * (P3[2] * V1[2] - P3[0] * V1[4]) + P1[2] * (P3[0] *
      V1[3] - P3[1] * V1[2]))) + (+1./3. * (TMP53 + TMP52))))));
}


void FFT1_1(complex<double> F2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP13; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + T3[0]; 
  F1[1] = +F2[1] + T3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (TMP13 * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) + F2[2] * M1)) + TMP15 * (F2[4] * (P1[0] + P1[3]) + (F2[5]
      * (P1[1] + cI * (P1[2])) + F2[2] * M1)));
  F1[3] = denom * cI * (TMP13 * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) - F2[3] * M1)) + TMP15 * (F2[4] * (+cI * (P1[2]) - P1[1])
      + (F2[5] * (P1[3] - P1[0]) - F2[3] * M1)));
  F1[4] = denom * - cI * (TMP13 * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] +
      cI * (P1[2])) - F2[4] * M1)) + TMP15 * (F2[2] * (P1[3] - P1[0]) + (F2[3]
      * (P1[1] + cI * (P1[2])) - F2[4] * M1)));
  F1[5] = denom * cI * (TMP13 * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) + F2[5] * M1)) + TMP15 * (F2[2] * (+cI * (P1[2]) - P1[1])
      + (F2[3] * (P1[0] + P1[3]) + F2[5] * M1)));
}

void FFT1_2_3_5_1(complex<double> F2[], complex<double> T3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P1[4]; 
//   double P2[4]; 
  int i; 
  complex<double> Ftmp[6]; 
  FFT1_1(F2, T3, COUP1, M1, W1, F1); 
  FFT2_1(F2, T3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFT3_1(F2, T3, COUP3, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFT5_1(F2, T3, COUP4, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFT1_2_4_5_1(complex<double> F2[], complex<double> T3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P1[4]; 
//   double P2[4]; 
  int i; 
  complex<double> Ftmp[6]; 
  FFT1_1(F2, T3, COUP1, M1, W1, F1); 
  FFT2_1(F2, T3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFT4_1(F2, T3, COUP3, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFT5_1(F2, T3, COUP4, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFT4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP9; 
  complex<double> TMP8; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP8 = -1. * (F1[2] * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI *
      (P1[2]))) + (F1[3] * (F2[4] * (P1[1] - cI * (P1[2])) + F2[5] * (P1[0] -
      P1[3])) + (F1[4] * (F2[2] * (P1[3] - P1[0]) + F2[3] * (P1[1] + cI *
      (P1[2]))) + F1[5] * (F2[2] * (P1[1] - cI * (P1[2])) - F2[3] * (P1[0] +
      P1[3])))));
  TMP11 = -1. * (F1[2] * (F2[4] * (P2[0] + P2[3]) + F2[5] * (P2[1] + cI *
      (P2[2]))) + (F1[3] * (F2[4] * (P2[1] - cI * (P2[2])) + F2[5] * (P2[0] -
      P2[3])) + (F1[4] * (F2[2] * (P2[3] - P2[0]) + F2[3] * (P2[1] + cI *
      (P2[2]))) + F1[5] * (F2[2] * (P2[1] - cI * (P2[2])) - F2[3] * (P2[0] +
      P2[3])))));
  TMP10 = -1. * (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI *
      (P3[2]))) + (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] -
      P3[3])) + (F1[4] * (F2[2] * (P3[3] - P3[0]) + F2[3] * (P3[1] + cI *
      (P3[2]))) + F1[5] * (F2[2] * (P3[1] - cI * (P3[2])) - F2[3] * (P3[0] +
      P3[3])))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * - 2. * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[2] * F2[4] +
      F1[3] * F2[5] + 2./3. * (P3[0] * OM3 * TMP10) - F1[4] * F2[2] - F1[5] *
      F2[3]) + (TMP9 * (F1[2] * F2[4] + F1[3] * F2[5] + 2./3. * (P3[0] * OM3 *
      TMP10) - F1[4] * F2[2] - F1[5] * F2[3]) + (P3[0] * 1./3. * (TMP8 - TMP11)
      + TMP10 * (P2[0] - P1[0])))) + 1./3. * (TMP10 * (TMP9 - TMP12))) + (P1[0]
      * (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) +
      (P2[0] * (F1[2] * F2[4] + F1[3] * F2[5] - F1[4] * F2[2] - F1[5] * F2[3])
      + (-1./3. * (TMP8) + 1./3. * (TMP11)))));
  T3[3] = denom * - cI * (OM3 * (P3[0] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2] - 4./3. * (P3[1] * OM3 * TMP10)) +
      (TMP9 * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] *
      F2[2] - 4./3. * (P3[1] * OM3 * TMP10)) + (P3[1] * 2./3. * (TMP8 - TMP11)
      + TMP10 * (P2[1] - P1[1])))) + P3[1] * (TMP12 * (F1[4] * F2[2] + F1[5] *
      F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] *
      F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) +
      (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P1[0] - P2[0])) + (F1[3] *
      (F2[4] * (P1[0] - P2[0]) + F2[5] * (P2[1] - P1[1])) + (F1[4] * (F2[2] *
      (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] * (P1[0] -
      P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[4] = denom * - cI * (OM3 * (P3[0] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5]
      + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2]
      * OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 *
      (F1[2] * F2[4] + F1[3] * F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 *
      (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P2[2] - P1[2]) + F2[5] * (-cI *
      (P2[0]) + cI * (P1[0]))) + (F1[3] * (F2[4] * (-cI * (P1[0]) + cI *
      (P2[0])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[0]) + cI * (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0])
      + cI * (P2[0])) + F2[3] * (P1[2] - P2[2]))))));
  T3[5] = denom * - cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2])
      + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10)
      - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) +
      TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[4] * F2[2] + F1[5] *
      F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] *
      F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) +
      (F1[2] * F2[4] * (P1[0] + P2[3] - P1[3] - P2[0]) + (F1[3] * F2[5] *
      (P2[3] + P2[0] - P1[3] - P1[0]) + (F1[4] * F2[2] * (P1[3] + P1[0] - P2[3]
      - P2[0]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[6] = denom * - cI * (OM3 * (P3[0] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2] - 4./3. * (P3[1] * OM3 * TMP10)) +
      (TMP9 * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] *
      F2[2] - 4./3. * (P3[1] * OM3 * TMP10)) + (P3[1] * 2./3. * (TMP8 - TMP11)
      + TMP10 * (P2[1] - P1[1])))) + P3[1] * (TMP12 * (F1[4] * F2[2] + F1[5] *
      F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] *
      F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) +
      (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P1[0] - P2[0])) + (F1[3] *
      (F2[4] * (P1[0] - P2[0]) + F2[5] * (P2[1] - P1[1])) + (F1[4] * (F2[2] *
      (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] * (P1[0] -
      P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[7] = denom * - 2. * cI * (OM3 * (P3[1] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2] - 2./3. * (P3[1] * OM3 * TMP10)) +
      (TMP9 * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] *
      F2[2] - 2./3. * (P3[1] * OM3 * TMP10)) + (P3[1] * 1./3. * (TMP8 - TMP11)
      + TMP10 * (P2[1] - P1[1])))) + 1./3. * (TMP10 * (TMP12 - TMP9))) + (P1[1]
      * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) +
      (P2[1] * - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] *
      F2[2]) + (-1./3. * (TMP11) + 1./3. * (TMP8)))));
  T3[8] = denom * - cI * (OM3 * (P3[1] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5]
      + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2]
      * OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 *
      - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) +
      TMP10 * (P2[1] - P1[1])))) + (F1[2] * F2[5] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (F1[3] * F2[4] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (F1[4] * F2[3] * (P1[2] - cI * (P2[1]) + cI * (P1[1])
      - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) -
      P2[2])))));
  T3[9] = denom * - cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[3] * F2[5] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] * F2[2])
      + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10)
      - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11) +
      TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 * - 1. * (F1[2] * F2[5] +
      F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + TMP10 * (P2[1] -
      P1[1])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P1[3] - P2[3]))
      + (F1[3] * (F2[4] * (P1[3] - P2[3]) + F2[5] * (P2[1] - P1[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] *
      (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[10] = denom * - cI * (OM3 * (P3[0] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5]
      + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2]
      * OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 *
      (F1[2] * F2[4] + F1[3] * F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 *
      (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P2[2] - P1[2]) + F2[5] * (-cI *
      (P2[0]) + cI * (P1[0]))) + (F1[3] * (F2[4] * (-cI * (P1[0]) + cI *
      (P2[0])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[0]) + cI * (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0])
      + cI * (P2[0])) + F2[3] * (P1[2] - P2[2]))))));
  T3[11] = denom * - cI * (OM3 * (P3[1] * (TMP12 * - 1. * (-cI * (F1[2] * F2[5]
      + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2]
      * OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) + cI *
      (F1[3] * F2[4] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP10)) + (P3[2]
      * 2./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 *
      - 1. * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) +
      TMP10 * (P2[1] - P1[1])))) + (F1[2] * F2[5] * (P1[2] - cI * (P2[1]) + cI
      * (P1[1]) - P2[2]) + (F1[3] * F2[4] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (F1[4] * F2[3] * (P1[2] - cI * (P2[1]) + cI * (P1[1])
      - P2[2]) + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) -
      P2[2])))));
  T3[12] = denom * - 2. * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (-cI * (F1[2] *
      F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 2./3. *
      (P3[2] * OM3 * TMP10)) + (TMP9 * (-cI * (F1[2] * F2[5] + F1[4] * F2[3]) +
      cI * (F1[3] * F2[4] + F1[5] * F2[2]) + 2./3. * (P3[2] * OM3 * TMP10)) +
      (P3[2] * 1./3. * (TMP8 - TMP11) + TMP10 * (P2[2] - P1[2])))) + 1./3. *
      (TMP10 * (TMP12 - TMP9))) + (P1[2] * (-cI * (F1[3] * F2[4] + F1[5] *
      F2[2]) + cI * (F1[2] * F2[5] + F1[4] * F2[3])) + (P2[2] * (-cI * (F1[2] *
      F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2])) + (-1./3.
      * (TMP11) + 1./3. * (TMP8)))));
  T3[13] = denom * - cI * (OM3 * (P3[2] * (TMP12 * - 1. * (F1[3] * F2[5] +
      F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] *
      F2[2]) + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 *
      TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11)
      + TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (-cI * (F1[3] * F2[4] +
      F1[5] * F2[2]) + cI * (F1[2] * F2[5] + F1[4] * F2[3])) + (TMP9 * (-cI *
      (F1[2] * F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2])) +
      TMP10 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P1[2] - P2[2]) + F2[5] *
      (-cI * (P2[3]) + cI * (P1[3]))) + (F1[3] * (F2[4] * (-cI * (P1[3]) + cI *
      (P2[3])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] * (F2[2] * (-cI * (P1[3])
      + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[14] = denom * - cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[3] * F2[5] +
      F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] *
      F2[2]) + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 *
      TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11)
      + TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[4] * F2[2] + F1[5] *
      F2[3] - F1[2] * F2[4] - F1[3] * F2[5]) + (TMP9 * (F1[2] * F2[4] + F1[3] *
      F2[5] - F1[4] * F2[2] - F1[5] * F2[3]) + TMP10 * (P2[0] - P1[0])))) +
      (F1[2] * F2[4] * (P1[0] + P2[3] - P1[3] - P2[0]) + (F1[3] * F2[5] *
      (P2[0] + P2[3] - P1[0] - P1[3]) + (F1[4] * F2[2] * (P1[0] + P1[3] - P2[0]
      - P2[3]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[15] = denom * - cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[3] * F2[5] +
      F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] *
      F2[2]) + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 *
      TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11)
      + TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[2] * F2[5] + F1[3] *
      F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + (TMP9 * - 1. * (F1[2] * F2[5] +
      F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]) + TMP10 * (P2[1] -
      P1[1])))) + (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P1[3] - P2[3]))
      + (F1[3] * (F2[4] * (P1[3] - P2[3]) + F2[5] * (P2[1] - P1[1])) + (F1[4] *
      (F2[2] * (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] *
      (P1[3] - P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[16] = denom * - cI * (OM3 * (P3[2] * (TMP12 * - 1. * (F1[3] * F2[5] +
      F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] *
      F2[2]) + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 *
      TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP8 - TMP11)
      + TMP10 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (-cI * (F1[3] * F2[4] +
      F1[5] * F2[2]) + cI * (F1[2] * F2[5] + F1[4] * F2[3])) + (TMP9 * (-cI *
      (F1[2] * F2[5] + F1[4] * F2[3]) + cI * (F1[3] * F2[4] + F1[5] * F2[2])) +
      TMP10 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P1[2] - P2[2]) + F2[5] *
      (-cI * (P2[3]) + cI * (P1[3]))) + (F1[3] * (F2[4] * (-cI * (P1[3]) + cI *
      (P2[3])) + F2[5] * (P2[2] - P1[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] * (F2[2] * (-cI * (P1[3])
      + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[17] = denom * - 2. * cI * (OM3 * (P3[3] * (TMP12 * - 1. * (F1[3] * F2[5] +
      F1[5] * F2[3] + 2./3. * (P3[3] * OM3 * TMP10) - F1[2] * F2[4] - F1[4] *
      F2[2]) + (TMP9 * (F1[3] * F2[5] + F1[5] * F2[3] + 2./3. * (P3[3] * OM3 *
      TMP10) - F1[2] * F2[4] - F1[4] * F2[2]) + (P3[3] * 1./3. * (TMP8 - TMP11)
      + TMP10 * (P2[3] - P1[3])))) + 1./3. * (TMP10 * (TMP12 - TMP9))) + (P1[3]
      * (F1[2] * F2[4] + F1[4] * F2[2] - F1[3] * F2[5] - F1[5] * F2[3]) +
      (P2[3] * (F1[3] * F2[5] + F1[5] * F2[3] - F1[2] * F2[4] - F1[4] * F2[2])
      + (-1./3. * (TMP11) + 1./3. * (TMP8)))));
}


void FFT3_2(complex<double> F1[], complex<double> T3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + T3[0]; 
  F2[1] = +F1[1] + T3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * (F1[3] * (P2[0] * (P2[3] * (T3[9] + T3[6] + T3[15] +
      T3[3] - cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P2[1] * - 1. *
      (T3[14] + T3[5] + cI * (T3[11] + T3[8]) - 2. * (T3[7] + T3[2])) + (P2[2]
      * (T3[8] + T3[11] - 2. * cI * (T3[12] + T3[2]) + cI * (T3[14] + T3[5])) +
      (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P1[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0] * (+cI * (T3[10] + T3[4]) -
      T3[6] - T3[3]) + (P1[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) - P1[2]
      * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[16] + T3[13] - cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6]) - T3[4] -
      T3[10]) + (P2[3] * (+2. * (T3[17]) + cI * (T3[11] + T3[8]) - 2. * (T3[7])
      - T3[5] - T3[14]) + (P1[1] * (T3[3] + T3[6] - T3[15] - T3[9]) + (P1[2] *
      (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] + T3[9] - T3[3] -
      T3[6]) + (P1[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P1[3] * (T3[5] +
      T3[14] - 2. * (T3[17])))))))) + (P2[2] * (P2[3] * (-2. * cI * (T3[17]) +
      cI * (T3[5] + T3[14]) + 2. * cI * (T3[12]) - T3[8] - T3[11]) + (P1[1] *
      (-cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9])) + (P1[2] * (-cI * (T3[4]
      + T3[10]) + cI * (T3[16] + T3[13])) + (P2[2] * (-cI * (T3[16] + T3[13]) +
      cI * (T3[4] + T3[10])) + (P1[0] * - 1. * (-2. * cI * (T3[2]) + cI *
      (T3[14] + T3[5])) - P1[3] * (-2. * cI * (T3[17]) + cI * (T3[5] +
      T3[14]))))))) + P2[3] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3])
      + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[3] * (+cI *
      (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. * (+cI * (T3[11] +
      T3[8]) - 2. * (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F1[2] * (P2[0] * (P2[1] * (T3[15] + T3[9] + 2. *
      (T3[3] + T3[6]) + cI * (T3[10] + T3[4])) + (P2[2] * (T3[16] + T3[13] + 2.
      * (T3[4] + T3[10]) - cI * (T3[6] + T3[3])) + (P1[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] +
      T3[13]) + (P2[3] * 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P1[0] *
      (T3[14] + T3[5] + 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] + 2. *
      (T3[17])) - P2[0] * (T3[14] + T3[5] + 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * - 1. * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) + cI * (T3[13] +
      T3[16])) + (P2[2] * - 1. * (T3[4] + T3[10] - cI * (T3[9] + T3[15]) + 2. *
      (T3[13] + T3[16])) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[2] *
      (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * - 1. * (T3[14] + T3[5] + 2.
      * (T3[2])) + (P1[3] * (T3[5] + T3[14] + 2. * (T3[17])) - P2[3] * (T3[5] +
      T3[14] + 2. * (T3[17])))))))) + (P2[1] * (P1[0] * - 1. * (T3[6] + T3[3] +
      cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[7]) + cI *
      (T3[12])) + (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P1[2] *
      (T3[8] + T3[11] + 2. * cI * (T3[12])) - P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8]))))))) + P2[2] * (P1[0] * (+cI * (T3[6] + T3[3]) - T3[10]
      - T3[4]) + (P1[3] * (T3[13] + T3[16] - cI * (T3[9] + T3[15])) + (P1[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) + (P1[2] * - 1. * (+cI * (T3[8] +
      T3[11]) - 2. * (T3[12])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M2 * (F1[4] * (P1[1] * (T3[3] + T3[6] - T3[15] -
      T3[9]) + (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P2[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (T3[14] * (P1[0] + P1[3] - P2[0] - P2[3]) + (T3[5] * (P1[3] + P1[0] -
      P2[3] - P2[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. * (T3[17] * (P2[3] -
      P1[3]))))))))) + F1[5] * (P1[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4]))
      + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[0] * (+cI *
      (T3[10] + T3[4]) - T3[6] - T3[3]) + (P2[3] * (T3[9] + T3[15] - cI *
      (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) -
      P1[2]) + (T3[8] * (P2[2] - cI * (P2[1]) + cI * (P1[1]) - P1[2]) + (T3[12]
      * 2. * (-cI * (P2[2]) + cI * (P1[2])) + 2. * (T3[7] * (P2[1] -
      P1[1]))))))))))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * (P2[3] * (T3[9] + T3[15] - cI *
      (T3[10] + T3[4]) + cI * (T3[13] + T3[16]) - T3[6] - T3[3]) + (P2[1] *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) + cI * (T3[11] + T3[8])) + (P2[2]
      * (T3[8] + T3[11] + cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) +
      (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI
      * (T3[10] + T3[4])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] +
      T3[8])) - P1[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P2[1] *
      (P2[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13] + cI * (T3[3] + T3[15]
      + T3[6] + T3[9])) + (P2[3] * (+2. * (T3[7]) + cI * (T3[11] + T3[8]) - 2.
      * (T3[17]) - T3[5] - T3[14]) + (P1[1] * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[0] * - 1. * (T3[14] + T3[5] + 2. * (T3[2]))
      + P1[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[2] * (P2[3] *
      (T3[8] + T3[11] - cI * (T3[5] + T3[14]) - 2. * cI * (T3[17]) + 2. * cI *
      (T3[12])) + (P1[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P1[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[2] * - 1. * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P1[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P1[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P2[3] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10] + T3[4]))
      + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P2[3] *
      (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P1[1] * - 1. * (+2. *
      (T3[7]) + cI * (T3[11] + T3[8])) - P1[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))))) + (F1[3] * (P2[0] * (P2[1] * - 1. * (T3[15] + T3[9] + cI
      * (T3[10] + T3[4]) - 2. * (T3[3] + T3[6])) + (P2[2] * (+2. * (T3[4] +
      T3[10]) + cI * (T3[6] + T3[3]) - T3[16] - T3[13]) + (P1[1] * (T3[15] +
      T3[9] - T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) +
      (P2[3] * 2. * (T3[5] + T3[14] - T3[17] - T3[2]) + (P1[0] * - 1. * (T3[14]
      + T3[5] - 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] - 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * (T3[3] + T3[6] + cI * (T3[13] + T3[16]) - 2. * (T3[9] + T3[15]))
      + (P2[2] * - 1. * (+2. * (T3[13] + T3[16]) + cI * (T3[9] + T3[15]) -
      T3[4] - T3[10]) + (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P1[2] *
      (T3[16] + T3[13] - T3[4] - T3[10]) + (P1[0] * - 1. * (T3[14] + T3[5] - 2.
      * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[1] * (P1[0] * (+cI * (T3[10]
      + T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] +
      T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[12]) + cI *
      (T3[7])) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      (P1[2] * (T3[8] + T3[11] - 2. * cI * (T3[12])) + P2[1] * (+cI * (T3[11] +
      T3[8]) - 2. * (T3[7]))))))) + P2[2] * (P1[0] * - 1. * (T3[10] + T3[4] +
      cI * (T3[6] + T3[3])) + (P1[3] * (T3[13] + T3[16] + cI * (T3[9] +
      T3[15])) + (P1[1] * (T3[11] + T3[8] + 2. * cI * (T3[7])) + (P1[2] * (+2.
      * (T3[12]) + cI * (T3[8] + T3[11])) - P2[2] * (+2. * (T3[12]) + cI *
      (T3[8] + T3[11]))))))))) + M2 * (F1[4] * (P1[0] * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] +
      T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) +
      (P2[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] -
      cI * (P1[1]) + cI * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) +
      cI * (P2[1]) - P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) +
      2. * (T3[7] * (P2[1] - P1[1]))))))))) + F1[5] * (P1[1] * (T3[3] + T3[15]
      + T3[6] + T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] *
      - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16]
      + T3[10] + T3[13]) + (T3[14] * (P2[0] + P1[3] - P1[0] - P2[3]) + (T3[5] *
      (P1[3] + P2[0] - P2[3] - P1[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. *
      (T3[17] * (P1[3] - P2[3]))))))))))));
  F2[4] = denom * cI * (F1[5] * (P2[0] * (P2[3] * (T3[6] + T3[3] - cI * (T3[10]
      + T3[4]) + cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P2[1] * - 1. *
      (T3[14] + T3[5] + 2. * (T3[7] + T3[2]) - cI * (T3[11] + T3[8])) + (P2[2]
      * (+cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2]) - T3[8] - T3[11])
      + (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] +
      T3[15] - cI * (T3[13] + T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10]
      + T3[4])) + (P1[1] * - 1. * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) +
      P1[2] * (T3[8] + T3[11] - 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[4] + T3[16] + T3[10] + T3[13] - cI * (T3[3] + T3[15] + T3[6] +
      T3[9])) + (P2[3] * (T3[5] + T3[14] + 2. * (T3[17]) + cI * (T3[11] +
      T3[8]) - 2. * (T3[7])) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9])
      + (P1[2] * - 1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * (T3[3] +
      T3[15] + T3[6] + T3[9]) + (P1[0] * (T3[14] + T3[5] + 2. * (T3[2])) -
      P1[3] * (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[2] * (P2[3] * - 1. *
      (T3[8] + T3[11] - 2. * cI * (T3[12]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17])) + (P1[1] * (+cI * (T3[3] + T3[15] + T3[6] + T3[9])) + (P1[2] *
      (+cI * (T3[4] + T3[16] + T3[10] + T3[13])) + (P2[2] * - 1. * (+cI *
      (T3[4] + T3[16] + T3[10] + T3[13])) + (P1[0] * - 1. * (+cI * (T3[14] +
      T3[5]) + 2. * cI * (T3[2])) + P1[3] * (+cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[17]))))))) + P2[3] * (P1[0] * (+cI * (T3[10] + T3[4]) - T3[6] -
      T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] + T3[16])) + (P2[3] *
      (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (P1[1] * - 1. * (+cI *
      (T3[11] + T3[8]) - 2. * (T3[7])) + P1[2] * (T3[8] + T3[11] - 2. * cI *
      (T3[12]))))))))) + (F1[4] * (P2[0] * (P2[1] * (+2. * (T3[3] + T3[6]) + cI
      * (T3[10] + T3[4]) - T3[15] - T3[9]) + (P2[2] * - 1. * (T3[16] + T3[13] +
      cI * (T3[6] + T3[3]) - 2. * (T3[4] + T3[10])) + (P1[1] * (T3[15] + T3[9]
      - T3[3] - T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[3] *
      2. * (T3[5] + T3[14] - T3[17] - T3[2]) + (P1[0] * - 1. * (T3[14] + T3[5]
      - 2. * (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) +
      P2[0] * (T3[14] + T3[5] - 2. * (T3[2]))))))))) + (P2[3] * (P2[1] * - 1. *
      (+2. * (T3[9] + T3[15]) + cI * (T3[13] + T3[16]) - T3[3] - T3[6]) +
      (P2[2] * (T3[4] + T3[10] + cI * (T3[9] + T3[15]) - 2. * (T3[13] +
      T3[16])) + (P1[1] * (T3[15] + T3[9] - T3[3] - T3[6]) + (P1[2] * (T3[16] +
      T3[13] - T3[4] - T3[10]) + (P1[0] * - 1. * (T3[14] + T3[5] - 2. *
      (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] - 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[1] * (P1[0] * - 1. * (T3[6]
      + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI *
      (T3[13] + T3[16])) + (P2[2] * - 2. * (T3[8] + T3[11] - cI * (T3[7]) + cI
      * (T3[12])) + (P1[1] * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) + (P1[2] *
      (T3[8] + T3[11] + 2. * cI * (T3[12])) - P2[1] * (+2. * (T3[7]) + cI *
      (T3[11] + T3[8]))))))) + P2[2] * (P1[0] * (+cI * (T3[6] + T3[3]) - T3[10]
      - T3[4]) + (P1[3] * (T3[13] + T3[16] - cI * (T3[9] + T3[15])) + (P1[1] *
      (T3[11] + T3[8] - 2. * cI * (T3[7])) + (P1[2] * - 1. * (+cI * (T3[8] +
      T3[11]) - 2. * (T3[12])) + P2[2] * (+cI * (T3[8] + T3[11]) - 2. *
      (T3[12]))))))))) + M2 * (F1[2] * (P1[1] * (T3[3] + T3[15] + T3[6] +
      T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) + (P2[1] * - 1. *
      (T3[3] + T3[15] + T3[6] + T3[9]) + (P2[2] * - 1. * (T3[4] + T3[16] +
      T3[10] + T3[13]) + (T3[14] * (P2[0] + P1[3] - P1[0] - P2[3]) + (T3[5] *
      (P1[3] + P2[0] - P2[3] - P1[0]) + (T3[2] * 2. * (P2[0] - P1[0]) + 2. *
      (T3[17] * (P1[3] - P2[3]))))))))) + F1[3] * (P1[0] * (+cI * (T3[10] +
      T3[4]) - T3[6] - T3[3]) + (P1[3] * (T3[9] + T3[15] - cI * (T3[13] +
      T3[16])) + (P2[0] * (T3[6] + T3[3] - cI * (T3[10] + T3[4])) + (P2[3] *
      (+cI * (T3[13] + T3[16]) - T3[9] - T3[15]) + (T3[11] * (P1[2] - cI *
      (P1[1]) + cI * (P2[1]) - P2[2]) + (T3[8] * (P1[2] - cI * (P1[1]) + cI *
      (P2[1]) - P2[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. *
      (T3[7] * (P1[1] - P2[1]))))))))))));
  F2[5] = denom * - cI * (F1[4] * (P2[0] * (P2[3] * (T3[9] + T3[6] + T3[15] +
      T3[3] + cI * (T3[13] + T3[10] + T3[16] + T3[4])) + (P2[1] * (+2. * (T3[7]
      + T3[2]) + cI * (T3[11] + T3[8]) - T3[14] - T3[5]) + (P2[2] * (T3[8] +
      T3[11] - cI * (T3[14] + T3[5]) + 2. * cI * (T3[12] + T3[2])) + (P1[0] *
      (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15]
      + cI * (T3[13] + T3[16])) + (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[1] * - 1. * (+2. * (T3[7]) + cI * (T3[11] + T3[8])) -
      P1[2] * (T3[8] + T3[11] + 2. * cI * (T3[12]))))))))) + (P2[1] * (P2[2] *
      (T3[16] + T3[13] - cI * (T3[3] + T3[6]) + cI * (T3[15] + T3[9]) - T3[4] -
      T3[10]) + (P2[3] * - 1. * (T3[5] + T3[14] + 2. * (T3[7]) + cI * (T3[11] +
      T3[8]) - 2. * (T3[17])) + (P1[1] * (T3[3] + T3[6] - T3[15] - T3[9]) +
      (P1[2] * (T3[4] + T3[10] - T3[16] - T3[13]) + (P2[1] * (T3[15] + T3[9] -
      T3[3] - T3[6]) + (P1[0] * (T3[14] + T3[5] - 2. * (T3[2])) + P1[3] *
      (T3[5] + T3[14] - 2. * (T3[17])))))))) + (P2[2] * (P2[3] * - 1. * (T3[8]
      + T3[11] - 2. * cI * (T3[17]) + cI * (T3[5] + T3[14]) + 2. * cI *
      (T3[12])) + (P1[1] * (-cI * (T3[15] + T3[9]) + cI * (T3[3] + T3[6])) +
      (P1[2] * (-cI * (T3[16] + T3[13]) + cI * (T3[4] + T3[10])) + (P2[2] *
      (-cI * (T3[4] + T3[10]) + cI * (T3[16] + T3[13])) + (P1[0] * (-2. * cI *
      (T3[2]) + cI * (T3[14] + T3[5])) + P1[3] * (-2. * cI * (T3[17]) + cI *
      (T3[5] + T3[14]))))))) + P2[3] * (P1[0] * - 1. * (T3[6] + T3[3] + cI *
      (T3[10] + T3[4])) + (P1[3] * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) + (P1[1] * (+2.
      * (T3[7]) + cI * (T3[11] + T3[8])) + P1[2] * (T3[8] + T3[11] + 2. * cI *
      (T3[12]))))))))) + (F1[5] * (P2[0] * (P2[1] * - 1. * (T3[15] + T3[9] + 2.
      * (T3[3] + T3[6]) - cI * (T3[10] + T3[4])) + (P2[2] * - 1. * (T3[16] +
      T3[13] + 2. * (T3[4] + T3[10]) + cI * (T3[6] + T3[3])) + (P1[1] * (T3[3]
      + T3[15] + T3[6] + T3[9]) + (P1[2] * (T3[4] + T3[16] + T3[10] + T3[13]) +
      (P2[3] * - 2. * (T3[5] + T3[17] + T3[2] + T3[14]) + (P1[0] * - 1. *
      (T3[14] + T3[5] + 2. * (T3[2])) + (P1[3] * (T3[5] + T3[14] + 2. *
      (T3[17])) + P2[0] * (T3[14] + T3[5] + 2. * (T3[2]))))))))) + (P2[3] *
      (P2[1] * (T3[3] + T3[6] + 2. * (T3[9] + T3[15]) - cI * (T3[13] + T3[16]))
      + (P2[2] * (T3[4] + T3[10] + 2. * (T3[13] + T3[16]) + cI * (T3[9] +
      T3[15])) + (P1[1] * - 1. * (T3[3] + T3[15] + T3[6] + T3[9]) + (P1[2] * -
      1. * (T3[4] + T3[16] + T3[10] + T3[13]) + (P1[0] * (T3[14] + T3[5] + 2. *
      (T3[2])) + (P1[3] * - 1. * (T3[5] + T3[14] + 2. * (T3[17])) + P2[3] *
      (T3[5] + T3[14] + 2. * (T3[17])))))))) + (P2[1] * (P1[0] * (T3[6] + T3[3]
      - cI * (T3[10] + T3[4])) + (P1[3] * (+cI * (T3[13] + T3[16]) - T3[9] -
      T3[15]) + (P2[2] * 2. * (T3[8] + T3[11] - cI * (T3[12]) + cI * (T3[7])) +
      (P1[1] * (+cI * (T3[11] + T3[8]) - 2. * (T3[7])) + (P1[2] * - 1. * (T3[8]
      + T3[11] - 2. * cI * (T3[12])) - P2[1] * (+cI * (T3[11] + T3[8]) - 2. *
      (T3[7]))))))) + P2[2] * (P1[0] * (T3[10] + T3[4] + cI * (T3[6] + T3[3]))
      + (P1[3] * - 1. * (T3[13] + T3[16] + cI * (T3[9] + T3[15])) + (P1[1] * -
      1. * (T3[11] + T3[8] + 2. * cI * (T3[7])) + (P1[2] * - 1. * (+2. *
      (T3[12]) + cI * (T3[8] + T3[11])) + P2[2] * (+2. * (T3[12]) + cI * (T3[8]
      + T3[11]))))))))) + M2 * (F1[2] * (P1[0] * (T3[6] + T3[3] + cI * (T3[10]
      + T3[4])) + (P1[3] * - 1. * (T3[9] + T3[15] + cI * (T3[13] + T3[16])) +
      (P2[0] * - 1. * (T3[6] + T3[3] + cI * (T3[10] + T3[4])) + (P2[3] * (T3[9]
      + T3[15] + cI * (T3[13] + T3[16])) + (T3[11] * (P2[2] - cI * (P1[1]) + cI
      * (P2[1]) - P1[2]) + (T3[8] * (P2[2] - cI * (P1[1]) + cI * (P2[1]) -
      P1[2]) + (T3[12] * 2. * (-cI * (P1[2]) + cI * (P2[2])) + 2. * (T3[7] *
      (P2[1] - P1[1]))))))))) + F1[3] * (P1[1] * (T3[15] + T3[9] - T3[3] -
      T3[6]) + (P1[2] * (T3[16] + T3[13] - T3[4] - T3[10]) + (P2[1] * (T3[3] +
      T3[6] - T3[15] - T3[9]) + (P2[2] * (T3[4] + T3[10] - T3[16] - T3[13]) +
      (T3[14] * (P2[0] + P2[3] - P1[0] - P1[3]) + (T3[5] * (P2[3] + P2[0] -
      P1[3] - P1[0]) + (T3[2] * 2. * (P1[0] - P2[0]) + 2. * (T3[17] * (P1[3] -
      P2[3]))))))))))));
}


void FFV6_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  complex<double> TMP0; 
  TMP1 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * (-cI * (TMP0) + 2. * cI * (TMP1)); 
}


void FFT5_2(complex<double> F1[], complex<double> T3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP13; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + T3[0]; 
  F2[1] = +F1[1] + T3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (TMP13 * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) + F1[2] * M2)) + TMP15 * (F1[4] * (P2[0] - P2[3]) +
      (F1[5] * (+cI * (P2[2]) - P2[1]) + F1[2] * M2)));
  F2[3] = denom * - cI * (TMP13 * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * -
      1. * (P2[0] + P2[3]) - F1[3] * M2)) + TMP15 * (F1[4] * (P2[1] + cI *
      (P2[2])) + (F1[5] * - 1. * (P2[0] + P2[3]) - F1[3] * M2)));
  F2[4] = denom * - cI * (TMP13 * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) - F1[4] * M2)) + TMP15 * (F1[2] * - 1. * (P2[0] +
      P2[3]) + (F1[3] * (+cI * (P2[2]) - P2[1]) - F1[4] * M2)));
  F2[5] = denom * cI * (TMP13 * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) + F1[5] * M2)) + TMP15 * (F1[2] * (P2[1] + cI * (P2[2]))
      + (F1[3] * (P2[0] - P2[3]) + F1[5] * M2)));
}


void VVT7_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP26; 
  double OM1; 
  complex<double> TMP13; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * TMP26 * (+cI * (TMP13 + TMP15)) - V2[2] * (+cI
      * (TMP13 + TMP15)));
  V1[3] = denom * (OM1 * P1[1] * TMP26 * (+cI * (TMP13 + TMP15)) - V2[3] * (+cI
      * (TMP13 + TMP15)));
  V1[4] = denom * (OM1 * P1[2] * TMP26 * (+cI * (TMP13 + TMP15)) - V2[4] * (+cI
      * (TMP13 + TMP15)));
  V1[5] = denom * (OM1 * P1[3] * TMP26 * (+cI * (TMP13 + TMP15)) - V2[5] * (+cI
      * (TMP13 + TMP15)));
}


void FFT2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP22; 
  double P2[4]; 
  complex<double> TMP23; 
  double P3[4]; 
  double OM3; 
  complex<double> denom; 
  complex<double> TMP24; 
  complex<double> TMP9; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  T3[0] = +F1[0] + F2[0]; 
  T3[1] = +F1[1] + F2[1]; 
  P3[0] = -T3[0].real(); 
  P3[1] = -T3[1].real(); 
  P3[2] = -T3[1].imag(); 
  P3[3] = -T3[0].imag(); 
  TMP24 = (F1[2] * (F2[4] * (P2[0] + P2[3]) + F2[5] * (P2[1] + cI * (P2[2]))) +
      (F1[3] * (F2[4] * (P2[1] - cI * (P2[2])) + F2[5] * (P2[0] - P2[3])) +
      (F1[4] * (F2[2] * (P2[0] - P2[3]) - F2[3] * (P2[1] + cI * (P2[2]))) +
      F1[5] * (F2[2] * (+cI * (P2[2]) - P2[1]) + F2[3] * (P2[0] + P2[3])))));
  TMP23 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
  TMP9 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP22 = (F1[2] * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI * (P1[2]))) +
      (F1[3] * (F2[4] * (P1[1] - cI * (P1[2])) + F2[5] * (P1[0] - P1[3])) +
      (F1[4] * (F2[2] * (P1[0] - P1[3]) - F2[3] * (P1[1] + cI * (P1[2]))) +
      F1[5] * (F2[2] * (+cI * (P1[2]) - P1[1]) + F2[3] * (P1[0] + P1[3])))));
  TMP12 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  T3[2] = denom * 2. * cI * (OM3 * (P3[0] * (TMP12 * (F1[2] * F2[4] + F1[3] *
      F2[5] + F1[4] * F2[2] + F1[5] * F2[3] - 2./3. * (P3[0] * OM3 * TMP23)) +
      (TMP9 * - 1. * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] *
      F2[3] - 2./3. * (P3[0] * OM3 * TMP23)) + (P3[0] * 1./3. * (TMP22 - TMP24)
      + TMP23 * (P2[0] - P1[0])))) + 1./3. * (TMP23 * (TMP9 - TMP12))) + (P1[0]
      * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) +
      (P2[0] * - 1. * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] *
      F2[3]) + (-1./3. * (TMP22) + 1./3. * (TMP24)))));
  T3[3] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[2] * F2[5] + F1[3] *
      F2[4] + 4./3. * (P3[1] * OM3 * TMP23) - F1[4] * F2[3] - F1[5] * F2[2]) +
      (TMP9 * (F1[2] * F2[5] + F1[3] * F2[4] + 4./3. * (P3[1] * OM3 * TMP23) -
      F1[4] * F2[3] - F1[5] * F2[2]) + (P3[1] * 2./3. * (TMP22 - TMP24) + TMP23
      * (P2[1] - P1[1])))) + P3[1] * (TMP12 * (F1[2] * F2[4] + F1[3] * F2[5] +
      F1[4] * F2[2] + F1[5] * F2[3]) + (TMP9 * - 1. * (F1[2] * F2[4] + F1[3] *
      F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + TMP23 * (P2[0] - P1[0])))) +
      (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P2[0] - P1[0])) + (F1[3] *
      (F2[4] * (P2[0] - P1[0]) + F2[5] * (P1[1] - P2[1])) + (F1[4] * (F2[2] *
      (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] * (P1[0] -
      P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[4] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (-cI * (F1[3] * F2[4] +
      F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP23)) + (TMP9 * (-cI * (F1[3] * F2[4] + F1[4] * F2[3]) + cI *
      (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP23)) + (P3[2]
      * 2./3. * (TMP22 - TMP24) + TMP23 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + (TMP9 *
      - 1. * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) +
      TMP23 * (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P1[2] - P2[2]) + F2[5] *
      (-cI * (P1[0]) + cI * (P2[0]))) + (F1[3] * (F2[4] * (-cI * (P2[0]) + cI *
      (P1[0])) + F2[5] * (P1[2] - P2[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[0]) + cI * (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0])
      + cI * (P2[0])) + F2[3] * (P1[2] - P2[2]))))));
  T3[5] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[2] * F2[4] + F1[5] *
      F2[3] + 4./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] * F2[2]) +
      (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP23) -
      F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP22 - TMP24) + TMP23
      * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[2] * F2[4] + F1[3] * F2[5] +
      F1[4] * F2[2] + F1[5] * F2[3]) + (TMP9 * - 1. * (F1[2] * F2[4] + F1[3] *
      F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + TMP23 * (P2[0] - P1[0])))) +
      (F1[2] * F2[4] * (P1[3] + P2[0] - P1[0] - P2[3]) + (F1[3] * F2[5] *
      (P1[3] + P1[0] - P2[3] - P2[0]) + (F1[4] * F2[2] * (P1[3] + P1[0] - P2[3]
      - P2[0]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[6] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[2] * F2[5] + F1[3] *
      F2[4] + 4./3. * (P3[1] * OM3 * TMP23) - F1[4] * F2[3] - F1[5] * F2[2]) +
      (TMP9 * (F1[2] * F2[5] + F1[3] * F2[4] + 4./3. * (P3[1] * OM3 * TMP23) -
      F1[4] * F2[3] - F1[5] * F2[2]) + (P3[1] * 2./3. * (TMP22 - TMP24) + TMP23
      * (P2[1] - P1[1])))) + P3[1] * (TMP12 * (F1[2] * F2[4] + F1[3] * F2[5] +
      F1[4] * F2[2] + F1[5] * F2[3]) + (TMP9 * - 1. * (F1[2] * F2[4] + F1[3] *
      F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + TMP23 * (P2[0] - P1[0])))) +
      (F1[2] * (F2[4] * (P1[1] - P2[1]) + F2[5] * (P2[0] - P1[0])) + (F1[3] *
      (F2[4] * (P2[0] - P1[0]) + F2[5] * (P1[1] - P2[1])) + (F1[4] * (F2[2] *
      (P1[1] - P2[1]) + F2[3] * (P1[0] - P2[0])) + F1[5] * (F2[2] * (P1[0] -
      P2[0]) + F2[3] * (P1[1] - P2[1]))))));
  T3[7] = denom * 2. * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[2] * F2[5] +
      F1[3] * F2[4] + 2./3. * (P3[1] * OM3 * TMP23) - F1[4] * F2[3] - F1[5] *
      F2[2]) + (TMP9 * (F1[2] * F2[5] + F1[3] * F2[4] + 2./3. * (P3[1] * OM3 *
      TMP23) - F1[4] * F2[3] - F1[5] * F2[2]) + (P3[1] * 1./3. * (TMP22 -
      TMP24) + TMP23 * (P2[1] - P1[1])))) + 1./3. * (TMP23 * (TMP12 - TMP9))) +
      (P1[1] * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3] * F2[4])
      + (P2[1] * (F1[2] * F2[5] + F1[3] * F2[4] - F1[4] * F2[3] - F1[5] *
      F2[2]) + (-1./3. * (TMP24) + 1./3. * (TMP22)))));
  T3[8] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (-cI * (F1[3] * F2[4] +
      F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP23)) + (TMP9 * (-cI * (F1[3] * F2[4] + F1[4] * F2[3]) + cI *
      (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP23)) + (P3[2]
      * 2./3. * (TMP22 - TMP24) + TMP23 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3] * F2[4]) + (TMP9 *
      (F1[2] * F2[5] + F1[3] * F2[4] - F1[4] * F2[3] - F1[5] * F2[2]) + TMP23 *
      (P2[1] - P1[1])))) + (F1[2] * F2[5] * (P2[2] - cI * (P1[1]) + cI *
      (P2[1]) - P1[2]) + (F1[3] * F2[4] * (P2[2] - cI * (P2[1]) + cI * (P1[1])
      - P1[2]) + (F1[4] * F2[3] * (P1[2] - cI * (P2[1]) + cI * (P1[1]) - P2[2])
      + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) - P2[2])))));
  T3[9] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[2] * F2[4] + F1[5] *
      F2[3] + 4./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] * F2[2]) +
      (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP23) -
      F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP22 - TMP24) + TMP23
      * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[4] * F2[3] + F1[5] * F2[2] -
      F1[2] * F2[5] - F1[3] * F2[4]) + (TMP9 * (F1[2] * F2[5] + F1[3] * F2[4] -
      F1[4] * F2[3] - F1[5] * F2[2]) + TMP23 * (P2[1] - P1[1])))) + (F1[2] *
      (F2[4] * (P2[1] - P1[1]) + F2[5] * (P2[3] - P1[3])) + (F1[3] * (F2[4] *
      (P2[3] - P1[3]) + F2[5] * (P1[1] - P2[1])) + (F1[4] * (F2[2] * (P1[1] -
      P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] * (P1[3] - P2[3]) +
      F2[3] * (P2[1] - P1[1]))))));
  T3[10] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (-cI * (F1[3] * F2[4] +
      F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP23)) + (TMP9 * (-cI * (F1[3] * F2[4] + F1[4] * F2[3]) + cI *
      (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP23)) + (P3[2]
      * 2./3. * (TMP22 - TMP24) + TMP23 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + (TMP9 *
      - 1. * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) +
      TMP23 * (P2[0] - P1[0])))) + (F1[2] * (F2[4] * (P1[2] - P2[2]) + F2[5] *
      (-cI * (P1[0]) + cI * (P2[0]))) + (F1[3] * (F2[4] * (-cI * (P2[0]) + cI *
      (P1[0])) + F2[5] * (P1[2] - P2[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[0]) + cI * (P1[0]))) + F1[5] * (F2[2] * (-cI * (P1[0])
      + cI * (P2[0])) + F2[3] * (P1[2] - P2[2]))))));
  T3[11] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (-cI * (F1[3] * F2[4] +
      F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] *
      OM3 * TMP23)) + (TMP9 * (-cI * (F1[3] * F2[4] + F1[4] * F2[3]) + cI *
      (F1[2] * F2[5] + F1[5] * F2[2]) + 4./3. * (P3[2] * OM3 * TMP23)) + (P3[2]
      * 2./3. * (TMP22 - TMP24) + TMP23 * (P2[2] - P1[2])))) + P3[2] * (TMP12 *
      (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3] * F2[4]) + (TMP9 *
      (F1[2] * F2[5] + F1[3] * F2[4] - F1[4] * F2[3] - F1[5] * F2[2]) + TMP23 *
      (P2[1] - P1[1])))) + (F1[2] * F2[5] * (P2[2] - cI * (P1[1]) + cI *
      (P2[1]) - P1[2]) + (F1[3] * F2[4] * (P2[2] - cI * (P2[1]) + cI * (P1[1])
      - P1[2]) + (F1[4] * F2[3] * (P1[2] - cI * (P2[1]) + cI * (P1[1]) - P2[2])
      + F1[5] * F2[2] * (P1[2] - cI * (P1[1]) + cI * (P2[1]) - P2[2])))));
  T3[12] = denom * 2. * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (-cI * (F1[3] *
      F2[4] + F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2]) + 2./3. *
      (P3[2] * OM3 * TMP23)) + (TMP9 * (-cI * (F1[3] * F2[4] + F1[4] * F2[3]) +
      cI * (F1[2] * F2[5] + F1[5] * F2[2]) + 2./3. * (P3[2] * OM3 * TMP23)) +
      (P3[2] * 1./3. * (TMP22 - TMP24) + TMP23 * (P2[2] - P1[2])))) + 1./3. *
      (TMP23 * (TMP12 - TMP9))) + (P1[2] * (-cI * (F1[2] * F2[5] + F1[5] *
      F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3])) + (P2[2] * (-cI * (F1[3] *
      F2[4] + F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2])) + (-1./3.
      * (TMP24) + 1./3. * (TMP22)))));
  T3[13] = denom * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (F1[2] * F2[4] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] * F2[2])
      + (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP23)
      - F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP22 - TMP24) +
      TMP23 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (-cI * (F1[2] * F2[5] +
      F1[5] * F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3])) + (TMP9 * (-cI *
      (F1[3] * F2[4] + F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2])) +
      TMP23 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P2[2] - P1[2]) + F2[5] *
      (-cI * (P1[3]) + cI * (P2[3]))) + (F1[3] * (F2[4] * (-cI * (P2[3]) + cI *
      (P1[3])) + F2[5] * (P1[2] - P2[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] * (F2[2] * (-cI * (P1[3])
      + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[14] = denom * cI * (OM3 * (P3[0] * (TMP12 * - 1. * (F1[2] * F2[4] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] * F2[2])
      + (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP23)
      - F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP22 - TMP24) +
      TMP23 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[2] * F2[4] + F1[3] *
      F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + (TMP9 * - 1. * (F1[2] * F2[4] +
      F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]) + TMP23 * (P2[0] -
      P1[0])))) + (F1[2] * F2[4] * (P1[3] + P2[0] - P1[0] - P2[3]) + (F1[3] *
      F2[5] * (P1[0] + P1[3] - P2[0] - P2[3]) + (F1[4] * F2[2] * (P1[0] + P1[3]
      - P2[0] - P2[3]) + F1[5] * F2[3] * (P1[3] + P2[0] - P1[0] - P2[3])))));
  T3[15] = denom * cI * (OM3 * (P3[1] * (TMP12 * - 1. * (F1[2] * F2[4] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] * F2[2])
      + (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP23)
      - F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP22 - TMP24) +
      TMP23 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (F1[4] * F2[3] + F1[5] *
      F2[2] - F1[2] * F2[5] - F1[3] * F2[4]) + (TMP9 * (F1[2] * F2[5] + F1[3] *
      F2[4] - F1[4] * F2[3] - F1[5] * F2[2]) + TMP23 * (P2[1] - P1[1])))) +
      (F1[2] * (F2[4] * (P2[1] - P1[1]) + F2[5] * (P2[3] - P1[3])) + (F1[3] *
      (F2[4] * (P2[3] - P1[3]) + F2[5] * (P1[1] - P2[1])) + (F1[4] * (F2[2] *
      (P1[1] - P2[1]) + F2[3] * (P1[3] - P2[3])) + F1[5] * (F2[2] * (P1[3] -
      P2[3]) + F2[3] * (P2[1] - P1[1]))))));
  T3[16] = denom * cI * (OM3 * (P3[2] * (TMP12 * - 1. * (F1[2] * F2[4] + F1[5]
      * F2[3] + 4./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] * F2[2])
      + (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 4./3. * (P3[3] * OM3 * TMP23)
      - F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 2./3. * (TMP22 - TMP24) +
      TMP23 * (P2[3] - P1[3])))) + P3[3] * (TMP12 * (-cI * (F1[2] * F2[5] +
      F1[5] * F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3])) + (TMP9 * (-cI *
      (F1[3] * F2[4] + F1[4] * F2[3]) + cI * (F1[2] * F2[5] + F1[5] * F2[2])) +
      TMP23 * (P2[2] - P1[2])))) + (F1[2] * (F2[4] * (P2[2] - P1[2]) + F2[5] *
      (-cI * (P1[3]) + cI * (P2[3]))) + (F1[3] * (F2[4] * (-cI * (P2[3]) + cI *
      (P1[3])) + F2[5] * (P1[2] - P2[2])) + (F1[4] * (F2[2] * (P1[2] - P2[2]) +
      F2[3] * (-cI * (P2[3]) + cI * (P1[3]))) + F1[5] * (F2[2] * (-cI * (P1[3])
      + cI * (P2[3])) + F2[3] * (P2[2] - P1[2]))))));
  T3[17] = denom * 2. * cI * (OM3 * (P3[3] * (TMP12 * - 1. * (F1[2] * F2[4] +
      F1[5] * F2[3] + 2./3. * (P3[3] * OM3 * TMP23) - F1[3] * F2[5] - F1[4] *
      F2[2]) + (TMP9 * (F1[2] * F2[4] + F1[5] * F2[3] + 2./3. * (P3[3] * OM3 *
      TMP23) - F1[3] * F2[5] - F1[4] * F2[2]) + (P3[3] * 1./3. * (TMP22 -
      TMP24) + TMP23 * (P2[3] - P1[3])))) + 1./3. * (TMP23 * (TMP12 - TMP9))) +
      (P1[3] * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5] * F2[3])
      + (P2[3] * (F1[2] * F2[4] + F1[5] * F2[3] - F1[3] * F2[5] - F1[4] *
      F2[2]) + (-1./3. * (TMP24) + 1./3. * (TMP22)))));
}


void VVT8_1(complex<double> V2[], complex<double> T3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP16; 
  complex<double> TMP15; 
  complex<double> TMP26; 
  complex<double> TMP13; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + T3[0]; 
  V1[1] = +V2[1] + T3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (P1[0] * - 1. * (P2[1] * T3[3] + P2[2] * T3[4] + P2[3] * T3[5] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[8] + P2[3] * T3[9]
      - P2[0] * T3[6]) + (P1[2] * (P2[1] * T3[11] + P2[2] * T3[12] + P2[3] *
      T3[13] - P2[0] * T3[10]) + P1[3] * (P2[1] * T3[15] + P2[2] * T3[16] +
      P2[3] * T3[17] - P2[0] * T3[14]))));
  TMP26 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP13 = (P1[0] * - 1. * (P2[1] * T3[6] + P2[2] * T3[10] + P2[3] * T3[14] -
      P2[0] * T3[2]) + (P1[1] * (P2[1] * T3[7] + P2[2] * T3[11] + P2[3] *
      T3[15] - P2[0] * T3[3]) + (P1[2] * (P2[1] * T3[8] + P2[2] * T3[12] +
      P2[3] * T3[16] - P2[0] * T3[4]) + P1[3] * (P2[1] * T3[9] + P2[2] * T3[13]
      + P2[3] * T3[17] - P2[0] * T3[5]))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (P2[0] * - TMP26 * (+cI * (TMP13 + TMP15)) + TMP16 * V2[2] *
      (+cI * (TMP13 + TMP15)));
  V1[3] = denom * (P2[1] * - TMP26 * (+cI * (TMP13 + TMP15)) + TMP16 * V2[3] *
      (+cI * (TMP13 + TMP15)));
  V1[4] = denom * (P2[2] * - TMP26 * (+cI * (TMP13 + TMP15)) + TMP16 * V2[4] *
      (+cI * (TMP13 + TMP15)));
  V1[5] = denom * (P2[3] * - TMP26 * (+cI * (TMP13 + TMP15)) + TMP16 * V2[5] *
      (+cI * (TMP13 + TMP15)));
}


void FFV7_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  complex<double> TMP0; 
  TMP1 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - 1. * (+cI * (TMP0) + 2. * cI * (TMP1)); 
}


void FFV3_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP8; 
  TMP8 = (F1[2] * (F2[4] * - 1. * (V3[2] + V3[5]) - F2[5] * (V3[3] + cI *
      (V3[4]))) + (F1[3] * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5] -
      V3[2])) + (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI *
      (V3[4]))) + F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5])))));
  vertex = COUP * - cI * TMP8; 
}


void VVV2_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP6; 
  complex<double> TMP5; 
  complex<double> TMP4; 
  complex<double> TMP3; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP5 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP4 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP7 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP3 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  vertex = COUP * (TMP3 * - 1. * (+cI * (TMP0 + TMP2)) + (+cI * (TMP4 * TMP5 +
      TMP6 * TMP7)));
}

void VVV1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP9 = -1. * (P1[0] * (V1[3] * (V3[5] * V2[4] - V3[4] * V2[5]) + (V1[4] *
      (V3[3] * V2[5] - V3[5] * V2[3]) + V1[5] * (V3[4] * V2[3] - V3[3] *
      V2[4]))) + (P1[1] * (V1[2] * (V3[4] * V2[5] - V3[5] * V2[4]) + (V1[4] *
      (V3[5] * V2[2] - V3[2] * V2[5]) + V1[5] * (V3[2] * V2[4] - V3[4] *
      V2[2]))) + (P1[2] * (V1[2] * (V3[5] * V2[3] - V3[3] * V2[5]) + (V1[3] *
      (V3[2] * V2[5] - V3[5] * V2[2]) + V1[5] * (V3[3] * V2[2] - V3[2] *
      V2[3]))) + P1[3] * (V1[2] * (V3[3] * V2[4] - V3[4] * V2[3]) + (V1[3] *
      (V3[4] * V2[2] - V3[2] * V2[4]) + V1[4] * (V3[2] * V2[3] - V3[3] *
      V2[2]))))));
  TMP10 = -1. * (P2[0] * (V1[3] * (V3[5] * V2[4] - V3[4] * V2[5]) + (V1[4] *
      (V3[3] * V2[5] - V3[5] * V2[3]) + V1[5] * (V3[4] * V2[3] - V3[3] *
      V2[4]))) + (P2[1] * (V1[2] * (V3[4] * V2[5] - V3[5] * V2[4]) + (V1[4] *
      (V3[5] * V2[2] - V3[2] * V2[5]) + V1[5] * (V3[2] * V2[4] - V3[4] *
      V2[2]))) + (P2[2] * (V1[2] * (V3[5] * V2[3] - V3[3] * V2[5]) + (V1[3] *
      (V3[2] * V2[5] - V3[5] * V2[2]) + V1[5] * (V3[3] * V2[2] - V3[2] *
      V2[3]))) + P2[3] * (V1[2] * (V3[3] * V2[4] - V3[4] * V2[3]) + (V1[3] *
      (V3[4] * V2[2] - V3[2] * V2[4]) + V1[4] * (V3[2] * V2[3] - V3[3] *
      V2[2]))))));
  vertex = COUP * (-cI * (TMP9) + cI * (TMP10)); 
}

void VVV1_2_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> tmp; 
  VVV1_0(V1, V2, V3, COUP1, vertex); 
  VVV2_0(V1, V2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV4_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP13; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP13 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  vertex = COUP * TMP13 * (-cI * (TMP0) + cI * (TMP2)); 
}


void FFV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  TMP1 = (F2[4] * F1[4] + F2[5] * F1[5] - F2[2] * F1[2] - F2[3] * F1[3]); 
  TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  vertex = COUP * TMP1 * (-cI * (TMP0) + cI * (TMP2)); 
}

void FFV1_2_3_4_0(complex<double> F1[], complex<double> F2[], complex<double>
    V3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> tmp; 
  FFV1_0(F1, F2, V3, COUP1, vertex); 
  FFV2_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
  FFV3_0(F1, F2, V3, COUP3, tmp); 
  vertex = vertex + tmp; 
  FFV4_0(F1, F2, V3, COUP4, tmp); 
  vertex = vertex + tmp; 
}


void FFV3_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
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
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] *
      (V3[3] + cI * (V3[4])))));
  F1[3] = denom * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[3] * (P1[0] * - 1. * (V3[2] + V3[5]) + (P1[1]
      * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] *
      (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * - 1. * (V3[2] + V3[5]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[2]
      + V3[5])))) + (F2[5] * (P1[0] * - 1. * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * - 1.
      * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] -
      cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] *
      (V3[2] + V3[5]))));
}


void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  complex<double> TMP9; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  TMP9 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP11 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (TMP11 * (F1[4] * (P2[3] - P2[0]) + (F1[5] * (P2[1] - cI
      * (P2[2])) + F1[2] * M2)) + TMP9 * (F1[4] * (P2[0] - P2[3]) + (F1[5] *
      (+cI * (P2[2]) - P2[1]) - F1[2] * M2)));
  F2[3] = denom * cI * (TMP11 * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * - 1.
      * (P2[0] + P2[3]) + F1[3] * M2)) + TMP9 * (F1[4] * - 1. * (P2[1] + cI *
      (P2[2])) + (F1[5] * (P2[0] + P2[3]) - F1[3] * M2)));
  F2[4] = denom * - cI * (TMP11 * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) + F1[4] * M2)) + TMP9 * (F1[2] * (P2[0] + P2[3])
      + (F1[3] * (P2[1] - cI * (P2[2])) - F1[4] * M2)));
  F2[5] = denom * - cI * (TMP11 * (F1[2] * - 1. * (P2[1] + cI * (P2[2])) +
      (F1[3] * (P2[3] - P2[0]) + F1[5] * M2)) + TMP9 * (F1[2] * (P2[1] + cI *
      (P2[2])) + (F1[3] * (P2[0] - P2[3]) - F1[5] * M2)));
}

void FFV1_2_3_4_2(complex<double> F1[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV1_2(F1, V3, COUP1, M2, W2, F2); 
  FFV2_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFV3_2(F1, V3, COUP3, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
  FFV4_2(F1, V3, COUP4, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}


void VVV2_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  complex<double> TMP17; 
  double P3[4]; 
  complex<double> TMP20; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP14; 
  complex<double> TMP19; 
  complex<double> TMP18; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +V1[0] + V2[0]; 
  V3[1] = +V1[1] + V2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP20 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP19 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP18 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP14 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP17 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP10 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP12 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * (OM3 * P3[0] * (TMP10 * (+cI * (TMP17 + TMP18)) + (-cI *
      (TMP12 * TMP19 + TMP14 * TMP20))) + (TMP10 * - 1. * (+cI * (P1[0] +
      P2[0])) + (+cI * (V1[2] * TMP12 + V2[2] * TMP14))));
  V3[3] = denom * (OM3 * P3[1] * (TMP10 * (+cI * (TMP17 + TMP18)) + (-cI *
      (TMP12 * TMP19 + TMP14 * TMP20))) + (TMP10 * - 1. * (+cI * (P1[1] +
      P2[1])) + (+cI * (V1[3] * TMP12 + V2[3] * TMP14))));
  V3[4] = denom * (OM3 * P3[2] * (TMP10 * (+cI * (TMP17 + TMP18)) + (-cI *
      (TMP12 * TMP19 + TMP14 * TMP20))) + (TMP10 * - 1. * (+cI * (P1[2] +
      P2[2])) + (+cI * (V1[4] * TMP12 + V2[4] * TMP14))));
  V3[5] = denom * (OM3 * P3[3] * (TMP10 * (+cI * (TMP17 + TMP18)) + (-cI *
      (TMP12 * TMP19 + TMP14 * TMP20))) + (TMP10 * - 1. * (+cI * (P1[3] +
      P2[3])) + (+cI * (V1[5] * TMP12 + V2[5] * TMP14))));
}

void VVV1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> denom; 
  double OM3; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +V1[0] + V2[0]; 
  V3[1] = +V1[1] + V2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP7 = -1. * (P1[0] * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P3[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP8 = -1. * (P2[0] * (P3[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P3[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P3[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P2[1] * (P3[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P3[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P3[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P2[2] * (P3[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P3[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P3[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P2[3] * (P3[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P3[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P3[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * cI * (V1[3] * (V2[4] * (P2[3] - P1[3]) + V2[5] * (P1[2] -
      P2[2])) + (V1[4] * (V2[3] * (P1[3] - P2[3]) + V2[5] * (P2[1] - P1[1])) +
      (V1[5] * (V2[3] * (P2[2] - P1[2]) + V2[4] * (P1[1] - P2[1])) + OM3 *
      P3[0] * (TMP7 - TMP8))));
  V3[3] = denom * - cI * (V1[2] * (V2[4] * (P1[3] - P2[3]) + V2[5] * (P2[2] -
      P1[2])) + (V1[4] * (V2[2] * (P2[3] - P1[3]) + V2[5] * (P1[0] - P2[0])) +
      (V1[5] * (V2[2] * (P1[2] - P2[2]) + V2[4] * (P2[0] - P1[0])) + OM3 *
      P3[1] * (TMP8 - TMP7))));
  V3[4] = denom * - cI * (V1[2] * (V2[3] * (P2[3] - P1[3]) + V2[5] * (P1[1] -
      P2[1])) + (V1[3] * (V2[2] * (P1[3] - P2[3]) + V2[5] * (P2[0] - P1[0])) +
      (V1[5] * (V2[2] * (P2[1] - P1[1]) + V2[3] * (P1[0] - P2[0])) + OM3 *
      P3[2] * (TMP8 - TMP7))));
  V3[5] = denom * - cI * (V1[2] * (V2[3] * (P1[2] - P2[2]) + V2[4] * (P2[1] -
      P1[1])) + (V1[3] * (V2[2] * (P2[2] - P1[2]) + V2[4] * (P1[0] - P2[0])) +
      (V1[4] * (V2[2] * (P1[1] - P2[1]) + V2[3] * (P2[0] - P1[0])) + OM3 *
      P3[3] * (TMP8 - TMP7))));
}

void VVV1_2_3(complex<double> V1[], complex<double> V2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
//   double P3[4]; 
  complex<double> denom; 
//   double OM3; 
  int i; 
  complex<double> Vtmp[6]; 
  VVV1_3(V1, V2, COUP1, M3, W3, V3); 
  VVV2_3(V1, V2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFV4_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  complex<double> TMP9; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  TMP9 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP11 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (TMP11 * (F2[4] * - 1. * (P1[0] + P1[3]) + (F2[5] * -
      1. * (P1[1] + cI * (P1[2])) + F2[2] * M1)) + TMP9 * (F2[4] * (P1[0] +
      P1[3]) + (F2[5] * (P1[1] + cI * (P1[2])) - F2[2] * M1)));
  F1[3] = denom * - cI * (TMP11 * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) + F2[3] * M1)) + TMP9 * (F2[4] * (P1[1] - cI * (P1[2])) +
      (F2[5] * (P1[0] - P1[3]) - F2[3] * M1)));
  F1[4] = denom * - cI * (TMP11 * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] +
      cI * (P1[2])) + F2[4] * M1)) + TMP9 * (F2[2] * (P1[0] - P1[3]) + (F2[3] *
      - 1. * (P1[1] + cI * (P1[2])) - F2[4] * M1)));
  F1[5] = denom * - cI * (TMP11 * (F2[2] * (P1[1] - cI * (P1[2])) + (F2[3] * -
      1. * (P1[0] + P1[3]) + F2[5] * M1)) + TMP9 * (F2[2] * (+cI * (P1[2]) -
      P1[1]) + (F2[3] * (P1[0] + P1[3]) - F2[5] * M1)));
}


void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  complex<double> TMP9; 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  TMP9 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP11 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * (TMP11 * (F2[4] * - 1. * (P1[0] + P1[3]) + (F2[5] * -
      1. * (P1[1] + cI * (P1[2])) - F2[2] * M1)) + TMP9 * (F2[4] * (P1[0] +
      P1[3]) + (F2[5] * (P1[1] + cI * (P1[2])) + F2[2] * M1)));
  F1[3] = denom * - cI * (TMP11 * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) - F2[3] * M1)) + TMP9 * (F2[4] * (P1[1] - cI * (P1[2])) +
      (F2[5] * (P1[0] - P1[3]) + F2[3] * M1)));
  F1[4] = denom * cI * (TMP11 * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) - F2[4] * M1)) + TMP9 * (F2[2] * (P1[0] - P1[3]) + (F2[3] * -
      1. * (P1[1] + cI * (P1[2])) + F2[4] * M1)));
  F1[5] = denom * cI * (TMP11 * (F2[2] * (P1[1] - cI * (P1[2])) + (F2[3] * - 1.
      * (P1[0] + P1[3]) - F2[5] * M1)) + TMP9 * (F2[2] * (+cI * (P1[2]) -
      P1[1]) + (F2[3] * (P1[0] + P1[3]) + F2[5] * M1)));
}

void FFV1_2_3_4_1(complex<double> F2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFV1_1(F2, V3, COUP1, M1, W1, F1); 
  FFV2_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFV3_1(F2, V3, COUP3, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
  FFV4_1(F2, V3, COUP4, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}


void FFV3_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  complex<double> TMP22; 
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
  TMP22 = -1. * (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI *
      (P3[2]))) + (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] -
      P3[3])) + (F1[4] * (F2[2] * (P3[3] - P3[0]) + F2[3] * (P3[1] + cI *
      (P3[2]))) + F1[5] * (F2[2] * (P3[1] - cI * (P3[2])) - F2[3] * (P3[0] +
      P3[3])))));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[4] * F2[2] + F1[5] * F2[3] - F1[2] * F2[4] - F1[3]
      * F2[5] - P3[0] * OM3 * TMP22);
  V3[3] = denom * - cI * (F1[2] * F2[5] + F1[3] * F2[4] + F1[4] * F2[3] + F1[5]
      * F2[2] - P3[1] * OM3 * TMP22);
  V3[4] = denom * - cI * (-cI * (F1[3] * F2[4] + F1[5] * F2[2]) + cI * (F1[2] *
      F2[5] + F1[4] * F2[3]) - P3[2] * OM3 * TMP22);
  V3[5] = denom * - cI * (F1[2] * F2[4] + F1[4] * F2[2] - F1[3] * F2[5] - F1[5]
      * F2[3] - P3[3] * OM3 * TMP22);
}


void FFV4_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  complex<double> TMP9; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  TMP9 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP11 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (TMP11 * (F1[4] * (P2[3] - P2[0]) + (F1[5] * (P2[1] - cI
      * (P2[2])) - F1[2] * M2)) + TMP9 * (F1[4] * (P2[0] - P2[3]) + (F1[5] *
      (+cI * (P2[2]) - P2[1]) + F1[2] * M2)));
  F2[3] = denom * cI * (TMP11 * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * - 1.
      * (P2[0] + P2[3]) - F1[3] * M2)) + TMP9 * (F1[4] * - 1. * (P2[1] + cI *
      (P2[2])) + (F1[5] * (P2[0] + P2[3]) + F1[3] * M2)));
  F2[4] = denom * cI * (TMP11 * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] * (+cI
      * (P2[2]) - P2[1]) - F1[4] * M2)) + TMP9 * (F1[2] * (P2[0] + P2[3]) +
      (F1[3] * (P2[1] - cI * (P2[2])) + F1[4] * M2)));
  F2[5] = denom * cI * (TMP11 * (F1[2] * - 1. * (P2[1] + cI * (P2[2])) + (F1[3]
      * (P2[3] - P2[0]) - F1[5] * M2)) + TMP9 * (F1[2] * (P2[1] + cI * (P2[2]))
      + (F1[3] * (P2[0] - P2[3]) + F1[5] * M2)));
}


void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P3[4]; 
  complex<double> TMP26; 
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
  TMP26 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3] - P3[0] * OM3 * TMP26);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4] - P3[1] * OM3 * TMP26);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]) - P3[2] * OM3 * TMP26);
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3] - P3[3] * OM3 * TMP26);
}


void FFV3_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
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
  F2[2] = denom * cI * (F1[2] * (P2[0] * - 1. * (V3[2] + V3[5]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[2]
      + V3[5])))) + (F1[3] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] * (V3[2]
      - V3[5]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI * (V3[4]) -
      V3[3]))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * -
      1. * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] *
      (V3[2] + V3[5]))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * (V3[2] - V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) + P2[3] *
      (V3[2] - V3[5])))) + (F1[5] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] *
      (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] *
      (+cI * (V3[4]) - V3[3]))));
  F2[5] = denom * - cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * - 1. * (V3[2] + V3[5]) +
      (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void FFV1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP23; 
  complex<double> TMP17; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP18; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP17 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP18 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (F1[4] * F2[4] + F1[5] * F2[5] - F1[2] * F2[2] - F1[3] * F2[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * TMP23 * (OM3 * P3[0] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[0]) + cI * (P2[0])));
  V3[3] = denom * TMP23 * (OM3 * P3[1] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[1]) + cI * (P2[1])));
  V3[4] = denom * TMP23 * (OM3 * P3[2] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[2]) + cI * (P2[2])));
  V3[5] = denom * TMP23 * (OM3 * P3[3] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[3]) + cI * (P2[3])));
}

void FFV1_2_3_4_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M3, double W3, complex<double> V3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
//   double P1[4]; 
//   double P2[4]; 
//   double P3[4]; 
//   double OM3; 
  complex<double> Vtmp[6]; 
  int i; 
  FFV1_3(F1, F2, COUP1, M3, W3, V3); 
  FFV2_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
  FFV3_3(F1, F2, COUP3, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
  FFV4_3(F1, F2, COUP4, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void VVV2_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP16; 
  complex<double> TMP15; 
  complex<double> denom; 
  double OM1; 
  complex<double> TMP9; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP9 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP16 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP11 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP12 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (-cI * (TMP15 * TMP16) + cI * (TMP11 * TMP12))
      + (V2[2] * - 1. * (+cI * (TMP9 + TMP11)) + (+cI * (V3[2] * TMP12 + P2[0]
      * TMP15))));
  V1[3] = denom * (OM1 * P1[1] * (-cI * (TMP15 * TMP16) + cI * (TMP11 * TMP12))
      + (V2[3] * - 1. * (+cI * (TMP9 + TMP11)) + (+cI * (V3[3] * TMP12 + P2[1]
      * TMP15))));
  V1[4] = denom * (OM1 * P1[2] * (-cI * (TMP15 * TMP16) + cI * (TMP11 * TMP12))
      + (V2[4] * - 1. * (+cI * (TMP9 + TMP11)) + (+cI * (V3[4] * TMP12 + P2[2]
      * TMP15))));
  V1[5] = denom * (OM1 * P1[3] * (-cI * (TMP15 * TMP16) + cI * (TMP11 * TMP12))
      + (V2[5] * - 1. * (+cI * (TMP9 + TMP11)) + (+cI * (V3[5] * TMP12 + P2[3]
      * TMP15))));
}


void VVV1_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  double OM1; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP6 = -1. * (P1[0] * (P2[1] * (V3[4] * V2[5] - V3[5] * V2[4]) + (P2[2] *
      (V3[5] * V2[3] - V3[3] * V2[5]) + P2[3] * (V3[3] * V2[4] - V3[4] *
      V2[3]))) + (P1[1] * (P2[0] * (V3[5] * V2[4] - V3[4] * V2[5]) + (P2[2] *
      (V3[2] * V2[5] - V3[5] * V2[2]) + P2[3] * (V3[4] * V2[2] - V3[2] *
      V2[4]))) + (P1[2] * (P2[0] * (V3[3] * V2[5] - V3[5] * V2[3]) + (P2[1] *
      (V3[5] * V2[2] - V3[2] * V2[5]) + P2[3] * (V3[2] * V2[3] - V3[3] *
      V2[2]))) + P1[3] * (P2[0] * (V3[4] * V2[3] - V3[3] * V2[4]) + (P2[1] *
      (V3[2] * V2[4] - V3[4] * V2[2]) + P2[2] * (V3[3] * V2[2] - V3[2] *
      V2[3]))))));
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * cI * (V2[3] * (V3[4] * (P2[3] - P1[3]) + V3[5] * (P1[2] -
      P2[2])) + (V2[4] * (V3[3] * (P1[3] - P2[3]) + V3[5] * (P2[1] - P1[1])) +
      (V2[5] * (V3[3] * (P2[2] - P1[2]) + V3[4] * (P1[1] - P2[1])) - P1[0] *
      OM1 * TMP6)));
  V1[3] = denom * - cI * (V2[2] * (V3[4] * (P1[3] - P2[3]) + V3[5] * (P2[2] -
      P1[2])) + (V2[4] * (V3[2] * (P2[3] - P1[3]) + V3[5] * (P1[0] - P2[0])) +
      (V2[5] * (V3[2] * (P1[2] - P2[2]) + V3[4] * (P2[0] - P1[0])) + P1[1] *
      OM1 * TMP6)));
  V1[4] = denom * - cI * (V2[2] * (V3[3] * (P2[3] - P1[3]) + V3[5] * (P1[1] -
      P2[1])) + (V2[3] * (V3[2] * (P1[3] - P2[3]) + V3[5] * (P2[0] - P1[0])) +
      (V2[5] * (V3[2] * (P2[1] - P1[1]) + V3[3] * (P1[0] - P2[0])) + P1[2] *
      OM1 * TMP6)));
  V1[5] = denom * - cI * (V2[2] * (V3[3] * (P1[2] - P2[2]) + V3[4] * (P2[1] -
      P1[1])) + (V2[3] * (V3[2] * (P2[2] - P1[2]) + V3[4] * (P1[0] - P2[0])) +
      (V2[4] * (V3[2] * (P1[1] - P2[1]) + V3[3] * (P2[0] - P1[0])) + P1[3] *
      OM1 * TMP6)));
}

void VVV1_2_1(complex<double> V2[], complex<double> V3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> V1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Vtmp[6]; 
//   double OM1; 
  VVV1_1(V2, V3, COUP1, M1, W1, V1); 
  VVV2_1(V2, V3, COUP2, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
}

void FFV4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP17; 
  double P3[4]; 
  double OM3; 
  complex<double> TMP24; 
  complex<double> TMP18; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P2[0] = F2[0].real(); 
  P2[1] = F2[1].real(); 
  P2[2] = F2[1].imag(); 
  P2[3] = F2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./pow(M3, 2); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP24 = (F1[2] * F2[2] + F1[3] * F2[3] + F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP17 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP18 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * TMP24 * (OM3 * P3[0] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[0]) + cI * (P2[0])));
  V3[3] = denom * TMP24 * (OM3 * P3[1] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[1]) + cI * (P2[1])));
  V3[4] = denom * TMP24 * (OM3 * P3[2] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[2]) + cI * (P2[2])));
  V3[5] = denom * TMP24 * (OM3 * P3[3] * (-cI * (TMP18) + cI * (TMP17)) + (-cI
      * (P1[3]) + cI * (P2[3])));
}

void VVS3_4_5_0(complex<double> V1[], complex<double> V2[], complex<double>
    S3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  VVS3_0(V1, V2, S3, COUP1, vertex); 
  VVS4_0(V1, V2, S3, COUP2, tmp); 
  vertex = vertex + tmp; 
  VVS5_0(V1, V2, S3, COUP3, tmp); 
  vertex = vertex + tmp; 
}

void VVS3_4_5_3(complex<double> V1[], complex<double> V2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, double M3, double W3,
    complex<double> S3[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Stmp[3]; 
//   double P3[4]; 
  complex<double> denom; 
  int i; 
  VVS3_3(V1, V2, COUP1, M3, W3, S3); 
  VVS4_3(V1, V2, COUP2, M3, W3, Stmp); 
  i = 2; 
  while (i < 3)
  {
    S3[i] = S3[i] + Stmp[i]; 
    i++; 
  }
  VVS5_3(V1, V2, COUP3, M3, W3, Stmp); 
  i = 2; 
  while (i < 3)
  {
    S3[i] = S3[i] + Stmp[i]; 
    i++; 
  }
}


void VVS1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP7 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP6 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * 4. * (-cI * (TMP6) + cI * (TMP7)); 
}


void VVS2_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP7; 
  complex<double> TMP6; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP7 = -1. * (P1[0] * (P2[1] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[3] * V1[5] - V2[5] * V1[3]) + P2[3] * (V2[4] * V1[3] - V2[3] *
      V1[4]))) + (P1[1] * (P2[0] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[4] - V2[4] *
      V1[2]))) + (P1[2] * (P2[0] * (V2[5] * V1[3] - V2[3] * V1[5]) + (P2[1] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[3] * V1[2] - V2[2] *
      V1[3]))) + P1[3] * (P2[0] * (V2[3] * V1[4] - V2[4] * V1[3]) + (P2[1] *
      (V2[4] * V1[2] - V2[2] * V1[4]) + P2[2] * (V2[2] * V1[3] - V2[3] *
      V1[2]))))));
  TMP6 = -1. * (P1[0] * (P2[1] * (V2[4] * V1[5] - V2[5] * V1[4]) + (P2[2] *
      (V2[5] * V1[3] - V2[3] * V1[5]) + P2[3] * (V2[3] * V1[4] - V2[4] *
      V1[3]))) + (P1[1] * (P2[0] * (V2[5] * V1[4] - V2[4] * V1[5]) + (P2[2] *
      (V2[2] * V1[5] - V2[5] * V1[2]) + P2[3] * (V2[4] * V1[2] - V2[2] *
      V1[4]))) + (P1[2] * (P2[0] * (V2[3] * V1[5] - V2[5] * V1[3]) + (P2[1] *
      (V2[5] * V1[2] - V2[2] * V1[5]) + P2[3] * (V2[2] * V1[3] - V2[3] *
      V1[2]))) + P1[3] * (P2[0] * (V2[4] * V1[3] - V2[3] * V1[4]) + (P2[1] *
      (V2[2] * V1[4] - V2[4] * V1[2]) + P2[2] * (V2[3] * V1[2] - V2[2] *
      V1[3]))))));
  vertex = COUP * S3[2] * (-cI * (TMP7) + cI * (TMP6)); 
}


void VVS2_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * - 2. * cI * S3[2] * (P1[1] * (P2[3] * V2[4] - P2[2] * V2[5])
      + (P1[2] * (P2[1] * V2[5] - P2[3] * V2[3]) + P1[3] * (P2[2] * V2[3] -
      P2[1] * V2[4])));
  V1[3] = denom * 2. * cI * S3[2] * (P1[0] * (P2[2] * V2[5] - P2[3] * V2[4]) +
      (P1[2] * (P2[3] * V2[2] - P2[0] * V2[5]) + P1[3] * (P2[0] * V2[4] - P2[2]
      * V2[2])));
  V1[4] = denom * 2. * cI * S3[2] * (P1[0] * (P2[3] * V2[3] - P2[1] * V2[5]) +
      (P1[1] * (P2[0] * V2[5] - P2[3] * V2[2]) + P1[3] * (P2[1] * V2[2] - P2[0]
      * V2[3])));
  V1[5] = denom * 2. * cI * S3[2] * (P1[0] * (P2[1] * V2[4] - P2[2] * V2[3]) +
      (P1[1] * (P2[2] * V2[2] - P2[0] * V2[4]) + P1[2] * (P2[0] * V2[3] - P2[1]
      * V2[2])));
}


void VVS3_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP8; 
  TMP8 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  vertex = COUP * - cI * TMP8 * S3[2]; 
}


void VVS3_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  double OM1; 
  complex<double> TMP9; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (-cI * (V2[2]) + cI * (P1[0] * OM1 * TMP9)); 
  V1[3] = denom * S3[2] * (-cI * (V2[3]) + cI * (P1[1] * OM1 * TMP9)); 
  V1[4] = denom * S3[2] * (-cI * (V2[4]) + cI * (P1[2] * OM1 * TMP9)); 
  V1[5] = denom * S3[2] * (-cI * (V2[5]) + cI * (P1[3] * OM1 * TMP9)); 
}


void VVS3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  complex<double> TMP8; 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP8 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP8; 
}


void VVS3_4_5_1(complex<double> V2[], complex<double> S3[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, double M1, double W1,
    complex<double> V1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
  complex<double> denom; 
  complex<double> Vtmp[6]; 
//   double OM1; 
  int i; 
  VVS3_1(V2, S3, COUP1, M1, W1, V1); 
  VVS4_1(V2, S3, COUP2, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
  VVS5_1(V2, S3, COUP3, M1, W1, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V1[i] = V1[i] + Vtmp[i]; 
    i++; 
  }
}


void VVS4_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP9; 
  complex<double> TMP13; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP15 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP8 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP13 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  vertex = COUP * S3[2] * (-cI * (TMP9 * TMP13) + cI * (TMP8 * TMP15)); 
}


void VVS4_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP15; 
  complex<double> TMP9; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (-cI * (P2[0] * TMP9) + cI * (V2[2] * TMP15)); 
  V1[3] = denom * S3[2] * (-cI * (P2[1] * TMP9) + cI * (V2[3] * TMP15)); 
  V1[4] = denom * S3[2] * (-cI * (P2[2] * TMP9) + cI * (V2[4] * TMP15)); 
  V1[5] = denom * S3[2] * (-cI * (P2[3] * TMP9) + cI * (V2[5] * TMP15)); 
}


void VVS4_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP15; 
  complex<double> TMP9; 
  complex<double> TMP13; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP15 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP8 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP13 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * (-cI * (TMP8 * TMP15) + cI * (TMP9 * TMP13)); 
}


void VVS5_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP15; 
  complex<double> TMP14; 
  complex<double> TMP9; 
  complex<double> TMP13; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP8 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP15 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP14 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP11 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP10 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP13 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP12 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  vertex = COUP * S3[2] * (TMP10 * (-cI * (TMP12 * TMP15) + cI * (TMP9 *
      TMP11)) + TMP14 * (-cI * (TMP8 * TMP11) + cI * (TMP12 * TMP13)));
}


void VVS5_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP15; 
  complex<double> TMP14; 
  complex<double> denom; 
  complex<double> TMP9; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP14 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP11 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP12 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (P3[0] * (-cI * (TMP12 * TMP15) + cI * (TMP9 *
      TMP11)) + TMP14 * (-cI * (V2[2] * TMP11) + cI * (P2[0] * TMP12)));
  V1[3] = denom * S3[2] * (P3[1] * (-cI * (TMP12 * TMP15) + cI * (TMP9 *
      TMP11)) + TMP14 * (-cI * (V2[3] * TMP11) + cI * (P2[1] * TMP12)));
  V1[4] = denom * S3[2] * (P3[2] * (-cI * (TMP12 * TMP15) + cI * (TMP9 *
      TMP11)) + TMP14 * (-cI * (V2[4] * TMP11) + cI * (P2[2] * TMP12)));
  V1[5] = denom * S3[2] * (P3[3] * (-cI * (TMP12 * TMP15) + cI * (TMP9 *
      TMP11)) + TMP14 * (-cI * (V2[5] * TMP11) + cI * (P2[3] * TMP12)));
}


void VVS5_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  complex<double> TMP11; 
  double P1[4]; 
  complex<double> TMP10; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP15; 
  complex<double> TMP14; 
  complex<double> denom; 
  complex<double> TMP9; 
  complex<double> TMP13; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP8 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP15 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP14 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP11 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP10 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP13 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP12 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * (TMP10 * (-cI * (TMP9 * TMP11) + cI * (TMP12 * TMP15)) +
      TMP14 * (-cI * (TMP12 * TMP13) + cI * (TMP8 * TMP11)));
}


void FFS1_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP0; 
  TMP0 = (F2[4] * F1[4] + F2[5] * F1[5] - F2[2] * F1[2] - F2[3] * F1[3]); 
  vertex = COUP * - cI * TMP0 * S3[2]; 
}


void FFS1_1(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) + F2[2] * M1));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) - F2[3] * M1));
  F1[4] = denom * - cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] +
      cI * (P1[2])) - F2[4] * M1));
  F1[5] = denom * cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) + F2[5] * M1));
}


void FFS1_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) - F1[2] * M2));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * -
      1. * (P2[0] + P2[3]) + F1[3] * M2));
  F2[4] = denom * cI * S3[2] * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] * (+cI
      * (P2[2]) - P2[1]) + F1[4] * M2));
  F2[5] = denom * - cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) - F1[5] * M2));
}


void FFS1_2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex)
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> tmp; 
  FFS1_0(F1, F2, S3, COUP1, vertex); 
  FFS2_0(F1, F2, S3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFS1_2_1(complex<double> F2[], complex<double> S3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[])
{
//   complex<double> cI = complex<double> (0., 1.); 
//   double P1[4]; 
  complex<double> denom; 
  int i; 
  complex<double> Ftmp[6]; 
  FFS1_1(F2, S3, COUP1, M1, W1, F1); 
  FFS2_1(F2, S3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFS1_2_2(complex<double> F1[], complex<double> S3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[])
{
//   complex<double> cI = complex<double> (0., 1.); 
  complex<double> Ftmp[6]; 
//   double P2[4]; 
  complex<double> denom; 
  int i; 
  FFS1_2(F1, S3, COUP1, M2, W2, F2); 
  FFS2_2(F1, S3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}


void FFS2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * - cI * TMP1 * S3[2]; 
}


void FFS2_1(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) - F2[2] * M1));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) + F2[3] * M1));
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) + F2[4] * M1));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) - F2[5] * M1));
}


void FFS2_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) + F1[2] * M2));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * -
      1. * (P2[0] + P2[3]) - F1[3] * M2));
  F2[4] = denom * - cI * S3[2] * (F1[2] * - 1. * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) - F1[4] * M2));
  F2[5] = denom * cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) + F1[5] * M2));
}


}  // end namespace $(namespace)s_HEF_MEKD2_1

