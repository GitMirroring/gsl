/* const/test.c
 * 
 * Copyright (C) 2003, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_test.h>

#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  gsl_ieee_env_setup ();

  /* Basic check to make sure the header files are functioning */

  {
    double c = GSL_CONST_MKS_SPEED_OF_LIGHT;
    double eps = GSL_CONST_MKS_VACUUM_PERMITTIVITY;
    double mu = GSL_CONST_MKS_VACUUM_PERMEABILITY;

    gsl_test_rel (c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mks)");
  }

  {
    double ly = GSL_CONST_CGS_LIGHT_YEAR;
    double c = GSL_CONST_CGS_SPEED_OF_LIGHT;
    double y = 365.25 * GSL_CONST_CGS_DAY;
    
    gsl_test_rel (ly, c * y, 1e-15, "light year (cgs)");
  }

  /* au = 149 597 870 700 m exactly, pc = 648000/pi au (bug #37408) */
  {
    gsl_test_rel (GSL_CONST_MKSA_ASTRONOMICAL_UNIT, 149597870700.0, 1e-16,
                  "astronomical unit (mksa)");
    gsl_test_rel (GSL_CONST_CGS_ASTRONOMICAL_UNIT, 14959787070000.0, 1e-16,
                  "astronomical unit (cgs)");
    gsl_test_rel (GSL_CONST_MKSA_PARSEC,
                  648000.0 / M_PI * GSL_CONST_MKSA_ASTRONOMICAL_UNIT, 1e-15,
                  "parsec (mksa)");
    gsl_test_rel (GSL_CONST_CGSM_PARSEC,
                  648000.0 / M_PI * GSL_CONST_CGSM_ASTRONOMICAL_UNIT, 1e-15,
                  "parsec (cgsm)");
  }

  {
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double eps = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
    double mu = GSL_CONST_MKSA_VACUUM_PERMEABILITY;

    gsl_test_rel (c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mksa)");
  }

  {
    double ly = GSL_CONST_CGSM_LIGHT_YEAR;
    double c = GSL_CONST_CGSM_SPEED_OF_LIGHT;
    double y = 365.25 * GSL_CONST_CGSM_DAY;
    
    gsl_test_rel (ly, c * y, 1e-15, "light year (cgsm)");
  }

  {
    double micro = GSL_CONST_NUM_MICRO;
    double mega = GSL_CONST_NUM_MEGA;
    double kilo = GSL_CONST_NUM_KILO;

    gsl_test_rel (mega/kilo, 1/(micro*kilo), 1e-10, "kilo (mega/kilo, 1/(micro*kilo))");
  }

  {
    double d = GSL_CONST_MKSA_DEBYE;
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double desu = d * c * 1000.0;
    
    gsl_test_rel (desu, 1e-18, 1e-10, "debye (esu)");
  }

  {
    double k = GSL_CONST_MKSA_BOLTZMANN;
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
    double s = 2 * pow(M_PI, 5.0) * pow(k, 4.0) / (15 * pow(c, 2.0) * pow(h, 3.0));
    double sigma = GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT;
    
    gsl_test_rel(s, sigma, 1e-10, "stefan boltzmann constant");
  }

  /* CODATA 2022: the 2019 SI defining constants are exact, and so are
     the constants derived only from them (bug #60635) */
  {
    double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
    double e = GSL_CONST_MKSA_ELECTRON_CHARGE;
    double k = GSL_CONST_MKSA_BOLTZMANN;
    double NA = GSL_CONST_NUM_AVOGADRO;

    gsl_test_rel (h, 6.62607015e-34, 1e-16, "planck constant (SI 2019)");
    gsl_test_rel (e, 1.602176634e-19, 1e-16, "electron charge (SI 2019)");
    gsl_test_rel (k, 1.380649e-23, 1e-16, "boltzmann constant (SI 2019)");
    gsl_test_rel (NA, 6.02214076e23, 1e-16, "avogadro constant (SI 2019)");

    gsl_test_rel (GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR, h / (2.0 * M_PI), 1e-15,
                  "hbar = h / 2 pi");
    gsl_test_rel (GSL_CONST_MKSA_ELECTRON_VOLT, e, 1e-16, "electron volt = e");
    gsl_test_rel (GSL_CONST_MKSA_FARADAY, NA * e, 1e-15, "faraday = N_A e");
    gsl_test_rel (GSL_CONST_MKSA_MOLAR_GAS, NA * k, 1e-15, "molar gas = N_A k");
    gsl_test_rel (GSL_CONST_MKSA_STANDARD_GAS_VOLUME, NA * k * 273.15 / 1e5, 1e-15,
                  "standard gas volume = R (273.15 K) / (100 kPa)");
    gsl_test_rel (GSL_CONST_CGSM_FARADAY, 0.1 * NA * e, 1e-15, "faraday (cgsm)");
  }

  /* measured CODATA 2022 values, and relations between them that hold
     to the published precision */
  {
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
    double me = GSL_CONST_MKSA_MASS_ELECTRON;
    double alpha = GSL_CONST_NUM_FINE_STRUCTURE;
    double mu0 = GSL_CONST_MKSA_VACUUM_PERMEABILITY;
    double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

    gsl_test_rel (GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT, 6.67430e-11, 1e-16,
                  "gravitational constant (CODATA 2022)");
    gsl_test_rel (alpha, 7.2973525643e-3, 1e-16, "fine structure constant (CODATA 2022)");
    gsl_test_rel (mu0 * eps0 * c * c, 1.0, 1e-10, "mu0 eps0 c^2 = 1");
    gsl_test_rel (mu0, 2.0 * alpha * h / (GSL_CONST_MKSA_ELECTRON_CHARGE *
                                          GSL_CONST_MKSA_ELECTRON_CHARGE * c),
                  1e-9, "mu0 = 2 alpha h / (e^2 c)");
    gsl_test_rel (GSL_CONST_MKSA_RYDBERG, alpha * alpha * me * c * c / 2.0, 1e-9,
                  "rydberg = alpha^2 m_e c^2 / 2");
    gsl_test_rel (GSL_CONST_MKSA_BOHR_RADIUS,
                  GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR / (alpha * me * c), 1e-9,
                  "bohr radius = hbar / (alpha m_e c)");
    gsl_test_rel (GSL_CONST_MKSA_SOLAR_MASS * GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT,
                  1.3271244e20, 1e-15, "solar mass = IAU nominal GM_sun / G");
  }


  exit (gsl_test_summary ());
}

