

#ifndef incl_StandardAtmosphere_h
#define incl_StandardAtmosphere_h


inline void standardAtmosphere(double h, double *rho, double *drhodh = NULL)
{
  // takes altitude h in metres and returns air density in kg/m^3

  double T, p, R, g, p0;

  g  = 9.81;   //   m / s^2
  R  = 287;    // J / (kg K)

  if (h < 11000.)

    p0 = 101325; //  N / m^2

    fact1 = 71.5/(288.2*11000);

    fact2 = g/R*11000/71.5;

    p  = p0 * pow((1.-fact1*h),fact2);
 
    T  = 288.2 - 71.5 * h / 11000;

    *rho = p / (R * T);

    if (drhodh != NULL)

      *drhodh = (p / (R*T) * 71.5/11000
               - p0 * pow((1.-fact1*h),fact2-1.) * fact2 * fact1) / (R * T);
  }
  else if (h < 20000.)
  {
    p0 = 22621; //  N / m^2

    fact1 = g/R*11000/216.7;

    fact2 = - g/R/216.7;

    p  = p0 * exp(fact1 + fact2*h);
 
    T  = 216.7;

    *rho = p / (R * T);

    if (drhodh != NULL) *drhodh = *rho * fact2;
  }
  else
  {
    std::cout << "\n  ERROR: StandardAtmosphere: not yet defined for h > 20,000 m!\n\n";
  }

  return;
}




#endif










