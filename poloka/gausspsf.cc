#include <poloka/gausspsf.h>

static double sqr(const double& x) { return x*x; }

GaussPsf::GaussPsf(const ReducedImage &Rim)
{
  double sigx, sigy, rho;
  Rim.GetPsfShapeParams(sigx, sigy, rho);
  beta = 1 - sqr(rho);
  alphax = -0.5 / sqr(sigx) / beta;
  alphay = -0.5 / sqr(sigy) / beta;
  rxy    = 2.* rho / (sigx * sigy);
  norm   = 1./ (2.* M_PI * sigx * sigy * sqrt(beta));  
}


double GaussPsf::Value(const int i, const int j, const double &Xc, const double &Yc, 
		     double &DpDx, double &DpDy) const
{
  static double dx = i-Xc;
  static double dy = j-Yc;
  static double val =  norm * exp(sqr(dx)*alphax + sqr(dy)*alphay + rxy*dx*dy);
  DpDx = (2.*dx*alphax + dy*rxy) * val;
  DpDy = (2.*dy*alphay + dx*rxy) * val;

  return val;
}

