#include <cmath>

/*
 * First-order null-geodesic equations in Kerr spacetime.
 *
 * Coordinates:
 *
 *     x^mu = (t, r, theta, phi)
 *
 * Units and signature:
 *
 *     G = c = M_BH = 1
 *     (-, +, +, +)
 *
 * Conserved photon quantities:
 *
 *     energy = E = -p_t
 *     angmom = Lz = p_phi
 *     carter = Q
 *
 * Dots represent differentiation with respect to an affine
 * parameter lambda:
 *
 *     tdot   = dt/dlambda
 *     rdot   = dr/dlambda
 *     thdot  = dtheta/dlambda
 *     phidot = dphi/dlambda.
 */


//Kerr function: Sigma = r^2 + a^2 cos^2(theta).
double sigma(double r, double phi, double theta, double time){
    double cosTerm = std::cos(theta);
    double sigmaOut = (r*r)+(a*a*cosTerm*cosTerm);
	return sigmaOut;
	}

//Kerr horizon function: Delta = r^2 - 2r + a^2.
double delta(double r, double phi, double theta, double time){
    double deltaOut = (r*r)+(-2.*r)+(a*a);
	  return deltaOut;
	}

//Time component of the photon tangent vector.
double tdot(double r, double phi, double theta, double time){
    double sinTerm = std::sin(theta);
    double energyTerm = ((sigma(r,phi,theta,time)*((r*r)+(a*a)))+(2.*r*a*a*sinTerm*sinTerm))*energy;
	double angmomTerm = -2.*a*r*angmom;
	double bottomTerm1 = (sigma(r,phi,theta,time) - (2.*r))*((r*r)+(a*a));
    double bottomTerm2 = 2.*r*a*a*sinTerm*sinTerm;
	return (energyTerm+angmomTerm)/(bottomTerm1+bottomTerm2);
	}

// Azimuthal component of the photon tangent vector.
double phidot(double r, double phi, double theta, double time){
    double sinTerm = std::sin(theta);
    double energyTerm = 2.*a*r*sinTerm*sinTerm*energy;
	double angmomTerm = (sigma(r,phi,theta,time) - (2.*r))*angmom;
    double bottomTerm1 = (sigma(r,phi,theta,time) - (2.*r))*((r*r)+(a*a))*sinTerm*sinTerm;
    double bottomTerm2 = 2.*r*a*a*sinTerm*sinTerm*sinTerm*sinTerm;
	return (energyTerm+angmomTerm)/(bottomTerm1+bottomTerm2);
	}

//Squared polar component:
double thdotsq(double r, double phi, double theta, double time){
    double cosTerm = std::cos(theta);
    double cotTerm = 1./std::tan(theta);
    double aTerm = a*a*cosTerm*cosTerm*energy*energy; // Should this be energy**2, except energy is normalized to 1
    double angmomTerm = -1.*angmom*angmom*cotTerm*cotTerm;
    double topTerm = carter+aTerm+angmomTerm;
    double bottomTerm = sigma(r,phi,theta,time)*sigma(r,phi,theta,time);
	return topTerm/bottomTerm;
	}

// Squared radial component.
double rdotsq(double r, double phi, double theta, double time){
	double energyTerm = energy*tdot(r, phi, theta, time);
	double angmomTerm = -1.*angmom*phidot(r, phi, theta, time);
    double carterTerm = -1.*sigma(r,phi,theta,time)*thdotsq(r,phi,theta,time);
	double frontTerm = delta(r, phi, theta, time) / sigma(r, phi, theta, time);
	return (frontTerm*(energyTerm+angmomTerm+carterTerm));
	}
