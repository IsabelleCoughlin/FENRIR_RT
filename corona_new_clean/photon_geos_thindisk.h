#include <cmath>

/*
 * Photon geodesic helper functions in Kerr spacetime.
 *
 * These functions compute derivatives of the photon coordinates with respect to
 * the affine parameter, lambda, along the photon geodesic :
 *
 *   dt/dlambda, dphi/dlambda, (dtheta/dlambda)^2, (dr/dlambda)^2
 *
 * They use the global photon constants of motion:
 *
 *   energy = E
 *   angmom = L_z
 *   carter = Q
 *
 * and the Kerr metric helper functions sigma() and delta().
 */


/* Return dt/dlambda for a photon with conserved E and L_z. */
double tdot(double r, double phi, double theta, double time){
    double sinTerm = std::sin(theta);
    double energyTerm = ((sigma(r,phi,theta,time)*((r*r)+(a*a)))+(2.*r*a*a*sinTerm*sinTerm))*energy;
	double angmomTerm = -2.*a*r*angmom;
	double bottomTerm1 = (sigma(r,phi,theta,time) - (2.*r))*((r*r)+(a*a));
    double bottomTerm2 = 2.*r*a*a*sinTerm*sinTerm;
	return (energyTerm+angmomTerm)/(bottomTerm1+bottomTerm2);
	}


/* Return dphi/dlambda for a photon with conserved E and L_z. */
double phidot(double r, double phi, double theta, double time){
    double sinTerm = std::sin(theta);
    double energyTerm = 2.*a*r*sinTerm*sinTerm*energy;
	double angmomTerm = (sigma(r,phi,theta,time) - (2.*r))*angmom;
    double bottomTerm1 = (sigma(r,phi,theta,time) - (2.*r))*((r*r)+(a*a))*sinTerm*sinTerm;
    double bottomTerm2 = 2.*r*a*a*sinTerm*sinTerm*sinTerm*sinTerm;
	return (energyTerm+angmomTerm)/(bottomTerm1+bottomTerm2);
	}


/*
 * Return (dtheta/dlambda)^2.
 *
 * The sign of dtheta/dlambda is handled elsewhere using thSqrtSwitch. This
 * function only returns the squared magnitude from the theta potential.
 */
double thdotsq(double r, double phi, double theta, double time){
    double cosTerm = std::cos(theta);
    double cotTerm = 1./std::tan(theta);
    double aTerm = a*a*cosTerm*cosTerm*energy;
    double angmomTerm = -1.*angmom*angmom*cotTerm*cotTerm;
    double topTerm = carter+aTerm+angmomTerm;
    double bottomTerm = sigma(r,phi,theta,time)*sigma(r,phi,theta,time);
	return topTerm/bottomTerm;
	}


/*
 * Return (dr/dlambda)^2.
 *
 * The sign of dr/dlambda is handled elsewhere using rSqrtSwitch. This function
 * only returns the squared magnitude from the radial potential.
 */
double rdotsq(double r, double phi, double theta, double time){
	double energyTerm = energy*tdot(r, phi, theta, time);
	double angmomTerm = -1.*angmom*phidot(r, phi, theta, time);
    double carterTerm = -1.*sigma(r,phi,theta,time)*thdotsq(r,phi,theta,time);
	double frontTerm = delta(r, phi, theta, time) / sigma(r, phi, theta, time);
	return (frontTerm*(energyTerm+angmomTerm+carterTerm));
	}