#include <cmath>

/*
 * Kerr metric helper functions in Boyer-Lindquist coordinates.
 *
 * Coordinates are ordered as:
 *   x^mu = (t, r, theta, phi)
 *
 * Units:
 *   G = c = M = 1
 *
 * The spin parameter a is a global variable defined elsewhere.
 * The metric is independent of phi and time, but those arguments are kept for
 * a consistent function interface throughout the code.
 */

 /* Sigma = r^2 + a^2 cos^2(theta). */
double sigma(double r, double phi, double theta, double time){
    double cosTerm = std::cos(theta);
    double sigmaOut = (r*r)+(a*a*cosTerm*cosTerm);
	return sigmaOut;
	}
	
/* Delta = r^2 - 2r + a^2. */
double delta(double r, double phi, double theta, double time){
    double deltaOut = (r*r)+(-2.*r)+(a*a);
	return deltaOut;
	}


    /*
 * Kerr A function:
 *
 *   A = (r^2 + a^2)^2 - Delta a^2 sin^2(theta)
 *
 * Used in LNRF/tetrad expressions.
 */
double aFunct(double r, double phi, double theta, double time){
    double sinTerm = std::sin(theta);
    double deltaTerm = delta(r,phi,theta,time)*a*a*sinTerm*sinTerm;
    double aFunctOut = (((r*r)+(a*a))*(((r*r)+(a*a))))-deltaTerm;
	return aFunctOut;
	}

/*
 * Energy-shift/helper factor used by some photon initialization routines.
 *
 * Depends on the photon emission angles imgAlpha and imgBeta in the local frame.
 *
 */
double epsFunct(double r, double phi, double theta, double time, double imgAlpha, double imgBeta){
    double sigmaTerm = sqrt(sigma(r,phi,theta,time));
    double deltaTerm = sqrt(delta(r,phi,theta,time));
    double aFunctTerm = sqrt(aFunct(r,phi,theta,time));
    double epsFunctTerm = ((sigmaTerm*deltaTerm)+((2.*a*r/sigmaTerm)*std::sin(theta)*std::sin(imgAlpha)*std::sin(imgBeta)))/aFunctTerm;
	return epsFunctTerm;
	}

/* g_tt component of Kerr metric. */
double gTT(double r, double phi, double theta, double time){
	return -1.*(1.-(2.*r/sigma(r,phi,theta,time)));
}

/* g_tphi component of Kerr metric. */
double gTPh(double r, double phi, double theta, double time){
    return (-2.*a*r*pow(std::sin(theta),2.))/sigma(r,phi,theta,time);
}

/* g_phit component. Equal to g_tphi because the metric is symmetric. */
double gPhT(double r, double phi, double theta, double time){
    return (-2.*a*r*pow(std::sin(theta),2.))/sigma(r,phi,theta,time);
}

/* g_rr component of Kerr metric. */
double gRR(double r, double phi, double theta, double time){
    return sigma(r,phi,theta,time)/delta(r,phi,theta,time);
}

/* g_thetatheta component of Kerr metric. */
double gThTh(double r, double phi, double theta, double time){
    return sigma(r,phi,theta,time);
}

/* g_phiphi component of Kerr metric. */
double gPhPh(double r, double phi, double theta, double time){
    return (pow(r,2.)+pow(a,2.)+((2.*r*pow(a*std::sin(theta),2.))/sigma(r,phi,theta,time)))*pow(std::sin(theta),2.);
}

/*
 * Lower the index of a contravariant vector using the Kerr metric:
 *
 *   oneForm_mu = g_mu_nu vec^nu
 *
 * Input:
 *   posVec  = (t, r, theta, phi)
 *   vec     = contravariant vector components
 *
 * Output:
 *   oneForm = covariant vector components
 *
 * Because Kerr has only one off-diagonal term, g_tphi = g_phit, the lowering is:
 *
 *   p_t     = g_tt p^t + g_tphi p^phi
 *   p_r     = g_rr p^r
 *   p_theta = g_thetatheta p^theta
 *   p_phi   = g_phit p^t + g_phiphi p^phi
 */
void vecToOneForm(double posVec[4], double vec[4], double oneForm[4]){
    //Defining local location variables to use
    double localTime = posVec[0];
    double localRadius = posVec[1];
    double localTheta = posVec[2];
    double localPhi = posVec[3];

    //unused?
    //Calculating sin(theta)^2 value
    //double sinVal = std::sin(localTheta);
    //double sinTerm = sinVal*sinVal;

    //Calculating Kerr components
    double gTTval = gTT(localRadius,localPhi,localTheta,localTime);
    double gPhTval = gPhT(localRadius,localPhi,localTheta,localTime);
    double gTPhVal = gPhTval;
    double gRRval = gRR(localRadius,localPhi,localTheta,localTime);
    double gThThVal = gThTh(localRadius,localPhi,localTheta,localTime);
    double gPhPhVal = gPhPh(localRadius,localPhi,localTheta,localTime);

    //Performing conversion
    oneForm[0] = gTTval*vec[0] + gTPhVal*vec[3];
    oneForm[1] = gRRval*vec[1];
    oneForm[2] = gThThVal*vec[2];
    oneForm[3] = gPhTval*vec[0] + gPhPhVal*vec[3];

}