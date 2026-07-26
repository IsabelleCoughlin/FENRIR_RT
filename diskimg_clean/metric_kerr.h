#include <cmath>
/*
 * Covariant Kerr metric in Boyer-Lindquist coordinates:
 *
 *     x^mu = (t, r, theta, phi)
 *
 * using geometrized units:
 *
 *     G = c = M_BH = 1
 *
 * and metric signature:
 *
 *     (-, +, +, +).
 *
 * The helper functions sigma() and delta() must implement:
 *
 *     Sigma = r^2 + a^2 cos^2(theta)
 *     Delta = r^2 - 2r + a^2.
 *
 * The metric is stationary and axisymmetric, so its components do
 * not depend on t or phi. Those arguments are retained for interface
 * consistency with the legacy code.
 */


 /*
 * Time-time component:
 *
 *     g_tt = -(1 - 2r/Sigma).
 */
double gTT(double r, double phi, double theta, double time){
	return -1.*(1.-(2.*r/sigma(r,phi,theta,time)));
}

/*
 * Mixed time-azimuth component:
 *
 *     g_tphi = -2ar sin^2(theta)/Sigma.
 */
double gTPh(double r, double phi, double theta, double time){
    return (-2.*a*r*pow(std::sin(theta),2.))/sigma(r,phi,theta,time);
}

/*
 * Metric symmetry requires:
 *
 *     g_phit = g_tphi.
 */
double gPhT(double r, double phi, double theta, double time){
    return (-2.*a*r*pow(std::sin(theta),2.))/sigma(r,phi,theta,time);
}

/*
 * Radial component:
 *
 *     g_rr = Sigma/Delta.
 */
double gRR(double r, double phi, double theta, double time){
    return sigma(r,phi,theta,time)/delta(r,phi,theta,time);
}

/*
 * Polar component:
 *
 *     g_thetatheta = Sigma.
 */
double gThTh(double r, double phi, double theta, double time){
    return sigma(r,phi,theta,time);
}

/*
 * Azimuthal component:
 *
 *     g_phiphi =
 *       [r^2 + a^2
 *        + 2ra^2 sin^2(theta)/Sigma] sin^2(theta).
 */
double gPhPh(double r, double phi, double theta, double time){
    return (pow(r,2.)+pow(a,2.)+((2.*r*pow(a*std::sin(theta),2.))/sigma(r,phi,theta,time)))*pow(std::sin(theta),2.);
}

/*
 * Lower the index of a contravariant vector:
 *
 *     V_mu = g_mn V^n.
 *
 * Inputs:
 *
 *     posVec = (t, r, theta, phi)
 *     vec    = (V^t, V^r, V^theta, V^phi)
 *
 * Output:
 *
 *     oneForm = (V_t, V_r, V_theta, V_phi).
 *
 * Because the Boyer-Lindquist Kerr metric has a nonzero t-phi
 * cross term:
 *
 *     V_t   = g_tt V^t + g_tphi V^phi
 *     V_phi = g_tphi V^t + g_phiphi V^phi.
 *
 * The radial and polar components have no cross terms.
 */
void vecToOneForm(double posVec[4], double vec[4], double oneForm[4]){
    //Defining local location variables to use
    double localTime = posVec[0];
    double localRadius = posVec[1];
    double localTheta = posVec[2];
    double localPhi = posVec[3];

    //Calculating sin(theta)^2 value
    double sinVal = std::sin(localTheta);
    double sinTerm = sinVal*sinVal;

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
