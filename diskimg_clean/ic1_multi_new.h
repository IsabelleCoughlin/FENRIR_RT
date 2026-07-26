
/*
 * Global simulation parameters
 *
 * Geometrized units are assumed:
 *
 *     G = c = M_BH = 1.
 *
 * Consequently, lengths and times are measured in units of
 *
 *     r_g = GM_BH/c^2
 *     t_g = GM_BH/c^3.
 *
 */

//Phystical parameters
double a;
double inclination;

//Accretion parameters
int numMdot; //the number of Mdot values being explored
double initMdot; //the first Mdot value being explored
double finalMdot; //the final Mdot being explored
double deltaMdot; //the change in Mdot between the different Mdot values
double Mdot; //the current mass accretion rate

//Position of observer
double obsRadius = 1000.; //radius of observer
double obsPhi = 0.; //phi angle of observer

//Integration variables
double imgSize; //length of image plane in rg
double maxStep = 100000000;  //maximal number of steps
double horizonStop = 0.01;  //minimum radial distance from event horizon before considered lost
double tolerance = 1000.; //Tolerance for automatic step refinement. The larger, the more finely sampled.

//The prefix for the different output files
std::string outFilePrefix;