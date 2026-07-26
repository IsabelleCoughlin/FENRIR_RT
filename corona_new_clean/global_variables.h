#include <cmath>

/*
 * Global ray-tracing state variables.
 *
 * These values are shared between the main tracing loop and helper functions.
 * They describe the current photon, the Kerr geometry, the corona frame, and
 * the disk impact point.
 */

/* Photon conserved quantities in Kerr spacetime. */
double energy;   // E = -p_t
double angmom;   // L_z = p_phi
double carter;   // Carter constant Q


//Defining the spacetime position vector of the photon, setting its components to the initial component
double posVec[4];

//Defining the momentum 4-vector of the photon
double momVec[4];

//Defining final momentum one-form for photon
double momOneForm[4];

//Defining energy of photon that is observed co-moving disk frame
double finalEnergy;

//Defining vector for disk velocity at the final position of the photon in B-L coordinates
double diskVelVec[4];

//Defining Lorentz factor of the disk element of the final position relative to the LNRF that is instantaneously at co-spatial with disk element
double lorentz;


/*
 * Disk-impact geometry.
 *
 * scaleHeightValue : disk surface height z(rho) at the impact point
 * rProjected       : cylindrical/projected radius rho = r sin(theta)
 */
//Defining two distances of interest to be calculated at end of integration
double scaleHeightValue;  //Defining the variable that will hold the scale height where the photon ends up
double rProjected;  //Defining the variable that will hold the radial distance from the final placement of the photon to the black hole, projected onto the plane

//Calculating radius of event horizon
double rEvent;// = 1. + pow(1.-(a*a),0.5);

//Calculating the inner-most stable orbital radius (r_isco)
double z1;// = 1. + (pow(1-a*a,(1./3.))*((pow(1+a,(1./3.)))+(pow(1-a,(1./3.)))));
double z2;// = sqrt((3.*a*a)+(z1*z1));
double rIsco;// = 3.+z2-sqrt((3-z1)*(3+z1+(2.*z2)));

//Efficiency and accretion rate (Mdot/Eddington) to calculate scale height
/*
 * Radiative efficiency eta = 1 - E(r_ISCO).
 *
 * efficiencyUpper and efficiencyLower are intermediate terms in the Kerr
 * circular-orbit energy evaluated at the ISCO.
 */
double efficiencyUpper;// = (rIsco*rIsco)-(2.*rIsco)+(a*sqrt(rIsco));
double efficiencyLower;// = (rIsco*rIsco)-(3.*rIsco)+(2.*a*sqrt(rIsco));
double efficiency;// = 1.-(efficiencyUpper/(rIsco*sqrt(efficiencyLower)));

//Calculating the front constant for scale height calculation
/*
 * Front factor in the disk scale-height equation:
 *
 *   heightFrontTerm = (3 / efficiency) * accretion
 *
 * so that:
 *
 *   z(rho) = heightFrontTerm * [1 - sqrt(r_ISCO / rho)]
 */
double heightFrontTerm;// = 2.*(3./(2.*efficiency))*accretion;

//Position of the corona in B-L coordinates
double initTime,initRadius,initTheta,initPhi;

//Corona rotational velocity (dphi/dt = rotOmega)
double rotOmega;

//Integration variables
double dStep;  //step size; changed sign
double rLimitLow;// = rEvent + horizonStop;  //minimum radial distance (currently r_horizon + 0.1) (changed to 0.5)

double rSqrtSwitch;
double thSqrtSwitch;
double hitDiskSwitch;

//Angles in corona co-moving frame
double imgAlpha,imgBeta;  //Defining the angle variables

//Basis vectors for orthonormal tetrad at corona
double eTvec[4],eRvec[4],eThVec[4],ePhVec[4];

//Value of metric components at corona position
double gTTval,gTPhVal,gRRval,gThThVal,gPhTval,gPhPhVal;

//Value of trig functions at position of corona given
double corCos,corSin,corCot;
