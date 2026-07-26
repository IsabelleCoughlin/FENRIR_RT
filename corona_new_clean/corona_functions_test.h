#include <cmath>

/*
 * Corona setup and photon initialization helpers.
 *
 * These functions place the corona, choose initial photon directions, build the
 * local orthonormal tetrad comoving with the corona, convert photon momentum
 * from the corona rest frame into Boyer-Lindquist coordinates, and compute the
 * photon constants of motion.
 *
 * Most routines use globals defined elsewhere: spin a, corona geometry settings,
 * metric functions, tetrad arrays, momentum arrays, and conserved quantities.
 */


//Functions used in calculating alpha and beta
/* Draw a random polar emission angle alpha for isotropic emission. */
double randAlpha(){
	double randVal = (double)rand()/RAND_MAX;
	return std::acos(-1. + (2.*randVal));
}
/* Draw a random azimuthal emission angle beta. */
double randBeta(){
	double randVal = (double)rand()/RAND_MAX;
	return 2.*M_PI*randVal;
}

/*
 * Set the initial corona position in Boyer-Lindquist coordinates.
 *
 * For a lamppost corona, place the source very close to the polar axis at
 * radius = height.
 *
 * For an off-axis corona, place it at cylindrical radius initCylRadius and at
 * heightAboveDisk above the finite-thickness disk surface.
 */
void findPosition(){
	if (coronaType == "lp"){
		initTime = 0.;
		initRadius = height;
		initTheta = (M_PI/180.)*0.01;
		initPhi = 0.;
	}else{
		if (isRisco == true){
			initCylRadius = initCylRadius*rIsco;
		}
		double thicknessAtCorona = heightFrontTerm*(1. - sqrt(rIsco/initCylRadius));
		initTime = 0.;
		initTheta = std::atan(initCylRadius/(thicknessAtCorona+heightAboveDisk));
		initRadius = initCylRadius/std::sin(initTheta);
		initPhi = coronaPhi;
	}
}

/*
 * Set the corona angular velocity Omega = dphi/dt.
 *
 * Options:
 *   lp                  : rotate with the LNRF/frame-dragging angular velocity
 *   offaxis + lnrf      : same LNRF choice
 *   offaxis + kep       : Keplerian angular velocity at initCylRadius
 *   otherwise           : use user-supplied rotOmegaInput
 */
void findOmega(){
	if ((coronaType == "lp") or ((coronaType == "offaxis") and (rotOmegaType == "lnrf"))){
		rotOmega = 2.*a*initRadius/aFunct(initRadius, initPhi, initTheta, initTime);
	}else if ((coronaType == "offaxis") and (rotOmegaType == "kep")){
		rotOmega = 1./(pow(initCylRadius,1.5) + a);
	}else{
		rotOmega = rotOmegaInput;
	}
}


/*
 * Time-basis component of the corona tetrad.
 * The tetrad is orthonormal and comoving with a source rotating at rotOmega.
 */

double eTT(double r, double phi, double theta, double time, double rotOmega){
	double bottomTerm1,bottomTerm2,bottomTerm3,eTTout;
	bottomTerm1 = gTT(r,phi,theta,time);
	bottomTerm2 = (gTPh(r,phi,theta,time) + gPhT(r,phi,theta,time))*rotOmega;
	bottomTerm3 = gPhPh(r,phi,theta,time)*rotOmega*rotOmega;
	eTTout = -1.*(bottomTerm1+bottomTerm2+bottomTerm3);
	eTTout = 1./sqrt(eTTout);
	return eTTout;
}


/* Phi component of the time-like tetrad basis vector. */
double eTPh(double r, double phi, double theta, double time, double rotOmega){
	return rotOmega*eTT(r,phi,theta,time,rotOmega);
}

/* Radial tetrad basis normalization. */
double eRR(double r, double phi, double theta, double time, double rotOmega){
	return 1./sqrt(gRR(r,phi,theta,time));
}

/* Polar tetrad basis normalization. */
double eThTh(double r, double phi, double theta, double time, double rotOmega){
	return 1./sqrt(gThTh(r,phi,theta,time));
}

/*
 * Coupling term used to construct the azimuthal tetrad basis vector while
 * keeping it orthogonal to the time-like basis vector.
 */
double cFunct(double r, double phi, double theta, double time, double rotOmega){
	double topTerm, bottomTerm;
	topTerm = -1.*(gPhT(r,phi,theta,time) + gPhPh(r,phi,theta,time)*rotOmega);
	bottomTerm = (gTT(r,phi,theta,time) + gTPh(r,phi,theta,time)*rotOmega);
	return topTerm/bottomTerm;
}

/* Phi component of the azimuthal tetrad basis vector. */
double ePhPh(double r, double phi, double theta, double time, double rotOmega){
	double cFunctVal = cFunct(r,phi,theta,time,rotOmega);
	double bottomTerm1,bottomTerm2,bottomTerm3,ePhPhOut;
	bottomTerm1 = gPhPh(r,phi,theta,time);
	bottomTerm2 = (gTPh(r,phi,theta,time) + gPhT(r,phi,theta,time))*cFunctVal;
	bottomTerm3 = gTT(r,phi,theta,time)*cFunctVal*cFunctVal;
	ePhPhOut = bottomTerm1 + bottomTerm2 + bottomTerm3;
	return 1./sqrt(ePhPhOut);
}

/* Time component of the azimuthal tetrad basis vector. */
double ePhT(double r, double phi, double theta, double time, double rotOmega){
	return cFunct(r,phi,theta,time,rotOmega)*ePhPh(r,phi,theta,time,rotOmega);
}

/* Photon energy in the corona rest frame. Arbitrarily normalized to 1. */
double restE(void){
	double restEout;
	restEout = 1.;
	return restEout;
}


/* Rest-frame time component of photon momentum. */
double restPt(void){
	double restPtOut;
	restPtOut = restE();
	return restPtOut;
}

/* Rest-frame radial momentum component. */
double restPr(double alpha, double beta){
	double restPrOut;
	restPrOut = restE()*std::cos(alpha);
	return restPrOut;
}


/* Rest-frame polar momentum component. */
double restPth(double alpha, double beta){
	double restPthOut;
	restPthOut = restE()*std::sin(alpha)*std::cos(beta);
	return restPthOut;
}

/* Rest-frame azimuthal momentum component. */
double restPph(double alpha, double beta){
	double restPphOut;
	restPphOut = restE()*std::sin(alpha)*std::sin(beta);
	return restPphOut;
}

/* Transform rest-frame photon momentum into B-L t component. */
double pT(double restMomVec[4],double eTvec[4], double eRvec[4], double eThVec[4], double ePhVec[4]){
	return (restMomVec[0]*eTvec[0]) + (restMomVec[3]*ePhVec[0]);
}
	

/* Transform rest-frame photon momentum into B-L r component. */
double pR(double restMomVec[4],double eTvec[4], double eRvec[4], double eThVec[4], double ePhVec[4]){
	return (restMomVec[1]*eRvec[1]);
}

/* Transform rest-frame photon momentum into B-L theta component. */
double pTh(double restMomVec[4],double eTvec[4], double eRvec[4], double eThVec[4], double ePhVec[4]){
	return (restMomVec[2]*eThVec[2]);
}

/* Transform rest-frame photon momentum into B-L phi component. */
double pPh(double restMomVec[4],double eTvec[4], double eRvec[4], double eThVec[4], double ePhVec[4]){
	return (restMomVec[0]*eTvec[3]) + (restMomVec[3]*ePhVec[3]);
}

/* Compute and store the corona tetrad vectors at the initial position. */
void findTetrad(){
	double eTTval,eTPhVal,eRRval,eThThVal,ePhTval,ePhPhVal;
	
	eTTval = eTT(initRadius,initPhi,initTheta,initTime,rotOmega);
	eTPhVal = eTPh(initRadius,initPhi,initTheta,initTime,rotOmega);
	eRRval = eRR(initRadius,initPhi,initTheta,initTime,rotOmega);
	eThThVal = eThTh(initRadius,initPhi,initTheta,initTime,rotOmega);
	ePhTval = ePhT(initRadius,initPhi,initTheta,initTime,rotOmega);
	ePhPhVal = ePhPh(initRadius,initPhi,initTheta,initTime,rotOmega);

    //Time basis vector
    eTvec[0] = eTTval;
    eTvec[1] = 0.;
    eTvec[2] = 0.;
    eTvec[3] = eTPhVal;
    //Radial basis vector
    eRvec[0] = 0.;
    eRvec[1] = eRRval;
    eRvec[2] = 0.;
    eRvec[3] = 0.;
    //Theta basis vector
    eThVec[0] = 0.;
    eThVec[1] = 0.;
    eThVec[2] = eThThVal;
    eThVec[3] = 0.;
    //Phi basis vector
    ePhVec[0] = ePhTval;
    ePhVec[1] = 0.;
    ePhVec[2] = 0.;
    ePhVec[3] = ePhPhVal;
    
}

/* Store metric components at the corona position for later reuse. */
void findComponents(){
    gTTval = gTT(initRadius,initPhi,initTheta,initTime);
    gTPhVal = gTPh(initRadius,initPhi,initTheta,initTime);
    gRRval = gRR(initRadius,initPhi,initTheta,initTime);
    gThThVal = gThTh(initRadius,initPhi,initTheta,initTime);
    gPhTval = gPhT(initRadius,initPhi,initTheta,initTime);
    gPhPhVal = gPhPh(initRadius,initPhi,initTheta,initTime);
}

/*
 * Calculate the initial photon momentum in B-L coordinates.
 *
 * The rest-frame photon direction is set by imgAlpha and imgBeta. The rest-frame
 * energy normalization is chosen so that the conserved photon energy at infinity
 * is normalized consistently through the tetrad transformation.
 */
void findMomentum(){
	double restPvec[4];
    restPvec[0] = -1./((gTTval*eTvec[0] + gTPhVal*eTvec[3]) + ((gTTval*ePhVec[0] + gTPhVal*ePhVec[3])*std::sin(imgAlpha)*std::sin(imgBeta)));//restPt();
    restPvec[1] = restPvec[0]*restPr(imgAlpha,imgBeta);
    restPvec[2] = restPvec[0]*restPth(imgAlpha,imgBeta);
    restPvec[3] = restPvec[0]*restPph(imgAlpha,imgBeta);
      
    //Calculating initial B-L photon momentum vector
    momVec[0] = pT(restPvec,eTvec,eRvec,eThVec,ePhVec);
    momVec[1] = pR(restPvec,eTvec,eRvec,eThVec,ePhVec);
    momVec[2] = pTh(restPvec,eTvec,eRvec,eThVec,ePhVec);
    momVec[3] = pPh(restPvec,eTvec,eRvec,eThVec,ePhVec);
}

/*
 * Compute conserved photon quantities in Kerr:
 *
 *   energy = -p_t
 *   angmom = p_phi
 *   carter = Carter constant Q
 */
void findConserved(){
	energy = -1.*((gTTval*momVec[0]) + (gTPhVal*momVec[3]));
    angmom = (gPhTval*momVec[0]) + (gPhPhVal*momVec[3]);
    carter = ((gThThVal*momVec[2])*(gThThVal*momVec[2])) - (corCos*corCos*a*a*energy*energy) + (angmom*angmom*corCot*corCot);
}
      