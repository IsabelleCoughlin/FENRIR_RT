#include <cmath>

/*
 * Disk geometry and velocity helper functions.
 *
 * This file provides:
 *
 *   1. The finite-thickness disk surface height:
 *
 *        z(rho) = heightFrontTerm * [1 - sqrt(r_ISCO / rho)]
 *
 *      where rho = r sin(theta) is the cylindrical/projected radius.
 *      In the Mdot-based version of the code,
 *
 *        heightFrontTerm = (3 / eta) * (Mdot / Mdot_Edd)
 *
 *      so this reproduces the Taylor & Reynolds / Fenrir disk-height
 *      prescription.
 *
 *   2. The disk 4-velocity at the photon impact point.
 *
 *      The disk material is assumed to move on circular azimuthal orbits:
 *
 *        u^r     = 0
 *        u^theta = 0
 *        Omega   = dphi/dt = 1 / (rho^(3/2) + a)
 *
 *      The time component u^t is computed by enforcing the 4-velocity
 *      normalization:
 *
 *        g_mu_nu u^mu u^nu = -1
 *
 *   3. The Lorentz factor of the disk material relative to the local
 *      non-rotating frame, LNRF. The LNRF basis functions below convert the
 *      coordinate-basis disk velocity into the locally measured frame.
 *
 * Assumptions:
 *   - Geometric units are used, G = c = M = 1.
 *   - The disk has no radial or polar velocity.
 *   - The azimuthal velocity uses a Keplerian Kerr-like form evaluated at the
 *     projected cylindrical radius rho, not the spherical radius r.
 *   - rIsco, heightFrontTerm, and a are global variables.
 *   - Metric helper functions such as gTT, gTPh, gPhPh, sigma, delta, and
 *     aFunct are defined elsewhere.
 */



/*
 * Return the vertical height of the disk surface at a given spacetime position.
 *
 * Inputs:
 *   radius : Boyer-Lindquist radial coordinate r
 *   theta  : Boyer-Lindquist polar angle theta
 *
 * The disk surface is written in terms of projected/cylindrical radius:
 *
 *   rho = r sin(theta)
 *
 * For rho < r_ISCO, the disk is treated as absent/truncated, so this returns
 * zero height. For rho >= r_ISCO, the height follows:
 *
 *   z(rho) = heightFrontTerm * [1 - sqrt(r_ISCO / rho)]
 *
 * The returned value is a vertical height z, not a polar angle theta.
 */
double scaleHeightFnct(double radius, double theta){
  double scaleHTprojR = radius*std::sin(theta);
  if (scaleHTprojR < rIsco){
    return 0.;
  }
  else{
    return heightFrontTerm*(1 - sqrt(rIsco/scaleHTprojR));
  }
}

/*
 * Return the disk angular velocity Omega = dphi/dt.
 *
 * The disk is assumed to rotate on circular, Keplerian-like orbits with:
 *
 *   Omega = 1 / (rho^(3/2) + a)
 *
 * where rho is the projected radius and a is the black-hole spin.
 *
 * The extra arguments are kept for interface consistency with the other disk
 * velocity functions, even though this expression only uses rProjected and a.
 */
double diskPhiVel(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
    return 1./(sqrt(rProjected*rProjected*rProjected)+(a));
}


/*
 * Return u^t, the time component of the disk 4-velocity.
 *
 * Given Omega = dphi/dt, the disk 4-velocity has the form:
 *
 *   u^mu = u^t (1, 0, 0, Omega)
 *
 * with u^r = u^theta = 0.
 *
 * u^t is determined by normalizing the 4-velocity:
 *
 *   g_mu_nu u^mu u^nu = -1
 *
 * Substituting the circular-orbit form gives:
 *
 *   u^t = 1 / sqrt[-(g_tt + 2 g_tphi Omega + g_phiphi Omega^2)]
 */
double diskTdot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
    double diskPhiVelValue = diskPhiVel(r,phi,theta,time,scaleHeightValue,rProjected);
    double gTTvalue = gTT(r,phi,theta,time);
    double gTPhValue = gTPh(r,phi,theta,time);
    double gPhPhValue = gPhPh(r,phi,theta,time);
    double diskTdotSqValue = -1.*(gTTvalue+(2.*gTPhValue*diskPhiVelValue)+(gPhPhValue*diskPhiVelValue*diskPhiVelValue));
    return sqrt(1./diskTdotSqValue);
}
/*
 * Return u^phi, the azimuthal component of the disk 4-velocity.
 *
 * Since Omega = dphi/dt = u^phi / u^t:
 *
 *   u^phi = Omega * u^t
 */
double diskPhiDot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
    double diskPhiVelValue = diskPhiVel(r,phi,theta,time,scaleHeightValue,rProjected);
    double diskTdotValue = diskTdot(r,phi,theta,time,scaleHeightValue,rProjected);
    return diskTdotValue*diskPhiVelValue;
}

/*
 * Return u^theta.
 *
 * The disk material is assumed to move only azimuthally, with no polar motion.
 */
double diskThDot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
	return 0.;
	}


    /*
 * Return u^r.
 *
 * The disk material is assumed to move only azimuthally, with no radial inflow
 * in this velocity prescription.
 */
double diskRdot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
	return 0.;
	}


/*
 * Fill diskVelVec with the disk 4-velocity at the photon impact position.
 *
 * posVec ordering:
 *   posVec[0] = t
 *   posVec[1] = r
 *   posVec[2] = theta
 *   posVec[3] = phi
 *
 * diskVelVec ordering:
 *   diskVelVec[0] = u^t
 *   diskVelVec[1] = u^r
 *   diskVelVec[2] = u^theta
 *   diskVelVec[3] = u^phi
 */
//Calculate the disk velocity at the position that the photon is intercepted
void diskVelocity(double posVec[4], double diskVelVec[4], double scaleHeightValue, double rProjected){
	diskVelVec[0] = diskTdot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //time
    diskVelVec[1] = diskRdot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //radius
    diskVelVec[2] = diskThDot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //theta
    diskVelVec[3] = diskPhiDot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //phi
    }

/*
 * LNRF lapse-like basis factor e^nu.
 *
 * These functions define the local non-rotating frame, or LNRF, basis used to
 * convert coordinate-basis velocity components into locally measured
 * orthonormal-frame components.
 */
double eNu(double r, double phi, double theta, double time){
	double eNuOut = (sigma(r,phi,theta,time)*delta(r,phi,theta,time))/aFunct(r,phi,theta,time);
	eNuOut = sqrt(eNuOut);
	return eNuOut;
}
/*
 * LNRF azimuthal basis factor e^psi.
 */
double ePsi(double r, double phi, double theta, double time){
	double sinVal = std::sin(theta);
	double ePsiOut = aFunct(r,phi,theta,time)/sigma(r,phi,theta,time);
	ePsiOut = sqrt((sinVal*sinVal)*ePsiOut);
	return ePsiOut;
}

/*
 * LNRF radial basis factor e^mu_1.
 */
double eMu1(double r, double phi, double theta, double time){
	return sqrt(sigma(r,phi,theta,time)/delta(r,phi,theta,time));
}

/*
 * LNRF polar basis factor e^mu_2.
 */
double eMu2(double r, double phi, double theta, double time){
	return sqrt(sigma(r,phi,theta,time));
}

/*
 * Frame-dragging angular velocity of the LNRF.
 *
 * This is the angular velocity omega of a zero-angular-momentum observer:
 *
 *   omega = -g_tphi / g_phiphi
 *
 * For the Kerr metric helper functions used here, this reduces to:
 *
 *   omega = 2 a r / A
 *
 * where A is returned by aFunct().
 */
double lnrfOmega(double r, double phi, double theta, double time){
	return (2.*a*r)/aFunct(r,phi,theta,time);
}

/*
 * Calculate the Lorentz factor gamma of the disk material as measured by the
 * co-spatial LNRF.
 *
 * Steps:
 *   1. Convert the coordinate-basis disk 4-velocity into the LNRF basis.
 *   2. Compute the local 3-speed squared:
 *
 *        v^2 = [(u^(r))^2 + (u^(theta))^2 + (u^(phi))^2] / (u^(t))^2
 *
 *   3. Return:
 *
 *        gamma = 1 / sqrt(1 - v^2)
 *
 * The phi component includes the subtraction of the LNRF frame-dragging
 * rotation, so it measures the disk's azimuthal motion relative to the local
 * non-rotating observer rather than relative to coordinates.
 */
double findLorentz(double posVec[4], double diskVelVec[4]){
	double lorentzVelVec[4];
	double lorentzVelSq;
	lorentzVelVec[0] = diskVelVec[0]*eNu(posVec[1],posVec[3],posVec[2],posVec[0]);
    lorentzVelVec[1] = diskVelVec[1]*eMu1(posVec[1],posVec[3],posVec[2],posVec[0]);
    lorentzVelVec[2] = diskVelVec[2]*eMu2(posVec[1],posVec[3],posVec[2],posVec[0]);
    lorentzVelVec[3] = (-1.*lnrfOmega(posVec[1],posVec[3],posVec[2],posVec[0])*ePsi(posVec[1],posVec[3],posVec[2],posVec[0])*diskVelVec[0])+(ePsi(posVec[1],posVec[3],posVec[2],posVec[0])*diskVelVec[3]);
    lorentzVelSq = ((lorentzVelVec[1]*lorentzVelVec[1]) + (lorentzVelVec[2]*lorentzVelVec[2]) + (lorentzVelVec[3]*lorentzVelVec[3]))/(lorentzVelVec[0]*lorentzVelVec[0]);
    return 1./sqrt(1.-lorentzVelSq);
}
