#include <cmath>
/*
 * Disk model and coordinate conventions
 * -------------------------------------
 *
 * Geometrized units are assumed: G = c = M_BH = 1.
 *
 * Boyer-Lindquist coordinates are ordered as: x^mu = (t, r, theta, phi).
 *
 * The pseudo-cylindrical radius is rho = r sin(theta).
 *
 * The upper disk photosphere is modeled as  z(rho) = heightFrontTerm [1 - sqrt(r_ISCO / rho)].
 *
 * For the original Fenrir prescription, heightFrontTerm = 3 mdot / eta,
 *
 * because the X-ray photosphere is placed at z = 2H, where
 *
 *     H = (3/2)(mdot/eta)
 *         [1 - sqrt(r_ISCO/rho)].
 *
 * This is a prescribed surface geometry rather than a dynamically
 * self-consistent thick-disk solution.
 */


// Return the height z of the upper disk photosphere above the equatorial plane.
double scaleHeightFnct(double radius, double theta){
  double scaleHTprojR = radius*std::sin(theta);
  if (scaleHTprojR < rIsco){
    return 0.;
  }
  else{
    return heightFrontTerm*(1 - sqrt(rIsco/scaleHTprojR));
  }
}

// Coordinate angular velocity Omega = dphi/dt.
double diskPhiVel(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
    return 1./(sqrt(rProjected*rProjected*rProjected)+(a));
}

//Calculate the time component u^t of the disk four-velocity.
double diskTdot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
    double diskPhiVelValue = diskPhiVel(r,phi,theta,time,scaleHeightValue,rProjected);
    double gTTvalue = gTT(r,phi,theta,time);
    double gTPhValue = gTPh(r,phi,theta,time);
    double gPhPhValue = gPhPh(r,phi,theta,time);
    double diskTdotSqValue = -1.*(gTTvalue+(2.*gTPhValue*diskPhiVelValue)+(gPhPhValue*diskPhiVelValue*diskPhiVelValue));
    return sqrt(1./diskTdotSqValue);
}

// Calculate the azimuthal four-velocity component:
double diskPhiDot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
    double diskPhiVelValue = diskPhiVel(r,phi,theta,time,scaleHeightValue,rProjected);
    double diskTdotValue = diskTdot(r,phi,theta,time,scaleHeightValue,rProjected);
    return diskTdotValue*diskPhiVelValue;
}


/*
 * The prescribed disk has no motion in the polar direction.
 *
 * This is an idealized stationary-surface approximation, not a
 * complete accretion-flow velocity field.
 */
double diskThDot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
	return 0.;
	}

/*
 * The prescribed disk has no radial inflow.
 *
 * A physical accretion disk has a small inward radial velocity, but
 * it is neglected in this ray-tracing model.
 */
double diskRdot(double r, double phi, double theta, double time, double scaleHeightValue, double rProjected){
	return 0.;
	}

/*
 * Construct the contravariant disk four-velocity
 *
 *     u^mu = (u^t, u^r, u^theta, u^phi)
 *
 * at the photon/disk intersection.
 *
 * posVec uses the ordering:
 *
 *     [t, r, theta, phi].
 */
void diskVelocity(double posVec[4], double diskVelVec[4], double scaleHeightValue, double rProjected){
	diskVelVec[0] = diskTdot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //time
    diskVelVec[1] = diskRdot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //radius
    diskVelVec[2] = diskThDot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //theta
    diskVelVec[3] = diskPhiDot(posVec[1],posVec[3],posVec[2],posVec[0],scaleHeightValue,rProjected); //phi
    }
