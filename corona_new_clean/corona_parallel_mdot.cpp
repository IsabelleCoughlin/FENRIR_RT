#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include "ic1_parallel.h"
#include "metric_kerr.h"
#include "global_variables.h"
#include "corona_functions_test.h"
#include "photon_geos_thindisk.h"
#include "disk_equations.h"
#include "propagate_rk4_thindisk.h"

/*
 * Ray-traces photons from a lamppost corona to a finite-thickness accretion disk in the Kerr metric.
 *
 * Goal:
 *   For a given black-hole spin and coronal height, this program traces photons
 *   emitted from the corona and records where/how they intersect disk surfaces
 *   corresponding to a grid of accretion rates. The disk height is computed using
 *   the Taylor & Reynolds / Fenrir prescription:
 *
 *       z(rho) = (3 / eta) * (Mdot / Mdot_Edd)
 *                * [1 - sqrt(r_ISCO / rho)]
 *
 *   where eta = 1 - E(r_ISCO) is the radiative efficiency (defined in Bardeen, J. M., Press, W. H., & Teukolsky, S. A. 1972, ApJ,178, 347), 
 *   r_ISCO is the Kerr ISCO radius, and rho = r sin(theta) is the cylindrical radius.
 *
 * Inputs from command line:
 *   argv[1]  n              Total number of photons to trace
 *   argv[2]  totProcs       Total number of parallel processes
 *   argv[3]  procNum        Index of this process, starting from 0
 *   argv[4]  a              Dimensionless black-hole spin
 *   argv[5]  height         Lamppost corona height
 *   argv[6]  numMdot        Number of accretion-rate / disk-height grid values
 *   argv[7]  initMdot       Initial Mdot / Mdot_Edd value
 *   argv[8]  finalMdot      Final Mdot / Mdot_Edd value
 *   argv[9]  outFilePrefix  Prefix for output files
 *
 * Output:
 *   Writes one text file per accretion-rate grid point and process:
 *
 *       <outFilePrefix>_<thicknessIndex>_<procNum>.txt
 *
 *   Each row corresponds to one photon and contains:
 *
 *       imgAlpha, imgBeta,
 *       finalEnergy / initialEnergy,
 *       final t, r, theta, phi,
 *       disk scale height,
 *       projected cylindrical radius,
 *       disk Lorentz factor,
 *       hitDiskSwitch
 *
 * Notes: In the current FENRIR github, the disk thickness is defined by just the frontHeightTerm, which is dependent on 
 * both efficiency and accretion rate. I have undone that change to revert the code back to being dependent on Mdot value.
 * This is more aligned with the original paper and results. 
 * 
 * NOTE: The Mdot grid must run from thickest to thinnest. We intentionally do not
 *  reset posVec inside this loop: after a photon reaches one disk surface, we
 *  continue the same geodesic to find intersections with lower/nested disk
 *  surfaces, avoiding repeated propagation from the corona.
*/


// Compatibility replacement for std::to_string on older cluster compilers.
#include <string>
#include <sstream>
namespace patch
{
    template < typename T > std::string to_string( const T& nInput )
    {
        std::ostringstream stm ;
        stm << nInput ;
        return stm.str() ;
    }
}

int main(int argc, char* argv[]){

	//This divides the job over multiple processes
	int n, totProcs, procNum, photonPerProc;
	n = atoi(argv[1]); //total number of photons traced from corona
	totProcs = atoi(argv[2]); //total number of process
	procNum = atoi(argv[3]); //index number of current process (starts at 0)
	photonPerProc = n/totProcs; //number of photons per process

	//Accessing user-input for spin (a) and corona height (height)
	a = atof(argv[4]);
	height = atof(argv[5]);

	//Accretion parameters (initAcc > finalAcc)
	numMdot = atoi(argv[6]); //Total number of Mdot values
	initMdot = atof(argv[7]); //Initial mdot that each process will consider
	finalMdot = atof(argv[8]); //Final mdot that each process will consider

	//Output file information
	outFilePrefix = argv[9]; //Prefix of input file (file name: outputFilePrefix_thicknessIndex_processIndex.txt)

	//Calculating the deltaAcc from the inputs
	if (numMdot > 1){
		deltaMdot = (finalMdot - initMdot)/(numMdot-1);
	}else{
		deltaMdot = 0.;
	}

	// Kerr horizon radius in units of r_g.
	rEvent = 1. + pow(1.-(a*a),0.5);

	// Bardeen et al. ISCO radius for circular equatorial orbits.
	z1 = 1. + (pow(1-a*a,(1./3.))*((pow(1+a,(1./3.)))+(pow(1-a,(1./3.))))); //
	z2 = sqrt((3.*a*a)+(z1*z1)); //
	if (a < 0.){
		rIsco = 3.+z2+sqrt((3-z1)*(3+z1+(2.*z2))); //
	}else{
		rIsco = 3.+z2-sqrt((3-z1)*(3+z1+(2.*z2)));
	}

	//Efficiency and accretion rate (Mdot/Eddington) to calculate scale height
	// Radiative efficiency eta = 1 - E(r_ISCO), used in the disk-height prescription.
	efficiencyUpper = (rIsco*rIsco)-(2.*rIsco)+(a*sqrt(rIsco)); 
	efficiencyLower = (rIsco*rIsco)-(3.*rIsco)+(2.*a*sqrt(rIsco)); 
	efficiency = 1.-(efficiencyUpper/(rIsco*sqrt(efficiencyLower))); 

	//Defining the minimum R cutoff (at which point, the photon will stop being traced)
	rLimitLow = rEvent + horizonStop;  //minimum radial distance (currently r_horizon + 0.01) //

	//Initializing the array of output streams
	std::ofstream outStreams[numMdot];

	//Opening up each of the output streams by looping over output stream array outStreams
	int fileIndex = 0;
	std::string fileName;
	while (fileIndex < numMdot){
		fileName = outFilePrefix + "_" + patch::to_string(fileIndex) + "_" + patch::to_string(procNum) + ".txt";
		outStreams[fileIndex].open(fileName.c_str());
		fileIndex += 1;
	}


  	//Finding initial positions
	findPosition(); //in corona_functions.h

  	//Finding rotational velocity Omega (dphi/dt)
  	findOmega(); //in corona_functions_test.h

  	//Calculating the corona orthonormal tetrad components
  	findTetrad(); //in corona_functions_test.h

  	//Calculating value of metric components at corona position
  	findComponents(); //in corona_functions_test.h

  	//Calculating trig values at position of corona (used for carter calculation)
  	corCos = std::cos(initTheta);
  	corSin = std::sin(initTheta);
  	corCot = corCos/corSin;

	//Angular griding for Dauser emissivity method
	double minAngle = ((a*a)*((initRadius*initRadius)-(2.*initRadius)+(a*a)))/(((initRadius*initRadius)+(a*a))*((initRadius*initRadius)+(a*a)));
	double maxAngle = M_PI-minAngle;
	double delAngle = (maxAngle - minAngle)/(n-1);

	//The angles aren't right for parallelization, but that's okay for now. Need to fix it.
	// TODO: CHECK IF The angular grid is split by photon index, which may not distribute angles evenly across processes.
	//Dummy variable for main photon tracing loop (each iteration is new photon)
 	int j = procNum*photonPerProc;

	//Accretion rate dummy variable (will change each time you hit a different disk level throughout a single j loop)
	int MdotIndex;

 	while (j < (procNum+1)*photonPerProc){
  		//Generating initial angle of photon path relative to polar azis
  		imgAlpha = M_PI - (minAngle + j*delAngle);//randAlpha();  //Calculating the image verticle angle (alpha)
  		imgBeta = 0.; //Due to azimuthal symmetry in LP case, I can choose one beta value.

		//Calculating corona rest frame photon momentum vector
		findMomentum(); //in corona_functions_test.h

		//Calculating conserved quantities from B-L momentum vector
      	findConserved(); //in corona_functions_test.h

      //Determining the sign of rdot based upon the verticle angle (alpha)
      if (imgAlpha < M_PI/2.){
        rSqrtSwitch = 1.;
      }else{
        rSqrtSwitch = -1.;
      }

      //Determining the sign of thetadot based upon the horizontal angle (beta)
      if ((imgBeta < M_PI/2.) || (1.5*M_PI < imgBeta)){
        thSqrtSwitch = 1.;
      }else{
        thSqrtSwitch = -1.;
      }

      //Giving the photon's 4-vector its initial conditions
      posVec[0] = initTime;
      posVec[1] = initRadius;
      posVec[2] = initTheta;
      posVec[3] = initPhi;

	  //Starting the propagation loop, where I will write to a different file for each Mdot value
	  MdotIndex = 0;
	  while (MdotIndex < numMdot){

		  //Calculating the accretion rate
		  accretion = initMdot + (MdotIndex*deltaMdot);

		  //Calculating the scale height normalization term from the accretion rate
		  heightFrontTerm = 2.*(3./(2.*efficiency))*accretion;

		  if (posVec[1]*cos(posVec[2]) > scaleHeightFnct(posVec[1],posVec[2])){ //scaleHeightFnct in disk_equations.h
			  //Resetting the hitDiskSwitch to zero. Changes to 1 if photon hits disk.
			  hitDiskSwitch = 0;

			  //Propagate the photon to the disk
		      propagate(posVec,momVec,dStep,tolerance,maxStep,rLimitLow,rLimitHigh,rEvent,scaleHeightValue,rProjected);

		  }

		  //Calculating the projected radius and the scale height of the disk
		  rProjected = posVec[1]*std::sin(posVec[2]);
		  scaleHeightValue = heightFrontTerm*(1 - sqrt(rIsco/rProjected));

		  //Calculating the one-form of the photon's momentum 4-vector
		  vecToOneForm(posVec, momVec, momOneForm);

		  //Calculating the disk's velocity 4-vector
		  diskVelocity(posVec, diskVelVec, scaleHeightValue, rProjected);

		  //Calculating the final energy of the photon (-1 x dot product of photon one-form momentum and the disk 4-velocity)
		  finalEnergy = (momOneForm[0]*diskVelVec[0]) + (momOneForm[1]*diskVelVec[1]) + (momOneForm[2]*diskVelVec[2]) + (momOneForm[3]*diskVelVec[3]);
		  finalEnergy = -1.*finalEnergy;

		  //Calculating the gamma (Lorentz factor) of the disk element as seen from LNRF.  This is done by multiplying the 4-vector by the basis one-forms (i.e. the dot product)
		  lorentz = findLorentz(posVec,diskVelVec);

		  //Printing to output file
		  outStreams[MdotIndex] << imgAlpha << " " << imgBeta << " " << (finalEnergy/energy) << " " << posVec[0] << " " << posVec[1] << " " << posVec[2] << " " << posVec[3] << " " << scaleHeightValue << " " << rProjected << " " << lorentz << " " << hitDiskSwitch << "\n";

		  //advancing the accretion index by 1
		  MdotIndex += 1;

	  	}
	//Incrementing j index (changing to next photon)
		j++;
	}
	//closing the various output files to be saved by looping over the stream array
	fileIndex = 0;
	while (fileIndex < numMdot){
		outStreams[fileIndex].close();
		fileIndex += 1;
	}

	return 0;  //returning integer value of 0 if executed to completion

} //end of main function
