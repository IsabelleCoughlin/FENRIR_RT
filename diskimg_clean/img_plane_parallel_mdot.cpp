#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "ic1_multi_new.h"
#include "variable_names.h"
#include "photon_geos_thindisk.h"
#include "metric_kerr.h"
#include "disk_equations_new.h"
#include "propagate_rk4_thindisk_new.h"

//Last update by Corbin Taylor at 11:28 am on 3/27/19


/*
 * Create one relativistic disk image for each requested value of
 * Mdot/Mdot_Edd. Each process traces a disjoint range of image columns.
 * Photons are propagated backward from a distant observer until they
 * reach the prescribed disk photosphere or the near-horizon cutoff.
 */


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
	
	//Getting total number of photons and rows/process from command line
	//This divides the vertical dimension into parts, where each process does a certain number of rows
	int n, totProcs, colPerProc, procNum;
	n = atoi(argv[1]); //total number of rows/columns for disk image (n x n photons)
	totProcs = atoi(argv[2]); //total number of processes
	procNum = atoi(argv[3]); //index number of process (starts at 0)
	colPerProc = n/totProcs;
	
    //Physical parameters
    a = atof(argv[4]);
	inclination = acos(atof(argv[5]));
	
	//Accretion parameters (initAcc > finalAcc)
	numMdot = atoi(argv[6]);
	initMdot = atof(argv[7]);
	finalMdot = atof(argv[8]);
	imgSize = atof(argv[9]);
	outFilePrefix = argv[10];
	
	//Calculating the deltaAcc from the inputs
	if (numMdot > 1){
		deltaMdot = (finalMdot - initMdot)/(numMdot-1);
	}else{
		deltaMdot = 0.;
	}
		
	
	//Calculating radius of event horizon
	rEvent = 1. + pow(1.-(a*a),0.5);

	//Calculating the inner-most stable orbital radius
	z1 = 1. + (pow(1-a*a,(1./3.))*((pow(1+a,(1./3.)))+(pow(1-a,(1./3.))))); //
	z2 = sqrt((3.*a*a)+(z1*z1)); //
	if (a < 0.){
		rIsco = 3.+z2+sqrt((3-z1)*(3+z1+(2.*z2))); 
	}else{
		rIsco = 3.+z2-sqrt((3-z1)*(3+z1+(2.*z2)));
	}

	//Efficiency and accretion rate (Mdot/Eddington) to calculate scale height
	efficiencyUpper = (rIsco*rIsco)-(2.*rIsco)+(a*sqrt(rIsco)); 
	efficiencyLower = (rIsco*rIsco)-(3.*rIsco)+(2.*a*sqrt(rIsco)); 
	efficiency = 1.-(efficiencyUpper/(rIsco*sqrt(efficiencyLower))); 

	//Defining the minimum R cutoff (at which point, the photon will stop being traced)
	rLimit = rEvent + horizonStop;  //minimum radial distance (currently r_horizon + 0.01) //
	
	//Setting up the initial photon position
	initRadius = obsRadius; //radial position of observer
	initPhi = obsPhi; //phi position of observer
	initTheta = inclination; //theta position of observer (should equal inclination angle, alpha)
	initTime = 0.; //initial time
	
	//Initializing the array of output streams
	std::ofstream outStreams[numMdot];
	
	//Filling stream name array
	int fileIndex = 0;
	
	//Opening up each of the output streams by looping over output stream array outStreams
	std::string fileName;
	while (fileIndex < numMdot){
		fileName = outFilePrefix + "_" + patch::to_string(fileIndex) + ".txt";
		outStreams[fileIndex].open(fileName.c_str());
		fileIndex += 1;
	}
	
    //Looping over NxN image plane
    int j,k; //defining looping dummy variables. j is x variable. k is y variable
    
    //Accretion rate dummy variable
    int MdotIndex;

    j = procNum*colPerProc; //setting dummy j to 0
    while (j < (procNum+1)*colPerProc){
        k = 0;
        while (k < n){

            imgX = ((j - (n/2.))+0.5)*(imgSize/n); //calculating x value for j value, could always do it outside of loop to not have to constantly do calculation
            imgY = ((k - (n/2))+0.5)*(imgSize/n); //calculating y value for k value
            imgB = pow(pow(imgX,2.)+pow(imgY,2.),0.5); //calculating Euclidean 2-D distance on image plane

            //calculating the phi angular momentum and carter constant of photon from image plane coordinates
            angmom = -1.*energy*imgB*std::sin(inclination)*(imgX/imgB); //phi angular momentum
			//NOTE: this is different than Corbin's original implementation
			// where most of the energy terms are omitted. Should not make a difference
			//since energy is fixed at 1/
            carter = energy*energy*pow(imgB*(imgY/imgB),2.) - energy*energy*pow(a*std::cos(initTheta),2.) + pow(angmom*std::cos(initTheta)/std::sin(initTheta), 2.); //carter constant

			//setting the switches for the sign of rdot and thetadot
            rSqrtSwitch = 1.;
            if (imgY < 0.){
                thSqrtSwitch = -1.;
            }else{
                thSqrtSwitch = 1.;
            }

			//setting the initial values of the photon's position 4-vector to be found in the ic file
            posVec[0] = initTime; //time
            posVec[1] = initRadius; //spherical radius
            posVec[2] = initTheta; //vertical angle theta
            posVec[3] = initPhi; //horizontal angle phi
            
			//Starting the propagation loop, where I will write to a different file for each Mdot value
			MdotIndex = 0;
			while (MdotIndex < numMdot){
				//Calculating the accretion rate
				Mdot = initMdot + (MdotIndex * deltaMdot);

				// Fenrir photosphere height z = 2H
				heightFrontTerm = (3.0 / efficiency) * Mdot;
				
				if (posVec[1]*cos(posVec[2]) > scaleHeightFnct(posVec[1],posVec[2])){
				
            		//Propagate the photon from the observer to the disk
            		propagate(posVec,momVec,dStep,tolerance,maxStep,rLimit,rEvent,scaleHeightValue,rProjected);
				
				}
            	//Calculating the pseudo-cylindrical radius and the vertical height of the disk above the mid-plane
            	rProjected = posVec[1]*std::sin(posVec[2]); //pseudo-cylindrical radius
            	scaleHeightValue = posVec[1]*std::cos(posVec[2]);

            	//Calculating the one-form of the photon's momentum 4-vector
            	vecToOneForm(posVec, momVec, momOneForm);

            	//Calculating the disk's velocity 4-vector
				diskVelocity(posVec, diskVelVec, scaleHeightValue, rProjected);
			
            	//Calculating the final energy of the photon (E = -p*U, dotting photon 4-momentum with disk 4-velocity)
            	finalEnergy = (momOneForm[0]*diskVelVec[0]) + (momOneForm[1]*diskVelVec[1]) + (momOneForm[2]*diskVelVec[2]) + (momOneForm[3]*diskVelVec[3]);
            	finalEnergy = -1.*finalEnergy;

				//Outputting data to 'myfile': x, y, g, final_t, final_r, final_theta, final_phi, disk_H, pseudo-cylindrical_r
            	outStreams[MdotIndex] << imgX << " " << imgY << " " << (energy/finalEnergy) << " " << posVec[0] << " " << posVec[1] << " " << posVec[2] << " " << posVec[3] << " " << scaleHeightValue << " " << rProjected << "\n";
				
				//advancing the accretion index by 1
				MdotIndex += 1;
			}
            //advancing k by 1
            k++;
        }
        //advancing j by 1
        j++;
    }
    //closing the various output files to be saved by looping over the stream array
    fileIndex = 0;
    while (fileIndex < numMdot){
		outStreams[fileIndex].close();
		fileIndex += 1;
	}
    return 0;  //returning integer value of 0 if executed to completion

}  //end of main function
