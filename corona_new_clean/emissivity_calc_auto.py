import numpy as np
import numpy.ma as mask
import matplotlib.pyplot as pl
import matplotlib
import csv
import sys
"""
Build a radial illumination profile from Fenrir ray-tracing output.

This script reads one NumPy ray-tracing output file for a given disk
thickness/Mdot index, selects photons that hit the disk between r_ISCO and an
outer cutoff, bins those hits by projected disk radius, and estimates the local
illumination flux in each radial bin.

Inputs:
    argv[1]  filePath        Directory containing input/output .npy files
    argv[2]  inFilePrefix    Prefix for the ray-tracing input file
    argv[3]  outFilePrefix   Prefix for the processed output file
    argv[4]  thicknessIndex  Disk thickness or Mdot grid index
    argv[5]  a               Dimensionless black-hole spin
    argv[6]  hCorona         Lamppost corona height
    argv[7]  nBins           Number of radial bins

Input file:
    <filePath><inFilePrefix><thicknessIndex>.npy

Output file:
    <filePath><outFilePrefix><thicknessIndex>.npy

Output contents:
    [[radial_bin_centers, flux_profile, mean_incident_angle], spin, hCorona]

Notes:
    The flux estimate weights each photon by the energy-shift factor g^Gamma,
    divides by the local disk Lorentz factor and approximate proper annular
    area, and sums photons within each radial bin.
"""

def rIsco(a):
	"""Calculates the inner-most stable circular orbit (or rIsco, in gravitational
	radii Rg) for a given dimensionless spin value. |a| <= 1

	Example: rIsco(0) => 6.0 """
	z1 = 1 + (((1-(a*a))**(1./3.))*(((1+a)**(1./3.))+((1-a)**(1./3.))))
	z2 = ((3*(a*a))+(z1*z1))**0.5
	if (a < 0) and (a >= -1):
		rOut = 3+z2+(((3-z1)*(3+z1+(2.*z2)))**0.5)
	elif (a >= 0) and (a <= 1):
		rOut = 3+z2-(((3-z1)*(3+z1+(2.*z2)))**0.5)
	else:
		print("ERROR: Spin parameter must be within range of -1 to 1, inclusive.")
		sys.exit()
	return rOut


# Kerr metric helper functions used for the proper disk area calculation.
def grSigma(r,a,theta):
	return (r**2.) + ((a*np.cos(theta))**2.)

def grDelta(r,a):
	return (r**2.) - (2.*r) + (a**2.)

def aFunct(r,a,theta):
	aFunctOut = (((r*r)+(a*a))**2.)-(a*a)*grDelta(r,a)*(np.sin(theta)**2.)
	return aFunctOut

def grArea(r,theta,dr,a,dtheta):
	rTerm = np.sqrt((grSigma(r,a,theta)/grDelta(r,a)) + (grSigma(r,a,theta)*((dtheta/dr)**2.)*(np.sin(theta)**2.)))
	phiTerm = np.sqrt((aFunct(r,a,theta)/grSigma(r,a,theta)))
	areaOut = (2.*np.pi)*phiTerm*rTerm*dr
	return areaOut

#Input argv: inFile, outFile, nBins, spin(a), corona_height(hCorona)
print(sys.argv[1:])
filePath = sys.argv[1]
inFilePrefix = sys.argv[2]
outFilePrefix = sys.argv[3]
thicknessIndex = int(sys.argv[4])

a = float(sys.argv[5])
hCorona = float(sys.argv[6])
nBins = int(sys.argv[7])

inFile = filePath + inFilePrefix + str(thicknessIndex) + '.npy'
outFile = filePath + outFilePrefix + str(thicknessIndex) + '.npy'

# Spectral index used in the illumination weighting g^Gamma.
specIndex = 2.

# Input array columns:
#   x, y              : photon launch/emission coordinates or angles
#   gRatio           : energy shift ratio
#   time,r,theta,phi : final photon position
#   scaleHeight      : disk height at hit point
#   projectedRadius  : cylindrical radius rho = r sin(theta)
#   gamma            : disk Lorentz factor
#   diskHitSwitch    : whether photon hit the disk
data = np.load(inFile)
x = data[0]
y = data[1]
gRatio = data[2]
time = data[3]
radius = data[4]
theta = data[5]
phi = data[6]
scaleHeight = data[7]
projectedRadius = data[8]
gamma = data[9]
diskHitSwitch = data[10]

rIn = rIsco(a)
rOut = 100 # Change if needs to be increased


# Keep only usable disk hits between ISCO and rOut, and remove extreme g-ratio outliers.
iRcut = np.where(np.logical_and((np.logical_and(np.logical_and(projectedRadius[:-1] > rIn,projectedRadius[:-1] < rOut), diskHitSwitch[:-1] > 0)),(gRatio[:-1] < 100.)))

specEnergy = (gRatio)[iRcut]
specRadius = (projectedRadius)[iRcut]
specTheta = (theta)[iRcut]
specPhi = (phi)[iRcut]
specGamma = (gamma)[iRcut]


# Delta is the photon emission angle measured relative to the disk direction. Radians
specDelta = np.pi - x[iRcut]

binLimArray = np.logspace(np.log10(rIn), np.log10(rOut), nBins + 1)
specBin = np.sqrt(binLimArray[:-1] * binLimArray[1:])
drArray = binLimArray[1:] - binLimArray[:-1]

dThetaArray = np.zeros(nBins)
numArray = np.zeros(nBins)
fluxArray = np.zeros(nBins)
dDelta = np.zeros(nBins)
binDelta = np.zeros(nBins)

l = 0
while (l < nBins):
	
	rMinBin = binLimArray[l]
	rMaxBin = binLimArray[l + 1]
	specBin[l] = (rMinBin+rMaxBin)/2.
	drArray[l] = rMaxBin-rMinBin
	
	# Photons landing in this radial bin.
	iInBin = np.where((specRadius >= rMinBin) & (specRadius < rMaxBin))
	radInBin = specRadius[iInBin]
	
	thetaInBin = specTheta[iInBin]
	gammaInBin = specGamma[iInBin]
	energyInBin = specEnergy[iInBin]
	deltaInBin = specDelta[iInBin]
	numArray[l] = len(radInBin)
	if (numArray[l] > 0):
		# Estimate how much the disk surface angle varies across this radial bin.
    	# This approximates dtheta in the finite-thickness annulus area formula.
		dThetaArray[l] = np.max(thetaInBin)-np.min(thetaInBin)
		
		# Store the mean incident angle diagnostic for this bin.
		dDelta[l] = np.max(deltaInBin) - np.min(deltaInBin)
		binDelta[l] = np.mean(deltaInBin)
		# Proper GR annular area for photons landing in this bin.
		grAreaInBin = grArea(radInBin,thetaInBin,drArray[l],a,dThetaArray[l])

		# Illumination contribution from photons in this bin.
        # The gamma and proper-area terms convert photon counts into local flux.
        # energyInBin**specIndex applies the redshift/blueshift weighting.
		fluxInBin = (1./(gammaInBin*grAreaInBin)) * (energyInBin**specIndex)
		fluxArray[l] = np.sum(fluxInBin)
	else:
		# Empty radial bins get zero flux and zero diagnostic angle.
		fluxInBin = 0.
		fluxArray[l] = 0.
		binDelta[l] = 0.
	l+=1

# Save radial bin centers, flux profile, mean delta angle, plus spin and height.
outArray = np.array([[specBin, fluxArray, binDelta], a, hCorona], dtype=object)
np.save(outFile, outArray)
