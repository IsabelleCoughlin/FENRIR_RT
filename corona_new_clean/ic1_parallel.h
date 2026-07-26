#include <cmath>

/*
 * Global input/configuration values for the parallel ray-tracing executable.
 *
 * This header stores run-level parameters that are shared across the ray-tracing
 * code: black-hole spin, corona geometry, disk-thickness/accretion-rate grid,
 * integration limits, and output naming.
 *
 * Most of these variables are set in main() from command-line arguments, then
 * read by helper functions in the corona, disk, metric, and propagation modules.
 *
 * NOTE:
 *   This header defines globals directly. That works for this code style, but it
 *   means the header should only be included in one translation unit, or the
 *   project must rely on a nonstandard single-file/include build pattern.
 */

/* Dimensionless black-hole spin. Set from argv[4]. */
double a;


/*
 * Current accretion-rate value, usually Mdot / Mdot_Edd.
 *
 * This is updated inside the Mdot/thickness loop and is used to compute the disk
 * scale-height normalization.
 */
double accretion = 0.0;


/*
 * Corona geometry selector.
 *
 * Supported values:
 *   "lp"      : lamppost corona on/near the spin axis
 *   "offaxis" : corona at finite cylindrical radius above the disk
 */
std::string coronaType = "lp";


/* Lamppost corona height. Used when coronaType == "lp". Set from argv[5]. */
double height;


/*
 * Disk-thickness / accretion-rate grid.
 *
 * Legacy names use "Thickness", but in the Mdot-based version these may represent
 * a grid in Mdot / Mdot_Edd. If so, consider renaming these to numMdot, initMdot,
 * finalMdot, deltaMdot, and mdot.
 */
int numMdot;          // Number of grid values to trace
double initMdot;      // First grid value
double finalMdot;     // Final grid value
double deltaMdot;     // Step between grid values
double Mdot;          // Current grid value


/*
 * Off-axis corona parameters. Used when coronaType == "offaxis".
 */
double initCylRadius = 1.;     // Corona cylindrical radius
double heightAboveDisk = 0.1;  // Vertical distance above local disk surface
double coronaPhi = 0.;         // Corona azimuthal position

/*
 * If true, initCylRadius is interpreted in units of r_ISCO.
 * If false, initCylRadius is interpreted directly in code length units.
 */
bool isRisco = true;


/*
 * Off-axis corona angular velocity prescription.
 *
 * Supported values:
 *   "lnrf" : rotate with the local non-rotating frame
 *   "kep"  : use Keplerian angular velocity
 *   other  : use rotOmegaInput directly
 */
std::string rotOmegaType = "kep";


/* User-supplied angular velocity, used when rotOmegaType is not "lnrf" or "kep". */
double rotOmegaInput;


/*
 * Random-number seed.
 *
 * Currently initialized from wall-clock time. If exact reproducibility is needed,
 * make this a command-line input instead.
 */
double seed = time(NULL);


/*
 * Photon integration controls.
 */
double maxStep = 100000000;  // Maximum number of integration steps
double horizonStop = 0.01;   // Stop if photon gets this close to event horizon
double rLimitHigh = 1000.;   // Treat photon as escaped beyond this radius

/*
 * Step-refinement tolerance.
 *
 * In this code's convention, larger values correspond to finer sampling.
 * Check propagate_rk4_* for the exact meaning before tuning.
 */
double tolerance = 1000.;


/*
 * Prefix used when writing output files.
 *
 * Output files are usually named:
 *   <outFilePrefix>_<gridIndex>_<procNum>.txt
 *
 * Set from argv[9].
 */
std::string outFilePrefix;