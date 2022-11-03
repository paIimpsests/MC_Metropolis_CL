/* Monte-Carlo implementation for simulation of a HS system using Metropolis
 * Scheme and Cell Lists
 * Antoine Castagnede | October 2022
 *
 *  TYPESETTING FOR ALL ASSOCIATED .SPH FILES :
 *  -------------------------------------------
 *  """
 *  &N
 *  x-length y-length z-length
 *  type x y z r
 *  """
 *  Where :
 *  - N is the number of particles
 *  - .-length is the size of the box in the . direction
 *  - type, x, y, z, and r are as defined for the Sphere structure below
 *  And the last line is repeated N times to describe all particles
 *  independantly
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>       // need to compile with `-lm` flag
#include "mt19937ar.c"

#define T 300           // absolute temperature
#define kB 1.380649e-23 // Boltzmann constant
#define eps kB*T        // unit of energy


// STRUCTURES
// ----------
typedef struct Sphere Sphere;
struct Sphere
{
	double x;	// reduced x coordinate 
	double y;	// reduced y coordinate
        double z;       // reduced z coordinate
	double r;	// reduced radius --- redundant with sigma ... might remove later on
        char type;      // particle type --- for visualization purposes only
        Sphere *next;   // pointer to the next particle in CL
        Sphere *prev;
        int CLindex;
};


// FUNCTIONS
// ---------
int overlapCheck1P(int n, int m);
int overlapCheckvar(Sphere one, Sphere two);
int overlapCheckNNP(void);
int overlapCheckCL(int n);
void sanityCheck();
void readInit(char *filename);
void writeCoords(char *filename);
int particleMove(double particleStepTune);
int volumeMove(double volumeStepTune);
double tuneStepSize(int nSuccess, int nCycles, double acceptanceRate);
void stepTuneWrapper(int MCCycle, double acceptanceRate);
double measurePF(void);
int initCL(void);
int resizeCL(void);
int updateCL(void);
int updateSingleCL(int n);
int retrieveIndex(int u, int v, int w);
int CLOverflow(void);

// GLOBAL VARIABLES
// ----------------
// System parameters
int N = 0;                      // [READ FROM IMPUT] number of particles in the system
double L = 0.0f;                // [READ FROM INPUT] reduced box size
double P = 5.0f;                // reduced pressure
double sigma = 1.0f;            // unit of length = 1 = particle size
Sphere *pSpheres = NULL;        // [READ FROM INPUT] pointer to the table of particles data
// CL parameters
Sphere **CLTable = NULL;
double sCell1D = 0.0f;
int nCell1D = 0;
int nCell = 0;


int overlapCheck1P(int n, int m)
{
        /*
         * Function:    overlapCheck1P
         * ---------------------------
         * Checks for overlap between one pair of particules by computing the
         * distance between their centers and comparing it to the sum of their
         * radii
         *
         * n:           index of the first particle in the data list
         * m:           index of the second particle in the data lsit
         *
         * return:      1 is there is overlap, 0 if not
         */

	int overlap = 0;
        double  dx = pSpheres[n].x - pSpheres[m].x,
                dy = pSpheres[n].y - pSpheres[m].y,
                dz = pSpheres[n].z - pSpheres[m].z;

                dx = dx - L * rint(dx / L);
                dy = dy - L * rint(dy / L);
                dz = dz - L * rint(dz / L);

	if ( dx * dx + dy * dy + dz * dz
           < ( (pSpheres[n].r + pSpheres[m].r) * (pSpheres[n].r + pSpheres[m].r) )
           )
	        overlap = 1;
	return overlap;
}

int overlapCheckvar(Sphere one, Sphere two)
{
        int overlap = 0;
        double dx = one.x - two.x,
               dy = one.y - two.y,
               dz = one.z - two.z;
        dx = dx - L * rint(dx / L);
        dy = dy - L * rint(dy / L);
        dz = dz - L * rint(dz / L);
        if ((dx * dx + dy * dy + dz * dz) < ((one.r + two.r) * (one.r + two.r)))
                overlap = 1;       
        return overlap;
}

int overlapCheckNNP(void)
{
        /* Function:     overlapCheckNNP
         * -----------------------------
         * Checks for overlap in the whole system, i.e. for N(N-1)/2 pairs
         *
         * return:     1 if there is overlap, 0 if not
         */
        int overlap = 0, n = 1, m = 0;
        while ((overlap == 0) && (n < N))
        {
                while ((overlap == 0) && (m < N))
                {
                        overlap = (n != m) * overlapCheck1P(n,m);
                        m++;
                }
                m = n;
                n++;
        }
        return overlap;
}


int overlapCheckCL(int n)
{
        /*
         * Function:    overlapCheckCL
         * ---------------------------
         * Checks for overlap between a given particle and all other particles
         * that can be found in neighbouring cells. Relies on the Cell Lists
         * implementation.
         *
         * n:           index of the particle to check overlap for
         *
         * return:      1 if there's overlap, 0 if not
         */
        int a = pSpheres[n].x / sCell1D;
        int b = pSpheres[n].y / sCell1D;
        int c = pSpheres[n].z / sCell1D;
        int index = 0;
        int overlap = 1;
        Sphere *current = NULL;
        for (int i = a-1; i < a+2; i++)
        {        
                for (int j = b-1; j < b+2; j++)
                {
                        for (int k = c-1; k < c+2; k++)
                        {
                                index = retrieveIndex(i,j,k);
                                current = CLTable[index];
                                while (current != NULL)
                                {
                                        if ((pSpheres+n != current) && overlapCheckvar(pSpheres[n], *current))
                                                overlap *= 0;
                                        current = current->next;
                                }
                        }
                }
        }
        return !(overlap);
}

void sanityCheck(void)
{
        /*
         * Function: sanityCheck
         * ---------------------
         * Ends program execution if there is overlap in the system or if
         * something is wrong with the cell lists
         */
        if (overlapCheckNNP())
        {
                printf("Something is VERY wrong... There's overlap...\n");
                exit(0);
        }
        if (CLOverflow())
        {
                printf("Something is VERY wrong... The Cell Lists are overflowing...\n");
                exit(0);
        }
}

void readInit(char *filename)
{
        /*
         * Function:    readInit
         * ---------------------
         * Initializes the table of data of all N particles by reading
         * it from a user supplied .sph file; particles coordinates
         * are written in reduced units; also retrieves **cubic** box size
         * NB: Reading is done only for files with typesetting indicated above
         * NB: It is assumed that the supplied init file is written in standard units so
         * conversion to reduced units is applied
         *
         * *filename:   pointer to the name of the initialization .txt file
         */

        FILE *initfile = NULL;
        int i = 0;
        double rTemp = 0.0f;
        initfile = fopen(filename, "r");
        if (initfile != NULL)
        {
                // Read value of N and allocate memory accordingly
                fscanf(initfile, "%*c%d%*c", &N);
                pSpheres = malloc(N * sizeof(Sphere));
                // Read value of box size
                fscanf(initfile, "%lf %*f %*f%*c", &L);
                // Populate table of particles data in standard units
                for (i = 0; i < N; i++)
                {
                        fscanf  ( initfile,
                                  "%c %lf %lf %lf %lf%*c",
                                  &pSpheres[i].type,
                                  &pSpheres[i].x,
                                  &pSpheres[i].y,
                                  &pSpheres[i].z,
                                  &pSpheres[i].r
                                );
                        // Find larger radius to rescale everything to obtain
                        // the expected sigma = 1 for the (larger) particle size
                        if (pSpheres[i].r > rTemp)
                                rTemp = pSpheres[i].r;
                }
                fclose(initfile);
                // Rescale everything w.r.t. larger particle size to obtain the
                // expected sigma = 1, also converts to reduced units
                L /= (2.0f * rTemp);
                for (i = 0; i < N; i++)
                {
                        pSpheres[i].x /= (2.0f * rTemp);
                        pSpheres[i].y /= (2.0f * rTemp);
                        pSpheres[i].z /= (2.0f * rTemp);
                        pSpheres[i].r /= (2.0f * rTemp);
                }
        }
}

void writeCoords(char *filename)
{
        /*
         * Function:    writeCoords
         * -------------------------
         * Writes the position, radius, and type data for all N
         * particles in a .sph file; standard units are used
         * NB: Writing is done using the typesetting indicated above
         *
         * *filename:   pointer to the name of the output .sph file
         *
         */
        FILE *outfile = NULL;
        int i = 0;
        double boxSize = L * sigma;
        outfile = fopen(filename, "a");
        if (outfile != NULL)
        {
                fprintf(outfile, "&%d\n", N);
                fprintf(outfile, "%lf %lf %lf\n", boxSize, boxSize, boxSize);
                for (i = 0; i < N; i++)
                        fprintf(outfile,
                                "%c %lf %lf %lf %lf\n",
                                pSpheres[i].type,
                                pSpheres[i].x * sigma,
                                pSpheres[i].y * sigma,
                                pSpheres[i].z * sigma,
                                pSpheres[i].r * sigma
                               );
                fclose(outfile);
        }
}

int particleMove(double particleStepTune)
{
        /*
         * Function: particleMove                
         * ---------------------
         * Tries to move a randomly selected particle in space by a small amount
         * and checks for overlap
         * If there is no overlap, overwrites the previous particle position
         * with the new one
         *
         * particleStepTune:    tuning factor for particle step size
         *
         * return:               0 if attempt failed, 1 if attempt succeeded
         */ 
        // Generate random displacements in x, y, z
        double delta[3] = {(genrand() - 0.5) / sigma * (particleStepTune),
                           (genrand() - 0.5) / sigma * (particleStepTune),
                           (genrand() - 0.5) / sigma * (particleStepTune)
                          };
        // Randomly select a particle in the system
        int n = (int) (genrand() * N);
        // Set up a backup of the moved particle in case of overlap
        Sphere  bufferSphere = {        pSpheres[n].x,
                                        pSpheres[n].y,
                                        pSpheres[n].z
                                };
        // Move randomly selected particle by abiding to PBC and nearest image
        // convention
        pSpheres[n].x = fmod((pSpheres[n].x + delta[0]) + 2 * L, L);
        pSpheres[n].y = fmod((pSpheres[n].y + delta[1]) + 2 * L, L);
        pSpheres[n].z = fmod((pSpheres[n].z + delta[2]) + 2 * L, L);
        // Update Cell List
        updateSingleCL(n);
        // Check for overlap and move back particle if need be
        if (overlapCheckCL(n))
        {
                pSpheres[n].x = bufferSphere.x;
                pSpheres[n].y = bufferSphere.y;
                pSpheres[n].z = bufferSphere.z;
                updateSingleCL(n);
                return 0;
        }
        else
                return 1;                
}

int volumeMove(double volumeStepTune)
{
        /*
         * Function:    volumeMove
         * -----------------------
         * Attempts to change the volume of the simulation box by a small
         * amount according to the known acceptance rule
         * 
         * volumeStepTune:      tuning factor for volume step size
         *
         * return:              0 if attempt failed, 1 if attempt succeeded
         */
        // Generate random volume change, define new adequate quantities
        double  vol = L * L * L,
                delta[2] = {(genrand() - 0.5) / (sigma) * (volumeStepTune),
                            genrand()
                           },
                newVol = vol + delta[0],                  
                newL = cbrt(newVol),                      
                ratio = newL / L,
                rule = exp(- P * (newVol - vol) + N * log(newVol / vol));
        // Scale the system according to the volume change
        for (int i = 0; i < N; i++)
        {
                pSpheres[i].x *= ratio;
                pSpheres[i].y *= ratio;
                pSpheres[i].z *= ratio;
        }
        L *= ratio;
        // Scale back the system to its previous state according to the
        // relevant acceptance rule, also update CLs accordingly
        if (delta[1] > rule) 
        { // Rule rejected
                for (int i = 0; i < N; i++)
                {
                        pSpheres[i].x /= ratio;
                        pSpheres[i].y /= ratio;
                        pSpheres[i].z /= ratio;
                }
                L /= ratio;
                return 0;
        }
        else if (ratio < 1)
        { // Box shrinks
                if (overlapCheckNNP())
                { // Overlap
                        for (int i = 0; i < N; i++)
                        {
                                pSpheres[i].x /= ratio;
                                pSpheres[i].y /= ratio;
                                pSpheres[i].z /= ratio;
                        }
                        L /= ratio;
                        return 0;
                }
                else
                { // No overlap
                        sCell1D *= ratio;
                        if (sCell1D < sigma)
                        // Cell size does not guarantee check for all possible
                        // overlaps
                        // Need to reduce cell size, thus number of cells
                                updateCL();
                        else
                        // Cell size is sufficient to check for all possible
                        // overlaps
                        // Only need to reduce cell size accordingly
                                resizeCL();
                        return 1;
                }
        }
        else
        { // Box expands = no overlap
                sCell1D *= ratio;
                if (sCell1D > (2.0f * sigma))
                // Cell size is such that we check for overlap where there
                // cannot be any
                // Need to reduce cell size, thus increase number of cells,
                // thus free existing CL table and re-allocate enough memory
                        initCL();
                else
                // Cell size is necessary (?) to check for all possible overlaps
                // Only need to increase cell size accordingly
                        resizeCL();
                return 1;
        }
}

double tuneStepSize(int nSuccess, int nCycles, double acceptanceRate)
{
        /*
         * Function:    tuneStepSize
         * -------------------------
         * Tunes the step size for volume or particle moves by +-2% depending
         * on whether the success rate is above or below the specified
         * acceptance rate
         *
         * NSuccess:    number of successes over the last nCycles cyles
         * acceptanceRate:      targeted acceptance rate
         * nCycles:     number of cycles to average the number of successes over
         *
         */
        if (nSuccess - (int) (acceptanceRate * nCycles) > 0)
                return 1.02;
        else
                return 0.98; 
}

void stepTuneWrapper(int MCCycle, double acceptanceRate)
{
        int i = 1,
            j = 0,
            sV = 0,                             // number of volumeMove() successes
            sP = 0,                             // number of particleMove() successes
            rPstep = 10000,                     // number of MC cycles before tuning particle step size
            rVstep = (N + 1) * rPstep;          // number of MC cycles before tuning volume step size
        double q = 0.0f,
               volumeStepTune = 1.0f,           // tune factor for volume step size
               particleStepTune = 1.0f;         // tune factor for particle step size
        for (i = 1; i < MCCycle + 1; i++)
        {
                for (j = 0; j < N; j++)
                {
                        q = genrand();
                        if (q < (1.0f / (float) (N + 1)))
                                sV += volumeMove(volumeStepTune);
                        else
                                sP += particleMove(particleStepTune);
                }
                if (i % rPstep == 0)
                {
                        sanityCheck();
                        particleStepTune *= tuneStepSize(sP, rPstep, acceptanceRate);
                        sP = 0;
                }
                if (i % rVstep == 0)
                {
                        volumeStepTune *= tuneStepSize(sV, rVstep, acceptanceRate);
                        sV = 0;
                }
        } 
}

double measurePF(void)
{
        /*
         * Function: measurePF
         * -------------------
         * Computes the reduced packing fraction for the given parameters
         *
         * return:     reduced packing fraction
         */
        double PF = (double) N * sigma * sigma * sigma * M_PI / (6.0f * L * L * L);
        return PF;
}

int initCL(void)
{
        /*
         * Function: initCL
         * ----------------
         * Initializes the table of cell lists according to the current state of
         * the system. Does it after allocating enough memory for the system to
         * be divided into cells of size close (>=) to the targeted size sigma. 
         *
         * USE-CASE:    At initialization; after volumeMove() if cell size gets
         *              > 2*sigma
         */
        // CL parameters that need to be updated with L
        sCell1D = L / ((int) (L / sigma));
        nCell1D = (int) (L / sCell1D);
        nCell = nCell1D * nCell1D * nCell1D;
        // Allocating memory for the nCell1D * nCell1D * nCell1D table of CLs
        free(CLTable);
        CLTable = malloc(nCell * sizeof(**CLTable));
        if (CLTable == NULL)
                exit(EXIT_FAILURE);
        // Initializing all CLs to a NULL pointer that also corresponds to the
        // last element of the CL
        for (int i = 0; i < nCell; i++)
                CLTable[i] = NULL;
        // Updating CLs according to the current state of the system
        for (int i = 0; i < N; i++)
        {
                // Retrieve index of relevant cell
                int u = pSpheres[i].x / sCell1D;
                int v = pSpheres[i].y / sCell1D;
                int w = pSpheres[i].z / sCell1D;
                int index = retrieveIndex(u, v, w);
                pSpheres[i].CLindex = index;
                // Place particle on top of relevant CL
                pSpheres[i].prev = NULL;
                pSpheres[i].next = CLTable[index];
                if (CLTable[index] != NULL)
                        CLTable[index]->prev = pSpheres+i;
                CLTable[index] = pSpheres+i;
        }
        return 0;
}

int resizeCL(void)
{
        /*
         * Function: resizeCL
         * ----------------
         * Re-initializes the table of cell lists according to the current state
         * of the system and the newly defined sCell1D. The values of nCell1D
         * and nCell are not changed so no extra memory is allocated. Only the
         * size of the cells has changed.
         *
         * USE-CASE:    After volumeMove() success, if the new cell size is
         *              neither < sigma nor > 2*sigma
         */
        // Initializing all CLs to a NULL pointer that also corresponds to the
        // last element of the CLs
        for (int i = 0; i < nCell; i++)
                CLTable[i] = NULL;
        // Updating CLs according to the current state of the system
        for (int i = 0; i < N; i++)
        {
                // Retrieve index of relevant cell
                int u = pSpheres[i].x / sCell1D;
                int v = pSpheres[i].y / sCell1D;
                int w = pSpheres[i].z / sCell1D;
                int index = retrieveIndex(u, v, w);
                pSpheres[i].CLindex = index;
                // Place particle on top of relevant CL
                pSpheres[i].prev = NULL;
                pSpheres[i].next = CLTable[index];
                if (CLTable[index] != NULL)
                        CLTable[index]->prev = pSpheres+i;
                CLTable[index] = pSpheres+i;
        }
        return 0;
}

int updateCL(void)
{
        /*
         * Function: updateCL
         * ------------------
         * Re-initializes the table of cell lists according to the current state
         * of the system. The values of sCell1D, nCell1D and nCell are changed
         * and the last two are reduced but no new memory is allocated, only
         * part of the previously allocated memory is now used.
         *
         * USE-CASE:    After volumeMove() success, if box shrinks and cell
         *              size gets < sigma
         */
        // Initializing all CLs to a NULL pointer that also corresponds to the
        // last element of the CL --- Done before updating nCell so the whole
        // table is updated
        for (int i = 0; i < nCell; i++)
                CLTable[i] = NULL;
        // CL parameters that need to be updated
        sCell1D = L / ((int) (L / sigma));
        nCell1D = (int) (L / sCell1D);
        nCell = nCell1D * nCell1D * nCell1D;
        // Updating CLs according to the current state of the system
        for (int i = 0; i < N; i++)
        {
                // Retrieve index of relevant cell
                int u = pSpheres[i].x / sCell1D;
                int v = pSpheres[i].y / sCell1D;
                int w = pSpheres[i].z / sCell1D;
                int index = retrieveIndex(u, v, w);
                pSpheres[i].CLindex = index;
                // Place particle on top of relevant CL
                pSpheres[i].prev = NULL;
                pSpheres[i].next = CLTable[index];
                if (CLTable[index] != NULL)
                        CLTable[index]->prev = pSpheres+i;
                CLTable[index] = pSpheres+i;
        }
        return 0;
}

int updateSingleCL(int n)
{
        int index = pSpheres[n].CLindex;
        int u = pSpheres[n].x / sCell1D,
            v = pSpheres[n].y / sCell1D,
            w = pSpheres[n].z / sCell1D;
        // Remove particle from old CL
        if (CLTable[index] == pSpheres+n)
                CLTable[index] = pSpheres[n].next;
        if (pSpheres[n].next != NULL)
                pSpheres[n].next->prev = pSpheres[n].prev;
        if (pSpheres[n].prev != NULL)
                pSpheres[n].prev->next = pSpheres[n].next;
        // Fetches the index of the new CL
        index = retrieveIndex(u, v, w);
        pSpheres[n].CLindex = index;
        // Place particle on top of new CL
        pSpheres[n].prev = NULL;
        pSpheres[n].next = CLTable[index];
        if (CLTable[index] != NULL)
                CLTable[index]->prev = pSpheres+n;
        CLTable[index] = pSpheres+n;
        return 0;
}

int retrieveIndex(int u, int v, int w)
{
        int index = 0;
        index = (w + 2 * nCell1D) % nCell1D * nCell1D * nCell1D
                + (v + 2 * nCell1D) % nCell1D * nCell1D
                + (u + 2 * nCell1D) % nCell1D;
        return index;

}

int CLOverflow(void)
{
        int count = 0;
        Sphere *current = NULL;
        for (int i = 0; i < nCell ; i++)
        {
                current = CLTable[i];
                while (current != NULL)
                {
                        count += 1;
                        current = current->next;
                }
        }
        return (count != N);
}

int main(int argc, char *argv[])
{
        // Measuring CPU running time
        struct timespec begin, end;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);

        // System parameters
        double PF = 0.0f;               // system packing fraction
        double rho = 0.0f;              // reduced number density

        // MC parameters
        init_genrand(486488);
        //init_genrand(2801);
        double proba = 0.0f,
               targetedAcceptanceRate = 0.30f,
               fixedStepSize = 1.0f;
        int i = 0,
            j = 0,
            MCCycle = 100000;

       

        // OUTPUT FILES CLEANING
        char *outfilename = "coords.sph";
        FILE *outfile = fopen(outfilename, "w+");
        if (outfile != NULL)
                fclose(outfile);
        char *densityfilename = "density.txt";
        FILE *densityoutfile = fopen(densityfilename, "w+");
        if (densityoutfile != NULL)
                fclose(densityoutfile);
        FILE *writefile = NULL;

        // SYSTEM INITIALIZATION        
        readInit("SC_init_6.sph");
        initCL();
        writeCoords(outfilename);



        // 1st half of simuation time
        // Tuning of step sizes to obtain targetedAcceptanceRate acceptance rate
        stepTuneWrapper(MCCycle / 2, targetedAcceptanceRate);

        // 2nd half of simulation time
        // Measurements of observables --- here: packing fraction to obtain number density
        for (i = 1; i < (MCCycle / 2 + 1); i++)
        {
                for (j = 0; j < N; j++)
                {
                        proba = genrand();
                        if (proba < (1.0f / (float) (N + 1)))
                                volumeMove(fixedStepSize);
                        else
                                particleMove(fixedStepSize);
                }
                PF += measurePF();
                if (i % 1000 == 0)
                {
                        writeCoords("coords.sph");
                        sanityCheck();
                        PF = PF / 1000;
                        rho = 6.0f * PF / (M_PI * sigma * sigma * sigma);
                        writefile = fopen("density.txt", "a");
                        fprintf(writefile, "%d\t%lf\n", i, rho);
                        fclose(writefile);
                        PF = 0.0f;
                }
        }

        free(CLTable);
        free(pSpheres);

        // Stop measuring CPU time
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        long seconds = end.tv_sec - begin.tv_sec;       
        long nanoseconds = end.tv_nsec - begin.tv_nsec;
        double elapsed = seconds + nanoseconds*1e-9;
        printf("\n\nProgram terminated correctly.\nElapsed CPU time:\t%.3fs.\n\n", elapsed);
        
        return 0;
}
