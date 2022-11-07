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





//    LIBRARIES
// ===============
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>       // need to compile with `-lm` flag
#include "mt19937ar.c"





//    PREPROCESSOR CONSTANTS
// ============================
#define T 300           // absolute temperature
#define kB 1.380649e-23 // Boltzmann constant
#define eps kB*T        // unit of energy





//    PROTOTYPES
// ================
// structures
// ----------------
typedef struct Particle Particle;
// functions
// ----------------
int overlapCheck(Particle* one, Particle* two);
int overlapCheckGlobal(void);
int overlapCheckCL(int n);
void sanityCheck();
void readInit(char *filename);
void writeCoords(char *filename);
int particleMove(double particleStepTune);
int volumeMove(double volumeStepTune);
double tuneStepSize(int nSuccess, int nCycles, double acceptanceRate);
void stepTuneWrapper(int MCCycle, double acceptanceRate);
double measurePF(void);
int buildCL(int usecase);
int updateSingleCL(Particle* one);
int retrieveIndex(int u, int v, int w);
int CLOverflow(void);
int snapshot(void);





//    GLOBAL VARIABLES
// ======================
// hard spheres system
// -----------------------
int N = 0;                      // [READ FROM IMPUT] number of particles in the system
double L = 0.0f;                // [READ FROM INPUT] reduced box size
double P = 5.0f;                // reduced pressure
double sigma = 1.0f;            // unit of length = 1 = particle size
Particle *particles = NULL;     // [READ FROM INPUT] pointer to the table of particles data
// cell lists
// -----------------------
Particle **CLTable = NULL;
double sCell1D = 0.0f;
int nCell1D = 0;
int nCell = 0;





//    STRUCTURES
// ================
struct Particle
{
	double x;	        // reduced x coordinate 
	double y;	        // reduced y coordinate
        double z;               // reduced z coordinate
	double r;	        // reduced radius --- redundant with sigma ... might remove later on
        char type;              // particle type --- for visualization purposes only
        Particle *next;         // pointer to the next particle in CL
        Particle *prev;         // pointer to the previous particle in CL
        int CLindex;            // index of particle in current CL
};





//    FUNCTIONS
// ===============
int overlapCheck(Particle* one, Particle* two)
{
        double dx = one->x - two->x,
               dy = one->y - two->y,
               dz = one->z - two->z;
        dx = dx - L * rint(dx / L);
        dy = dy - L * rint(dy / L);
        dz = dz - L * rint(dz / L);
        if ((dx * dx + dy * dy + dz * dz) < ((one->r + two->r) * (one->r + two->r)))
                return 1;       
        return 0;
}



int overlapCheckGlobal(void)
{
        /* Function:     overlapCheckGlobal
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
                        overlap = (n != m) * overlapCheck(particles+n, particles+m);
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
        int a = particles[n].x / sCell1D;
        int b = particles[n].y / sCell1D;
        int c = particles[n].z / sCell1D;
        int index;
        Particle *current = NULL;
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
                                        if ((particles+n != current) && overlapCheck(particles+n, current))
                                                return 1;
                                        current = current->next;
                                }
                        }
                }
        }
        return 0;
}



void sanityCheck(void)
{
        /*
         * Function: sanityCheck
         * ---------------------
         * Ends program execution if there is overlap in the system
         */
        if (overlapCheckGlobal())
        {
                printf("Something is VERY wrong... There's overlap...\n");
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
                if (fscanf(initfile, "%*c%d%*c", &N) != 1)
                {
                        printf("\nERROR --- Wrong input file\n\n");
                        exit(0);
                }
                particles = malloc(N * sizeof(Particle));
                // Read value of box size
                if (fscanf(initfile, "%lf %*f %*f%*c", &L) != 1)
                {
                        printf("\nERROR --- Wrong input file\n\n");
                        exit(0);
                }
                // Populate table of particles data in standard units
                for (i = 0; i < N; i++)
                {
                        if (fscanf  ( initfile,
                                  "%c %lf %lf %lf %lf%*c",
                                  &particles[i].type,
                                  &particles[i].x,
                                  &particles[i].y,
                                  &particles[i].z,
                                  &particles[i].r
                                )
                            != 5
                           )
                        {
                                printf("\nERROR --- Wrong input file\n\n");
                                exit(0);
                        }
                        // Find larger radius to rescale everything to obtain
                        // the expected sigma = 1 for the (larger) particle size
                        if (particles[i].r > rTemp)
                                rTemp = particles[i].r;
                }
                fclose(initfile);
                // Rescale everything w.r.t. larger particle size to obtain the
                // expected sigma = 1, also converts to reduced units
                L /= (2.0f * rTemp);
                for (i = 0; i < N; i++)
                {
                        particles[i].x /= (2.0f * rTemp);
                        particles[i].y /= (2.0f * rTemp);
                        particles[i].z /= (2.0f * rTemp);
                        particles[i].r /= (2.0f * rTemp);
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
                                particles[i].type,
                                particles[i].x * sigma,
                                particles[i].y * sigma,
                                particles[i].z * sigma,
                                particles[i].r * sigma
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
        Particle *selected = particles+n;
        // Set up a backup of the moved particle in case of overlap
        Particle bufferParticle = {     selected->x,
                                        selected->y,
                                        selected->z
                                  };
        // Move randomly selected particle by abiding to PBC and nearest image
        // convention
        selected->x = fmod((selected->x + delta[0]) + 2 * L, L);
        selected->y = fmod((selected->y + delta[1]) + 2 * L, L);
        selected->z = fmod((selected->z + delta[2]) + 2 * L, L);
        // Update Cell List
        updateSingleCL(particles+n);
        // Check for overlap and move back particle if need be
        if (overlapCheckCL(n))
        {
                selected->x = bufferParticle.x;
                selected->y = bufferParticle.y;
                selected->z = bufferParticle.z;
                updateSingleCL(particles+n);
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
        // Rejects move if rule is not respected
        if (delta[1] > rule)
                return 0;
        // Scale the system according to the volume change
        for (int i = 0; i < N; i++)
        {
                particles[i].x *= ratio;
                particles[i].y *= ratio;
                particles[i].z *= ratio;
        }
        L *= ratio;
        // Checks for overlap in case of box shrinking and rescale cells
        // accordingly in each case
        if (ratio < 1)
        { // Box shrinks
                int overlap = 0;
                sCell1D *= ratio;
                if (sCell1D < sigma)
                // Cell size does not guarantee check for all possible
                // overlaps
                // Need to reduce cell size, thus number of cells
                        buildCL(3);
                else
                // Cell size is sufficient to check for all possible
                // overlaps
                // Only need to reduce cell size accordingly
                        buildCL(2);

                for (int i = 0; i < N; i++)
                {
                        if (overlapCheckCL(i))
                        {
                                overlap = 1;
                                break;
                        }
                }

                if (overlap)
                { // Overlap
                        for (int i = 0; i < N; i++)
                        {
                                particles[i].x /= ratio;
                                particles[i].y /= ratio;
                                particles[i].z /= ratio;
                        }
                        L /= ratio;
                        sCell1D /= ratio;
                        buildCL(3);
                        return 0;
                }
                else
                { // No overlap
                        return 1;
                }
        }
        else
        { // Box expands = no overlap
                //sCell1D *= ratio;
                //if (sCell1D > (2.0f * sigma))
                // Cell size is such that we check for overlap where there
                // cannot be any
                // Need to reduce cell size, thus increase number of cells,
                // thus free existing CL table and re-allocate enough memory
                //        buildCL(1);
                //else
                // Cell size is necessary (?) to check for all possible overlaps
                // Only need to increase cell size accordingly
                //        buildCL(2);
                buildCL(3);
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



int updateSingleCL(Particle *one)
{
        int index = one->CLindex;
        int u = one->x / sCell1D,
            v = one->y / sCell1D,
            w = one->z / sCell1D;
        // Remove particle from old CL
        if (CLTable[index] == one)
                CLTable[index] = one->next;
        if (one->next != NULL)
                one->next->prev = one->prev;
        if (one->prev != NULL)
                one->prev->next = one->next;
        // Fetches the index of the new CL
        index = retrieveIndex(u, v, w);
        one->CLindex = index;
        // Place particle on top of new CL
        one->prev = NULL;
        one->next = CLTable[index];
        if (CLTable[index] != NULL)
                CLTable[index]->prev = one;
        CLTable[index] = one;
        return 0;
}



int buildCL(int usecase)
{
        switch (usecase)
        {
                case 1:         // former initCL()
                // sCell goes back to ~sigma
                // nCell goes up, i.e. more cells
                // more memory is allocated to accomodate for increase in number of cells
                {
                        sCell1D = L / ((int) (L / sigma));
                        nCell1D = (int) (L / sCell1D);
                        nCell = nCell1D * nCell1D * nCell1D;
                        free(CLTable);
                        CLTable = malloc(nCell * sizeof(**CLTable));
                        if (CLTable == NULL)
                                exit(EXIT_FAILURE);
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        break;
                }
                case 2:         // former resizeCL() // get rid of this altogether
                // sigma < sCell < 2*sigma
                // same number of cells
                // no more memory allocated
                {
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        break;
                }
                case 3:         // former updateCL()
                // sCell goes back to ~sigma
                // nCell goes down, i.e. less cells
                // no more memory allocated, we juste use less of it
                {
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        sCell1D = L / ((int) (L / sigma));
                        nCell1D = (int) (L / sCell1D);
                        nCell = nCell1D * nCell1D * nCell1D;
                        break;
                }
                default:
                        return 1;
        }
        for (int i = 0; i < N; i++)
        {
                Particle *one = particles+i;
                int u = one->x / sCell1D;
                int v = one->y / sCell1D;
                int w = one->z / sCell1D;
                int index = retrieveIndex(u, v, w);
                one->CLindex = index;
                one->prev = NULL;
                one->next = CLTable[index];
                if (CLTable[index] != NULL)
                        CLTable[index]->prev = one;
                CLTable[index] = one;
        }
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
        Particle *current = NULL;
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



int snapshot(void)
{
        // output files cleansing
        char *snapcoords = "snap_coords.sph";
        FILE *snap1 = fopen(snapcoords, "w+");
        if (snap1 != NULL)
                fclose(snap1);
        char *snaprho = "snap_density.txt";
        FILE *snap2 = fopen(snaprho, "w+");
        if (snap2 != NULL)
                fclose(snap2);


        // writing snapshots
        writeCoords(snapcoords);
        snap2 = fopen(snaprho, "a");
        fprintf(snap2, "%lf\n", 6.0f * measurePF() / (M_PI * sigma * sigma * sigma));
        fclose(snap2);


        return 0;
}





//    MAIN
// ==========
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
        buildCL(1);
        writeCoords(outfilename);



        // 1st half of simuation time
        // Tuning of step sizes to obtain targetedAcceptanceRate acceptance rate
        stepTuneWrapper(MCCycle / 2, targetedAcceptanceRate);

        // 2nd half of simulation time
        // Measurements of observables --- here: packing fraction to obtain number density
        for (i = 0; i < (MCCycle / 2); i++)
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
                if (i % 999 == 0)
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
                if (i == 34729)
                        snapshot();
        }

        free(CLTable);
        free(particles);

        // Stop measuring CPU time
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        long seconds = end.tv_sec - begin.tv_sec;       
        long nanoseconds = end.tv_nsec - begin.tv_nsec;
        double elapsed = seconds + nanoseconds*1e-9;
        printf("\n\nProgram terminated correctly.\nElapsed CPU time:\t%.3fs.\n\n", elapsed);
        
        return 0;
}
