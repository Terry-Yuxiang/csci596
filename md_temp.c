/*******************************************************************************
Molecular dynamics (MD) simulation with the Lennard-Jones potential.

USAGE

%cc -o md md.c -lm
%md < md.in (see md.h for the input-file format)
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "md_temp.h"
#include <time.h>

#define GRID_SIZE 10

// Initialize output
void InitOutput() {
 	FILE *totalFile = fopen("total_res.txt", "w");
    FILE *layersFile = fopen("layers_res.txt", "w");

    // Add column headers for total_res.txt
    fprintf(totalFile, "Time Temperature Pressure PotentialEnergy TotalEnergy\n");

    // Dynamically generate column headers for layers_res.txt
    fprintf(layersFile, "Time");
    for (int i = 0; i < GRID_SIZE; i++) {
        fprintf(layersFile, " Layer%d", i + 1);
    }
    fprintf(layersFile, "\n");

    fclose(totalFile);
    fclose(layersFile);
}

int main(int argc, char **argv) {
	double cpu, cpu1, cpu2;
	InitParams();  
	InitConf(); 
	ComputeAccel(); /* Computes initial accelerations */ 
	InitOutput();
	cpu1 = ((double) clock())/CLOCKS_PER_SEC;
	for (stepCount=1; stepCount<=StepLimit; stepCount++) {
		SingleStep(); 
		if (stepCount%StepAvg == 0) EvalProps();
	}
	cpu2 = ((double) clock())/CLOCKS_PER_SEC;
	cpu = cpu2-cpu1;
	printf("Execution time (s) = %le\n",cpu);
	return 0;
}

/*----------------------------------------------------------------------------*/
void InitParams() {
/*------------------------------------------------------------------------------
	Initializes parameters.
------------------------------------------------------------------------------*/
	int k;
	double rr,ri2,ri6,r1;

	/* Reads control parameters */
	scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
	scanf("%le",&Density);
	scanf("%le",&InitTemp);
	scanf("%le",&DeltaT);
	scanf("%d",&StepLimit);
	scanf("%d",&StepAvg);

	/* Computes basic parameters */
	DeltaTH = 0.5*DeltaT;
	for (k=0; k<3; k++) {
		Region[k] = InitUcell[k]/pow(Density/4.0,1.0/3.0);
		RegionH[k] = 0.5*Region[k];
	}

	/* Constants for potential truncation */
	rr = RCUT*RCUT; ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1=sqrt(rr);
	Uc = 4.0*ri6*(ri6 - 1.0);
	Duc = -48.0*ri6*(ri6 - 0.5)/r1;
}

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to face-centered cubic (fcc) lattice positions.  
	rv are initialized with a random velocity corresponding to Temperature.  
------------------------------------------------------------------------------*/
	double c[3],gap[3],e[3],vSum[3],vMag;
	int j,n,k,nX,nY,nZ;
	double seed;
	/* FCC atoms in the original unit cell */
	double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
	                         {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}; 

	/* Sets up a face-centered cubic (fcc) lattice */
	for (k=0; k<3; k++) gap[k] = Region[k]/InitUcell[k];
	nAtom = 0;
	for (nZ=0; nZ<InitUcell[2]; nZ++) {
		c[2] = nZ*gap[2];
		for (nY=0; nY<InitUcell[1]; nY++) {
			c[1] = nY*gap[1];
			for (nX=0; nX<InitUcell[0]; nX++) {
				c[0] = nX*gap[0];
				for (j=0; j<4; j++) {
					for (k=0; k<3; k++)
						r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
					++nAtom;
				}
			}
		}
	}

	/* Generates random velocities */
	seed = 13597.0;
	vMag = sqrt(3*InitTemp);
	for(k=0; k<3; k++) vSum[k] = 0.0;
	for(n=0; n<nAtom; n++) {

		// Determine the layer index of the particle
		int zIdx = (int)(r[n][2] / (Region[2] / InitUcell[2])); // Layer index based on z-coordinate

		// Adjust the velocity magnitude based on the layer index
		if (zIdx == 0) { 
			// For the first layer, use higher temperature
			vMag = sqrt(3 * 10.0 * InitTemp); // Higher temperature is twice the initial temperature
		} else {
			// For other layers, use the default temperature
			vMag = sqrt(3 * InitTemp);
		}

		RandVec3(e,&seed);
		for (k=0; k<3; k++) {
			rv[n][k] = vMag*e[k];
			vSum[k] = vSum[k] + rv[n][k];
		}
	}
	/* Makes the total momentum zero */
	for (k=0; k<3; k++) vSum[k] = vSum[k]/nAtom;
	for (n=0; n<nAtom; n++) for(k=0; k<3; k++) rv[n][k] = rv[n][k] - vSum[k];
}

/*----------------------------------------------------------------------------*/
void ComputeAccel() {
/*------------------------------------------------------------------------------
	Acceleration, ra, are computed as a function of atomic coordinates, r,
	using the Lennard-Jones potential.  The sum of atomic potential energies,
	potEnergy, is also computed.   
------------------------------------------------------------------------------*/
	double dr[3],f,fcVal,rrCut,rr,ri2,ri6,r1;
	int j1,j2,n,k;

	rrCut = RCUT*RCUT;
	for (n=0; n<nAtom; n++) for (k=0; k<3; k++) ra[n][k] = 0.0;
	potEnergy = 0.0;

	/* Doubly-nested loop over atomic pairs */
	for (j1=0; j1<nAtom-1; j1++) {
		for (j2=j1+1; j2<nAtom; j2++) {
			/* Computes the squared atomic distance */
			for (rr=0.0, k=0; k<3; k++) {
				dr[k] = r[j1][k] - r[j2][k];
				/* Chooses the nearest image */
				dr[k] = dr[k] - SignR(RegionH[k],dr[k]-RegionH[k])
				              - SignR(RegionH[k],dr[k]+RegionH[k]);
				rr = rr + dr[k]*dr[k];
			}
			/* Computes acceleration & potential within the cut-off distance */
			if (rr < rrCut) {
				ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1 = sqrt(rr);
				fcVal = 48.0*ri2*ri6*(ri6-0.5) + Duc/r1;
				for (k=0; k<3; k++) {
					f = fcVal*dr[k];
					ra[j1][k] = ra[j1][k] + f;
					ra[j2][k] = ra[j2][k] - f;
				}
				potEnergy = potEnergy + 4.0*ri6*(ri6-1.0) - Uc - Duc*(r1-RCUT);
			} 
		} 
	}
}

/*----------------------------------------------------------------------------*/
void SingleStep() {
/*------------------------------------------------------------------------------
	r & rv are propagated by DeltaT in time using the velocity-Verlet method.
------------------------------------------------------------------------------*/
	int n,k;

	HalfKick(); /* First half kick to obtain v(t+Dt/2) */
	for (n=0; n<nAtom; n++) /* Update atomic coordinates to r(t+Dt) */
		for (k=0; k<3; k++) r[n][k] = r[n][k] + DeltaT*rv[n][k];
	ApplyBoundaryCond();
	ComputeAccel(); /* Computes new accelerations, a(t+Dt) */
	HalfKick(); /* Second half kick to obtain v(t+Dt) */
}

/*----------------------------------------------------------------------------*/
void HalfKick() {
/*------------------------------------------------------------------------------
	Accelerates atomic velocities, rv, by half the time step.
------------------------------------------------------------------------------*/
	int n,k;
	for (n=0; n<nAtom; n++)
		for (k=0; k<3; k++) rv[n][k] = rv[n][k] + DeltaTH*ra[n][k];
}

/*----------------------------------------------------------------------------*/
void ApplyBoundaryCond() {
/*------------------------------------------------------------------------------
	Applies periodic boundary conditions to atomic coordinates.
------------------------------------------------------------------------------*/
	int n,k;
	for (n=0; n<nAtom; n++) 
		for (k=0; k<3; k++) 
			r[n][k] = r[n][k] - SignR(RegionH[k],r[n][k])
			                  - SignR(RegionH[k],r[n][k]-Region[k]);
}

void EvalProps() {
/*------------------------------------------------------------------------------
    The layers evaluate physical properties: kinetic energy, temperature, etc., at each layer.
------------------------------------------------------------------------------*/
    double layerKinEnergy[GRID_SIZE] = {0.0};  // Sum of kinetic energy in each layer
    double layerCount[GRID_SIZE] = {0.0};      // The number of particles in each layer
    double vv;
    int n, k;

    // Open files
    FILE *totalFile = fopen("total_res.txt", "a");
    FILE *layersFile = fopen("layers_res.txt", "a");

    // Calculate the kinetic energy of each particle and categorize it into the corresponding layer
    for (n = 0; n < nAtom; n++) {
        vv = 0.0;
        for (k = 0; k < 3; k++)
            vv += rv[n][k] * rv[n][k]; // Kinetic energy is proportional to the square of velocity

        double kineticEnergy = 0.5 * vv; // Kinetic energy of a single particle
        int zIdx = (int)(r[n][2] / (Region[2] / GRID_SIZE)); // Determine the layer based on the z-coordinate

        if (zIdx >= 0 && zIdx < GRID_SIZE) {
            layerKinEnergy[zIdx] += kineticEnergy; // Accumulate kinetic energy for the corresponding layer
            layerCount[zIdx]++; // Count the number of particles in this layer
        }
    }

    // Output the average kinetic energy and temperature of each layer to layers_res.txt
    fprintf(layersFile, "%9.6f", stepCount * DeltaT);
    for (int i = 0; i < GRID_SIZE; i++) {
        double avgKinEnergy = layerCount[i] > 0 ? layerKinEnergy[i] / layerCount[i] : 0.0;
        double temperature = avgKinEnergy * 2.0 / 3.0; // Calculate temperature based on kinetic energy
        fprintf(layersFile, " %9.6f", temperature); // Save the temperature of each layer
    }
    fprintf(layersFile, "\n"); // Newline

    // Calculate total kinetic energy and temperature
    kinEnergy = 0.0;
    for (int i = 0; i < GRID_SIZE; i++) {
        kinEnergy += layerKinEnergy[i];
    }
    kinEnergy /= nAtom; // Average kinetic energy
    potEnergy /= nAtom; // Average potential energy
    totEnergy = kinEnergy + potEnergy;
    temperature = kinEnergy * 2.0 / 3.0;

    // Save total properties to total_res.txt
    fprintf(totalFile, "%9.6f %9.6f %9.6f %9.6f\n",
            stepCount * DeltaT, temperature, potEnergy, totEnergy);

    // Close files
    fclose(totalFile);
    fclose(layersFile);
}