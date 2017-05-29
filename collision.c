#include <math.h>
#include "collision.h"
#include "LBDefinitions.h"

#define FLUID 0

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    // Apply BGK update rule from Eq. 14
    for(int q = 0; q < Q; q++)
        currentCell[q] = currentCell[q] - (currentCell[q] - feq[q]) / (*tau);
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
    //const int num_cells = (int)pow(totalLength, D);
    for(int z = 1; z <= xlength; z++){
        for(int y = 1; y <= xlength; y++ ){
            for(int x = 1; x <= xlength; x++){
                double density = 0;
                double velocity[3] = {0};
                double feq[19] = {0};
                double *currentCell = &collideField[Q*(z*totalLengthSq + y*totalLength + x)];
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
}

