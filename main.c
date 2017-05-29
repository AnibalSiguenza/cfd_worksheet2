#ifndef _MAIN_C_
#define _MAIN_C_

#include <math.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){
    double *collideField=NULL;
    double *streamField=NULL;
    int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    int timesteps;
    int timestepsPerPlotting;

    if(readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv) != 0){
        printf("Error during the lecture of parameters (mayble no file indicated?)\n");
        return -1 ;
     }

    totalLength=xlength+2;
    totalLengthSq=totalLength*totalLength;
    const int num_cells = (int)pow(totalLength, D);
    const int size_field = Q * num_cells;
    collideField = (double*)calloc(size_field,sizeof(double));
    streamField = (double*)calloc(size_field,sizeof(double));
    flagField = (int*)calloc(num_cells,sizeof(int));

    initialiseFields(collideField, streamField, flagField, xlength);
    treatBoundary(collideField, flagField, velocityWall, xlength);
    writeVtkOutput(collideField, flagField, "simulation_data/init_values", 0, xlength);
    writeVtkOutputFlags(collideField, flagField, "simulation_data/init_values", 0, xlength);

    for(int t = 0; t <= timesteps; t++){
        double *swap=NULL;
        doStreaming(collideField, streamField, flagField, xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField, flagField, &tau, xlength);
        treatBoundary(collideField, flagField, velocityWall, xlength);
        if (t % timestepsPerPlotting == 0){
            writeVtkOutput(collideField, flagField, argv[1], t, xlength);
        }
    }

    return 0;
}

#endif
