#ifndef _MAIN_C_
#define _MAIN_C_

#include <math.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"
#include <string.h>

void sanity_test(int t, int xlength, double const * const collideField, double const * const streamField, char * message){
    int x=xlength/2;
    int y=xlength/2;
    int z[5]={0,1,xlength/2,xlength,xlength+1};
    double density=0;
    double velocity[3];


    for(int i=0;i<5;i++){
        printf("[%s]Sanity test at t=%i x=%i y =%i z=%i\n",message,t,x,y,z[i]);
        const double * const collcurrentCell = &collideField[Q*(z[i]*totalLengthSq + y*totalLength + x)];
        const double * const stremcurrentCell = &streamField[Q*(z[i]*totalLengthSq + y*totalLength + x)];
        computeDensity(collcurrentCell , &density);
        computeVelocity(collcurrentCell , &density, velocity);
        printf("collideField density=%f velocity=[%f,%f,%f]\n",density,velocity[0],velocity[1],velocity[2]);
        printf("f=");
        for(int j=0;j<Q;j++){
            printf("%f,",collcurrentCell[j]);
        }
        computeDensity(stremcurrentCell, &density);
        computeVelocity(stremcurrentCell, &density, velocity);
        printf("\n");
        printf("streamField density=%f velocity=[%f,%f,%f]\n",density,velocity[0],velocity[1],velocity[2]);
        printf("f=");
        for(int i=0;i<Q;i++){
            printf("%f,",stremcurrentCell [i]);
        }
        printf("\n");
        if(density<0.8||density>1.3){
            ERROR("Bad density ");
        }

    }

    printf("\n");
}


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
    //writeVtkOutput(collideField, flagField, "simulation_data/init_values", 0, xlength);
    //writeVtkOutputFlags(collideField, flagField, "simulation_data/init_values", 0, xlength);
    //sanity_test(0, xlength,collideField, streamField,"before for");
    for(int t = 0; t <= timesteps; t++){
        double *swap=NULL;
        doStreaming(collideField, streamField, flagField, xlength);
        //sanity_test(t, xlength,collideField, streamField,"after streaming");
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField, flagField, &tau, xlength);
        //sanity_test(t, xlength,collideField, streamField,"after collision");
        treatBoundary(collideField, flagField, velocityWall, xlength);
        if (t % timestepsPerPlotting == 0){
            writeVtkOutput(collideField, flagField, argv[1], t, xlength);
        }
    }

    return 0;
}

#endif
