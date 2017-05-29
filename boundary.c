#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, const int* const flagField, const double * const wallVelocity, int xlength){
    int xMax=xlength+1;
    for (int x=0;x<=xMax;x++){
        for (int y=0;y<=xMax;y++){
            for (int z=0;z<=xMax;z++){
                if(flagField[z*totalLengthSq+y*totalLength+x]!=0){
                    setBoundaryCell(collideField, flagField, wallVelocity, xlength, xMax,x,y,z);
                }
            }
        }
    }
}

void setBoundaryCell(double *collideField, const int* const flagField, const double* const wallVelocity, int xlength, int xMax, int x, int y, int z){
    int currentCell=z*totalLengthSq+y*totalLength+x;
    double density;
    for(int i=0;i<=8;i++){ //this part would need to be modify for a difrent Q
        int x_neighbor=x+LATTICEVELOCITIES[i][0];
        int y_neighbor=y+LATTICEVELOCITIES[i][1];
        int z_neighbor=z+LATTICEVELOCITIES[i][2];
        int neighborCell=z_neighbor*totalLengthSq+y_neighbor*totalLength+x_neighbor;
        if(x_neighbor >= 0 && x_neighbor <= xMax && y_neighbor >= 0 && y_neighbor <= xMax && z_neighbor >= 0 && z_neighbor <= xMax && flagField[neighborCell] == 0){
            if(flagField[currentCell] == 1){
                collideField[Q*currentCell + i] = collideField[Q*currentCell + 18 - i];
            } else if(flagField[currentCell == 2]) {
                int c_dot_u = 0;
                for(int k = 0; k < D; k++)
                    c_dot_u += wallVelocity[k]*LATTICEVELOCITIES[i][k];
                computeDensity(&collideField[Q*neighborCell], &density);
                collideField[Q*currentCell + i] = collideField[Q*currentCell + 18 - i] + 2*LATTICEWEIGHTS[i]*density*c_dot_u / C_S_INVSQR; 
            }
        }
    }
    for(int i=10;i<=18;i++){ //this part would need to be modify for a difrent Q
        int x_neighbor=x+LATTICEVELOCITIES[i][0];
        int y_neighbor=y+LATTICEVELOCITIES[i][1];
        int z_neighbor=z+LATTICEVELOCITIES[i][2];
        int neighborCell=z_neighbor*totalLengthSq+y_neighbor*totalLength+x_neighbor;
        if(x_neighbor >= 0 && x_neighbor <= xMax && y_neighbor >= 0 && y_neighbor <= xMax && z_neighbor >= 0 && z_neighbor <= xMax && flagField[neighborCell] == 0){
            if(flagField[currentCell] == 1){
                collideField[Q*currentCell + i] = collideField[Q*currentCell + 18 - i];
            } else if(flagField[currentCell == 2]) {
                int c_dot_u = 0;
                for(int k = 0; k < D; k++)
                    c_dot_u += wallVelocity[k]*LATTICEVELOCITIES[i][k];
                computeDensity(&collideField[Q*neighborCell], &density);
                collideField[Q*currentCell + i] = collideField[Q*currentCell + 18 - i] + 2*LATTICEWEIGHTS[i]*density*c_dot_u / C_S_INVSQR; 
            }
        }
    }
}

int fluidNeighbor(int x, int y, int z, const int * const flagField, int xMax, int xlength, int totalLengthSq){
    if((x<0)||(x>xMax)||(y<0)||(y>xMax)||(z<0)||(z>xMax)){
        return 0;
    }
    if(flagField[z*totalLengthSq+y*totalLength+x]==0){
        return 1;
    }else{
        return 0;
    }
}
