#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
    for(int z = 1; z < xlength + 1; z++){
        for(int y = 1; y < xlength + 1; y++){
            for(int x = 1; x < xlength + 1; x++){
                int current_position = Q*(z*totalLengthSq + y*totalLength + x);
                for(int q = 0; q < Q; q++){
                    int x_neighbor = x - LATTICEVELOCITIES[q][0];
                    int y_neighbor = y - LATTICEVELOCITIES[q][1];
                    int z_neighbor = z - LATTICEVELOCITIES[q][2];
                    int neighbor_position = Q*(z_neighbor*totalLengthSq + y_neighbor*totalLength + x_neighbor);
                    streamField[current_position + q] = collideField[neighbor_position + q];
                }
            }
        }
    }
}
