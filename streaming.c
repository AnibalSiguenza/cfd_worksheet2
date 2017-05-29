#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
    for(int x=1;x<=xlength;x++){
        for(int y=1;y<=xlength;y++){
            for(int z=1;z<=xlength;z++){
                for (int i=0;i<Q;i++){
                    streamField[Q*(z*totalLengthSq+y*totalLength+x)+i]=collideField[Q*((z-LATTICEVELOCITIES[i][2])*totalLengthSq\
                            +(y-LATTICEVELOCITIES[i][1])*totalLength+(x-LATTICEVELOCITIES[i][0]))+i]; //Corresponding to equation 14
                }
            }
        }
    }
}
