#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
    if(argc!=2){
        printf("No file name given");
        return -1;
    }
    double velocityWallx;
    double velocityWally;
    double velocityWallz;
    
    printf("Reading parameters from file %s\n",argv[1]);
    READ_INT(argv[1], *xlength);
    READ_DOUBLE(argv[1], *tau);
    READ_DOUBLE(argv[1], velocityWallx);
    READ_DOUBLE(argv[1], velocityWally);
    READ_DOUBLE(argv[1], velocityWallz);
    READ_INT(argv[1], *timesteps);
    READ_INT(argv[1], *timestepsPerPlotting);

    velocityWall[0]=velocityWallx;
    velocityWall[1]=velocityWally;
    velocityWall[2]=velocityWallz;
    return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength)
{
     //defining the flags
     //const int FLUID = 0;
     const int NO_SLIP = 1;
     const int MOVING_WALL = 2;
     
/***********************************************************************/   
/***********************************************************************/  
    //setting all fields on z=0 and z = xlength+2 (top and bottom walls lying on x-y plane)
    for(int x=0; x<totalLength; ++x)
    {
        for(int y=0; y<totalLength; ++y)
        {    
            //initialiing NO_SLIP flagField at plane z=0 
            flagField[0*totalLengthSq + y*totalLength + x] = NO_SLIP;

            //initalizing MOVING_WALL flagField at plane z=xlength+2
            flagField[(xlength+1)*totalLengthSq + y*totalLength + x] = MOVING_WALL;
            
        } //end of y loop
        
    }//end of x loop

/***********************************************************************/   
/***********************************************************************/   


/***********************************************************************/       
/***********************************************************************/   
    //setting all fields on y=0 and y=xlength+1 (side walls lying in x-z plane)
    for(int x=0; x<totalLength; ++x)
    {
        for(int z=1; z<=xlength; ++z)
        {    
            //initialiing flagField at y=0
            flagField[z*totalLengthSq + 0*totalLength + x] = NO_SLIP;

            //initalizing flagField at y=xlength+1
            flagField[z*totalLengthSq + (xlength+1)*totalLength + x] = NO_SLIP;
       
        } //end of z loop

    } //end of x loop

/***********************************************************************/
/***********************************************************************/     




/***********************************************************************/   
/***********************************************************************/
    //setting all fields on x=0 and x=xlength+1 (side walls lying in y-z plane)
    for(int y=1; y<=xlength; ++y)
    {
        for(int z=1; z<=xlength; ++z)
        {    
            //initialiing flagField at x=0
            flagField[z*totalLengthSq + y*totalLength + 0] = NO_SLIP;

            //initalizing flagField at x=xlength+1
            flagField[z*totalLengthSq + y*totalLength+ xlength+1] = NO_SLIP;
       
        } //end of z loop

    } //end of x loop

/***********************************************************************/  
/***********************************************************************/ 

//initializing streamField and collideField array
/***********************************************************************/  
/***********************************************************************/ 
for(int x=0; x<=xlength+1; ++x)
{
    for(int y=0; y<=xlength+1; ++y)
    {
        for(int z=0; z<=xlength+1; ++z)
        {
            for(int j=0; j<Q; ++j)
            {
                collideField[Q*(z*totalLengthSq + y*totalLength + x) + j] = LATTICEWEIGHTS[j];
                streamField[Q*(z*totalLengthSq + y*totalLength + x) + j] =  LATTICEWEIGHTS[j];
            } //end of j loop
        }//end of z loop
    }//end of y loop
}//end of x loop
/***********************************************************************/  
/***********************************************************************/    

}
