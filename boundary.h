#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, const int* const flagField, const double * const wallVelocity, int xlength);

void setBoundaryCell(double *collideField, const int* const flagField, const double* const wallVelocity, int xlength, int xMax, int x, int y, int z);

int fluidNeighbor(int x, int y, int z, const int * const flagField, int xMax, int xlength, int xlength_square);
#endif

