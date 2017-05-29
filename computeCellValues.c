#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
    *density = 0;
    for(int i = 0; i < Q; i++)
        *density += currentCell[i];
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    double rhs[3] = {0}; // Right-hand side of eq. 10
    double inv_density = 1 / (*density);

    for(int q = 0; q < Q; q++)
        for(int d = 0; d < D; d++)
            rhs[d] += currentCell[q]*LATTICEVELOCITIES[q][d];

    for(int d = 0; d < D; d++)
        velocity[d] = inv_density*rhs[d];
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
    double velocity_dot_velocity = dot_product(velocity, velocity, D);
    double vel_sound_square = C_S * C_S;
    double vel_sound_fourth = vel_sound_square * vel_sound_square;
    double constant = velocity_dot_velocity / (2 * vel_sound_square);
    double direction_dot_velocity = 0;
    for(int q = 0; q < Q; q++){
        for(int d = 0; d < D; d++)
            direction_dot_velocity += velocity[d]*LATTICEVELOCITIES[q][d];
        // Eq. 11 Worksheet
        feq[q] = LATTICEWEIGHTS[q] * (*density) * (1 + direction_dot_velocity / vel_sound_square +
                                 (direction_dot_velocity*direction_dot_velocity) / (2 * vel_sound_fourth) -
                                 constant);
    }
}

double dot_product(const double * const u, const double * const v, int size){
    double result = 0;
    for(int i = 0; i < size; i++)
        result += u[i]*v[i];
    return result;
}

double dot_product_int(const int * const u, const double * const v, int size){
    double result = 0;
    for(int i = 0; i < size; i++)
        result += u[i]*v[i];
    return result;
}
