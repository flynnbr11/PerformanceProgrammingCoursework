/* $Id: coord.h,v 1.2 2002/01/08 12:32:48 spb Exp spb $
 * This file defines static arrays that contains the primary coordinates
 * of the particles,
 *
 *  Nbody	Number of particles
 *  Npair	Number of particle pairs
 *  pos		Position of the particles
 *  r           distance of partice from central mass 
 * vel		velocity of the particles
 *  f		Forces acting on each particle
 *  visc		viscosity coefficient for each particle
 *  mass		mass of each particle
 *  delta_pos	seperation vector for each particle pair
 *  delta_r		seperation for each particle pair
 */

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif
#define Nbody 4*1024
#define  Npair ((Nbody*(Nbody-1))/2)

enum{ Xcoord=0, Ycoord, Zcoord, Ndim };

DEF double pos[Nbody][Ndim], vel[Nbody][Ndim];
DEF double f[Nbody][Ndim], visc[Nbody], mass[Nbody], radius[Nbody];
DEF double delta_pos[Nbody*Nbody][3];
DEF double r[Nbody];
DEF double delta_r[Nbody*Nbody];
DEF double wind[Ndim];
DEF int collisions;


#define G 2.0
#define M_central 1000.0

void evolve(int Nstep, double dt);