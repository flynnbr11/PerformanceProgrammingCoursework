/*
 *  Simple molecular dynamics code.
 *  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
 *
 * This program implements:
 *     long range inverse square forces between particles. F = G * m1*m2 / r**2
 *     viscosity term     F = -u V
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * Coordinates are relative to a large central mass and the entire system is moving relative to the
 * viscous media.
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * This program was developed as part of a code optimisation course
 * and is therefore deliberately inefficient.
 *
 */
#include <math.h>
#include "coord.h"

void visc_force(int N,double *f, double *visc, double *vel);
void add_norm(int N,double *r, double *delta);
double force(double W, double delta, double r);
void wind_force(int N,double *f, double *visc, double vel);





void evolve(int count,double dt){
	int  step;
	int i,j,k,l;
	int collided;
	double Size;
/*
 * Loop over timesteps.
 */
  for(step = 1;step<=count;step++){
    printf("timestep %d\n",step);
    printf("collisions %d\n",collisions);

/* set the viscosity term in the force calculation */
		for(i=0;i<Nbody;i++){
			for(j=0;j<Ndim;j++){
				f[i][j] = -visc[i] * vel[i][j];
				f[i][j] = f[i][j] - visc[i]*wind[j];
			}
		}

/* calculate distance from central mass */
    for(k=0;k<Nbody;k++){
      r[k] = 0.0;
    }
    
		//  add_norm(Nbody,r,pos[i]);
	 for(k=0;k<Nbody;k++){
      for(i=0;i<Ndim;i++){
			  r[k] += (pos[k][i] * pos[k][i]);
			}
    }
    for(k=0;k<Nbody;k++){
      r[k] = sqrt(r[k]);
    }
   /* calculate central force */
    for(i=0;i<Nbody;i++){
		  for(l=0;l<Ndim;l++){
            f[i][l] = f[i][l] - 
            		((G*mass[i]*M_central*pos[i][l])/pow(r[i],3.0));
			}
		}
/* calculate pairwise separation of particles */
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
	      for(l=0;l<Ndim;l++){
          delta_pos[k][l] = pos[i][l] - pos[j][l];
        }
        k = k + 1;
      }
    }

/* calculate norm of seperation vector */
    for(k=0;k<Npair;k++){
      delta_r[k] = 0.0;
    }

//			  add_norm(Npair,delta_r,delta_pos[i]);
	  for(k=0;k<Npair;k++){
      for(i=0;i<Ndim;i++){
        delta_r[k] += (delta_pos[k][i] * delta_pos[k][i]);
      }
    }
    for(k=0;k<Npair;k++){
      delta_r[k] = sqrt(delta_r[k]);
    }

/*
* add pairwise forces.
*/
		k = 0;
		for(i=0;i<Nbody;i++){
		  for(j=i+1;j<Nbody;j++){
		    Size = radius[i] + radius[j];
		    collided=0;
		    for(l=0;l<Ndim;l++){
		      
		      if( delta_r[k] >= Size ){
		        f[i][l] = f[i][l] - G*mass[i]*mass[j]*delta_pos[k][l]/(pow(delta_r[k],3.0));
		        f[j][l] = f[j][l] + G*mass[i]*mass[j]*delta_pos[k][l]/(pow(delta_r[k],3.0));
		      }
		      else{
		        f[i][l] = f[i][l] + G*mass[i]*mass[j]*delta_pos[k][l]/(pow(delta_r[k],3.0));
		        f[j][l] = f[j][l] - G*mass[i]*mass[j]*delta_pos[k][l]/(pow(delta_r[k],3.0));
						collided=1;
		      }

				  if( collided == 1 && l==0){
					  collisions++;
					}
		  	}
		    k = k + 1;
		  }
		}

/* update positions */
    for(i=0;i<Nbody;i++){
      for(j=0;j<Ndim;j++){
        pos[i][j] = pos[i][j] + dt * vel[i][j];
        vel[i][j] = vel[i][j] + dt * (f[i][j]/mass[i]);
      }
    }

  }
}




