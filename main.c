/*
FILE  main.c

DOC   This program simulates a three-dimensional Lennard-Jones fluid.

LANG  C plus MPI message passing.

HIS   1) Parallel version originally written by University of
HIS      Southampton.
HIS   2) Adopted to run on the Intel iPSC/2 computer at Daresbury
HIS      Laboratory, and placed in the public domain through the
HIS      CCP5 programme.
HIS   3) Adopted to run on the Intel Paragon and enhanced by David W.
HIS      Walker at Oak Ridge National Laboratory, Tennessess, USA
HIS   4) Converted to use MPI message passing library by David W. Walker
HIS      at University of Wales Cardiff in July 1997.
HIS   5) Converted from Fortran to C in November 1997.

*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
//#include <sys/types.h>
#include <sys/time.h>
#include "moldyn.h"

double get_tick()
{
   struct timeval tim;
   gettimeofday(&tim, NULL);
   double sec = (double)tim.tv_sec;
   return (double)tim.tv_sec+((double)tim.tv_usec/1000000.0);
}

int main ( int argc, char *argv[])
{
   float sigma, rcut, dt, eqtemp, dens, boxlx, boxly, boxlz, sfx, sfy, sfz, sr6, vrcut, dvrcut, dvrc12, freex; 
   int nstep, nequil, iscale, nc, mx, my, mz, iprint;
   float *rx, *ry, *rz, *vx, *vy, *vz, *fx, *fy, *fz;
   float ace, acv, ack, acp, acesq, acvsq, acksq, acpsq, vg, kg, wg;
   int   *head, *list;
   int   natoms=0;
   int ierror;
   int jstart, step, itemp;
   float potential, virial, kinetic;
   float tmpx;
   int i, icell;
//   float tsum1, tsum2, psum, *vpArray;
 //  int j;
   double elapsedTime = (float)0;
   double startTime;
   double endTime;


   ierror = input_parameters (&sigma, &rcut, &dt, &eqtemp, &dens, &boxlx, &boxly, &boxlz, &sfx, &sfy, &sfz, &sr6, &vrcut, &dvrcut, &dvrc12, &freex, &nstep, &nequil, &iscale, &nc, &natoms, &mx, &my, &mz, &iprint);
   //printf ("\nReturned from input_parameters, natoms = %d\n", natoms);

   //printf("\nmx = %d, my = %d, mz = %d\n",mx,my,mz);
   rx = (float *)malloc(2*natoms*sizeof(float));
   ry = (float *)malloc(2*natoms*sizeof(float));
   rz = (float *)malloc(2*natoms*sizeof(float));
   vx = (float *)malloc(natoms*sizeof(float));
   vy = (float *)malloc(natoms*sizeof(float));
   vz = (float *)malloc(natoms*sizeof(float));
   fx = (float *)malloc(natoms*sizeof(float));
   fy = (float *)malloc(natoms*sizeof(float));
   fz = (float *)malloc(natoms*sizeof(float));
   list = (int *)malloc(2*natoms*sizeof(int));
   head= (int *)malloc((mx+2)*(my+2)*(mz+2)*sizeof(int));
//   vpArray = (float *)malloc(2*natoms*sizeof(float));
   assert(rx);
   assert(ry);
   assert(rz);
   assert(vx);
   assert(vy);
   assert(vz);
   assert(fx);
   assert(fy);
   assert(fz);
   assert(list);
   assert(head);
//   assert(vpArray);
  // printf ("\nFinished allocating memory\n");

   initialise_particles (rx, ry, rz, vx, vy, vz, nc);
 //  printf ("\nReturned from initialise_particles\n");

   loop_initialise(&ace, &acv, &ack, &acp, &acesq, &acvsq, &acksq, &acpsq, sigma, rcut, dt);
//   printf ("\nReturned from loop_initialise\n");

//   output_particles(rx,ry,rz,vx,vy,vz,fx,fy,fz,0);
      movout (rx, ry, rz, vx, vy, vz, sfx, sfy, sfz, head, list, mx, my, mz, natoms);
//   printf ("\nReturned from movout\n");
   //   check_cells(rx, ry, rz, head, list, mx, my, mz, natoms,0,0);
      startTime = get_tick();
      force (&potential, &virial, rx, ry, rz, fx, fy, fz, sigma, rcut, vrcut, dvrc12, dvrcut, head, list, mx, my, mz, natoms,0, sfx, sfy, sfz);
     endTime = get_tick();
     elapsedTime += (endTime - startTime);
//   printf ("\nReturned from force: potential = %f, virial = %f, kinetic = %f\n",potential, virial, kinetic);
//   output_particles(rx,ry,rz,vx,vy,vz,fx,fy,fz,10);
//   exit(0);


   for(step=1;step<=nstep;step++){
     // if(step>=85)printf ("\nStarted step %d\n",step);
      movea (rx, ry, rz, vx, vy, vz, fx, fy, fz, dt, natoms);
//      check_cells(rx, ry, rz, head, list, mx, my, mz, natoms,step,step);
    //  if(step>85)printf ("\nReturned from movea\n");
      movout (rx, ry, rz, vx, vy, vz, sfx, sfy, sfz, head, list, mx, my, mz, natoms);
   //  if(step>85) printf ("\nReturned from movout\n");
  //    check_cells(rx, ry, rz, head, list, mx, my, mz, natoms,step,step);
      startTime = get_tick();
      force (&potential, &virial, rx, ry, rz, fx, fy, fz, sigma, rcut, vrcut, dvrc12, dvrcut, head, list, mx, my, mz, natoms, step, sfx, sfy, sfz);
     endTime = get_tick();
      elapsedTime += (endTime - startTime);
//      printf ("\nReturned from force: potential = %f, virial = %f, kinetic = %f\n",potential, virial, kinetic);
 //     fflush(stdout);
      moveb (&kinetic, vx, vy, vz, fx, fy, fz, dt, natoms);
 //     check_cells(rx, ry, rz, head, list, mx, my, mz, natoms,step,step);
  // if(step>85)   printf ("\nReturned from moveb: potential = %f, virial = %f, kinetic = %f\n",potential, virial, kinetic);
//      potential = 0.0;
//      virial = 0.0;
//      for (i=0;i<natoms;i+=512){
//          tsum1 = 0.0;
//	  tsum2    = 0.0;
//	  for(j=0;j<512;j++) tsum1 += vpArray[i+j];
//	  for(j=0;j<512;j++) tsum2 += vpArray[natoms+i+j];
//	  potential += tsum1;
//	  virial += tsum2;
//      }
      sum_energies (potential, kinetic, virial, &vg, &wg, &kg);
      hloop (kinetic, step, vg, wg, kg, freex, dens, sigma, eqtemp, &tmpx, &ace, &acv, &ack, &acp, &acesq, &acvsq, &acksq, &acpsq, vx, vy, vz, iscale, iprint, nequil, natoms);
   }

   tidyup (ace, ack, acv, acp, acesq, acksq, acvsq, acpsq, nstep, nequil);
   printf("\n%f seconds have elapsed\n", elapsedTime);

   
   return 0;
}
