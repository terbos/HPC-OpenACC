#include <openacc.h>
#include <math.h>
#include <stdio.h>

#include "moldyn.h"

void force (float *pval, float *vval, float *rx, float *ry, float *rz, float *fx, float *fy, float *fz, float sigma, float rcut, float vrcut, float dvrc12, float dvrcut, int *head, int *list, int mx, int my, int mz, int natoms, int step, float sfx, float sfy, float sfz)
{
   float sigsq  = sigma*sigma;
   float rcutsq = rcut*rcut;
   float potential = 0.0;
   float virial    = 0.0;
#pragma acc kernels loop copyin(rx[0:2*natoms],ry[0:2*natoms],rz[0:2*natoms],list[0:2*natoms],head[0:(mx+2)*(my+2)*(mz+2)]) copyout(fx[0:natoms],fy[0:natoms],fz[0:natoms]) reduction(+:potential,virial)
   for(int i=0;i<natoms;++i){ // parallelize this loop
	 float rxi = rx[i];
	 float ryi = ry[i];
	 float rzi = rz[i];
	 float fxi = 0.0;
	 float fyi = 0.0;
	 float fzi = 0.0;
	 int xi = (int)((rxi+0.5)/sfx) + 1;
	 int yi = (int)((ryi+0.5)/sfy) + 1;
	 int zi = (int)((rzi+0.5)/sfz) + 1;
         if(xi > mx) xi = mx;
         if(yi > my) yi = my;
         if(zi > mz) zi = mz;
	 int icell = xi + (mx+2)*(yi+zi*(my+2));
         for (int ix=-1;ix<=1;ix++) // these inner loops not parallelized
             for (int jx=-1;jx<=1;jx++)
                 for (int kx=-1;kx<=1;kx++){
		     int xcell = ix+xi;
		     int ycell = jx+yi;
		     int zcell = kx+zi;
                     int jcell = xcell + (mx+2)*(ycell+(my+2)*zcell);
                     int j = head[jcell];
                     while (j>=0) {
                         if (j!=i) {
                             float rxij = rxi - rx[j];
                             float ryij = ryi - ry[j];
                             float rzij = rzi - rz[j];
                             float rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
                             if (rijsq < rcutsq) {
				 float rij = (float) sqrt ((double)rijsq);
                                 float sr2 = sigsq/rijsq;
                                 float sr6 = sr2*sr2*sr2;
                                 float vij = sr6*(sr6-1.0) - vrcut - dvrc12*(rij-rcut);
                                 float wij = sr6*(sr6-0.5) + dvrcut*rij;
                                 float fij = wij/rijsq;
                                 float fxij = fij*rxij;
                                 float fyij = fij*ryij;
                                 float fzij = fij*rzij;
                                 wij *= 0.5;
                                 vij *= 0.5;
                                 potential += vij;
                                 virial    += wij;
                                 fxi       += fxij;
                                 fyi       += fyij;
                                 fzi       += fzij;
                             }
			 }
                         j = list[j];
                      }
	         }
         *(fx+i) = 48.0*fxi;
         *(fy+i) = 48.0*fyi;
         *(fz+i) = 48.0*fzi;
   } // end of parallelized loop
   potential *= 4.0;
   virial    *= 48.0/3.0;
   *pval = potential;
   *vval = virial;
}
