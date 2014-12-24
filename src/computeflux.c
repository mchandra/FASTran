#include "timestepper.h"

void computeFlux(struct timeStepper ts[ARRAY_ARGS 1],
                 REAL ***primPetsc, REAL ***divFlux)
{

  Vec primPetscVecOldLocal;
  DMGetLocalVector(ts->dmda, &primPetscVecOldLocal);

    /* Compute Div(flux) at t=n */
  DMGlobalToLocalBegin(ts->dmda, 
                       ts->primPetscVecOld,
                       INSERT_VALUES,
                       primPetscVecOldLocal);
  DMGlobalToLocalEnd(ts->dmda,
                     ts->primPetscVecOld,
                     INSERT_VALUES,
                     primPetscVecOldLocal);

  REAL ***primPetscOld;
  DMDAVecGetArrayDOF(ts->dmda, primPetscVecOldLocal, &primPetscOld); 

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  for(int j=ts->X2Start; j<ts->X2Start+ts->X2Size; j++)
  {
    for (int i=ts->X1Start; i<ts->X1Start+ts->X1Size; i++)
    {
      REAL X1RightEdge  = X1_ZONE_RIGHT_EDGE(i);
      REAL X1LeftEdge   = X1_ZONE_LEFT_EDGE(i);
      REAL X1Center     = X1_ZONE_CENTER(i);
      
      REAL X2Center     = X1_ZONE_CENTER(j);
      REAL X2BottomEdge = X2_ZONE_BOTTOM_EDGE(j);
      REAL X2TopEdge    = X2_ZONE_TOP_EDGE(j);

      REAL T_i_j              = primPetsc[j][i][0];
      REAL T_iPlus1_j         = primPetsc[j][i+1][0];
      REAL T_iMinus1_j        = primPetsc[j][i-1][0];
      REAL T_i_jPlus1         = primPetsc[j+1][i][0];
      REAL T_i_jMinus1        = primPetsc[j-1][i][0];
      REAL T_iPlus1_jPlus1    = primPetsc[j+1][i+1][0];
      REAL T_iMinus1_jPlus1   = primPetsc[j+1][i-1][0];
      REAL T_iPlus1_jMinus1   = primPetsc[j-1][i+1][0];
      REAL T_iMinus1_jMinus1  = primPetsc[j-1][i-1][0];

      REAL T_old_i_j              = primPetscOld[j][i][0];
      REAL T_old_iPlus1_j         = primPetscOld[j][i+1][0];
      REAL T_old_iMinus1_j        = primPetscOld[j][i-1][0];
      REAL T_old_i_jPlus1         = primPetscOld[j+1][i][0];
      REAL T_old_i_jMinus1        = primPetscOld[j-1][i][0];
      REAL T_old_iPlus1_jPlus1    = primPetscOld[j+1][i+1][0];
      REAL T_old_iMinus1_jPlus1   = primPetscOld[j+1][i-1][0];
      REAL T_old_iPlus1_jMinus1   = primPetscOld[j-1][i+1][0];
      REAL T_old_iMinus1_jMinus1  = primPetscOld[j-1][i-1][0];
      

      #if (TIME_STEPPING==FULLY_IMPLICIT)
      divFlux[j][i][0] = 
     -(
        (1./DX1)
      * (   D(X1RightEdge, X2Center)*bx(X1RightEdge, X2Center)
          * (
               (bx(X1RightEdge, X2Center) * (T_iPlus1_j - T_i_j)/DX1)
             + (  by(X1RightEdge, X2Center)/DX2 
                * limiter4( T_i_jPlus1 - T_i_j, T_iPlus1_jPlus1 - T_iPlus1_j,
                            T_i_j - T_i_jMinus1, T_iPlus1_j - T_iPlus1_jMinus1
                          )
               )
            )
         -  D(X1LeftEdge, X2Center)*bx(X1LeftEdge, X2Center)
          * (
               (bx(X1LeftEdge, X2Center) * (T_i_j - T_iMinus1_j)/DX1)
             + (  by(X1LeftEdge, X2Center)/DX2
                * limiter4(T_iMinus1_jPlus1 - T_iMinus1_j, T_i_jPlus1 - T_i_j,
                           T_iMinus1_j - T_iMinus1_jMinus1, T_i_j - T_i_jMinus1
                          )
                   
               )   
            )
        )
      + 
        (1./DX2)
      * (   D(X1Center, X2TopEdge)*by(X1Center, X2TopEdge)
          * (
               (by(X1Center, X2TopEdge) * (T_i_jPlus1 - T_i_j)/DX2)
             + (  bx(X1Center, X2TopEdge)/DX1
                * limiter4(T_iPlus1_j - T_i_j, T_iPlus1_jPlus1 - T_i_jPlus1,
                           T_i_j - T_iMinus1_j, T_i_jPlus1 - T_iMinus1_jPlus1
                          )
               )
            )
        -   D(X1Center, X2BottomEdge)*by(X1Center, X2BottomEdge)
          * (
               (by(X1Center, X2BottomEdge) * (T_i_j - T_i_jMinus1)/DX2)
             + (  bx(X1Center, X2BottomEdge)/DX1
                * limiter4(T_iPlus1_jMinus1 - T_i_jMinus1, T_iPlus1_j - T_i_j,
                           T_i_jMinus1 - T_iMinus1_jMinus1, T_i_j - T_iMinus1_j
                          )
               )
            )
        )
      );       
      #elif (TIME_STEPPING==SEMI_IMPLICIT)
      divFlux[j][i][0] = 
     -(
        (1./DX1)
      * (   D(X1RightEdge, X2Center)*bx(X1RightEdge, X2Center)
          * (
               (bx(X1RightEdge, X2Center) * (T_iPlus1_j - T_i_j)/DX1)
             + (  by(X1RightEdge, X2Center)/DX2 
                * limiter4( T_old_i_jPlus1 - T_old_i_j, T_old_iPlus1_jPlus1 - T_old_iPlus1_j,
                            T_old_i_j - T_old_i_jMinus1, T_old_iPlus1_j - T_old_iPlus1_jMinus1
                          )
               )
            )
         -  D(X1LeftEdge, X2Center)*bx(X1LeftEdge, X2Center)
          * (
               (bx(X1LeftEdge, X2Center) * (T_i_j - T_iMinus1_j)/DX1)
             + (  by(X1LeftEdge, X2Center)/DX2
                * limiter4(T_old_iMinus1_jPlus1 - T_old_iMinus1_j, T_old_i_jPlus1 - T_old_i_j,
                           T_old_iMinus1_j - T_old_iMinus1_jMinus1, T_old_i_j - T_old_i_jMinus1
                          )
                   
               )   
            )
        )
      + 
        (1./DX2)
      * (   D(X1Center, X2TopEdge)*by(X1Center, X2TopEdge)
          * (
               (by(X1Center, X2TopEdge) * (T_i_jPlus1 - T_i_j)/DX2)
             + (  bx(X1Center, X2TopEdge)/DX1
                * limiter4(T_old_iPlus1_j - T_old_i_j, T_old_iPlus1_jPlus1 - T_old_i_jPlus1,
                           T_old_i_j - T_old_iMinus1_j, T_old_i_jPlus1 - T_old_iMinus1_jPlus1
                          )
               )
            )
        -   D(X1Center, X2BottomEdge)*by(X1Center, X2BottomEdge)
          * (
               (by(X1Center, X2BottomEdge) * (T_i_j - T_i_jMinus1)/DX2)
             + (  bx(X1Center, X2BottomEdge)/DX1
                * limiter4(T_old_iPlus1_jMinus1 - T_old_i_jMinus1, T_old_iPlus1_j - T_old_i_j,
                           T_old_i_jMinus1 - T_old_iMinus1_jMinus1, T_old_i_j - T_old_iMinus1_j
                          )
               )
            )
        )
      );       
      #endif

    }
  }

  DMDAVecRestoreArrayDOF(ts->dmda, primPetscVecOldLocal, &primPetscOld); 
  DMRestoreLocalVector(ts->dmda, &primPetscVecOldLocal);
}
