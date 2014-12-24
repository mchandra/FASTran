#include "timestepper.h"


PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residualPetscVec,
                               void *ptr)
{
  struct timeStepper *ts = (struct timeStepper*)ptr;

  int X1Start, X2Start, X3Start;
  int X1Size, X2Size, X3Size;

  DMDAGetCorners(ts->dmda,
                 &X1Start, &X2Start, &X3Start,
                 &X1Size, &X2Size, &X3Size);

#if (TIME_STEPPING == CRANK_NICHOLSON) 
  if (ts->computeOldSourceTermsAndOldDivOfFluxes)
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

    REAL ***primOldLocal;
    REAL ***divFluxOldGlobal;
    
    DMDAVecGetArrayDOF(ts->dmda, primPetscVecOldLocal, &primOldLocal);
    DMDAVecGetArrayDOF(ts->dmda, ts->divFluxPetscVecOld, &divFluxOldGlobal);

    computeFlux(ts, primOldLocal, divFluxOldGlobal);

    DMDAVecRestoreArrayDOF(ts->dmda, primPetscVecOldLocal, &primOldLocal);
    DMDAVecRestoreArrayDOF(ts->dmda, ts->divFluxPetscVecOld, &divFluxOldGlobal);

    DMRestoreLocalVector(ts->dmda, &primPetscVecOldLocal);

    /* All old sources and divFluxes have now been computed */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  }
#endif

  Vec primPetscVecLocal;
  DMGetLocalVector(ts->dmda, &primPetscVecLocal);

  /* Exchange ghost zone data. */
  DMGlobalToLocalBegin(ts->dmda, 
                       primPetscVec,
                       INSERT_VALUES,
                       primPetscVecLocal);
  DMGlobalToLocalEnd(ts->dmda,
                     primPetscVec,
                     INSERT_VALUES,
                     primPetscVecLocal);

  REAL ***primLocal;
  REAL ***residualGlobal;
  REAL ***divFluxOldGlobal;
  REAL ***divFluxGlobal;
  REAL ***primOldGlobal;

  DMDAVecGetArrayDOF(ts->dmda, primPetscVecLocal, &primLocal);
  DMDAVecGetArrayDOF(ts->dmda, residualPetscVec, &residualGlobal);
  DMDAVecGetArrayDOF(ts->dmda, ts->divFluxPetscVecOld, &divFluxOldGlobal);
  DMDAVecGetArrayDOF(ts->dmda, ts->divFluxPetscVec, &divFluxGlobal);
  DMDAVecGetArrayDOF(ts->dmda, ts->primPetscVecOld, &primOldGlobal); 

  computeFlux(ts, primLocal, divFluxGlobal);

  for (int j=ts->X2Start; j<ts->X2Start+ts->X2Size; j++)
  {
    for (int i=ts->X1Start; i<ts->X1Start+ts->X1Size; i++)
    {
      residualGlobal[j][i][0] = 
           (primLocal[j][i][0] - primOldGlobal[j][i][0])/ts->dt
         + divFluxGlobal[j][i][0];
//#if (TIME_STEPPING == BACKWARD_EULER)
//      residualGlobal[j][i][0] = 
//           (primLocal[j][i][0] - primOldGlobal[j][i][0])/ts->dt
//         + divFluxGlobal[j][i][0];
//#elif (TIME_STEPPING == CRANK_NICHOLSON)
//      residualGlobal[j][i][0] = 
//           (primLocal[j][i][0] - primOldGlobal[j][i][0])/ts->dt
//     + 0.5*(divFluxOldGlobal[j][i][0] + divFluxGlobal[j][i][0];
//#endif
    }
  }

  DMRestoreLocalVector(ts->dmda, &primPetscVecLocal);

  DMDAVecRestoreArrayDOF(ts->dmda, primPetscVecLocal, &primLocal);
  DMDAVecRestoreArrayDOF(ts->dmda, residualPetscVec, &residualGlobal);
  DMDAVecRestoreArrayDOF(ts->dmda, ts->divFluxPetscVecOld, &divFluxOldGlobal);
  DMDAVecRestoreArrayDOF(ts->dmda, ts->divFluxPetscVec, &divFluxGlobal);
  DMDAVecRestoreArrayDOF(ts->dmda, ts->primPetscVecOld, &primOldGlobal); 

  return(0);
}
