#include "timestepper.h"

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1])
{
  SNESCreate(PETSC_COMM_WORLD, &ts->snes);

  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
               DMDA_STENCIL_BOX,
               N1, N2,
               PETSC_DECIDE, PETSC_DECIDE,
               1, NG, PETSC_NULL, PETSC_NULL, &ts->dmda);
  
  SNESSetDM(ts->snes, ts->dmda);

  DMDAGetCorners(ts->dmda,
                 &ts->X1Start, &ts->X2Start, &ts->X3Start,
                 &ts->X1Size, &ts->X2Size, &ts->X3Size);
  
  DMCreateGlobalVector(ts->dmda, &ts->primPetscVecOld);
  DMCreateGlobalVector(ts->dmda, &ts->divFluxPetscVecOld);
  DMCreateGlobalVector(ts->dmda, &ts->divFluxPetscVec);
  DMCreateGlobalVector(ts->dmda, &ts->residualPetscVec);
  DMCreateGlobalVector(ts->dmda, &ts->primPetscVec);

  VecSet(ts->primPetscVecOld, 0.);
  VecSet(ts->divFluxPetscVecOld, 0.);
  VecSet(ts->divFluxPetscVec, 0.);
  VecSet(ts->residualPetscVec, 0.);
  VecSet(ts->primPetscVec, 0.);

  SNESSetFunction(ts->snes, ts->residualPetscVec,
                  computeResidual, ts);

  SNESSetFromOptions(ts->snes);

  ts->dt = COURANT * (DX1*DX1 + DX2*DX2)/(4.*D0);
  ts->t = START_TIME;
  ts->tDump = START_TIME;

  ts->timeStepCounter = 0;
  ts->dumpCounter = 0;

  ts->computeOldSourceTermsAndOldDivOfFluxes = 0;

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "           Memory allocation complete\n\n");

  int numProcs;
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  PetscPrintf(PETSC_COMM_WORLD,
              " Number of MPI procs being used       = %d\n",
              numProcs);
  PetscPrintf(PETSC_COMM_WORLD,
              " Grid size                            = %d x %d\n",
              N1, N2);
  PetscPrintf(PETSC_COMM_WORLD,
              " Grid points in each MPI process      = %d x %d\n",
              ts->X1Size, ts->X2Size);

  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  /* Set the initial conditions */
  initialConditions(ts);

  /* Output the initial conditions */
  VecCopy(ts->primPetscVecOld, ts->primPetscVec);
  diagnostics(ts);

}


/* -----------------------
 * Implicit time stepping:
 * -----------------------
 * Go directly from t=n to t=n+1 using
 * 
 *    (U^(n+1) - U^n)/dt + 0.5*(grad(F)^n + grad(F)^(n+1))
 *                       + 0.5*(sources^n + sources^(n+1)) = 0
 */
void timeStep(struct timeStepper ts[ARRAY_ARGS 1])
{
  PetscErrorCode errCode;


  ts->computeOldSourceTermsAndOldDivOfFluxes = 1;

  VecCopy(ts->primPetscVecOld, ts->primPetscVec);
  errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
  CHKERRQ(errCode);
  VecCopy(ts->primPetscVec, ts->primPetscVecOld);

  ts->t = ts->t + ts->dt;
  ts->timeStepCounter++;

  PetscPrintf(PETSC_COMM_WORLD,
              "\nCompleted step %d, current time = %.5f, dt = %.5f\n\n",
              ts->timeStepCounter, ts->t, ts->dt);

  diagnostics(ts);
}

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecDestroy(&ts->primPetscVecOld);
  VecDestroy(&ts->divFluxPetscVecOld);
  VecDestroy(&ts->divFluxPetscVec);
  VecDestroy(&ts->residualPetscVec);
  VecDestroy(&ts->primPetscVec);

  DMDestroy(&ts->dmda);

  SNESDestroy(&ts->snes);

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "# Memory deallocation complete #\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");
}
