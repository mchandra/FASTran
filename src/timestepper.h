#ifndef TIMESTEPPER_H_
#define TIMESTEPPER_H_

#include <petsc.h>
#include <petscviewerhdf5.h>
#include <unistd.h>
#include "inputs.h"
#include "limiters.h"

#define BACKWARD_EULER  (0)
#define CRANK_NICHOLSON (1)
#define FULLY_IMPLICIT  (2)
#define SEMI_IMPLICIT   (3)

struct timeStepper
{
  REAL t, dt, tDump;
  int timeStepCounter;
  int dumpCounter;

  SNES snes;
  DM dmda;

  Vec primPetscVec;
  Vec residualPetscVec;
  Vec primPetscVecOld;
  Vec divFluxPetscVecOld;
  Vec divFluxPetscVec;

  int computeOldSourceTermsAndOldDivOfFluxes;

  int X1Start, X1Size;
  int X2Start, X2Size;
  int X3Start, X3Size;
};

/* User functions */

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStep(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1]);

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);
/* Internal functions */

PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residalPetscVec,
                               void *ptr);

void diagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

void computeFlux(struct timeStepper ts[ARRAY_ARGS 1],
                 REAL ***primPetsc, REAL ***divFlux);
#endif
