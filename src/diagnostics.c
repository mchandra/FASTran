#include "timestepper.h"

void diagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  if (ts->t > ts->tDump || fabs(ts->t - ts->tDump) < 1e-10
                        || fabs(ts->t - FINAL_TIME) < 1e-10 
     )
  {
    char primVarsFileName[50];
    sprintf(primVarsFileName, "%s%06d.h5", DUMP_FILE_PREFIX, ts->dumpCounter);

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, primVarsFileName,
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) ts->primPetscVec, "primVars");
    VecView(ts->primPetscVec, viewer);
    PetscViewerDestroy(&viewer);

    PetscPrintf(PETSC_COMM_WORLD,
                "\nDumped primitive variables at t = %f in %s\n\n",
                ts->t, primVarsFileName);
    
    ts->tDump += DT_DUMP;
    ts->dumpCounter++;
  }
  
}
