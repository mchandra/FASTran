#include "inputs.h"
#include "timestepper.h"

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  REAL ***primOldGlobal;
  DMDAVecGetArrayDOF(ts->dmda, ts->primPetscVecOld, &primOldGlobal);

  for (int j=ts->X2Start; j<ts->X2Start+ts->X2Size; j++)
  {
    for (int i=ts->X1Start; i<ts->X1Start+ts->X1Size; i++)
    {
      REAL X1Center = (X1_A + X1_B)/2.;
      REAL X2Center = (X2_A + X2_B)/2.;

      REAL X1 = X1_ZONE_CENTER(i);
      REAL X2 = X2_ZONE_CENTER(j);

      REAL r = sqrt( pow(X1 - X1Center, 2.) + pow(X2 - X2Center, 2) );

      REAL theta = atan(X2/X1);

//      primOldGlobal[j][i][0] = 1 + 0.2*exp(-r*r/0.005);

//      if (r < INITIAL_RADIUS)
//      {
//        primOldGlobal[j][i][0] = TEMP_INSIDE;
//      } else
//      {
//        primOldGlobal[j][i][0] = TEMP_OUTSIDE;
//      }

      if (r > 0.5 && r < 0.7 && theta > -M_PI/12. && theta < M_PI/12. && X1 < 0.)
      {
        primOldGlobal[j][i][0] = TEMP_INSIDE;
      } else
      {
        primOldGlobal[j][i][0] = TEMP_OUTSIDE;
      }
    }
  }

  DMDAVecRestoreArrayDOF(ts->dmda, ts->primPetscVecOld, &primOldGlobal);

}

REAL bx(REAL X1, REAL X2)
{
  return -X2/sqrt(X1*X1 + X2*X2);
}

REAL by(REAL X1, REAL X2)
{
  return X1/sqrt(X1*X1 + X2*X2);
}

/* Diffusion coefficient */
REAL D(REAL X1, REAL X2)
{
  return D0;
}
