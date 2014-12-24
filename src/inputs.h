#ifndef GRIM_INPUT_H_
#define GRIM_INPUT_H_

/* Macros */
#define DX1 ( (X1_B - X1_A)/((REAL)N1) )
#define DX2 ( (X2_B - X2_A)/((REAL)N2) )

#define X1_ZONE_CENTER(i)     (X1_A + (i + 0.5)*DX1)
#define X2_ZONE_CENTER(j)     (X1_A + (j + 0.5)*DX2)

#define X1_ZONE_RIGHT_EDGE(i) (X1_A + (i + 1)*DX1)
#define X1_ZONE_LEFT_EDGE(i)  (X1_A + (i)*DX1)

#define X2_ZONE_BOTTOM_EDGE(j) (X2_A + (j)*DX2)
#define X2_ZONE_TOP_EDGE(j)    (X2_A + (j+1)*DX2)


/* Immutable constants. Should not be changed */
#define M_PI  (3.141592653589793238462643383279)
#define NG    (2)

#define MINMOD    (0)
#define VAN_LEER  (1)
#define MC        (2)

/* End of immutable constants */
  
#define    DUMP_FILE_PREFIX "data"
#define    DT              (0.0001)
#define    DT_DUMP         (0.0001)
#define    START_TIME      (0.)
#define    FINAL_TIME      (4.)
#define    COURANT         (100000.)
#define    TIME_STEPPING   (SEMI_IMPLICIT)
#define    REAL             double
#define    ARRAY_ARGS       const restrict static

/* Domain inputs */
#define    N1               (200)
#define    N2               (200)
#define  USE_OPENMP 0


/* Initial condition parameters */
#define    INITIAL_RADIUS    (0.05)
#define    TEMP_INSIDE       (10)
#define    TEMP_OUTSIDE      (0.01)

/* Diffusion coefficient */
#define    D0    (1.)

/* Domain */
#define    X1_A  (-1.)
#define    X1_B  (1.)
#define    X2_A  (-1.)
#define    X2_B  (1.)

/* Reconstruction */
#define    LIMITER  (MC)

/* Magnetic field unit vectors */
REAL bx(REAL X1, REAL X2);
REAL by(REAL X1, REAL X2);

/* Diffusion coefficient */
REAL D(REAL X1, REAL X2);

#endif /* GRIM_INPUT_H_ */
