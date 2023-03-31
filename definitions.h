#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   EIGHT_WAVES
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  P_OUT                          0
#define  P_IN                           1
#define  GAMMA                          2

/* [Beg] user-defined constants (do not change this line) */

#define  RING_AVERAGE                   NO
#define  INTERNAL_BOUNDARY              YES
#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  SHOCK_FLATTENING               NO

#define  UNIT_VELOCITY                  1.e7
#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    CONST_Rsun

/* [End] user-defined constants (do not change this line) */
