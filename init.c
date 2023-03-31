/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#define PI (3.141592653589793)

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
    
// ############ Parameters Solar Wind ############
  double R0 = 215.032;  // [solar radii] Solar Wind reference distance in solar radii (1 au)
    
  double V_SW = 4;      // [100 km s-1] Solar Wind velocity at 1 au
  double rho_SW = 5;    // [Mp cm-3] Solar Wind density at 1 au
  //double* T_SW = (double*) malloc(10*sizeof(double));
  double T_SW;
  double T_SW0 = 1.e5;   // [K] Solar Wind temperature at 1 au
  // Constants for the temperature-with-distance model
  double T_alpha = 1.037070592;
  double T_beta = -0.26564;
  double T_gamma = 2.64085;
  double T_delta = -0.65746;
  double B0_SW;
  B0_SW = (2.5)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY); // [code units] Solar Wind magnetic field at 1 solar radius(2.5 G)
    
  double mu;
  mu = MeanMolecularWeight(v);  // [no units] Solar Wind mean molecular mass ~0.6
  double W_SW = 0.018792;       // [rad s-1] Sun's angular velocity (at equator)
    
    
  T_SW = 0.6e5;  //T_SW0 * (T_alpha * pow(x1/R0,T_beta) + T_gamma * pow(x1/R0,T_delta))/2;
    
  /* Temperature in KELVIN */
  v[RHO] = rho_SW * (R0/x1 * R0/x1);
  v[VX1] = V_SW;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  //#if HAVE_ENERGY
  v[PRS] = (v[RHO]*T_SW)/(KELVIN*mu);
  //#endif
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = B0_SW * (1/x1)*(1/x1);
  v[BX2] = 0.0;
  v[BX3] = -B0_SW* (W_SW*1/V_SW) *(1/x1)*sin(x2);
    
  //printf("Double value magnetics 01= %lf %lf %lf \n", v[BX3]);
  //printf("Double value magnetics 03= %lf %lf %lf \n", v[BX3]);
  
  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{/*
    double Bx_ISM, By_ISM, Bz_ISM;
    Bx_ISM = (1.588e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY);
    By_ISM = (-1.112e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY);
    Bz_ISM = (-2.186e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY);
    
    B0[0] = Bx_ISM*sin(x2)*cos(x3) + By_ISM*sin(x2)*sin(x3) + Bz_ISM*cos(x2);
    B0[1] = Bx_ISM*cos(x2)*cos(x3) + By_ISM*cos(x2)*sin(x3) - Bz_ISM*sin(x2);
    B0[2] = -Bx_ISM*sin(x3) + By_ISM*cos(x3);
  */
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 
 
 *********************************************************************** */


{
  int i, j, k;
  double  *x1, *x2, *x3, xc, yc, zc;
    
  // ############ Parameters Solar Wind ############
  double R0 = 215.032;  // [solar radii] Solar Wind reference distance in solar radii (1 au)
        
  double V_SW = 4;      // [100 km s-1] Solar Wind velocity at 1 au
  double rho_SW = 5;    // [Mp cm-3] Solar Wind density at 1 au
  double T_SW;
  double T_SW0 = 1.e5;   // [K] Solar Wind temperature at 1 au
  // Constants for the temperature-with-distance model
  double T_alpha = 1.037070592;
  double T_beta = -0.26564;
  double T_gamma = 2.64085;
  double T_delta = -0.65746;
  double B0_SW;
  B0_SW = (2.5)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY); // [code units] Solar Wind magnetic field at 1 solar radius (2.5 G)

  double mu;
  mu = MeanMolecularWeight(d->Vc);  // [no units] Solar Wind mean molecular mass ~0.6
  double W_SW = 0.018792;           // [rad s-1] Sun's angular velocity (at equator)
  

  // ############ Parameters Interstellar Medium ############
  double V_ISM = 0.26;
  double rho_ISM = 0.1;
  double T_ISM = 8000;

  double Bx_ISM, By_ISM, Bz_ISM;
  //Bx_ISM = (1.5883e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY); // 1.5883 // 0.5294
  //By_ISM = (2.2684e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY); // 2.2684 // 0.7561
  //Bz_ISM = (1.1540e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY); // 1.1540 // 0.3847
    
  By_ISM = (3.e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY);  // Set ISM magnetic field parallel to y-axis
  //Bx_ISM = (3.e-6)/sqrt(4*3.141592653589793*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY);
    
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    // -- check solution inside domain --
    TOT_LOOP(k,j,i){ // Try DOM_LOOP as well
        xc = x1[i]*sin(x2[j])*cos(x3[k]);
        yc = x1[i]*cos(x2[j])*cos(x3[k]);
        zc = -x1[i]*sin(x3[k]);
        
        if (xc<-3000 && x1[i]>104500){ // Change x1[i] depending on resolution
            //if (g_time>=250000){
                d->Vc[RHO][k][j][i] = rho_ISM;
                d->Vc[VX1][k][j][i] = V_ISM * (sin(x2[j]))*cos(x3[k]);
                d->Vc[VX2][k][j][i] = V_ISM * (cos(x2[j]))*cos(x3[k]);
                d->Vc[VX3][k][j][i] = -V_ISM * sin(x3[k]);
                d->Vc[PRS][k][j][i] = ((d->Vc[RHO][k][j][i])*T_ISM)/(mu*KELVIN);
                
                
                d->Vc[BX1][k][j][i] = By_ISM*sin(x2[j])*sin(x3[k]);
                d->Vc[BX2][k][j][i] = By_ISM*cos(x2[j])*sin(x3[k]);
                d->Vc[BX3][k][j][i] = By_ISM*cos(x3[k]);
                // For perpendicular magnetic field
                
            
                /*
                d->Vc[BX1][k][j][i] = Bx_ISM*sin(x2[j])*cos(x3[k]);
                d->Vc[BX2][k][j][i] = Bx_ISM*cos(x2[j])*cos(x3[k]);
                d->Vc[BX3][k][j][i] = -Bx_ISM*sin(x3[k]);
                // For parallel magnetic field
                */
            
                /*
                 d->Vc[BX1][k][j][i] = Bx_ISM*sin(x2[j])*cos(x3[k]) + By_ISM*sin(x2[j])*sin(x3[k]) + Bz_ISM*cos(x2[j]);
                 d->Vc[BX2][k][j][i] = Bx_ISM*cos(x2[j])*cos(x3[k]) + By_ISM*cos(x2[j])*sin(x3[k]) - Bz_ISM*sin(x2[j]);
                 d->Vc[BX3][k][j][i] = -Bx_ISM*sin(x3[k]) + By_ISM*cos(x3[k]);
                 // For accurate magnetic field
                 */
            //}
        }
        /*
        if (x1[i] > 104500 && x2[j] > -0.942 && x2[j] < 2.199 && x3[k] > 0.960 && x3[k] < 4.102){
            d->Vc[BX1][k][j][i] = Bx_ISM*sin(x2[j])*cos(x3[k]) + By_ISM*sin(x2[j])*sin(x3[k]) + Bz_ISM*cos(x2[j]);
            d->Vc[BX2][k][j][i] = Bx_ISM*cos(x2[j])*cos(x3[k]) + By_ISM*cos(x2[j])*sin(x3[k]) - Bz_ISM*sin(x2[j]);
            d->Vc[BX3][k][j][i] = -Bx_ISM*sin(x3[k]) + By_ISM*cos(x3[k]);
        }
        */  // Used to have a separate internal boundary for ISM magnetic field. Overlapping and separate internal boundaries do not seem to work very well so this was not used.
        
    }
  }
  
  /*printf("Double value = %lf \n", mu);*/
  if (side == X1_BEG){
  /* -- select the boundary side -- */
    if (box->vpos == CENTER){
    /* -- select the variable position -- */
       /*printf("Double value pres = %lf \n", d->Vc[PRS][k][j][i]);*/
       BOX_LOOP(box,k,j,i){
      /* -- Loop over boundary zones -- */
               //T_SW = T_SW0 * T_alpha * pow(x1[i],T_beta);
               T_SW = 0.6e5;//T_SW0 * (T_alpha * pow(x1[i]/R0,T_beta) + T_gamma * pow(x1[i]/R0,T_delta))/2;
           
               d->Vc[RHO][k][j][i] = rho_SW * (R0/x1[i] * R0/x1[i]);
               d->Vc[VX1][k][j][i] = V_SW;
               d->Vc[VX2][k][j][i] = 0.0;
               d->Vc[VX3][k][j][i] = 0.0;
               d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*T_SW/(mu*KELVIN);
               d->Vc[BX1][k][j][i] = B0_SW * (1/x1[i])*(1/x1[i]);
               d->Vc[BX2][k][j][i] = 0.0;
               d->Vc[BX3][k][j][i] = -B0_SW * (W_SW * 1/V_SW) *(1/x1[i])*sin(x2[i]);
               
               //printf("Double value pres = %lf %ld %ld %ld\n", d->Vc[PRS][k][j][i], i, j, k);
               //printf("Double value rho = %lf %ld %ld %ld\n", d->Vc[RHO][k][j][i], i, j, k);
               //printf("Double value temp = %lf \n", T);
               //printf("Double value mu = %lf \n", mu);
               //printf("Double value KELVIN = %lf %ld %ld %ld\n", KELVIN);
               //printf("Double value prescal = %lf \n", d->Vc[RHO][k][j][i]*T/(mu*KELVIN));
               //if (d->Vc[PRS][k][j][i]<0){
               //printf("Double value pres negative ryeye = %lf \n", d->Vc[PRS][k][j][i]);
               //printf("Double value magnetics 1= %lf %d %d %d %lf %lf\n", d->Vc[BX3][k][j][i], k, j, i, x1[i], x3[k]);
    
     }
    }
   }
  if (side == X1_END){
  /* -- select the boundary side -- */
    if (box->vpos == CENTER){
    /* -- select the variable position -- */
      BOX_LOOP(box,k,j,i){
        /*
        xc = x1[i]*sin(x2[j])*cos(x3[k]);
        yc = x1[i]*cos(x2[j])*cos(x3[k]);
        zc = -x1[i]*sin(x3[k]);
         */
        /*printf("Double value rh0 = %lf %lf %lf %lf\n", x1[i], xc, x3[k], x2[j]-1.57);*/
          /*
        if (xc<-3000){
            
            d->Vc[RHO][k][j][i] = rho_ISM;
            d->Vc[VX1][k][j][i] = V_ISM * (sin(x2[j]))*cos(x3[k]);
            d->Vc[VX2][k][j][i] = V_ISM * (cos(x2[j]))*cos(x3[k]);
            d->Vc[VX3][k][j][i] = -V_ISM * sin(x3[k]);
            d->Vc[PRS][k][j][i] = ((d->Vc[RHO][k][j][i])*T_ISM)/(mu*KELVIN);
            d->Vc[BX1][k][j][i] = 0; 
            d->Vc[BX2][k][j][i] = 0.0;
            d->Vc[BX3][k][j][i] = 0;
            
            
            
            //printf("Double value rh0 = %lf %lf %lf %lf %lf\n", x1[i], xc, x3[k], x2[j]-1.57, d->Vc[VX1][k][j][i]);
            //printf("Double value pres = %lf \n", d->Vc[PRS][k][j][i]);
            //printf("Double value temp = %lf \n", T);
            //printf("Double value mu = %lf \n", mu);
            //printf("Double value prescal = %lf \n", d->Vc[RHO][k][j][i]*T/(mu*KELVIN));
            //printf("Double value rh0 = %lf %lf %lf \n", d->Vc[VX1][k][j][i], sin(x2[k]), x2[k]);
            //printf("Double value magnetics 1= %lf %d %d %d %lf %lf\n", d->Vc[BX3][k][j][i], k, j, i, x1[i], x3[k]);
         
       }else{
         
         d->Vc[RHO][k][j][i] = 2.e-5;
         d->Vc[VX1][k][j][i] = 0; //sin(x2[j])*cos(x3[k]);
         d->Vc[VX2][k][j][i] = 0; //cos(x2[j])*cos(x3[k]);
         d->Vc[VX3][k][j][i] = 0; //sin(x3[k]);
         d->Vc[PRS][k][j][i] = ((d->Vc[RHO][k][j][i])*T_ISM)/(mu*KELVIN);
         d->Vc[BX1][k][j][i] = 0; 
         d->Vc[BX2][k][j][i] = 0.0;
         d->Vc[BX3][k][j][i] = 0;
         
      }*/
    }
   }
  } //encloawa all
   
  
  else if (box->vpos == X2FACE){ /* -- y staggered field -- */
#ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = g_inputParam[B_PRE];
#endif
  }
  else if (box->vpos == X3FACE){ /* -- z staggered field -- */
  #ifdef STAGGERED_MHD
  BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = g_inputParam[B_PRE];
  #endif
  }
    
 }
 



#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
