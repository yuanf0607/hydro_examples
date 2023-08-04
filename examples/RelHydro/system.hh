/*
 * =====================================================================================
 *
 *       Filename:  system.hh
 *
 *    Description:  Systems of advection type equations
 *
 *        Version:  1.0
 *        Created:  11.06.2018 17:58:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once
//ensures that the contents of this header file are included only once during compilation 
//to prevent duplicate definitions

#include "./eos.hh"
//contains the equation of state class SimpleGammaLaw

template <int Tndim = 1 , typename T = double, typename EOS_t=SimpleGammaLaw<T>> class SRHD {
//as a template class, with template parameters Tndim, T, and EOS_t. 
//These template parameters are used to customize the class based on 
//the number of dimensions, data type, and the specific equation of state implementation

//Tndim: The template parameter representing the number of dimensions in the simulation. 
//It has a default value of 1.

//T: The template parameter representing the data type used for floating-point calculations, 
//such as float or double. 
//It also has a default value of double.

//EOS_t: The template parameter representing the equation of state class to be used in the simulation. 
//It defaults to SimpleGammaLaw<T>, 
//indicating that the default equation of state is the SimpleGammaLaw class with the specified data type.

public:
  // PRIMITVE VARIABLES ARE:
  // \rho, \varepsilon, u_i
  
  enum { RHOB = 0, EPS, WVX, WVY, WVZ, NUM_VARS };
  //enum (short for enumeration) is a user-defined data type that allows you to define a set of named constants, 
  //typically representing a set of related integral values.
  enum { RHOSTAR = 0, TAUENERGY, STX, STY, STZ };
//two enumerations are defined to provide named constants for accessing variables in the primitive and conservative variable arrays.

//RHOB, EPS, WVX, WVY, and WVZ are named constants representing the indices of the primitive variables. 
//They stand for density, specific internal energy, and the three components of fluid velocity, respectively. 
//NUM_VARS is the total number of primitive variables (5 in this case).

//RHOSTAR, TAUENERGY, STX, STY, and STZ are named constants representing the indices of the conservative variables. 
//They are used to access different components of the conserved variable array.


/******
Conserved Variables:
For relativistic hydrodynamics (SRHD) simulations, the typical conserved variables are as follows:
U[RHOB]: Rest mass density (mass per unit volume in the fluid's rest frame).
U[RHOSTAR]: Relativistic mass density times Lorentz factor (used in relativistic hydrodynamics).
U[EPS]: Specific internal energy (energy per unit mass).
U[STX], U[STY], U[STZ]: Momentum components in the x, y, and z directions, respectively.
U[TAUENERGY]: Total energy density (sum of rest mass energy and internal energy).

These conserved variables obey conservation laws and are particularly useful for numerical simulations 
because they help maintain stability and accuracy during calculations. 
However, they are not directly related to observable quantities or easily interpretable physical properties.

Primitive Variables:
Primitive variables include quantities such as rest mass density (rho), pressure (P), velocity components (v_x, v_y, v_z), and specific internal energy (epsilon). 

The conversion from conserved variables to primitive variables (and vice versa) often involves solving a system of equations 
that incorporates the equation of state (EOS), which relates pressure to other thermodynamic quantities.

******/
  static constexpr bool needs_c2p = true;
  // A boolean value indicating whether the class requires a conversion from primitive to conservative variables. 
  //It is set to true in this case, meaning that the conversion is needed
  static constexpr int ndim = Tndim;
  //A constant integer representing the number of dimensions. It takes the value of the template parameter Tndim.
  static constexpr int num_vars = NUM_VARS;
  //A constant integer representing the total number of primitive variables. It is set to NUM_VARS.
  static constexpr int num_aux = 0;
  //A constant integer representing the number of auxiliary variables. It is set to 0 in this case, indicating that there are no auxiliary variables.
  static T rho_atmo;
  static T eps_atmo;
  //Static member variables of type T, representing the atmospheric density and specific internal energy, respectively. 
  //They are initialized to a small value (1.e-10).
  static constexpr double c2p_tol = 1.e-10;
  //A constant double representing a tolerance value used in the conversion from primitive to conservative variables.



  template <int dir, typename Td>
  //template <int dir, typename Td>: This line defines the template parameters for the function. dir is an integer template parameter, and Td is a type template parameter. 
  //The dir parameter represents the direction in which the flux is computed, and Td is the data type used for the primitive variable array U.


      //int dir: This template parameter represents the direction in which the flux is computed. 
      //It allows the caller of the function to specify the direction as an integer value.

      //typename Td: This template parameter represents the data type of the primitive variable array U. It is used to allow the caller to choose the data type of the array elements. 
      //For example, if Td is double, the std::array<double, num_vars> will be used.

  
  static inline decltype(auto) compute_flux(std::array<Td, num_vars> &U) {
  //static inline decltype(auto): This line indicates that the function is a static member function of the class SRHD, and decltype(auto) is used as the return type. 
  //(std::array<Td, num_vars> &U): This is the function parameter list. 
      //It takes a single parameter U, which is a reference to an std::array of type Td with num_vars elements. 
  //The std::array type represents an array with a fixed size, and Td represents the data type of the elements in the array.

  //decltype(auto): This is the return type of the function. 
  //The use of decltype(auto) allows the return type to be automatically deduced based on the return statement inside the function. 
  //In this case, the return type will be std::array<Td, num_vars>.

  //****Explain the funtion:
    //The compute_flux function is defined here. It is a static member function of the class.
    //It takes two template parameters: dir, which represents the direction in which the flux is computed, and Td, the data type for the primitive variable array U.
    //The function returns a std::array containing the fluxes of the hydrodynamic variables.
    //The purpose of this function is to compute the fluxes for all components of the hydrodynamic variables based on the primitive variable array U. 
    //The fluxes are computed using the SRHD equations and the equation of state (EOS_t). 
    //It uses the SimpleGammaLaw<T>::press_eps__rho_eps_th function from the eos.hh file to obtain the pressure and total specific energy from the primitive variables.

    //The function implementation would be specific to the actual calculations performed to compute the fluxes for the hydrodynamic variables based on the provided primitive variable array U. 
    //The function would utilize the template parameters dir and Td to perform the appropriate calculations and access the elements of the U array with the specified data type.


    // Compute fluxes for all components, here we have only one!

    auto const z2 = U[WVX] * U[WVX] + U[WVY] * U[WVY] + U[WVZ] * U[WVZ];
    //U[WVX] retrieves the x-component of velocity from the U array.
    auto const lorentz = std::sqrt(1. + z2);
    //lorentz is the Lorentz factor, representing the relativistic correction factor due to the fluid's velocity.


    // Need EOS call here!!
    T eps_tot;
    auto press = EOS_t::press_eps__rho_eps_th(eps_tot, U[RHOB], U[EPS]); 
    //The pressure represents the thermodynamic pressure of the fluid at the specific location.
    //the index RHOB would indicate the position of the density value.
    //the index EPS would indicate the position of the specific internal energy value.

    //eps_tot represents the total energy per unit mass of the fluid, including both its internal energy and kinetic energy.

    //In a hydrodynamics simulation, the U array typically stores the primitive variables of the fluid at each grid cell or spatial location. 
    //The primitive variables include quantities like density, specific internal energy, and velocity components.

    auto const rhoh = U[RHOB] + (press + eps_tot);
    //rhoh = the relativistic enthalpy density of the fluid
    //It combines the rest mass energy density (U[RHOB]) with the pressure and total specific energy (eps_tot) to account for the relativistic effects in the fluid.
    auto const stx = rhoh * U[WVX];
    //stx = momentum in the x-direction
    auto const sty = rhoh * U[WVY];
    auto const stz = rhoh * U[WVZ];

//    auto const taud = rhoh*lorentz - U[RHOB];
//    This should also work without double precision for small velocities
    auto const taud = (press + eps_tot) * lorentz + U[RHOB]*z2/(lorentz+1.);
    //taud = conserved energy density = covariant energy momentum tensor diagonal term
    //This quantity is essential for maintaining energy conservation and accurately describing the fluid's behavior in relativistic hydrodynamics simulations.

    auto result = std::array<Td, num_vars>{ //a nested initializer list
        {U[RHOB] * U[WVX + dir], (taud)*U[WVX + dir], stx * U[WVX + dir],
         sty * U[WVX + dir], stz * U[WVX + dir]}};
    //The computed values inside the nested initializer list represent the flux components in the x-direction 
    //(due to the dir parameter of the compute_flux function).
    result[STX + dir] += press;
    //This line updates the component of the flux that corresponds to the x-direction.
    //The variable press contains the pressure of the fluid at the specific grid cell, which is added to the x-component of the flux.
    //This update is necessary to compute the net flow of conserved quantities (mass, momentum, and energy) in the x-direction across a control surface 
    return result;
  };

  template <typename Tstorage>
  static inline void add_source(Tstorage &storage, Tstorage &source) {
  //The add_source function is defined here. It is a static member function of the class.
  //The function takes two template parameters: Tstorage, representing the type of the storage arrays storage and source.
  //This function is used to add source terms to the storage array, but in this case, the function implementation is empty, 
  //indicating that this system does not have any source terms.

    // This system has no source
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

  template <typename Tstorage> static inline void fill_aux(Tstorage &UL) {
  //The function takes a template parameter Tstorage, representing the type of the storage array UL.
  //Similar to the add_source function, this function is also empty, indicating that this system does not use any auxiliary variables.
    return;
  };

  template <typename Td>
  static inline void switch_to_cons_single(std::array<Td, num_vars> &U) {
  //The function takes a template parameter Td, representing the data type for the primitive variable array U.
  //This function is responsible for converting a single set of primitive variables to conservative variables in the array U. 
  //It uses the SRHD equations and the equation of state (EOS_t) to perform the conversion.

    auto const z2 = U[WVX] * U[WVX] + U[WVY] * U[WVY] + U[WVZ] * U[WVZ];
    auto const lorentz2 = 1. + z2;
    auto const lorentz = std::sqrt(lorentz2);

    // Need EOS call here!!
    T eps_tot;
    auto press = EOS_t::press_eps__rho_eps_th(eps_tot, U[RHOB], U[EPS]);

    U[RHOSTAR] *= lorentz;
    auto const rhohW = U[RHOSTAR] + (press + eps_tot) * lorentz;

    U[STX] *= rhohW;
    U[STY] *= rhohW;
    U[STZ] *= rhohW;

//    U[TAUENERGY] = rhohW * lorentz - press -  U[RHOSTAR];
//    This should also work without double precision for small velocities
   U[TAUENERGY]  = press*z2 + eps_tot * lorentz2 + U[RHOSTAR]*z2/(lorentz+1.);

    return;
  };

  template <typename Td>
  static inline void switch_to_prims_single(std::array<Td, num_vars> &U) {
  //This function is responsible for converting a single set of conservative variables to primitive variables in the array U.
  //It uses the SRHD equations and the equation of state (EOS_t) to perform the conversion.

//    for(auto & x : U) std::cout << x << " ;  ";
//    std::cout << std::endl;

    if (U[RHOSTAR] < rho_atmo) { //checking if below threshold value
    //If so, it means that the current cell is in a region of rarefaction waves (a region of lower density and pressure compared to the surroundings)
      U[RHOB] = rho_atmo;
      U[EPS] = eps_atmo;
      U[WVX] = 0;
      U[WVY] = 0;
      U[WVZ] = 0;
      //Setting them to specific values representing a low-density and low-pressure region
      return;
    }

//    T press = 0.;
    T epst; //A variable to store the specific total energy (internal energy + rest mass energy) during the iteration.
    T press = EOS_t::press_eps__rho_eps_th(epst, U[RHOB], U[EPS]); 
        //A variable to store the current pressure value calculated from the equation of state (EOS) based on U[RHOB] and U[EPS].
    T press_prev = 1.e99;
    T press_pp = 1.e109;
    //Variables to store the previous and the previous-to-previous pressure values, respectively. 

    std::array<Td, num_vars> Uout;//An array to store the updated values of primitive variables during the iteration.

    auto const rhostari = 1./U[RHOSTAR];
    //The reciprocal of the conserved variable U[RHOSTAR], calculated once to avoid repeated calculations.

    auto const Snorm = std::sqrt(U[STX] * U[STX] + U[STY] * U[STY]+ U[STZ] * U[STZ])*rhostari;
    //The magnitude of the spatial velocity vector, calculated using the conserved variables U[STX], U[STY], and U[STZ].

    int nn = 0; //A counter variable to keep track of the number of iterations.
    T lorentz = 1.; //variable to store Lorentz factor 
    T hWi =  1.; //variable to store the reciprocal of the enthalpy factor,

    bool print = false;
    while ((std::abs(press - press_prev) > c2p_tol * press) && (nn < 1000)) {
      //This line starts a loop that iteratively calculates and updates the primitive variables (U[RHOB] and U[EPS]) until convergence is achieved. 

//      if(nn>500 && print== false){  print=true; nn=0; press= 0; press_prev = 1.e99;};


      auto const hW = (U[TAUENERGY] + press) * rhostari + 1.; // enthalpy factor
      hWi = 1. / hW; //inverse of enthalpy factor


      auto const v = std::min( std::abs(Snorm*hWi), 0.9999499987499375); //magnitude of the fluid velocity
      auto const v2 = v*v;
      auto const lorentz2i = std::abs(1.-(v2));
      auto const lorentzi = std::sqrt(lorentz2i); //the Lorentz factor
      lorentz = 1./lorentzi;


      Uout[RHOB] = U[RHOSTAR] * lorentzi;
      //updates the value of the rest mass density 

//      auto eps = Uout[RHOB]*(hW * lorentzi - 1.) - press;
      auto eps = U[TAUENERGY] * lorentz2i - (U[RHOSTAR] /(lorentz+1.) + press)*v2; //calculate the specific internal energy
      eps =  std::max(eps, Uout[RHOB]*eps_atmo);

      Uout[EPS] = EOS_t::eps_th__eps_rho(eps, Uout[RHOB]); 

      press_pp = press_prev; 
      press_prev = press;
      //store the previous pressure values before performing the EOS call in the next iteration.


      // Need EOS call here!!
      press = EOS_t::press_eps__rho_eps_th(eps, Uout[RHOB], Uout[EPS]);
      press = std::abs(press); //std::max(-(U[RHOSTAR]+U[TAUENERGY]), press);
      //Calculate the pressure press based on the current guess of primitive variables (Uout[RHOB] and Uout[EPS]) 
      //using the equation of state (EOS).
      
      
      //the goal is to find the correct primitive variables (e.g., pressure) from the conserved variables using an iterative approach
      //using the equation of state (EOS) to calculate pressure based on the guessed primitive variables 
      //and then refining the guess iteratively until the pressure converges to a solution.

      //Aitken accelerator
      //a numerical technique used to accelerate the convergence of the iterative process
      // Instead of using the last two pressure values directly, 
      //the Aitken method uses the last three pressure values to estimate the true pressure value more efficiently.

      auto const R = (press - press_prev)/(press_prev - press_pp);
      auto const Paitken = std::abs(press_prev + (press-press_prev)/(1.-R));

     if(std::abs(R) < 1. && nn > 2){
	      press_pp = press_prev;
	      press_prev = press;
	      press=Paitken;
     };


      nn++;
    };



    assert(nn < 1000); // Abort if C2P fails!
    //check that the loop didn't exceed a maximum number of iterations

     auto const conv =  hWi * rhostari * lorentz;

      U[RHOB] =Uout[RHOB];
      U[EPS] = Uout[EPS] ;
      U[WVX] = U[STX] * conv;
      U[WVY] = U[STY] * conv;
      U[WVZ] = U[STZ] * conv;
    //These updated primitive variables now represent the fluid state at the current cell.

    return; 
    //The iterative process ensures that the EOS and conservation laws are satisfied, 
    //and the final result represents the physical state of the fluid at the given cell.
  };


  template <typename Tstorage> static inline void switch_to_cons(Tstorage &UL) {
    //This function is responsible for converting all primitive variables in the storage array UL to conservative variables. 
    //It does this for all points in 1D, 2D, or 3D depending on the dimensionality specified by Tndim. 
    //It utilizes the switch_to_cons_single function to perform the conversion.

    UL.is_primitive = false;
    //This flag indicates that the UL array now holds conservative variables after the conversion.
    if (ndim == 1) {
      //If ndim is equal to 1 (1D simulation), the function iterates over each grid point in the 1D domain.
      //#pragma omp simd
      for (int i = 0; i < UL.grid.extent[0]; ++i) {

        std::array<typename Tstorage::data_t,
                   num_vars + num_aux>
            Ul;
        //A temporary std::array named Ul is created with size num_vars + num_aux. 
        //This array will hold the primitive variables for the current grid point.
        for (int nv = 0; nv < num_vars; ++nv) {
          Ul[nv] = UL[i + UL.grid.extent[0] * (nv)];
        }

        for (int nv = 0; nv < num_aux; ++nv) {
          Ul[nv + num_vars] = UL.aux[i + UL.grid.extent[0] * (nv)];
        }

        switch_to_cons_single(Ul);
        //The function switch_to_cons_single is called, passing the Ul array as an argument. 
        //This function converts the primitive variables to conservative variables using the hydrodynamics equations and the equation of state.
        for (int nv = 0; nv < num_vars; ++nv) {
          UL[i + UL.grid.extent[0] * (nv)] = Ul[nv];
        }
      } // for i
    }

    if (ndim == 2) {
#pragma omp parallel for
      for (int j = 0; j < UL.grid.extent[1]; ++j)
        #pragma omp simd
        for (int i = 0; i < UL.grid.extent[0]; ++i) {

          std::array<typename Tstorage::data_t,
                     num_vars + num_aux>
              Ul;

          for (int nv = 0; nv < num_vars; ++nv) {
            Ul[nv] = UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          for (int nv = 0; nv < num_aux; ++nv) {
            Ul[num_vars + nv] =
                UL.aux[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          switch_to_cons_single(Ul);

          for (int nv = 0; nv < num_vars; ++nv) {
            UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)] = Ul[nv];
          }

        } // for i
    };

    if (ndim == 3) {
      assert(!"Not implemented yet");
    };
    return; // No C2P needed
  };

  template <typename Tstorage>
  static inline void switch_to_prims(Tstorage &UL) {
  //This function is responsible for converting all conservative variables in the storage array UL to primitive variables. 
  //It does this for all points in 1D, 2D, or 3D depending on the dimensionality specified by Tndim. 
  //It utilizes the switch_to_prims_single function to perform the conversion.  

    UL.is_primitive = true;
    if (ndim == 1) {
      //#pragma omp simd
      for (int i = 0; i < UL.grid.extent[0]; ++i) {

        std::array<typename Tstorage::data_t,
                   num_vars + num_aux>
            Ul;

        for (int nv = 0; nv < num_vars; ++nv) {
          Ul[nv] = UL[i + UL.grid.extent[0] * (nv)];
        }

        for (int nv = 0; nv < num_aux; ++nv) {
          Ul[nv + num_vars] = UL.aux[i + UL.grid.extent[0] * (nv)];
        }

        switch_to_prims_single(Ul);
        //The function switch_to_cons_single is called, passing the Ul array as an argument. 
        //This function converts the primitive variables to conservative variables using the hydrodynamics equations and the equation of state.

        for (int nv = 0; nv < num_vars; ++nv) {
          UL[i + UL.grid.extent[0] * (nv)] = Ul[nv];
        }
      } // for i
    }

    if (ndim == 2) {
      //If ndim is equal to 2 (2D simulation), the function iterates over each grid point in the 2D domain using nested loops (2D grid).
#pragma omp parallel for
      for (int j = 0; j < UL.grid.extent[1]; ++j)
        //#pragma omp simd
        for (int i = 0; i < UL.grid.extent[0]; ++i) {

          std::array<typename Tstorage::data_t,
                     num_vars + num_aux>
              Ul;

          for (int nv = 0; nv < num_vars; ++nv) {
            Ul[nv] = UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          for (int nv = 0; nv < num_aux; ++nv) {
            Ul[num_vars + nv] =
                UL.aux[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          switch_to_prims_single(Ul);

          for (int nv = 0; nv < num_vars; ++nv) {
            UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)] = Ul[nv];
          }

        } // for i
    };

    if (ndim == 3) {
      assert(!"Not implemented yet");
    };
    return; // No C2P needed
  };

  template <int dir, typename Td>
  static inline decltype(auto)
  compute_max_characteristics(std::array<Td, num_vars> &U) {
  //It takes two template parameters: dir, representing the direction in which the characteristics are computed, 
  //and Td, the data type for the primitive variable array U.

  //The function returns a std::array containing the maximum characteristics (speed of sound) for the given primitive variable array U.

  //The purpose of this function is to compute the maximum characteristics of the hydrodynamic variables based on the primitive variable array U. 
  //It uses the SRHD equations and the equation of state (EOS_t) to perform the calculations.

    auto const z2 = U[WVX] * U[WVX] + U[WVY] * U[WVY] + U[WVZ] * U[WVZ]; //the square of the magnitude of the velocity vector (v_x^2 + v_y^2 + v_z^2).
    auto const lorentz2 = 1. + z2; //the square of the Lorentz factor (gamma) 
    auto const lorentz = std::sqrt(lorentz2); // Lorentz factor
    auto const lorentzi = 1./lorentz; // inverse of the Lorentz factor 
    auto const lorentzi2 = 1./lorentz2; //the square of the inverse of the Lorentz factor

    auto const v2 = z2 * lorentzi2;

    // Need EOS call here!!
    auto const cs2 = EOS_t::cs2__rho_eps_th(U[RHOB], U[EPS]); //compute the square of the speed of sound

    auto const z2_par = U[WVX + dir] * U[WVX + dir]; //the square of the component of velocity in the direction specified by the template parameter dir.

    auto const tmp = std::sqrt( 
	cs2 * lorentzi2 * (1. - (z2_par + (z2 - z2_par) * cs2) *lorentzi2)
	);

//    auto const tmp = lorentzi2*std::sqrt( 
//	cs2 * (1. +(z2 - z2_par) * (1.-cs2)));


    auto const p1 = U[WVX + dir] * lorentzi * (1. - cs2);

    auto const invden =1./ (1. - v2 * cs2);

    return std::array<T, 2>{{(p1 + tmp) * invden, ((p1 - tmp) * invden)}};
  };
};

template <int Tndim, typename T, typename EOS_t> T SRHD<Tndim,T,EOS_t>::rho_atmo = 1.e-10;

template <int Tndim, typename T,typename EOS_t> T SRHD<Tndim,T,EOS_t>::eps_atmo = 1.e-10;


//The class SRHD defined in this file provides functionalities for a relativistic hydrodynamics simulation. 
//It includes functions for computing fluxes, converting between primitive and conservative variables, and calculating the maximum characteristics of the system. 
//The simulation behavior can be customized using the template parameters Tndim, T, and EOS_t, 
//and it relies on the equation of state defined in the eos.hh file (SimpleGammaLaw in this case).