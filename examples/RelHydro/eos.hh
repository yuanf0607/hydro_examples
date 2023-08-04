#pragma once

template <typename T = double> class SimpleGammaLaw {
  //The template parameter T represents the data type used for floating-point calculations, such as float or double. 
  //It has a default value of double.

  //The class SimpleGammaLaw contains a static member variable Gamma, representing the adiabatic index or the specific heat ratio (gamma). 
  //It is a static member variable, which means it is shared across all instances of the class. 
  //The specific value of Gamma can be set by the user for their specific simulation.
  public:

  static T Gamma;
  template <typename Td>
  //In summary, T is a class-level type template parameter that determines the overall precision of the SimpleGammaLaw class, 
  //and Td is a local type template parameter used within the member functions to specify the data type of function arguments and local variables. 
  static inline T press_eps__rho_eps_th(T& eps_tot, Td const rho, Td const eps_th) {
    //calculate pressure from the total specific energy and density
          //the eps_tot parameter is passed by reference so that any modifications made inside the press_eps__rho_eps_th 
          //function will be reflected in the original variable used as an argument.
          //rho and eps_th are passed by value, and the function cannot modify them directly but can work with local copies of these values.
    eps_tot = eps_th;
    //eps_th = specific internal energy
    //eps_tot = specific total energy  
    return eps_tot * (Gamma - 1.);
    
  };

  template <typename Td>
  static inline T eps_th__eps_rho(Td const eps_tot, Td const rho) {
    //convert between total specific energy and specific internal energy
    return eps_tot;
  };

  template <typename Td> static inline T cs2__rho_eps_th(Td const rho, Td const eps_th) {
    //calculate the square of the speed of sound
    auto const eps_tot = eps_th;
    auto press = eps_tot * (Gamma - 1.);

    auto rhoh = eps_tot*Gamma + rho;
    return Gamma * press / rhoh;
  };
};
//Each function is templated on the data type T and takes a template parameter Td, representing the data type for the function arguments.

template <typename T> T SimpleGammaLaw<T>::Gamma = 2.;
//This line defines the static member variable Gamma of the SimpleGammaLaw class. 
//The value 2.0 is assigned to Gamma, which serves as the default value for the adiabatic index.

//In summary, the eos.hh file contains the class template SimpleGammaLaw, representing a simple gamma-law equation of state used in simulations to model the behavior of fluids. 
//The class provides functions to calculate pressure, convert between specific internal energy and total specific energy, and calculate the square of the speed of sound. 

//The values of Gamma and the specific calculations for these functions can be customized based on the data type T specified by the user. 

//The class provides a basic implementation for a specific equation of state but can be extended 
  //to support other equations of state by defining additional member functions and providing appropriate calculations.