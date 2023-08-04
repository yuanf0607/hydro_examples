/*
 * =====================================================================================
 *
 *       Filename:  initial_data.hh
 *
 *    Description:  Initial data
 *
 *        Version:  1.0
 *        Created:  19.06.2018 22:57:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */


#pragma once

struct Advected_Wave1D{
//sets up the initial condition of a smooth wave function
//It iterates over all degrees of freedom (ijk) and sets the initial state (U[ijk]) to the value of the wave function at each spatial coordinate (x)
  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    //Implement u = exp (-2 cos(2 pi x))
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for simd
    for(int ijk = 0; ijk < U.ndof; ++ijk){
      auto const x = U.grid.get_coords(0,ijk);

      // Implement wave
      // U[ijk] = ...
      U[ijk] = std::exp(-2 * std::cos(2 * M_PI * x));//Yuan
      
      //U[ijk] = std::cos(2 * M_PI * x); //Yuan

    }
  };
};

struct Advected_Step1D{
//sets up the initial condition of a discontinuous step function
//It iterates over all degrees of freedom (ijk) and sets the initial state (U[ijk]) to 1 if the spatial coordinate (x) falls within the specified range, or 0 otherwise
  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    //Implement step function
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for
    for(int ijk = 0; ijk < U.ndof; ++ijk){
      auto const x = U.grid.get_coords(0,ijk);

        // Implement step function
	// U[ijk] = ...
  if(std::abs(x - 0.5) < 0.25) {
          U[ijk] = 1;
      } else {
          U[ijk] = 0;
      }


    }
  };
};

