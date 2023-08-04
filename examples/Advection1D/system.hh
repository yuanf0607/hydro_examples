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


//defines the SimpleAdvectionSystem class, which encompasses the details of the advection system to be solved
template <int Tndim = 1, typename T = double> class SimpleAdvectionSystem {

public:
// define the static attributes of the system
  static constexpr bool needs_c2p = false;
  static constexpr int ndim = Tndim;
  static constexpr int num_vars = 1;
  static constexpr int num_aux = 0;

  //static constexpr T uL = 1.; // CONSTANT ADVECTION SPEED (original line)
  static constexpr T advection_velocity = 1.; //Yuan


//computing the fluxes based on the system's current state U.
  template <int dir,typename Td>
  static inline decltype(auto) compute_flux(std::array<Td, num_vars> &U) {

    // Compute fluxes for all components, here we have only one!


    Td u = U[0]; //Yuan
    // Compute ∂tu + u∂xu
    Td term = u * u; // Non-linear term //Yuan

    return std::array<Td, num_vars>{{
    	 // IMPLEMENT FLUX FOR ALL COMPONENTS (just one)
    //U[0]*advection_velocity //Yuan

    term // Flux for the non-linear advection equation //Yuan
    }};
  };

//add source terms to the system, but in this case, the system has no source terms so it's empty
  template <typename Tstorage>
  static inline void add_source(Tstorage &storage, Tstorage &source) {

    // This system has no source
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

//fill auxiliary data fields, but there's none in this simple advection system, so it's empty
  template <typename Tstorage>
  static inline void fill_aux(Tstorage &storage) {

    // This system has auxiliary data
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

//switch between conserved and primitive variables. In this simple system, no conversion is needed
  template <typename Tstorage>
  static inline void switch_to_cons(Tstorage &storage) {
    storage.is_primitive = false;
    return; // No C2P needed
  };

  template <typename Td>
  static inline void switch_to_cons_single(std::array<Td, num_vars> &U) {
    return; // No C2P needed
  };

  template <typename Td>
  static inline void switch_to_prims_single(std::array<Td, num_vars> &U) {
    return; // No C2P needed
  };

  template <typename Tstorage>
  static inline void switch_to_prims(Tstorage &storage) {
    storage.is_primitive = false;
    return; // No C2P needed
  };


//computes the maximum speed of the system's characteristics, which is necessary for the CFL condition in time-stepping
  template <int dir,typename Td>
  static inline decltype(auto)
  compute_max_characteristics(std::array<Td, num_vars> &U) {
    // IMPLEMENT CHARACTERISTICS:
    // +-c 
    return std::array<T, 2>{{
       // +-c ..
       advection_velocity, - advection_velocity //Yuan
    }};
  };
};
