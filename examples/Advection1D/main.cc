/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  1D Advection test problem
 *
 *        Version:  1.0
 *        Created:  19.06.2018 13:51:38
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */


#include "../../advance.hh"
#include "../../advect.hh"
#include "../../boundaries.hh"
#include "../../output.hh"
#include "../../reconstruct.hh"
#include "../../riemann.hh"
#include "../../storage.hh"
#include "./system.hh"
#include "./initial_data.hh"
#include <vector> 
//#include <Eigen/Dense>

static constexpr int this_ndim = 1; //the number of dimensions
const double global_cfl = 0.2; //the Courant–Friedrichs–Lewy (CFL) condition

const double global_final_time = 10; //the final simulation time 

const int global_ngz = 4; //ghost zones
const int global_ext = 1024; //number of discrete grid points
const double global_dx = 1./global_ext; //grid extent and resolution
const double global_origin_x = 0.; //grid origin
const int output_iteration = 4000000; //output iteration

static constexpr bool HO = true; // High order scheme (4th) //sets the high order scheme flag


//for convenience, provide shorter names for some types used in the code
using this_system_t = SimpleAdvectionSystem<1,double>;

using this_grid_t = SimpleCartesianGrid<this_ndim>;
using this_storage = SimpleStorage<this_grid_t,this_system_t,double>;

using this_reconstruct_t = WenoZ_Reconstruct<false,double>;

using this_riemann_t = HLL_RiemannSolver<this_system_t>;
using this_hrsc_t = McCorquodale_FV<HO,this_system_t, double,this_grid_t,
      this_reconstruct_t, this_riemann_t>;




//Select 1D initial data
using initial_data_t = Advected_Wave1D; //for the wave initial condition  
//using initial_data_t = Advected_Step1D; //for the step initial condition




static this_grid_t global_grid;



int main(){



std::vector<double> errors; // to store the errors
std::vector<double> dx_values; // to store the dx values
  // loop over a range of dx values
  for(int n = 5; n <= 10; n++) { // modify the range as necessary
  int global_ext = pow(2, n);
  double global_dx = 1./global_ext;


//set up the numerical grid and the storage for the system state
  std::array<int,this_ndim> ngz {{global_ngz}};
  std::array<int,this_ndim> extent {{global_ext + 2*global_ngz}};
  std::array<double,this_ndim> dx{{global_dx}};
  std::array<double,this_ndim> origin{{global_origin_x - global_dx*global_ngz}};







  global_grid = this_grid_t(std::move(extent), std::move(ngz), std::move(dx), std::move(origin));


  this_storage U(global_grid);
    this_storage U0(global_grid); // to store the initial state

    // save the initial state
    initial_data_t::initial_data(U);
    U0 = U; // save the initial state in U0


//set up the numerical scheme, which includes setting up the reconstruction, Riemann solver, and a Runge-Kutta 4th order time-stepping scheme
  this_hrsc_t HRSC(U.grid);

  //Select a time stepper:

  auto tstepper = TimeStepperRK4<this_system_t,this_storage>(U.grid);
  tstepper.set_timestep_from_CFL(global_cfl);

//The evolution function is set up with boundary conditions and advection calculations
//It takes the current system state, updates the boundaries, computes the fluxes, and outputs the updated state
  auto evolve_func = [&] (auto & U, auto &Uout) {

    //First boundaries
    
    PeriodicBC<this_ndim>(U);

    //Second fluxes

    HRSC.advect(U, Uout);

  };



//The initial conditions are set and output using the chosen initial condition function
//The state is then switched to volume representation and the boundary conditions are updated
  initial_data_t::initial_data(U);

  //Output initial data
  PeriodicBC<this_ndim>(U);
  output(U,tstepper.get_iteration());

  //Initial data is pointwise, now switch to volumes
  HRSC.switch_to_volume(U);
  PeriodicBC<this_ndim>(U);

//The system is evolved over time with the Runge-Kutta scheme until the final time is reached. 
//The system state is output at each output_iteration.
  while(tstepper.get_current_time() < global_final_time){
    tstepper.advance(U,evolve_func);

    if(tstepper.get_iteration() % output_iteration == 0){
      PeriodicBC<this_ndim>(U);
      output(U,tstepper.get_iteration());
    }
  };
//If the final time was not exactly reached, an extra time step is performed to reach it. 
//Then the final system state is output.
  //If not yet at right time, adjust!
  auto t_remain = global_final_time - tstepper.get_current_time();

  if(t_remain>0){
    tstepper.set_timestep(t_remain);
    tstepper.advance(U,evolve_func);
  };


  //Output final data
  PeriodicBC<this_ndim>(U);
  output(U,tstepper.get_iteration());



// Calculate the error
    double error = 0.0;
    for (int i = 0; i < U.grid.ndof; ++i) {
        double u_nt = U[i];
        double u_0 = U0[i];
        error += std::abs(u_nt - u_0);
    }
    error /= U.grid.ndof; // average error per grid cell
    //std::cout << "Error: " << error << std::endl;

    // Store the error and the dx value
    errors.push_back(error);
    dx_values.push_back(global_dx);
  }


  std::cout << "Errors: ";
for(const auto& error : errors) {
    std::cout << error << " ";
}
std::cout << std::endl;

std::cout << "dx_values: ";
for(const auto& dx : dx_values) {
    std::cout << dx << " ";
}
std::cout << std::endl;
  return 0;
/*
  Eigen::VectorXd x(dx_values.size());
Eigen::VectorXd y(errors.size());

for (size_t i = 0; i < dx_values.size(); ++i) {
    x(i) = log(dx_values[i]);
    y(i) = log(errors[i]);
}

// Perform the linear fit

Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
Eigen::MatrixXd A(x.size(), 2);
A << ones, x;

Eigen::VectorXd result = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);

double A_val = exp(result(0)); // A is the exponent of the first coefficient
double k = result(1); // k is the second coefficient

std::cout << "A = " << A_val << ", k = " << k << std::endl;

*/

};




