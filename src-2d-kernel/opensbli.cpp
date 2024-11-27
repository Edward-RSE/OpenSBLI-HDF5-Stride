#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#define OPS_2D
#define OPS_API 2
#include "ops_seq.h"
#include "opensbliblock00_kernels.h"
#include "io.h"
#include "io_stride.h"
int main(int argc, char **argv) {
  // Initializing OPS
  ops_init(argc, argv, 1);
  // Set restart to 1 to restart the simulation from HDF5 file
  restart = 0;
  // User defined constant values
  niter = 100;
  double rkold[] = {(1.0 / 4.0), (3.0 / 20.0), (3.0 / 5.0)};
  double rknew[] = {(2.0 / 3.0), (5.0 / 12.0), (3.0 / 5.0)};
  dt = 0.001;
  block0np0 = 100;
  block0np1 = 100;
  Delta0block0 = 2.0 / (block0np0);
  Delta1block0 = 2.0 / (block0np1);
  write_output_file = 500;
  HDF5_timing = 1;
  gama = 1.4;
  shock_filter_control = 1.00000000000000;
  gamma_m1 = -1 + gama;
  invDelta0block0 = 1.0 / (Delta0block0);
  invDelta1block0 = 1.0 / (Delta1block0);
  inv_gamma_m1 = 1.0 / ((-1 + gama));
  invgamma_m1 = 1.0 / (gamma_m1);
  ops_decl_const("Delta0block0", 1, "double", &Delta0block0);
  ops_decl_const("Delta1block0", 1, "double", &Delta1block0);
  ops_decl_const("HDF5_timing", 1, "int", &HDF5_timing);
  ops_decl_const("block0np0", 1, "int", &block0np0);
  ops_decl_const("block0np1", 1, "int", &block0np1);
  ops_decl_const("dt", 1, "double", &dt);
  ops_decl_const("gama", 1, "double", &gama);
  ops_decl_const("gamma_m1", 1, "double", &gamma_m1);
  ops_decl_const("invDelta0block0", 1, "double", &invDelta0block0);
  ops_decl_const("invDelta1block0", 1, "double", &invDelta1block0);
  ops_decl_const("inv_gamma_m1", 1, "double", &inv_gamma_m1);
  ops_decl_const("invgamma_m1", 1, "double", &invgamma_m1);
  ops_decl_const("niter", 1, "int", &niter);
  ops_decl_const("shock_filter_control", 1, "double", &shock_filter_control);
  ops_decl_const("simulation_time", 1, "double", &simulation_time);
  ops_decl_const("start_iter", 1, "int", &start_iter);
  ops_decl_const("write_output_file", 1, "int", &write_output_file);
  // Define and Declare OPS Block
  ops_block opensbliblock00 = ops_decl_block(2, "opensbliblock00");
#include "defdec_data_set.h"

  int dat_stride[] = {2, 2};
  HDF5_IO_Init_0_opensbliblock00_strided(opensbliblock00, dat_stride);

// Define and declare stencils
#include "stencils.h"
#include "bc_exchanges.h"
  // Init OPS partition
  double partition_start0, elapsed_partition_start0, partition_end0, elapsed_partition_end0;
  ops_timers(&partition_start0, &elapsed_partition_start0);
  ops_partition("");
  ops_timers(&partition_end0, &elapsed_partition_end0);
  ops_printf("-----------------------------------------\n MPI partition and reading input file time: %lf\n "
             "-----------------------------------------\n",
             elapsed_partition_end0 - elapsed_partition_start0);
  // Restart procedure
  ops_printf("\033[1;32m");
  if (restart == 1) {
    ops_printf("OpenSBLI is restarting from the input file: restart.h5\n");
  } else {
    ops_printf("OpenSBLI is starting from the initial condition.\n");
  }
  ops_printf("\033[0m");
  // Constants from HDF5 restart file
  if (restart == 1) {
    ops_get_const_hdf5("simulation_time", 1, "double", (char *)&simulation_time, "restart.h5");
    ops_get_const_hdf5("iter", 1, "int", (char *)&start_iter, "restart.h5");
  } else {
    simulation_time = 0.0;
    start_iter = 0;
  }
  tstart = simulation_time;

  if (restart == 0) {
    int iteration_range_15_block0[] = {-5, block0np0 + 5, -5, block0np1 + 5};
    ops_par_loop(opensbliblock00Kernel015, "Grid_based_initialisation0", opensbliblock00, 2, iteration_range_15_block0,
                 ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(x0_B0, 1, stencil_0_00_00_2, "double", OPS_RW),
                 ops_arg_dat(x1_B0, 1, stencil_0_00_00_2, "double", OPS_RW), ops_arg_idx());
  }

  // Initialize loop timers
  double cpu_start0, elapsed_start0, cpu_end0, elapsed_end0;
  ops_timers(&cpu_start0, &elapsed_start0);
  double inner_start, elapsed_inner_start;
  double inner_end, elapsed_inner_end;
  ops_timers(&inner_start, &elapsed_inner_start);
  for (iter = start_iter; iter <= start_iter + niter - 1; iter++) {
    simulation_time = tstart + dt * ((iter - start_iter) + 1);
    if (fmod(iter + 1, 100) == 0) {
      ops_timers(&inner_end, &elapsed_inner_end);
      ops_printf("Iteration: %d. Time-step: %.3e. Simulation time: %.5f. Time/iteration: %lf.\n", iter + 1, dt,
                 simulation_time, (elapsed_inner_end - elapsed_inner_start) / 100);
      ops_NaNcheck(rho_B0);
      ops_timers(&inner_start, &elapsed_inner_start);
    }

    ops_halo_transfer(periodicBC_direction0_side0_11_block0);
    ops_halo_transfer(periodicBC_direction0_side1_12_block0);
    ops_halo_transfer(periodicBC_direction1_side0_13_block0);
    ops_halo_transfer(periodicBC_direction1_side1_14_block0);
    int iteration_range_16_block0[] = {0, block0np0, 0, block0np1};
    ops_par_loop(opensbliblock00Kernel016, "Save equations", opensbliblock00, 2, iteration_range_16_block0,
                 ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                 ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                 ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                 ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                 ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                 ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

    for (stage = 0; stage <= 2; stage++) {
      int iteration_range_4_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4};
      ops_par_loop(opensbliblock00Kernel004, "CRu1", opensbliblock00, 2, iteration_range_4_block0,
                   ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(u1_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_5_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4};
      ops_par_loop(opensbliblock00Kernel005, "CRu0", opensbliblock00, 2, iteration_range_5_block0,
                   ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(u0_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_7_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4};
      ops_par_loop(opensbliblock00Kernel007, "CRp", opensbliblock00, 2, iteration_range_7_block0,
                   ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(u0_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(u1_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(p_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_3_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4};
      ops_par_loop(opensbliblock00Kernel003, "CRa", opensbliblock00, 2, iteration_range_3_block0,
                   ops_arg_dat(p_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(a_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_0_block0[] = {-1, block0np0 + 1, 0, block0np1};
      ops_par_loop(opensbliblock00Kernel000, "LFWeno_reconstruction_0_direction", opensbliblock00, 2,
                   iteration_range_0_block0, ops_arg_dat(a_B0, 1, stencil_0_01_00_3, "double", OPS_READ),
                   ops_arg_dat(p_B0, 1, stencil_0_23_00_7, "double", OPS_READ),
                   ops_arg_dat(rhoE_B0, 1, stencil_0_23_00_7, "double", OPS_READ),
                   ops_arg_dat(rho_B0, 1, stencil_0_23_00_7, "double", OPS_READ),
                   ops_arg_dat(rhou0_B0, 1, stencil_0_23_00_7, "double", OPS_READ),
                   ops_arg_dat(rhou1_B0, 1, stencil_0_23_00_7, "double", OPS_READ),
                   ops_arg_dat(u0_B0, 1, stencil_0_23_00_7, "double", OPS_READ),
                   ops_arg_dat(u1_B0, 1, stencil_0_01_00_3, "double", OPS_READ),
                   ops_arg_dat(wk0_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(wk1_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(wk2_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(wk3_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_1_block0[] = {0, block0np0, -1, block0np1 + 1};
      ops_par_loop(opensbliblock00Kernel001, "LFWeno_reconstruction_1_direction", opensbliblock00, 2,
                   iteration_range_1_block0, ops_arg_dat(a_B0, 1, stencil_0_00_01_3, "double", OPS_READ),
                   ops_arg_dat(p_B0, 1, stencil_0_00_23_7, "double", OPS_READ),
                   ops_arg_dat(rhoE_B0, 1, stencil_0_00_23_7, "double", OPS_READ),
                   ops_arg_dat(rho_B0, 1, stencil_0_00_23_7, "double", OPS_READ),
                   ops_arg_dat(rhou0_B0, 1, stencil_0_00_23_7, "double", OPS_READ),
                   ops_arg_dat(rhou1_B0, 1, stencil_0_00_23_7, "double", OPS_READ),
                   ops_arg_dat(u0_B0, 1, stencil_0_00_01_3, "double", OPS_READ),
                   ops_arg_dat(u1_B0, 1, stencil_0_00_23_7, "double", OPS_READ),
                   ops_arg_dat(wk4_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(wk5_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(wk6_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(wk7_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_10_block0[] = {0, block0np0, 0, block0np1};
      ops_par_loop(opensbliblock00Kernel010, "LFWeno Residual", opensbliblock00, 2, iteration_range_10_block0,
                   ops_arg_dat(wk0_B0, 1, stencil_0_10_00_3, "double", OPS_READ),
                   ops_arg_dat(wk1_B0, 1, stencil_0_10_00_3, "double", OPS_READ),
                   ops_arg_dat(wk2_B0, 1, stencil_0_10_00_3, "double", OPS_READ),
                   ops_arg_dat(wk3_B0, 1, stencil_0_10_00_3, "double", OPS_READ),
                   ops_arg_dat(wk4_B0, 1, stencil_0_00_10_3, "double", OPS_READ),
                   ops_arg_dat(wk5_B0, 1, stencil_0_00_10_3, "double", OPS_READ),
                   ops_arg_dat(wk6_B0, 1, stencil_0_00_10_3, "double", OPS_READ),
                   ops_arg_dat(wk7_B0, 1, stencil_0_00_10_3, "double", OPS_READ),
                   ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE));

      int iteration_range_18_block0[] = {0, block0np0, 0, block0np1};
      ops_par_loop(opensbliblock00Kernel018, "Sub stage advancement", opensbliblock00, 2, iteration_range_18_block0,
                   ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(rho_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_2, "double", OPS_WRITE),
                   ops_arg_gbl(&rknew[stage], 1, "double", OPS_READ));

      int iteration_range_17_block0[] = {0, block0np0, 0, block0np1};
      ops_par_loop(opensbliblock00Kernel017, "Temporal solution advancement", opensbliblock00, 2,
                   iteration_range_17_block0, ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_2, "double", OPS_READ),
                   ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_RW),
                   ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_RW),
                   ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_RW),
                   ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_2, "double", OPS_RW),
                   ops_arg_gbl(&rkold[stage], 1, "double", OPS_READ));

      ops_halo_transfer(periodicBC_direction0_side0_11_block0);
      ops_halo_transfer(periodicBC_direction0_side1_12_block0);
      ops_halo_transfer(periodicBC_direction1_side0_13_block0);
      ops_halo_transfer(periodicBC_direction1_side1_14_block0);
    }
    // if (fmod(1 + iter, write_output_file) == 0 || iter == 0) {
    //   HDF5_IO_Write_0_opensbliblock00_dynamic(opensbliblock00, iter, rho_B0, rhou0_B0, rhou1_B0, rhoE_B0, x0_B0,
    //   x1_B0,
    //                                           HDF5_timing);
    // }
  }
  ops_timers(&cpu_end0, &elapsed_end0);
  ops_printf("\nTimings are:\n");
  ops_printf("-----------------------------------------\n");
  ops_printf("Total Wall time %lf\n", elapsed_end0 - elapsed_start0);

  HDF5_IO_Write_0_opensbliblock00(opensbliblock00, rho_B0, rhou0_B0, rhou1_B0, rhoE_B0, x0_B0, x1_B0, HDF5_timing);
  HDF5_IO_Write_0_opensbliblock00_strided(opensbliblock00, dat_stride, rho_B0);
  ops_exit();
  // Main program end
}