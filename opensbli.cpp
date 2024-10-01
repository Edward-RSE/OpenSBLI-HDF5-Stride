#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#define OPS_1D
#define OPS_API 2
#include "ops_seq.h"
#include "opensbliblock00_kernels.h"
#include "io.h"
int main(int argc, char **argv)
{
// Initializing OPS
ops_init(argc,argv,1);
// Set restart to 1 to restart the simulation from HDF5 file
restart = 0;
// User defined constant values
block0np0 = 200;
Delta0block0 = 1.0/block0np0;
HDF5_timing = 0;
c0 = 0.5;
niter = 1.0/0.001;
double rkold[] = {(1.0/4.0), (3.0/20.0), (3.0/5.0)};
double rknew[] = {(2.0/3.0), (5.0/12.0), (3.0/5.0)};
dt = 0.001;
invDelta0block0 = 1.0/(Delta0block0);
ops_decl_const("Delta0block0" , 1, "double", &Delta0block0);
ops_decl_const("HDF5_timing" , 1, "int", &HDF5_timing);
ops_decl_const("block0np0" , 1, "int", &block0np0);
ops_decl_const("c0" , 1, "double", &c0);
ops_decl_const("dt" , 1, "double", &dt);
ops_decl_const("invDelta0block0" , 1, "double", &invDelta0block0);
ops_decl_const("niter" , 1, "int", &niter);
ops_decl_const("simulation_time" , 1, "double", &simulation_time);
ops_decl_const("start_iter" , 1, "int", &start_iter);
// Define and Declare OPS Block
ops_block opensbliblock00 = ops_decl_block(1, "opensbliblock00");
#include "defdec_data_set.h"
// Define and declare stencils
#include "stencils.h"
#include "bc_exchanges.h"
// Init OPS partition
double partition_start0, elapsed_partition_start0, partition_end0, elapsed_partition_end0;
ops_timers(&partition_start0, &elapsed_partition_start0);
ops_partition("");
ops_timers(&partition_end0, &elapsed_partition_end0);
ops_printf("-----------------------------------------\n MPI partition and reading input file time: %lf\n -----------------------------------------\n", elapsed_partition_end0-elapsed_partition_start0);
// Restart procedure
ops_printf("\033[1;32m");
if (restart == 1){
ops_printf("OpenSBLI is restarting from the input file: restart.h5\n");
}
else {
ops_printf("OpenSBLI is starting from the initial condition.\n");
}
ops_printf("\033[0m");
// Constants from HDF5 restart file
if (restart == 1){
ops_get_const_hdf5("simulation_time", 1, "double", (char*)&simulation_time, "restart.h5");
ops_get_const_hdf5("iter", 1, "int", (char*)&start_iter, "restart.h5");
}
else {
simulation_time = 0.0;
start_iter = 0;
}
tstart = simulation_time;

if (restart == 0){
int iteration_range_0_block0[] = {-5, block0np0 + 5};
ops_par_loop(opensbliblock00Kernel000, "Grid_based_initialisation0", opensbliblock00, 1, iteration_range_0_block0,
ops_arg_dat(phi_B0, 1, stencil_0_00_1, "double", OPS_WRITE),
ops_arg_dat(x0_B0, 1, stencil_0_00_1, "double", OPS_RW),
ops_arg_idx());
}

// Initialize loop timers
double cpu_start0, elapsed_start0, cpu_end0, elapsed_end0;
ops_timers(&cpu_start0, &elapsed_start0);
double inner_start, elapsed_inner_start;
double inner_end, elapsed_inner_end;
ops_timers(&inner_start, &elapsed_inner_start);
for(iter=start_iter; iter<=start_iter+niter - 1; iter++)
{
simulation_time = tstart + dt*((iter - start_iter)+1);
if(fmod(iter+1, 100) == 0){
        ops_timers(&inner_end, &elapsed_inner_end);
        ops_printf("Iteration: %d. Time-step: %.3e. Simulation time: %.5f. Time/iteration: %lf.\n", iter+1, dt, simulation_time, (elapsed_inner_end - elapsed_inner_start)/100);
        ops_NaNcheck(phi_B0);
        ops_timers(&inner_start, &elapsed_inner_start);
}

ops_halo_transfer(periodicBC_direction0_side0_5_block0);
ops_halo_transfer(periodicBC_direction0_side1_6_block0);
int iteration_range_7_block0[] = {0, block0np0};
ops_par_loop(opensbliblock00Kernel007, "Save equations", opensbliblock00, 1, iteration_range_7_block0,
ops_arg_dat(phi_B0, 1, stencil_0_00_1, "double", OPS_READ),
ops_arg_dat(phi_RKold_B0, 1, stencil_0_00_1, "double", OPS_WRITE));

for(stage=0; stage<=2; stage++)
{
int iteration_range_3_block0[] = {0, block0np0};
ops_par_loop(opensbliblock00Kernel003, "Convective CD phi_B0 x0 ", opensbliblock00, 1, iteration_range_3_block0,
ops_arg_dat(phi_B0, 1, stencil_0_22_0, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_1, "double", OPS_WRITE));

int iteration_range_4_block0[] = {0, block0np0};
ops_par_loop(opensbliblock00Kernel004, "Convective residual ", opensbliblock00, 1, iteration_range_4_block0,
ops_arg_dat(wk0_B0, 1, stencil_0_00_1, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00_1, "double", OPS_WRITE));

int iteration_range_9_block0[] = {0, block0np0};
ops_par_loop(opensbliblock00Kernel009, "Sub stage advancement", opensbliblock00, 1, iteration_range_9_block0,
ops_arg_dat(Residual0_B0, 1, stencil_0_00_1, "double", OPS_READ),
ops_arg_dat(phi_RKold_B0, 1, stencil_0_00_1, "double", OPS_READ),
ops_arg_dat(phi_B0, 1, stencil_0_00_1, "double", OPS_WRITE),
ops_arg_gbl(&rknew[stage], 1, "double", OPS_READ));

int iteration_range_8_block0[] = {0, block0np0};
ops_par_loop(opensbliblock00Kernel008, "Temporal solution advancement", opensbliblock00, 1, iteration_range_8_block0,
ops_arg_dat(Residual0_B0, 1, stencil_0_00_1, "double", OPS_READ),
ops_arg_dat(phi_RKold_B0, 1, stencil_0_00_1, "double", OPS_RW),
ops_arg_gbl(&rkold[stage], 1, "double", OPS_READ));

ops_halo_transfer(periodicBC_direction0_side0_5_block0);
ops_halo_transfer(periodicBC_direction0_side1_6_block0);
}
}
ops_timers(&cpu_end0, &elapsed_end0);
ops_printf("\nTimings are:\n");
ops_printf("-----------------------------------------\n");
ops_printf("Total Wall time %lf\n",elapsed_end0-elapsed_start0);

HDF5_IO_Write_0_opensbliblock00(opensbliblock00, phi_B0, x0_B0, HDF5_timing);
hdf5_strided(opensbliblock00, phi_B0, x0_B0, 2, HDF5_timing);

ops_exit();
//Main program end
}