#include <stdlib.h> 
#include <string.h> 
#include <math.h> 
#include "constants.h"
#define OPS_3D
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
block0np0 = 64;
block0np1 = 64;
block0np2 = 64;
Delta0block0 = 6.0/block0np0;
Delta1block0 = 2.0/(block0np1-1);
Delta2block0 = 3.0/block0np2;
niter = 10;
double rkB[] = {0.149659021999300, 0.379210312999900, 0.822955029386900, 0.699450455948800, 0.153057247968100};
double rkA[] = {0, -0.417890474500000, -1.19215169464300, -1.69778469247100, -1.51418344425700};
dt = 0.0001;
Twall = 1.0;
Minf = 0.4;
gama = 1.4;
aCF = 100.0;
c1 = 0;
c2 = 0;
Re = 400.0;
Pr = 0.7;
c0 = -1;
stretch = 1.7;
lx0 = 6.0;
lx2 = 3.0;
shock_filter_control = 1.00000000000000;
gamma_m1 = -1 + gama;
Ducros_check = 0.0500000000000000;
Ducros_select = 0.0500000000000000;
inv_rfact0_block0 = 1.0/Delta0block0;
inv_rfact2_block0 = 1.0/Delta2block0;
inv_rfact1_block0 = 1.0/Delta1block0;
filter_frequency = 25;
DRP_filt = 0.100000000000000;
HDF5_timing = 0;
write_output_file = 5;
write_slices = 5;
inv2Delta0block0 = 1.0/(Delta0block0*Delta0block0);
inv2Delta1block0 = 1.0/(Delta1block0*Delta1block0);
inv2Delta2block0 = 1.0/(Delta2block0*Delta2block0);
inv2Minf = 1.0/(Minf*Minf);
invDelta0block0 = 1.0/(Delta0block0);
invDelta1block0 = 1.0/(Delta1block0);
invDelta2block0 = 1.0/(Delta2block0);
invPr = 1.0/(Pr);
invRe = 1.0/(Re);
inv_gamma_m1 = 1.0/((-1 + gama));
invgama = 1.0/(gama);
invgamma_m1 = 1.0/(gamma_m1);
invlx0 = 1.0/(lx0);
invlx2 = 1.0/(lx2);
invniter = 1.0/(niter);
ops_decl_const("DRP_filt" , 1, "double", &DRP_filt);
ops_decl_const("Delta0block0" , 1, "double", &Delta0block0);
ops_decl_const("Delta1block0" , 1, "double", &Delta1block0);
ops_decl_const("Delta2block0" , 1, "double", &Delta2block0);
ops_decl_const("Ducros_check" , 1, "double", &Ducros_check);
ops_decl_const("Ducros_select" , 1, "double", &Ducros_select);
ops_decl_const("HDF5_timing" , 1, "int", &HDF5_timing);
ops_decl_const("Minf" , 1, "double", &Minf);
ops_decl_const("Pr" , 1, "double", &Pr);
ops_decl_const("Re" , 1, "double", &Re);
ops_decl_const("Twall" , 1, "double", &Twall);
ops_decl_const("aCF" , 1, "double", &aCF);
ops_decl_const("block0np0" , 1, "int", &block0np0);
ops_decl_const("block0np1" , 1, "int", &block0np1);
ops_decl_const("block0np2" , 1, "int", &block0np2);
ops_decl_const("c0" , 1, "double", &c0);
ops_decl_const("c1" , 1, "double", &c1);
ops_decl_const("c2" , 1, "double", &c2);
ops_decl_const("dt" , 1, "double", &dt);
ops_decl_const("filter_frequency" , 1, "int", &filter_frequency);
ops_decl_const("gama" , 1, "double", &gama);
ops_decl_const("gamma_m1" , 1, "double", &gamma_m1);
ops_decl_const("inv2Delta0block0" , 1, "double", &inv2Delta0block0);
ops_decl_const("inv2Delta1block0" , 1, "double", &inv2Delta1block0);
ops_decl_const("inv2Delta2block0" , 1, "double", &inv2Delta2block0);
ops_decl_const("inv2Minf" , 1, "double", &inv2Minf);
ops_decl_const("invDelta0block0" , 1, "double", &invDelta0block0);
ops_decl_const("invDelta1block0" , 1, "double", &invDelta1block0);
ops_decl_const("invDelta2block0" , 1, "double", &invDelta2block0);
ops_decl_const("invPr" , 1, "double", &invPr);
ops_decl_const("invRe" , 1, "double", &invRe);
ops_decl_const("inv_gamma_m1" , 1, "double", &inv_gamma_m1);
ops_decl_const("inv_rfact0_block0" , 1, "double", &inv_rfact0_block0);
ops_decl_const("inv_rfact1_block0" , 1, "double", &inv_rfact1_block0);
ops_decl_const("inv_rfact2_block0" , 1, "double", &inv_rfact2_block0);
ops_decl_const("invgama" , 1, "double", &invgama);
ops_decl_const("invgamma_m1" , 1, "double", &invgamma_m1);
ops_decl_const("invlx0" , 1, "double", &invlx0);
ops_decl_const("invlx2" , 1, "double", &invlx2);
ops_decl_const("invniter" , 1, "double", &invniter);
ops_decl_const("lx0" , 1, "double", &lx0);
ops_decl_const("lx2" , 1, "double", &lx2);
ops_decl_const("niter" , 1, "int", &niter);
ops_decl_const("shock_filter_control" , 1, "double", &shock_filter_control);
ops_decl_const("simulation_time" , 1, "double", &simulation_time);
ops_decl_const("start_iter" , 1, "int", &start_iter);
ops_decl_const("stretch" , 1, "double", &stretch);
ops_decl_const("write_output_file" , 1, "int", &write_output_file);
ops_decl_const("write_slices" , 1, "int", &write_slices);
// Define and Declare OPS Block
ops_block opensbliblock00 = ops_decl_block(3, "opensbliblock00");
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
int iteration_range_41_block0[] = {-5, block0np0 + 5, -5, block0np1 + 5, -5, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel041, "Grid_based_initialisation0", opensbliblock00, 3, iteration_range_41_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(x0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(x1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(x2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_idx());
}

int iteration_range_43_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel043, "MetricsEquation evaluation", opensbliblock00, 3, iteration_range_43_block0,
ops_arg_dat(x1_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(D11_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(detJ_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

ops_halo_transfer(periodicBC_direction0_side0_44_block0);
ops_halo_transfer(periodicBC_direction0_side1_45_block0);
int iteration_range_46_block0[] = {-4, block0np0 + 5, 0, 1, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel046, "Metric_copy_block0 boundary dir1 side0", opensbliblock00, 3, iteration_range_46_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_55_00_20, "double", OPS_RW),
ops_arg_dat(detJ_B0, 1, stencil_0_00_55_00_20, "double", OPS_RW));

int iteration_range_47_block0[] = {-4, block0np0 + 5, block0np1 - 1, block0np1, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel047, "Metric_copy_block0 boundary dir1 side1", opensbliblock00, 3, iteration_range_47_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_55_00_20, "double", OPS_RW),
ops_arg_dat(detJ_B0, 1, stencil_0_00_55_00_20, "double", OPS_RW));

ops_halo_transfer(periodicBC_direction2_side0_48_block0);
ops_halo_transfer(periodicBC_direction2_side1_49_block0);
int iteration_range_51_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel051, "MetricsEquation evaluation", opensbliblock00, 3, iteration_range_51_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(SD111_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

char slice_name0[80];
sprintf(slice_name0, "0");
ops_write_plane_group_hdf5({{2, block0np2/2}}, slice_name0, {{x0_B0, x1_B0}});
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
        ops_NaNcheck(rho_B0);
        ops_timers(&inner_start, &elapsed_inner_start);
}

ops_halo_transfer(periodicBC_direction0_side0_35_block0);
ops_halo_transfer(periodicBC_direction0_side1_36_block0);
int iteration_range_37_block0[] = {-4, block0np0 + 5, 0, 1, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel037, "IsothermalWall boundary dir1 side0", opensbliblock00, 3, iteration_range_37_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_51_00_15, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW));

int iteration_range_38_block0[] = {-4, block0np0 + 5, block0np1 - 1, block0np1, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel038, "IsothermalWall boundary dir1 side1", opensbliblock00, 3, iteration_range_38_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_15_00_15, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW));

ops_halo_transfer(periodicBC_direction2_side0_39_block0);
ops_halo_transfer(periodicBC_direction2_side1_40_block0);
for(stage=0; stage<=4; stage++)
{
int iteration_range_5_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel005, "CRu0_B0", opensbliblock00, 3, iteration_range_5_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_7_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel007, "CRu1_B0", opensbliblock00, 3, iteration_range_7_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_9_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel009, "CRu2_B0", opensbliblock00, 3, iteration_range_9_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_34_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel034, "CRphi_B0", opensbliblock00, 3, iteration_range_34_block0,
ops_arg_dat(x1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(phi_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_23_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel023, "CRp_B0", opensbliblock00, 3, iteration_range_23_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_11_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel011, "CRT_B0", opensbliblock00, 3, iteration_range_11_block0,
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(T_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_20_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel020, "CRH_B0", opensbliblock00, 3, iteration_range_20_block0,
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(H_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_26_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel026, "CRmu_B0", opensbliblock00, 3, iteration_range_26_block0,
ops_arg_dat(T_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(mu_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_4_block0[] = {0, block0np0, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel004, "Derivative evaluation CD u0_B0 xi0 ", opensbliblock00, 3, iteration_range_4_block0,
ops_arg_dat(u0_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_6_block0[] = {0, block0np0, -2, block0np1 + 2, 0, block0np2};
ops_par_loop(opensbliblock00Kernel006, "Derivative evaluation CD u1_B0 xi0 ", opensbliblock00, 3, iteration_range_6_block0,
ops_arg_dat(u1_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_8_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel008, "Derivative evaluation CD u2_B0 xi0 ", opensbliblock00, 3, iteration_range_8_block0,
ops_arg_dat(u2_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_10_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel010, "Derivative evaluation CD T_B0 xi0 ", opensbliblock00, 3, iteration_range_10_block0,
ops_arg_dat(T_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_12_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel012, "Derivative evaluation CD u0_B0 xi1 ", opensbliblock00, 3, iteration_range_12_block0,
ops_arg_dat(u0_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

int iteration_range_13_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel013, "Derivative evaluation CD u1_B0 xi1 ", opensbliblock00, 3, iteration_range_13_block0,
ops_arg_dat(u1_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

int iteration_range_14_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel014, "Derivative evaluation CD u2_B0 xi1 ", opensbliblock00, 3, iteration_range_14_block0,
ops_arg_dat(u2_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

int iteration_range_15_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel015, "Derivative evaluation CD T_B0 xi1 ", opensbliblock00, 3, iteration_range_15_block0,
ops_arg_dat(T_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

int iteration_range_16_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel016, "Derivative evaluation CD u0_B0 xi2 ", opensbliblock00, 3, iteration_range_16_block0,
ops_arg_dat(u0_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_17_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel017, "Derivative evaluation CD u1_B0 xi2 ", opensbliblock00, 3, iteration_range_17_block0,
ops_arg_dat(u1_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_18_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel018, "Derivative evaluation CD u2_B0 xi2 ", opensbliblock00, 3, iteration_range_18_block0,
ops_arg_dat(u2_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_19_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel019, "Derivative evaluation CD T_B0 xi2 ", opensbliblock00, 3, iteration_range_19_block0,
ops_arg_dat(T_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk11_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_32_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel032, "Convective terms", opensbliblock00, 3, iteration_range_32_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(H_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(phi_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

int iteration_range_33_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel033, "Viscous terms", opensbliblock00, 3, iteration_range_33_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(SD111_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(T_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(mu_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_44_22_27, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk11_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_22_11, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_00_22_11, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00_00_22_11, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_idx());

int iteration_range_71_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel071, "Temporal solution advancement", opensbliblock00, 3, iteration_range_71_block0,
ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_gbl(&rkA[stage], 1, "double", OPS_READ),
ops_arg_gbl(&rkB[stage], 1, "double", OPS_READ));

ops_halo_transfer(periodicBC_direction0_side0_35_block0);
ops_halo_transfer(periodicBC_direction0_side1_36_block0);
int iteration_range_37_block0[] = {-4, block0np0 + 5, 0, 1, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel037, "IsothermalWall boundary dir1 side0", opensbliblock00, 3, iteration_range_37_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_51_00_15, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW));

int iteration_range_38_block0[] = {-4, block0np0 + 5, block0np1 - 1, block0np1, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel038, "IsothermalWall boundary dir1 side1", opensbliblock00, 3, iteration_range_38_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_15_00_15, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_55_00_23, "double", OPS_RW));

ops_halo_transfer(periodicBC_direction2_side0_39_block0);
ops_halo_transfer(periodicBC_direction2_side1_40_block0);
}
if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_62_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel062, "User kernel: Block 0: Zero the filter array", opensbliblock00, 3, iteration_range_62_block0,
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

ops_halo_transfer(periodicBC_direction0_side0_76_block0);
ops_halo_transfer(periodicBC_direction0_side1_77_block0);
ops_halo_transfer(periodicBC_direction2_side0_78_block0);
ops_halo_transfer(periodicBC_direction2_side1_79_block0);
}

int iteration_range_69_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel069, "User kernel: None", opensbliblock00, 3, iteration_range_69_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(E_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(M_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(TT_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(T_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(a_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(mu_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(p_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(pp_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhomean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2u2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u0u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u1u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u1u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2u2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));

if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_63_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel063, "User kernel: Block 0: DRP filter calculation direction x", opensbliblock00, 3, iteration_range_63_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_44_00_00_19, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_44_00_00_19, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_44_00_00_19, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_44_00_00_19, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_44_00_00_19, "double", OPS_READ),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));
}

if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_64_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel064, "User kernel: Block 0: DRP filter update direction x", opensbliblock00, 3, iteration_range_64_block0,
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));
}

if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_65_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel065, "User kernel: Block 0: DRP filter calculation direction y", opensbliblock00, 3, iteration_range_65_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_44_00_19, "double", OPS_READ),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());
}

if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_66_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel066, "User kernel: Block 0: DRP filter update direction y", opensbliblock00, 3, iteration_range_66_block0,
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));
}

if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_67_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel067, "User kernel: Block 0: DRP filter calculation direction z", opensbliblock00, 3, iteration_range_67_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_44_19, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_44_19, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_44_19, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_44_19, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_44_19, "double", OPS_READ),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));
}

if (fmod(1 + iter,filter_frequency) == 0){
int iteration_range_68_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel068, "User kernel: Block 0: DRP filter update direction z", opensbliblock00, 3, iteration_range_68_block0,
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));
}

ops_halo_transfer(periodicBC_direction0_side0_72_block0);
ops_halo_transfer(periodicBC_direction0_side1_73_block0);
ops_halo_transfer(periodicBC_direction2_side0_74_block0);
ops_halo_transfer(periodicBC_direction2_side1_75_block0);
int iteration_range_52_block0[] = {-4, block0np0 + 5, -4, block0np1 + 5, -4, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel052, "User kernel: Constituent Relations evaluation", opensbliblock00, 3, iteration_range_52_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(a_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));

int iteration_range_56_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel056, "UserDefinedEquations evaluation", opensbliblock00, 3, iteration_range_56_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_22_44_22_35, "double", OPS_READ),
ops_arg_dat(kappa_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());

int iteration_range_57_block0[] = {-5, block0np0 + 5, -5, block0np1 + 5, -5, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel057, "User kernel: Zero the work arrays", opensbliblock00, 3, iteration_range_57_block0,
ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_58_block0[] = {-1, block0np0 + 1, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel058, "User kernel: WENO reconstruction direction 0", opensbliblock00, 3, iteration_range_58_block0,
ops_arg_dat(a_B0, 1, stencil_0_01_00_00_5, "double", OPS_READ),
ops_arg_dat(kappa_B0, 1, stencil_0_32_00_00_13, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_34_00_00_17, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_01_00_00_5, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_01_00_00_5, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_59_block0[] = {0, block0np0, -1, block0np1 + 1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel059, "User kernel: WENO reconstruction direction 1", opensbliblock00, 3, iteration_range_59_block0,
ops_arg_dat(a_B0, 1, stencil_0_00_01_00_5, "double", OPS_READ),
ops_arg_dat(kappa_B0, 1, stencil_0_00_32_00_13, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_01_00_5, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_34_00_17, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_01_00_5, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_60_block0[] = {0, block0np0, 0, block0np1, -1, block0np2 + 1};
ops_par_loop(opensbliblock00Kernel060, "User kernel: WENO reconstruction direction 2", opensbliblock00, 3, iteration_range_60_block0,
ops_arg_dat(a_B0, 1, stencil_0_00_00_01_5, "double", OPS_READ),
ops_arg_dat(kappa_B0, 1, stencil_0_00_00_32_13, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_01_5, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_01_5, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_34_17, "double", OPS_READ),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_61_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel061, "User kernel: Non-linear filter application", opensbliblock00, 3, iteration_range_61_block0,
ops_arg_dat(D11_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(kappa_B0, 1, stencil_0_12_12_12_21, "double", OPS_READ),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(WENO_filter_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_idx());

if (fmod(1 + iter,write_output_file) == 0 || iter == 0){
HDF5_IO_Write_1_opensbliblock00_dynamic(opensbliblock00, iter, rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0, x0_B0, x1_B0, x2_B0, D11_B0, HDF5_timing);
}

if (fmod(1 + iter,write_slices) == 0){
char slice_name0[80];
sprintf(slice_name0, "%d", iter + 1);
ops_write_plane_group_hdf5({{2, block0np2/2}}, slice_name0, {{rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0, WENO_filter_B0}});
}

if (fmod(1 + iter,write_slices) == 0){
char slice_name0[80];
sprintf(slice_name0, "%d", iter + 1);
ops_write_plane_group_hdf5({{1, 20}}, slice_name0, {{rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0, WENO_filter_B0}});
ops_write_plane_group_hdf5({{1, block0np1 - 20}}, slice_name0, {{rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0, WENO_filter_B0}});
}

}
ops_timers(&cpu_end0, &elapsed_end0);
ops_printf("\nTimings are:\n");
ops_printf("-----------------------------------------\n");
ops_printf("Total Wall time %lf\n",elapsed_end0-elapsed_start0);

int iteration_range_70_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel070, "User kernel: None", opensbliblock00, 3, iteration_range_70_block0,
ops_arg_dat(E_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(M_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(TT_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(T_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(a_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(mu_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(p_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(pp_mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhomean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2u2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u0u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u1u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u1u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2u0mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2u1mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(u2u2mean_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));

HDF5_IO_Write_1_opensbliblock00(opensbliblock00, rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0, x0_B0, x1_B0, x2_B0, D11_B0, HDF5_timing);
HDF5_IO_Write_0_opensbliblock00(opensbliblock00, u2u2mean_B0, rhou2mean_B0, rhou2u1mean_B0, u0u0mean_B0, rhou1u0mean_B0, E_mean_B0, u1u0mean_B0, u1u1mean_B0, rhou2u0mean_B0, rhou0mean_B0, rhou1mean_B0, pp_mean_B0, rhou2u2mean_B0, u2mean_B0, M_mean_B0, u2u0mean_B0, p_mean_B0, a_mean_B0, T_mean_B0, rhou0u0mean_B0, rhomean_B0, mu_mean_B0, u2u1mean_B0, TT_mean_B0, rhou1u1mean_B0, u0mean_B0, u1mean_B0, D11_B0, HDF5_timing);
ops_exit();
//Main program end 
}