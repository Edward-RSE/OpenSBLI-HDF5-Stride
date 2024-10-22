#ifndef OPENSBLIBLOCK00_KERNEL_H
#define OPENSBLIBLOCK00_KERNEL_H
 void opensbliblock00Kernel015(ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0,
ACC<double> &x0_B0, ACC<double> &x1_B0, const int *idx)
{
   double d = 0.0;
   double p = 0.0;
   double u0 = 0.0;
   double u1 = 0.0;
   x0_B0(0,0) = Delta0block0*idx[0];

   x1_B0(0,0) = Delta1block0*idx[1];

   d = 1 + 0.2*sin(M_PI*(x0_B0(0,0) + x1_B0(0,0)));

   u0 = 1.00000000000000;

   u1 = -0.500000000000000;

   p = 1.00000000000000;

   rho_B0(0,0) = d;

   rhou0_B0(0,0) = d*u0;

   rhou1_B0(0,0) = d*u1;

   rhoE_B0(0,0) = p*inv_gamma_m1 + 0.5*d*((u0*u0) + (u1*u1));

}

 void opensbliblock00Kernel016(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double> &rhou0_RKold_B0, ACC<double>
&rhou1_RKold_B0)
{
   rho_RKold_B0(0,0) = rho_B0(0,0);

   rhou0_RKold_B0(0,0) = rhou0_B0(0,0);

   rhou1_RKold_B0(0,0) = rhou1_B0(0,0);

   rhoE_RKold_B0(0,0) = rhoE_B0(0,0);

}

void opensbliblock00Kernel004(const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, ACC<double> &u0_B0)
{
   u0_B0(0,0) = rhou0_B0(0,0)/rho_B0(0,0);

}

void opensbliblock00Kernel008(const ACC<double> &rho_B0, const ACC<double> &rhou1_B0, ACC<double> &u1_B0)
{
   u1_B0(0,0) = rhou1_B0(0,0)/rho_B0(0,0);

}

 void opensbliblock00Kernel006(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, ACC<double> &p_B0)
{
    p_B0(0,0) = (-1 + gama)*(-(1.0/2.0)*(u0_B0(0,0)*u0_B0(0,0))*rho_B0(0,0) -
      (1.0/2.0)*(u1_B0(0,0)*u1_B0(0,0))*rho_B0(0,0) + rhoE_B0(0,0));

}

void opensbliblock00Kernel003(const ACC<double> &p_B0, const ACC<double> &rho_B0, ACC<double> &a_B0)
{
   a_B0(0,0) = sqrt((gama*p_B0(0,0)/rho_B0(0,0)));

}

 void opensbliblock00Kernel000(const ACC<double> &a_B0, const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const
ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, ACC<double> &wk0_B0, ACC<double> &wk1_B0, ACC<double> &wk2_B0, ACC<double> &wk3_B0)
{
   double AVG_0_0_LEV_00 = 0.0;
   double AVG_0_0_LEV_01 = 0.0;
   double AVG_0_0_LEV_02 = 0.0;
   double AVG_0_0_LEV_03 = 0.0;
   double AVG_0_0_LEV_10 = 0.0;
   double AVG_0_0_LEV_12 = 0.0;
   double AVG_0_0_LEV_20 = 0.0;
   double AVG_0_0_LEV_21 = 0.0;
   double AVG_0_0_LEV_22 = 0.0;
   double AVG_0_0_LEV_23 = 0.0;
   double AVG_0_0_LEV_30 = 0.0;
   double AVG_0_0_LEV_31 = 0.0;
   double AVG_0_0_LEV_32 = 0.0;
   double AVG_0_0_LEV_33 = 0.0;
   double AVG_0_a = 0.0;
   double AVG_0_rho = 0.0;
   double AVG_0_u0 = 0.0;
   double AVG_0_u1 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double omega_2 = 0.0;
   AVG_0_u0 = ((1.0/2.0))*(u0_B0(0,0) + u0_B0(1,0));

   AVG_0_a = ((1.0/2.0))*(a_B0(0,0) + a_B0(1,0));

   AVG_0_rho = ((1.0/2.0))*(rho_B0(0,0) + rho_B0(1,0));

   AVG_0_u1 = ((1.0/2.0))*(u1_B0(0,0) + u1_B0(1,0));

   inv_AVG_a = 1.0/(AVG_0_a);

   inv_AVG_rho = 1.0/(AVG_0_rho);

    AVG_0_0_LEV_00 = -(1.0/2.0)*(-2 - (AVG_0_u0*AVG_0_u0)*(inv_AVG_a*inv_AVG_a) -
      (AVG_0_u1*AVG_0_u1)*(inv_AVG_a*inv_AVG_a) + (AVG_0_u0*AVG_0_u0)*(inv_AVG_a*inv_AVG_a)*gama +
      (AVG_0_u1*AVG_0_u1)*(inv_AVG_a*inv_AVG_a)*gama);

   AVG_0_0_LEV_01 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u0;

   AVG_0_0_LEV_02 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u1;

   AVG_0_0_LEV_03 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

   AVG_0_0_LEV_10 = AVG_0_u1*inv_AVG_rho;

   AVG_0_0_LEV_12 = -inv_AVG_rho;

    AVG_0_0_LEV_20 = -0.353553390593274*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) - (AVG_0_u0*AVG_0_u0)*gama -
      (AVG_0_u1*AVG_0_u1)*gama + 2*AVG_0_a*AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

   AVG_0_0_LEV_21 = 0.707106781186547*(-gama*AVG_0_u0 + AVG_0_a + AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

   AVG_0_0_LEV_22 = -0.707106781186547*gamma_m1*AVG_0_u1*inv_AVG_a*inv_AVG_rho;

   AVG_0_0_LEV_23 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

    AVG_0_0_LEV_30 = 0.353553390593274*(-(AVG_0_u0*AVG_0_u0) - (AVG_0_u1*AVG_0_u1) + (AVG_0_u0*AVG_0_u0)*gama +
      (AVG_0_u1*AVG_0_u1)*gama + 2*AVG_0_a*AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

   AVG_0_0_LEV_31 = -0.707106781186547*(-AVG_0_u0 + gama*AVG_0_u0 + AVG_0_a)*inv_AVG_a*inv_AVG_rho;

   AVG_0_0_LEV_32 = -0.707106781186547*gamma_m1*AVG_0_u1*inv_AVG_a*inv_AVG_rho;

   AVG_0_0_LEV_33 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

    CF_00 = p_B0(-2,0)*AVG_0_0_LEV_01 + rhou0_B0(-2,0)*AVG_0_0_LEV_00 + p_B0(-2,0)*u0_B0(-2,0)*AVG_0_0_LEV_03 +
      u0_B0(-2,0)*rhoE_B0(-2,0)*AVG_0_0_LEV_03 + u0_B0(-2,0)*rhou0_B0(-2,0)*AVG_0_0_LEV_01 +
      u0_B0(-2,0)*rhou1_B0(-2,0)*AVG_0_0_LEV_02;

   CF_10 = rhou0_B0(-2,0)*AVG_0_0_LEV_10 + u0_B0(-2,0)*rhou1_B0(-2,0)*AVG_0_0_LEV_12;

    CF_20 = p_B0(-2,0)*AVG_0_0_LEV_21 + rhou0_B0(-2,0)*AVG_0_0_LEV_20 + p_B0(-2,0)*u0_B0(-2,0)*AVG_0_0_LEV_23 +
      u0_B0(-2,0)*rhoE_B0(-2,0)*AVG_0_0_LEV_23 + u0_B0(-2,0)*rhou0_B0(-2,0)*AVG_0_0_LEV_21 +
      u0_B0(-2,0)*rhou1_B0(-2,0)*AVG_0_0_LEV_22;

    CF_30 = p_B0(-2,0)*AVG_0_0_LEV_31 + rhou0_B0(-2,0)*AVG_0_0_LEV_30 + p_B0(-2,0)*u0_B0(-2,0)*AVG_0_0_LEV_33 +
      u0_B0(-2,0)*rhoE_B0(-2,0)*AVG_0_0_LEV_33 + u0_B0(-2,0)*rhou0_B0(-2,0)*AVG_0_0_LEV_31 +
      u0_B0(-2,0)*rhou1_B0(-2,0)*AVG_0_0_LEV_32;

    CS_00 = rho_B0(-2,0)*AVG_0_0_LEV_00 + rhoE_B0(-2,0)*AVG_0_0_LEV_03 + rhou0_B0(-2,0)*AVG_0_0_LEV_01 +
      rhou1_B0(-2,0)*AVG_0_0_LEV_02;

   CS_10 = rho_B0(-2,0)*AVG_0_0_LEV_10 + rhou1_B0(-2,0)*AVG_0_0_LEV_12;

    CS_20 = rho_B0(-2,0)*AVG_0_0_LEV_20 + rhoE_B0(-2,0)*AVG_0_0_LEV_23 + rhou0_B0(-2,0)*AVG_0_0_LEV_21 +
      rhou1_B0(-2,0)*AVG_0_0_LEV_22;

    CS_30 = rho_B0(-2,0)*AVG_0_0_LEV_30 + rhoE_B0(-2,0)*AVG_0_0_LEV_33 + rhou0_B0(-2,0)*AVG_0_0_LEV_31 +
      rhou1_B0(-2,0)*AVG_0_0_LEV_32;

    CF_01 = p_B0(-1,0)*AVG_0_0_LEV_01 + rhou0_B0(-1,0)*AVG_0_0_LEV_00 + p_B0(-1,0)*u0_B0(-1,0)*AVG_0_0_LEV_03 +
      u0_B0(-1,0)*rhoE_B0(-1,0)*AVG_0_0_LEV_03 + u0_B0(-1,0)*rhou0_B0(-1,0)*AVG_0_0_LEV_01 +
      u0_B0(-1,0)*rhou1_B0(-1,0)*AVG_0_0_LEV_02;

   CF_11 = rhou0_B0(-1,0)*AVG_0_0_LEV_10 + u0_B0(-1,0)*rhou1_B0(-1,0)*AVG_0_0_LEV_12;

    CF_21 = p_B0(-1,0)*AVG_0_0_LEV_21 + rhou0_B0(-1,0)*AVG_0_0_LEV_20 + p_B0(-1,0)*u0_B0(-1,0)*AVG_0_0_LEV_23 +
      u0_B0(-1,0)*rhoE_B0(-1,0)*AVG_0_0_LEV_23 + u0_B0(-1,0)*rhou0_B0(-1,0)*AVG_0_0_LEV_21 +
      u0_B0(-1,0)*rhou1_B0(-1,0)*AVG_0_0_LEV_22;

    CF_31 = p_B0(-1,0)*AVG_0_0_LEV_31 + rhou0_B0(-1,0)*AVG_0_0_LEV_30 + p_B0(-1,0)*u0_B0(-1,0)*AVG_0_0_LEV_33 +
      u0_B0(-1,0)*rhoE_B0(-1,0)*AVG_0_0_LEV_33 + u0_B0(-1,0)*rhou0_B0(-1,0)*AVG_0_0_LEV_31 +
      u0_B0(-1,0)*rhou1_B0(-1,0)*AVG_0_0_LEV_32;

    CS_01 = rho_B0(-1,0)*AVG_0_0_LEV_00 + rhoE_B0(-1,0)*AVG_0_0_LEV_03 + rhou0_B0(-1,0)*AVG_0_0_LEV_01 +
      rhou1_B0(-1,0)*AVG_0_0_LEV_02;

   CS_11 = rho_B0(-1,0)*AVG_0_0_LEV_10 + rhou1_B0(-1,0)*AVG_0_0_LEV_12;

    CS_21 = rho_B0(-1,0)*AVG_0_0_LEV_20 + rhoE_B0(-1,0)*AVG_0_0_LEV_23 + rhou0_B0(-1,0)*AVG_0_0_LEV_21 +
      rhou1_B0(-1,0)*AVG_0_0_LEV_22;

    CS_31 = rho_B0(-1,0)*AVG_0_0_LEV_30 + rhoE_B0(-1,0)*AVG_0_0_LEV_33 + rhou0_B0(-1,0)*AVG_0_0_LEV_31 +
      rhou1_B0(-1,0)*AVG_0_0_LEV_32;

    CF_02 = p_B0(0,0)*AVG_0_0_LEV_01 + rhou0_B0(0,0)*AVG_0_0_LEV_00 + p_B0(0,0)*u0_B0(0,0)*AVG_0_0_LEV_03 +
      u0_B0(0,0)*rhoE_B0(0,0)*AVG_0_0_LEV_03 + u0_B0(0,0)*rhou0_B0(0,0)*AVG_0_0_LEV_01 +
      u0_B0(0,0)*rhou1_B0(0,0)*AVG_0_0_LEV_02;

   CF_12 = rhou0_B0(0,0)*AVG_0_0_LEV_10 + u0_B0(0,0)*rhou1_B0(0,0)*AVG_0_0_LEV_12;

    CF_22 = p_B0(0,0)*AVG_0_0_LEV_21 + rhou0_B0(0,0)*AVG_0_0_LEV_20 + p_B0(0,0)*u0_B0(0,0)*AVG_0_0_LEV_23 +
      u0_B0(0,0)*rhoE_B0(0,0)*AVG_0_0_LEV_23 + u0_B0(0,0)*rhou0_B0(0,0)*AVG_0_0_LEV_21 +
      u0_B0(0,0)*rhou1_B0(0,0)*AVG_0_0_LEV_22;

    CF_32 = p_B0(0,0)*AVG_0_0_LEV_31 + rhou0_B0(0,0)*AVG_0_0_LEV_30 + p_B0(0,0)*u0_B0(0,0)*AVG_0_0_LEV_33 +
      u0_B0(0,0)*rhoE_B0(0,0)*AVG_0_0_LEV_33 + u0_B0(0,0)*rhou0_B0(0,0)*AVG_0_0_LEV_31 +
      u0_B0(0,0)*rhou1_B0(0,0)*AVG_0_0_LEV_32;

    CS_02 = rho_B0(0,0)*AVG_0_0_LEV_00 + rhoE_B0(0,0)*AVG_0_0_LEV_03 + rhou0_B0(0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(0,0)*AVG_0_0_LEV_02;

   CS_12 = rho_B0(0,0)*AVG_0_0_LEV_10 + rhou1_B0(0,0)*AVG_0_0_LEV_12;

    CS_22 = rho_B0(0,0)*AVG_0_0_LEV_20 + rhoE_B0(0,0)*AVG_0_0_LEV_23 + rhou0_B0(0,0)*AVG_0_0_LEV_21 +
      rhou1_B0(0,0)*AVG_0_0_LEV_22;

    CS_32 = rho_B0(0,0)*AVG_0_0_LEV_30 + rhoE_B0(0,0)*AVG_0_0_LEV_33 + rhou0_B0(0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(0,0)*AVG_0_0_LEV_32;

    CF_03 = p_B0(1,0)*AVG_0_0_LEV_01 + rhou0_B0(1,0)*AVG_0_0_LEV_00 + p_B0(1,0)*u0_B0(1,0)*AVG_0_0_LEV_03 +
      u0_B0(1,0)*rhoE_B0(1,0)*AVG_0_0_LEV_03 + u0_B0(1,0)*rhou0_B0(1,0)*AVG_0_0_LEV_01 +
      u0_B0(1,0)*rhou1_B0(1,0)*AVG_0_0_LEV_02;

   CF_13 = rhou0_B0(1,0)*AVG_0_0_LEV_10 + u0_B0(1,0)*rhou1_B0(1,0)*AVG_0_0_LEV_12;

    CF_23 = p_B0(1,0)*AVG_0_0_LEV_21 + rhou0_B0(1,0)*AVG_0_0_LEV_20 + p_B0(1,0)*u0_B0(1,0)*AVG_0_0_LEV_23 +
      u0_B0(1,0)*rhoE_B0(1,0)*AVG_0_0_LEV_23 + u0_B0(1,0)*rhou0_B0(1,0)*AVG_0_0_LEV_21 +
      u0_B0(1,0)*rhou1_B0(1,0)*AVG_0_0_LEV_22;

    CF_33 = p_B0(1,0)*AVG_0_0_LEV_31 + rhou0_B0(1,0)*AVG_0_0_LEV_30 + p_B0(1,0)*u0_B0(1,0)*AVG_0_0_LEV_33 +
      u0_B0(1,0)*rhoE_B0(1,0)*AVG_0_0_LEV_33 + u0_B0(1,0)*rhou0_B0(1,0)*AVG_0_0_LEV_31 +
      u0_B0(1,0)*rhou1_B0(1,0)*AVG_0_0_LEV_32;

    CS_03 = rho_B0(1,0)*AVG_0_0_LEV_00 + rhoE_B0(1,0)*AVG_0_0_LEV_03 + rhou0_B0(1,0)*AVG_0_0_LEV_01 +
      rhou1_B0(1,0)*AVG_0_0_LEV_02;

   CS_13 = rho_B0(1,0)*AVG_0_0_LEV_10 + rhou1_B0(1,0)*AVG_0_0_LEV_12;

    CS_23 = rho_B0(1,0)*AVG_0_0_LEV_20 + rhoE_B0(1,0)*AVG_0_0_LEV_23 + rhou0_B0(1,0)*AVG_0_0_LEV_21 +
      rhou1_B0(1,0)*AVG_0_0_LEV_22;

    CS_33 = rho_B0(1,0)*AVG_0_0_LEV_30 + rhoE_B0(1,0)*AVG_0_0_LEV_33 + rhou0_B0(1,0)*AVG_0_0_LEV_31 +
      rhou1_B0(1,0)*AVG_0_0_LEV_32;

    CF_04 = p_B0(2,0)*AVG_0_0_LEV_01 + rhou0_B0(2,0)*AVG_0_0_LEV_00 + p_B0(2,0)*u0_B0(2,0)*AVG_0_0_LEV_03 +
      u0_B0(2,0)*rhoE_B0(2,0)*AVG_0_0_LEV_03 + u0_B0(2,0)*rhou0_B0(2,0)*AVG_0_0_LEV_01 +
      u0_B0(2,0)*rhou1_B0(2,0)*AVG_0_0_LEV_02;

   CF_14 = rhou0_B0(2,0)*AVG_0_0_LEV_10 + u0_B0(2,0)*rhou1_B0(2,0)*AVG_0_0_LEV_12;

    CF_24 = p_B0(2,0)*AVG_0_0_LEV_21 + rhou0_B0(2,0)*AVG_0_0_LEV_20 + p_B0(2,0)*u0_B0(2,0)*AVG_0_0_LEV_23 +
      u0_B0(2,0)*rhoE_B0(2,0)*AVG_0_0_LEV_23 + u0_B0(2,0)*rhou0_B0(2,0)*AVG_0_0_LEV_21 +
      u0_B0(2,0)*rhou1_B0(2,0)*AVG_0_0_LEV_22;

    CF_34 = p_B0(2,0)*AVG_0_0_LEV_31 + rhou0_B0(2,0)*AVG_0_0_LEV_30 + p_B0(2,0)*u0_B0(2,0)*AVG_0_0_LEV_33 +
      u0_B0(2,0)*rhoE_B0(2,0)*AVG_0_0_LEV_33 + u0_B0(2,0)*rhou0_B0(2,0)*AVG_0_0_LEV_31 +
      u0_B0(2,0)*rhou1_B0(2,0)*AVG_0_0_LEV_32;

    CS_04 = rho_B0(2,0)*AVG_0_0_LEV_00 + rhoE_B0(2,0)*AVG_0_0_LEV_03 + rhou0_B0(2,0)*AVG_0_0_LEV_01 +
      rhou1_B0(2,0)*AVG_0_0_LEV_02;

   CS_14 = rho_B0(2,0)*AVG_0_0_LEV_10 + rhou1_B0(2,0)*AVG_0_0_LEV_12;

    CS_24 = rho_B0(2,0)*AVG_0_0_LEV_20 + rhoE_B0(2,0)*AVG_0_0_LEV_23 + rhou0_B0(2,0)*AVG_0_0_LEV_21 +
      rhou1_B0(2,0)*AVG_0_0_LEV_22;

    CS_34 = rho_B0(2,0)*AVG_0_0_LEV_30 + rhoE_B0(2,0)*AVG_0_0_LEV_33 + rhou0_B0(2,0)*AVG_0_0_LEV_31 +
      rhou1_B0(2,0)*AVG_0_0_LEV_32;

    CF_05 = p_B0(3,0)*AVG_0_0_LEV_01 + rhou0_B0(3,0)*AVG_0_0_LEV_00 + p_B0(3,0)*u0_B0(3,0)*AVG_0_0_LEV_03 +
      u0_B0(3,0)*rhoE_B0(3,0)*AVG_0_0_LEV_03 + u0_B0(3,0)*rhou0_B0(3,0)*AVG_0_0_LEV_01 +
      u0_B0(3,0)*rhou1_B0(3,0)*AVG_0_0_LEV_02;

   CF_15 = rhou0_B0(3,0)*AVG_0_0_LEV_10 + u0_B0(3,0)*rhou1_B0(3,0)*AVG_0_0_LEV_12;

    CF_25 = p_B0(3,0)*AVG_0_0_LEV_21 + rhou0_B0(3,0)*AVG_0_0_LEV_20 + p_B0(3,0)*u0_B0(3,0)*AVG_0_0_LEV_23 +
      u0_B0(3,0)*rhoE_B0(3,0)*AVG_0_0_LEV_23 + u0_B0(3,0)*rhou0_B0(3,0)*AVG_0_0_LEV_21 +
      u0_B0(3,0)*rhou1_B0(3,0)*AVG_0_0_LEV_22;

    CF_35 = p_B0(3,0)*AVG_0_0_LEV_31 + rhou0_B0(3,0)*AVG_0_0_LEV_30 + p_B0(3,0)*u0_B0(3,0)*AVG_0_0_LEV_33 +
      u0_B0(3,0)*rhoE_B0(3,0)*AVG_0_0_LEV_33 + u0_B0(3,0)*rhou0_B0(3,0)*AVG_0_0_LEV_31 +
      u0_B0(3,0)*rhou1_B0(3,0)*AVG_0_0_LEV_32;

    CS_05 = rho_B0(3,0)*AVG_0_0_LEV_00 + rhoE_B0(3,0)*AVG_0_0_LEV_03 + rhou0_B0(3,0)*AVG_0_0_LEV_01 +
      rhou1_B0(3,0)*AVG_0_0_LEV_02;

   CS_15 = rho_B0(3,0)*AVG_0_0_LEV_10 + rhou1_B0(3,0)*AVG_0_0_LEV_12;

    CS_25 = rho_B0(3,0)*AVG_0_0_LEV_20 + rhoE_B0(3,0)*AVG_0_0_LEV_23 + rhou0_B0(3,0)*AVG_0_0_LEV_21 +
      rhou1_B0(3,0)*AVG_0_0_LEV_22;

    CS_35 = rho_B0(3,0)*AVG_0_0_LEV_30 + rhoE_B0(3,0)*AVG_0_0_LEV_33 + rhou0_B0(3,0)*AVG_0_0_LEV_31 +
      rhou1_B0(3,0)*AVG_0_0_LEV_32;

   max_lambda_00 = shock_filter_control*fmax(fabs(u0_B0(1,0)), fabs(u0_B0(0,0)));

   max_lambda_11 = max_lambda_00;

   max_lambda_22 = shock_filter_control*fmax(fabs(a_B0(0,0) + u0_B0(0,0)), fabs(a_B0(1,0) + u0_B0(1,0)));

   max_lambda_33 = shock_filter_control*fmax(fabs(-u0_B0(1,0) + a_B0(1,0)), fabs(-u0_B0(0,0) + a_B0(0,0)));

    beta_0 = ((1.0/3.0))*((CS_04*max_lambda_00 + CF_04)*(CS_04*max_lambda_00 + CF_04)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_04*max_lambda_00 + CF_04) + ((25.0/6.0))*(CS_03*max_lambda_00 +
      CF_03))*(CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(CS_02*max_lambda_00 +
      CF_02)*(-(31.0/6.0)*(CS_03*max_lambda_00 + CF_03) + ((5.0/3.0))*(CS_02*max_lambda_00 + CF_02) +
      ((11.0/6.0))*(CS_04*max_lambda_00 + CF_04));

    beta_1 = ((1.0/3.0))*((CS_03*max_lambda_00 + CF_03)*(CS_03*max_lambda_00 + CF_03)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_03*max_lambda_00 + CF_03) + ((13.0/6.0))*(CS_02*max_lambda_00 +
      CF_02))*(CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(CS_01*max_lambda_00 +
      CF_01)*(-(13.0/6.0)*(CS_02*max_lambda_00 + CF_02) + ((2.0/3.0))*(CS_01*max_lambda_00 + CF_01) +
      ((5.0/6.0))*(CS_03*max_lambda_00 + CF_03));

    beta_2 = ((5.0/6.0))*((CS_02*max_lambda_00 + CF_02)*(CS_02*max_lambda_00 + CF_02)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_00*max_lambda_00 + CF_00) + ((11.0/6.0))*(CS_02*max_lambda_00 +
      CF_02))*(CS_00*max_lambda_00 + CF_00) + ((1.0/2.0))*(CS_01*max_lambda_00 +
      CF_01)*(-(31.0/6.0)*(CS_02*max_lambda_00 + CF_02) - (19.0/6.0)*(CS_00*max_lambda_00 + CF_00) +
      ((25.0/6.0))*(CS_01*max_lambda_00 + CF_01));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_0 = (-(7.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_00*max_lambda_00 + CF_00) +
      ((11.0/12.0))*(CS_02*max_lambda_00 + CF_02))*omega_2 + (-(1.0/12.0)*(CS_01*max_lambda_00 + CF_01) +
      ((1.0/6.0))*(CS_03*max_lambda_00 + CF_03) + ((5.0/12.0))*(CS_02*max_lambda_00 + CF_02))*omega_1 +
      (-(1.0/12.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(CS_03*max_lambda_00 + CF_03))*omega_0 + Recon_0;

    beta_0 = ((1.0/3.0))*((-CS_05*max_lambda_00 + CF_05)*(-CS_05*max_lambda_00 + CF_05)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_05*max_lambda_00 + CF_05) + ((25.0/6.0))*(-CS_04*max_lambda_00 +
      CF_04))*(-CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(-CS_03*max_lambda_00 +
      CF_03)*(-(31.0/6.0)*(-CS_04*max_lambda_00 + CF_04) + ((5.0/3.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((11.0/6.0))*(-CS_05*max_lambda_00 + CF_05));

    beta_1 = ((1.0/3.0))*((-CS_04*max_lambda_00 + CF_04)*(-CS_04*max_lambda_00 + CF_04)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_04*max_lambda_00 + CF_04) + ((13.0/6.0))*(-CS_03*max_lambda_00 +
      CF_03))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
      CF_02)*(-(13.0/6.0)*(-CS_03*max_lambda_00 + CF_03) + ((2.0/3.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/6.0))*(-CS_04*max_lambda_00 + CF_04));

    beta_2 = ((5.0/6.0))*((-CS_03*max_lambda_00 + CF_03)*(-CS_03*max_lambda_00 + CF_03)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_03*max_lambda_00 + CF_03) + ((25.0/6.0))*(-CS_02*max_lambda_00 +
      CF_02))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01)*(-(19.0/6.0)*(-CS_02*max_lambda_00 + CF_02) + ((2.0/3.0))*(-CS_01*max_lambda_00 + CF_01) +
      ((11.0/6.0))*(-CS_03*max_lambda_00 + CF_03));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_0 = (-(7.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_05*max_lambda_00 + CF_05) +
      ((11.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*omega_0 + (-(1.0/12.0)*(-CS_01*max_lambda_00 + CF_01) +
      ((1.0/6.0))*(-CS_03*max_lambda_00 + CF_03) + ((5.0/12.0))*(-CS_02*max_lambda_00 + CF_02))*omega_2 +
      (-(1.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*omega_1 + Recon_0;

    beta_0 = ((1.0/3.0))*((CS_14*max_lambda_11 + CF_14)*(CS_14*max_lambda_11 + CF_14)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_14*max_lambda_11 + CF_14) + ((25.0/6.0))*(CS_13*max_lambda_11 +
      CF_13))*(CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(CS_12*max_lambda_11 +
      CF_12)*(-(31.0/6.0)*(CS_13*max_lambda_11 + CF_13) + ((5.0/3.0))*(CS_12*max_lambda_11 + CF_12) +
      ((11.0/6.0))*(CS_14*max_lambda_11 + CF_14));

    beta_1 = ((1.0/3.0))*((CS_13*max_lambda_11 + CF_13)*(CS_13*max_lambda_11 + CF_13)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_13*max_lambda_11 + CF_13) + ((13.0/6.0))*(CS_12*max_lambda_11 +
      CF_12))*(CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(CS_11*max_lambda_11 +
      CF_11)*(-(13.0/6.0)*(CS_12*max_lambda_11 + CF_12) + ((2.0/3.0))*(CS_11*max_lambda_11 + CF_11) +
      ((5.0/6.0))*(CS_13*max_lambda_11 + CF_13));

    beta_2 = ((5.0/6.0))*((CS_12*max_lambda_11 + CF_12)*(CS_12*max_lambda_11 + CF_12)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_10*max_lambda_11 + CF_10) + ((11.0/6.0))*(CS_12*max_lambda_11 +
      CF_12))*(CS_10*max_lambda_11 + CF_10) + ((1.0/2.0))*(CS_11*max_lambda_11 +
      CF_11)*(-(31.0/6.0)*(CS_12*max_lambda_11 + CF_12) - (19.0/6.0)*(CS_10*max_lambda_11 + CF_10) +
      ((25.0/6.0))*(CS_11*max_lambda_11 + CF_11));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_1 = (-(7.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_10*max_lambda_11 + CF_10) +
      ((11.0/12.0))*(CS_12*max_lambda_11 + CF_12))*omega_2 + (-(1.0/12.0)*(CS_11*max_lambda_11 + CF_11) +
      ((1.0/6.0))*(CS_13*max_lambda_11 + CF_13) + ((5.0/12.0))*(CS_12*max_lambda_11 + CF_12))*omega_1 +
      (-(1.0/12.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(CS_13*max_lambda_11 + CF_13))*omega_0 + Recon_1;

    beta_0 = ((1.0/3.0))*((-CS_15*max_lambda_11 + CF_15)*(-CS_15*max_lambda_11 + CF_15)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_15*max_lambda_11 + CF_15) + ((25.0/6.0))*(-CS_14*max_lambda_11 +
      CF_14))*(-CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(-CS_13*max_lambda_11 +
      CF_13)*(-(31.0/6.0)*(-CS_14*max_lambda_11 + CF_14) + ((5.0/3.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((11.0/6.0))*(-CS_15*max_lambda_11 + CF_15));

    beta_1 = ((1.0/3.0))*((-CS_14*max_lambda_11 + CF_14)*(-CS_14*max_lambda_11 + CF_14)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_14*max_lambda_11 + CF_14) + ((13.0/6.0))*(-CS_13*max_lambda_11 +
      CF_13))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
      CF_12)*(-(13.0/6.0)*(-CS_13*max_lambda_11 + CF_13) + ((2.0/3.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/6.0))*(-CS_14*max_lambda_11 + CF_14));

    beta_2 = ((5.0/6.0))*((-CS_13*max_lambda_11 + CF_13)*(-CS_13*max_lambda_11 + CF_13)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_13*max_lambda_11 + CF_13) + ((25.0/6.0))*(-CS_12*max_lambda_11 +
      CF_12))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11)*(-(19.0/6.0)*(-CS_12*max_lambda_11 + CF_12) + ((2.0/3.0))*(-CS_11*max_lambda_11 + CF_11) +
      ((11.0/6.0))*(-CS_13*max_lambda_11 + CF_13));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_1 = (-(7.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_15*max_lambda_11 + CF_15) +
      ((11.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*omega_0 + (-(1.0/12.0)*(-CS_11*max_lambda_11 + CF_11) +
      ((1.0/6.0))*(-CS_13*max_lambda_11 + CF_13) + ((5.0/12.0))*(-CS_12*max_lambda_11 + CF_12))*omega_2 +
      (-(1.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*omega_1 + Recon_1;

    beta_0 = ((1.0/3.0))*((CS_24*max_lambda_22 + CF_24)*(CS_24*max_lambda_22 + CF_24)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_24*max_lambda_22 + CF_24) + ((25.0/6.0))*(CS_23*max_lambda_22 +
      CF_23))*(CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(CS_22*max_lambda_22 +
      CF_22)*(-(31.0/6.0)*(CS_23*max_lambda_22 + CF_23) + ((5.0/3.0))*(CS_22*max_lambda_22 + CF_22) +
      ((11.0/6.0))*(CS_24*max_lambda_22 + CF_24));

    beta_1 = ((1.0/3.0))*((CS_23*max_lambda_22 + CF_23)*(CS_23*max_lambda_22 + CF_23)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_23*max_lambda_22 + CF_23) + ((13.0/6.0))*(CS_22*max_lambda_22 +
      CF_22))*(CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(CS_21*max_lambda_22 +
      CF_21)*(-(13.0/6.0)*(CS_22*max_lambda_22 + CF_22) + ((2.0/3.0))*(CS_21*max_lambda_22 + CF_21) +
      ((5.0/6.0))*(CS_23*max_lambda_22 + CF_23));

    beta_2 = ((5.0/6.0))*((CS_22*max_lambda_22 + CF_22)*(CS_22*max_lambda_22 + CF_22)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_20*max_lambda_22 + CF_20) + ((11.0/6.0))*(CS_22*max_lambda_22 +
      CF_22))*(CS_20*max_lambda_22 + CF_20) + ((1.0/2.0))*(CS_21*max_lambda_22 +
      CF_21)*(-(31.0/6.0)*(CS_22*max_lambda_22 + CF_22) - (19.0/6.0)*(CS_20*max_lambda_22 + CF_20) +
      ((25.0/6.0))*(CS_21*max_lambda_22 + CF_21));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_2 = (-(7.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_20*max_lambda_22 + CF_20) +
      ((11.0/12.0))*(CS_22*max_lambda_22 + CF_22))*omega_2 + (-(1.0/12.0)*(CS_21*max_lambda_22 + CF_21) +
      ((1.0/6.0))*(CS_23*max_lambda_22 + CF_23) + ((5.0/12.0))*(CS_22*max_lambda_22 + CF_22))*omega_1 +
      (-(1.0/12.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(CS_23*max_lambda_22 + CF_23))*omega_0 + Recon_2;

    beta_0 = ((1.0/3.0))*((-CS_25*max_lambda_22 + CF_25)*(-CS_25*max_lambda_22 + CF_25)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_25*max_lambda_22 + CF_25) + ((25.0/6.0))*(-CS_24*max_lambda_22 +
      CF_24))*(-CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(-CS_23*max_lambda_22 +
      CF_23)*(-(31.0/6.0)*(-CS_24*max_lambda_22 + CF_24) + ((5.0/3.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((11.0/6.0))*(-CS_25*max_lambda_22 + CF_25));

    beta_1 = ((1.0/3.0))*((-CS_24*max_lambda_22 + CF_24)*(-CS_24*max_lambda_22 + CF_24)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_24*max_lambda_22 + CF_24) + ((13.0/6.0))*(-CS_23*max_lambda_22 +
      CF_23))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
      CF_22)*(-(13.0/6.0)*(-CS_23*max_lambda_22 + CF_23) + ((2.0/3.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/6.0))*(-CS_24*max_lambda_22 + CF_24));

    beta_2 = ((5.0/6.0))*((-CS_23*max_lambda_22 + CF_23)*(-CS_23*max_lambda_22 + CF_23)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_23*max_lambda_22 + CF_23) + ((25.0/6.0))*(-CS_22*max_lambda_22 +
      CF_22))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21)*(-(19.0/6.0)*(-CS_22*max_lambda_22 + CF_22) + ((2.0/3.0))*(-CS_21*max_lambda_22 + CF_21) +
      ((11.0/6.0))*(-CS_23*max_lambda_22 + CF_23));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_2 = (-(7.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_25*max_lambda_22 + CF_25) +
      ((11.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*omega_0 + (-(1.0/12.0)*(-CS_21*max_lambda_22 + CF_21) +
      ((1.0/6.0))*(-CS_23*max_lambda_22 + CF_23) + ((5.0/12.0))*(-CS_22*max_lambda_22 + CF_22))*omega_2 +
      (-(1.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*omega_1 + Recon_2;

    beta_0 = ((1.0/3.0))*((CS_34*max_lambda_33 + CF_34)*(CS_34*max_lambda_33 + CF_34)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_34*max_lambda_33 + CF_34) + ((25.0/6.0))*(CS_33*max_lambda_33 +
      CF_33))*(CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(CS_32*max_lambda_33 +
      CF_32)*(-(31.0/6.0)*(CS_33*max_lambda_33 + CF_33) + ((5.0/3.0))*(CS_32*max_lambda_33 + CF_32) +
      ((11.0/6.0))*(CS_34*max_lambda_33 + CF_34));

    beta_1 = ((1.0/3.0))*((CS_33*max_lambda_33 + CF_33)*(CS_33*max_lambda_33 + CF_33)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_33*max_lambda_33 + CF_33) + ((13.0/6.0))*(CS_32*max_lambda_33 +
      CF_32))*(CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(CS_31*max_lambda_33 +
      CF_31)*(-(13.0/6.0)*(CS_32*max_lambda_33 + CF_32) + ((2.0/3.0))*(CS_31*max_lambda_33 + CF_31) +
      ((5.0/6.0))*(CS_33*max_lambda_33 + CF_33));

    beta_2 = ((5.0/6.0))*((CS_32*max_lambda_33 + CF_32)*(CS_32*max_lambda_33 + CF_32)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_30*max_lambda_33 + CF_30) + ((11.0/6.0))*(CS_32*max_lambda_33 +
      CF_32))*(CS_30*max_lambda_33 + CF_30) + ((1.0/2.0))*(CS_31*max_lambda_33 +
      CF_31)*(-(31.0/6.0)*(CS_32*max_lambda_33 + CF_32) - (19.0/6.0)*(CS_30*max_lambda_33 + CF_30) +
      ((25.0/6.0))*(CS_31*max_lambda_33 + CF_31));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_3 = (-(7.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_30*max_lambda_33 + CF_30) +
      ((11.0/12.0))*(CS_32*max_lambda_33 + CF_32))*omega_2 + (-(1.0/12.0)*(CS_31*max_lambda_33 + CF_31) +
      ((1.0/6.0))*(CS_33*max_lambda_33 + CF_33) + ((5.0/12.0))*(CS_32*max_lambda_33 + CF_32))*omega_1 +
      (-(1.0/12.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(CS_33*max_lambda_33 + CF_33))*omega_0 + Recon_3;

    beta_0 = ((1.0/3.0))*((-CS_35*max_lambda_33 + CF_35)*(-CS_35*max_lambda_33 + CF_35)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_35*max_lambda_33 + CF_35) + ((25.0/6.0))*(-CS_34*max_lambda_33 +
      CF_34))*(-CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(-CS_33*max_lambda_33 +
      CF_33)*(-(31.0/6.0)*(-CS_34*max_lambda_33 + CF_34) + ((5.0/3.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((11.0/6.0))*(-CS_35*max_lambda_33 + CF_35));

    beta_1 = ((1.0/3.0))*((-CS_34*max_lambda_33 + CF_34)*(-CS_34*max_lambda_33 + CF_34)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_34*max_lambda_33 + CF_34) + ((13.0/6.0))*(-CS_33*max_lambda_33 +
      CF_33))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
      CF_32)*(-(13.0/6.0)*(-CS_33*max_lambda_33 + CF_33) + ((2.0/3.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/6.0))*(-CS_34*max_lambda_33 + CF_34));

    beta_2 = ((5.0/6.0))*((-CS_33*max_lambda_33 + CF_33)*(-CS_33*max_lambda_33 + CF_33)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_33*max_lambda_33 + CF_33) + ((25.0/6.0))*(-CS_32*max_lambda_33 +
      CF_32))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31)*(-(19.0/6.0)*(-CS_32*max_lambda_33 + CF_32) + ((2.0/3.0))*(-CS_31*max_lambda_33 + CF_31) +
      ((11.0/6.0))*(-CS_33*max_lambda_33 + CF_33));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_3 = (-(7.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_35*max_lambda_33 + CF_35) +
      ((11.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*omega_0 + (-(1.0/12.0)*(-CS_31*max_lambda_33 + CF_31) +
      ((1.0/6.0))*(-CS_33*max_lambda_33 + CF_33) + ((5.0/12.0))*(-CS_32*max_lambda_33 + CF_32))*omega_2 +
      (-(1.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*omega_1 + Recon_3;

    wk0_B0(0,0) = 0.707106781186547*AVG_0_rho*Recon_2*inv_AVG_a + 0.707106781186547*AVG_0_rho*Recon_3*inv_AVG_a +
      Recon_0;

    wk1_B0(0,0) = AVG_0_u0*Recon_0 + 0.707106781186547*(-AVG_0_a + AVG_0_u0)*AVG_0_rho*Recon_3*inv_AVG_a +
      0.707106781186547*(AVG_0_a + AVG_0_u0)*AVG_0_rho*Recon_2*inv_AVG_a;

    wk2_B0(0,0) = AVG_0_u1*Recon_0 - AVG_0_rho*Recon_1 + 0.707106781186547*AVG_0_rho*AVG_0_u1*Recon_2*inv_AVG_a +
      0.707106781186547*AVG_0_rho*AVG_0_u1*Recon_3*inv_AVG_a;

    wk3_B0(0,0) = (((1.0/2.0))*(AVG_0_u0*AVG_0_u0) + ((1.0/2.0))*(AVG_0_u1*AVG_0_u1))*Recon_0 -
      AVG_0_rho*AVG_0_u1*Recon_1 + 0.707106781186547*(((AVG_0_a*AVG_0_a) + ((1.0/2.0))*((AVG_0_u0*AVG_0_u0) +
      (AVG_0_u1*AVG_0_u1))*gamma_m1)*invgamma_m1 + AVG_0_a*AVG_0_u0)*AVG_0_rho*Recon_2*inv_AVG_a +
      0.707106781186547*(((AVG_0_a*AVG_0_a) + ((1.0/2.0))*((AVG_0_u0*AVG_0_u0) +
      (AVG_0_u1*AVG_0_u1))*gamma_m1)*invgamma_m1 - AVG_0_a*AVG_0_u0)*AVG_0_rho*Recon_3*inv_AVG_a;

}

 void opensbliblock00Kernel001(const ACC<double> &a_B0, const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const
ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, ACC<double> &wk4_B0, ACC<double> &wk5_B0, ACC<double> &wk6_B0, ACC<double> &wk7_B0)
{
   double AVG_1_1_LEV_00 = 0.0;
   double AVG_1_1_LEV_01 = 0.0;
   double AVG_1_1_LEV_02 = 0.0;
   double AVG_1_1_LEV_03 = 0.0;
   double AVG_1_1_LEV_10 = 0.0;
   double AVG_1_1_LEV_11 = 0.0;
   double AVG_1_1_LEV_20 = 0.0;
   double AVG_1_1_LEV_21 = 0.0;
   double AVG_1_1_LEV_22 = 0.0;
   double AVG_1_1_LEV_23 = 0.0;
   double AVG_1_1_LEV_30 = 0.0;
   double AVG_1_1_LEV_31 = 0.0;
   double AVG_1_1_LEV_32 = 0.0;
   double AVG_1_1_LEV_33 = 0.0;
   double AVG_1_a = 0.0;
   double AVG_1_rho = 0.0;
   double AVG_1_u0 = 0.0;
   double AVG_1_u1 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double omega_2 = 0.0;
   AVG_1_u0 = ((1.0/2.0))*(u0_B0(0,0) + u0_B0(0,1));

   AVG_1_a = ((1.0/2.0))*(a_B0(0,0) + a_B0(0,1));

   AVG_1_rho = ((1.0/2.0))*(rho_B0(0,0) + rho_B0(0,1));

   AVG_1_u1 = ((1.0/2.0))*(u1_B0(0,0) + u1_B0(0,1));

   inv_AVG_a = 1.0/(AVG_1_a);

   inv_AVG_rho = 1.0/(AVG_1_rho);

    AVG_1_1_LEV_00 = -(1.0/2.0)*(-2 - (AVG_1_u0*AVG_1_u0)*(inv_AVG_a*inv_AVG_a) -
      (AVG_1_u1*AVG_1_u1)*(inv_AVG_a*inv_AVG_a) + (AVG_1_u0*AVG_1_u0)*(inv_AVG_a*inv_AVG_a)*gama +
      (AVG_1_u1*AVG_1_u1)*(inv_AVG_a*inv_AVG_a)*gama);

   AVG_1_1_LEV_01 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u0;

   AVG_1_1_LEV_02 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u1;

   AVG_1_1_LEV_03 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

   AVG_1_1_LEV_10 = -AVG_1_u0*inv_AVG_rho;

   AVG_1_1_LEV_11 = inv_AVG_rho;

    AVG_1_1_LEV_20 = -0.353553390593274*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) - (AVG_1_u0*AVG_1_u0)*gama -
      (AVG_1_u1*AVG_1_u1)*gama + 2*AVG_1_a*AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

   AVG_1_1_LEV_21 = -0.707106781186547*gamma_m1*AVG_1_u0*inv_AVG_a*inv_AVG_rho;

   AVG_1_1_LEV_22 = 0.707106781186547*(-gama*AVG_1_u1 + AVG_1_a + AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

   AVG_1_1_LEV_23 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

    AVG_1_1_LEV_30 = 0.353553390593274*(-(AVG_1_u0*AVG_1_u0) - (AVG_1_u1*AVG_1_u1) + (AVG_1_u0*AVG_1_u0)*gama +
      (AVG_1_u1*AVG_1_u1)*gama + 2*AVG_1_a*AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

   AVG_1_1_LEV_31 = -0.707106781186547*gamma_m1*AVG_1_u0*inv_AVG_a*inv_AVG_rho;

   AVG_1_1_LEV_32 = -0.707106781186547*(-AVG_1_u1 + gama*AVG_1_u1 + AVG_1_a)*inv_AVG_a*inv_AVG_rho;

   AVG_1_1_LEV_33 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

    CF_00 = p_B0(0,-2)*AVG_1_1_LEV_02 + rhou1_B0(0,-2)*AVG_1_1_LEV_00 + p_B0(0,-2)*u1_B0(0,-2)*AVG_1_1_LEV_03 +
      u1_B0(0,-2)*rhoE_B0(0,-2)*AVG_1_1_LEV_03 + u1_B0(0,-2)*rhou0_B0(0,-2)*AVG_1_1_LEV_01 +
      u1_B0(0,-2)*rhou1_B0(0,-2)*AVG_1_1_LEV_02;

   CF_10 = rhou1_B0(0,-2)*AVG_1_1_LEV_10 + u1_B0(0,-2)*rhou0_B0(0,-2)*AVG_1_1_LEV_11;

    CF_20 = p_B0(0,-2)*AVG_1_1_LEV_22 + rhou1_B0(0,-2)*AVG_1_1_LEV_20 + p_B0(0,-2)*u1_B0(0,-2)*AVG_1_1_LEV_23 +
      u1_B0(0,-2)*rhoE_B0(0,-2)*AVG_1_1_LEV_23 + u1_B0(0,-2)*rhou0_B0(0,-2)*AVG_1_1_LEV_21 +
      u1_B0(0,-2)*rhou1_B0(0,-2)*AVG_1_1_LEV_22;

    CF_30 = p_B0(0,-2)*AVG_1_1_LEV_32 + rhou1_B0(0,-2)*AVG_1_1_LEV_30 + p_B0(0,-2)*u1_B0(0,-2)*AVG_1_1_LEV_33 +
      u1_B0(0,-2)*rhoE_B0(0,-2)*AVG_1_1_LEV_33 + u1_B0(0,-2)*rhou0_B0(0,-2)*AVG_1_1_LEV_31 +
      u1_B0(0,-2)*rhou1_B0(0,-2)*AVG_1_1_LEV_32;

    CS_00 = rho_B0(0,-2)*AVG_1_1_LEV_00 + rhoE_B0(0,-2)*AVG_1_1_LEV_03 + rhou0_B0(0,-2)*AVG_1_1_LEV_01 +
      rhou1_B0(0,-2)*AVG_1_1_LEV_02;

   CS_10 = rho_B0(0,-2)*AVG_1_1_LEV_10 + rhou0_B0(0,-2)*AVG_1_1_LEV_11;

    CS_20 = rho_B0(0,-2)*AVG_1_1_LEV_20 + rhoE_B0(0,-2)*AVG_1_1_LEV_23 + rhou0_B0(0,-2)*AVG_1_1_LEV_21 +
      rhou1_B0(0,-2)*AVG_1_1_LEV_22;

    CS_30 = rho_B0(0,-2)*AVG_1_1_LEV_30 + rhoE_B0(0,-2)*AVG_1_1_LEV_33 + rhou0_B0(0,-2)*AVG_1_1_LEV_31 +
      rhou1_B0(0,-2)*AVG_1_1_LEV_32;

    CF_01 = p_B0(0,-1)*AVG_1_1_LEV_02 + rhou1_B0(0,-1)*AVG_1_1_LEV_00 + p_B0(0,-1)*u1_B0(0,-1)*AVG_1_1_LEV_03 +
      u1_B0(0,-1)*rhoE_B0(0,-1)*AVG_1_1_LEV_03 + u1_B0(0,-1)*rhou0_B0(0,-1)*AVG_1_1_LEV_01 +
      u1_B0(0,-1)*rhou1_B0(0,-1)*AVG_1_1_LEV_02;

   CF_11 = rhou1_B0(0,-1)*AVG_1_1_LEV_10 + u1_B0(0,-1)*rhou0_B0(0,-1)*AVG_1_1_LEV_11;

    CF_21 = p_B0(0,-1)*AVG_1_1_LEV_22 + rhou1_B0(0,-1)*AVG_1_1_LEV_20 + p_B0(0,-1)*u1_B0(0,-1)*AVG_1_1_LEV_23 +
      u1_B0(0,-1)*rhoE_B0(0,-1)*AVG_1_1_LEV_23 + u1_B0(0,-1)*rhou0_B0(0,-1)*AVG_1_1_LEV_21 +
      u1_B0(0,-1)*rhou1_B0(0,-1)*AVG_1_1_LEV_22;

    CF_31 = p_B0(0,-1)*AVG_1_1_LEV_32 + rhou1_B0(0,-1)*AVG_1_1_LEV_30 + p_B0(0,-1)*u1_B0(0,-1)*AVG_1_1_LEV_33 +
      u1_B0(0,-1)*rhoE_B0(0,-1)*AVG_1_1_LEV_33 + u1_B0(0,-1)*rhou0_B0(0,-1)*AVG_1_1_LEV_31 +
      u1_B0(0,-1)*rhou1_B0(0,-1)*AVG_1_1_LEV_32;

    CS_01 = rho_B0(0,-1)*AVG_1_1_LEV_00 + rhoE_B0(0,-1)*AVG_1_1_LEV_03 + rhou0_B0(0,-1)*AVG_1_1_LEV_01 +
      rhou1_B0(0,-1)*AVG_1_1_LEV_02;

   CS_11 = rho_B0(0,-1)*AVG_1_1_LEV_10 + rhou0_B0(0,-1)*AVG_1_1_LEV_11;

    CS_21 = rho_B0(0,-1)*AVG_1_1_LEV_20 + rhoE_B0(0,-1)*AVG_1_1_LEV_23 + rhou0_B0(0,-1)*AVG_1_1_LEV_21 +
      rhou1_B0(0,-1)*AVG_1_1_LEV_22;

    CS_31 = rho_B0(0,-1)*AVG_1_1_LEV_30 + rhoE_B0(0,-1)*AVG_1_1_LEV_33 + rhou0_B0(0,-1)*AVG_1_1_LEV_31 +
      rhou1_B0(0,-1)*AVG_1_1_LEV_32;

    CF_02 = p_B0(0,0)*AVG_1_1_LEV_02 + rhou1_B0(0,0)*AVG_1_1_LEV_00 + p_B0(0,0)*u1_B0(0,0)*AVG_1_1_LEV_03 +
      u1_B0(0,0)*rhoE_B0(0,0)*AVG_1_1_LEV_03 + u1_B0(0,0)*rhou0_B0(0,0)*AVG_1_1_LEV_01 +
      u1_B0(0,0)*rhou1_B0(0,0)*AVG_1_1_LEV_02;

   CF_12 = rhou1_B0(0,0)*AVG_1_1_LEV_10 + u1_B0(0,0)*rhou0_B0(0,0)*AVG_1_1_LEV_11;

    CF_22 = p_B0(0,0)*AVG_1_1_LEV_22 + rhou1_B0(0,0)*AVG_1_1_LEV_20 + p_B0(0,0)*u1_B0(0,0)*AVG_1_1_LEV_23 +
      u1_B0(0,0)*rhoE_B0(0,0)*AVG_1_1_LEV_23 + u1_B0(0,0)*rhou0_B0(0,0)*AVG_1_1_LEV_21 +
      u1_B0(0,0)*rhou1_B0(0,0)*AVG_1_1_LEV_22;

    CF_32 = p_B0(0,0)*AVG_1_1_LEV_32 + rhou1_B0(0,0)*AVG_1_1_LEV_30 + p_B0(0,0)*u1_B0(0,0)*AVG_1_1_LEV_33 +
      u1_B0(0,0)*rhoE_B0(0,0)*AVG_1_1_LEV_33 + u1_B0(0,0)*rhou0_B0(0,0)*AVG_1_1_LEV_31 +
      u1_B0(0,0)*rhou1_B0(0,0)*AVG_1_1_LEV_32;

    CS_02 = rho_B0(0,0)*AVG_1_1_LEV_00 + rhoE_B0(0,0)*AVG_1_1_LEV_03 + rhou0_B0(0,0)*AVG_1_1_LEV_01 +
      rhou1_B0(0,0)*AVG_1_1_LEV_02;

   CS_12 = rho_B0(0,0)*AVG_1_1_LEV_10 + rhou0_B0(0,0)*AVG_1_1_LEV_11;

    CS_22 = rho_B0(0,0)*AVG_1_1_LEV_20 + rhoE_B0(0,0)*AVG_1_1_LEV_23 + rhou0_B0(0,0)*AVG_1_1_LEV_21 +
      rhou1_B0(0,0)*AVG_1_1_LEV_22;

    CS_32 = rho_B0(0,0)*AVG_1_1_LEV_30 + rhoE_B0(0,0)*AVG_1_1_LEV_33 + rhou0_B0(0,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,0)*AVG_1_1_LEV_32;

    CF_03 = p_B0(0,1)*AVG_1_1_LEV_02 + rhou1_B0(0,1)*AVG_1_1_LEV_00 + p_B0(0,1)*u1_B0(0,1)*AVG_1_1_LEV_03 +
      u1_B0(0,1)*rhoE_B0(0,1)*AVG_1_1_LEV_03 + u1_B0(0,1)*rhou0_B0(0,1)*AVG_1_1_LEV_01 +
      u1_B0(0,1)*rhou1_B0(0,1)*AVG_1_1_LEV_02;

   CF_13 = rhou1_B0(0,1)*AVG_1_1_LEV_10 + u1_B0(0,1)*rhou0_B0(0,1)*AVG_1_1_LEV_11;

    CF_23 = p_B0(0,1)*AVG_1_1_LEV_22 + rhou1_B0(0,1)*AVG_1_1_LEV_20 + p_B0(0,1)*u1_B0(0,1)*AVG_1_1_LEV_23 +
      u1_B0(0,1)*rhoE_B0(0,1)*AVG_1_1_LEV_23 + u1_B0(0,1)*rhou0_B0(0,1)*AVG_1_1_LEV_21 +
      u1_B0(0,1)*rhou1_B0(0,1)*AVG_1_1_LEV_22;

    CF_33 = p_B0(0,1)*AVG_1_1_LEV_32 + rhou1_B0(0,1)*AVG_1_1_LEV_30 + p_B0(0,1)*u1_B0(0,1)*AVG_1_1_LEV_33 +
      u1_B0(0,1)*rhoE_B0(0,1)*AVG_1_1_LEV_33 + u1_B0(0,1)*rhou0_B0(0,1)*AVG_1_1_LEV_31 +
      u1_B0(0,1)*rhou1_B0(0,1)*AVG_1_1_LEV_32;

    CS_03 = rho_B0(0,1)*AVG_1_1_LEV_00 + rhoE_B0(0,1)*AVG_1_1_LEV_03 + rhou0_B0(0,1)*AVG_1_1_LEV_01 +
      rhou1_B0(0,1)*AVG_1_1_LEV_02;

   CS_13 = rho_B0(0,1)*AVG_1_1_LEV_10 + rhou0_B0(0,1)*AVG_1_1_LEV_11;

    CS_23 = rho_B0(0,1)*AVG_1_1_LEV_20 + rhoE_B0(0,1)*AVG_1_1_LEV_23 + rhou0_B0(0,1)*AVG_1_1_LEV_21 +
      rhou1_B0(0,1)*AVG_1_1_LEV_22;

    CS_33 = rho_B0(0,1)*AVG_1_1_LEV_30 + rhoE_B0(0,1)*AVG_1_1_LEV_33 + rhou0_B0(0,1)*AVG_1_1_LEV_31 +
      rhou1_B0(0,1)*AVG_1_1_LEV_32;

    CF_04 = p_B0(0,2)*AVG_1_1_LEV_02 + rhou1_B0(0,2)*AVG_1_1_LEV_00 + p_B0(0,2)*u1_B0(0,2)*AVG_1_1_LEV_03 +
      u1_B0(0,2)*rhoE_B0(0,2)*AVG_1_1_LEV_03 + u1_B0(0,2)*rhou0_B0(0,2)*AVG_1_1_LEV_01 +
      u1_B0(0,2)*rhou1_B0(0,2)*AVG_1_1_LEV_02;

   CF_14 = rhou1_B0(0,2)*AVG_1_1_LEV_10 + u1_B0(0,2)*rhou0_B0(0,2)*AVG_1_1_LEV_11;

    CF_24 = p_B0(0,2)*AVG_1_1_LEV_22 + rhou1_B0(0,2)*AVG_1_1_LEV_20 + p_B0(0,2)*u1_B0(0,2)*AVG_1_1_LEV_23 +
      u1_B0(0,2)*rhoE_B0(0,2)*AVG_1_1_LEV_23 + u1_B0(0,2)*rhou0_B0(0,2)*AVG_1_1_LEV_21 +
      u1_B0(0,2)*rhou1_B0(0,2)*AVG_1_1_LEV_22;

    CF_34 = p_B0(0,2)*AVG_1_1_LEV_32 + rhou1_B0(0,2)*AVG_1_1_LEV_30 + p_B0(0,2)*u1_B0(0,2)*AVG_1_1_LEV_33 +
      u1_B0(0,2)*rhoE_B0(0,2)*AVG_1_1_LEV_33 + u1_B0(0,2)*rhou0_B0(0,2)*AVG_1_1_LEV_31 +
      u1_B0(0,2)*rhou1_B0(0,2)*AVG_1_1_LEV_32;

    CS_04 = rho_B0(0,2)*AVG_1_1_LEV_00 + rhoE_B0(0,2)*AVG_1_1_LEV_03 + rhou0_B0(0,2)*AVG_1_1_LEV_01 +
      rhou1_B0(0,2)*AVG_1_1_LEV_02;

   CS_14 = rho_B0(0,2)*AVG_1_1_LEV_10 + rhou0_B0(0,2)*AVG_1_1_LEV_11;

    CS_24 = rho_B0(0,2)*AVG_1_1_LEV_20 + rhoE_B0(0,2)*AVG_1_1_LEV_23 + rhou0_B0(0,2)*AVG_1_1_LEV_21 +
      rhou1_B0(0,2)*AVG_1_1_LEV_22;

    CS_34 = rho_B0(0,2)*AVG_1_1_LEV_30 + rhoE_B0(0,2)*AVG_1_1_LEV_33 + rhou0_B0(0,2)*AVG_1_1_LEV_31 +
      rhou1_B0(0,2)*AVG_1_1_LEV_32;

    CF_05 = p_B0(0,3)*AVG_1_1_LEV_02 + rhou1_B0(0,3)*AVG_1_1_LEV_00 + p_B0(0,3)*u1_B0(0,3)*AVG_1_1_LEV_03 +
      u1_B0(0,3)*rhoE_B0(0,3)*AVG_1_1_LEV_03 + u1_B0(0,3)*rhou0_B0(0,3)*AVG_1_1_LEV_01 +
      u1_B0(0,3)*rhou1_B0(0,3)*AVG_1_1_LEV_02;

   CF_15 = rhou1_B0(0,3)*AVG_1_1_LEV_10 + u1_B0(0,3)*rhou0_B0(0,3)*AVG_1_1_LEV_11;

    CF_25 = p_B0(0,3)*AVG_1_1_LEV_22 + rhou1_B0(0,3)*AVG_1_1_LEV_20 + p_B0(0,3)*u1_B0(0,3)*AVG_1_1_LEV_23 +
      u1_B0(0,3)*rhoE_B0(0,3)*AVG_1_1_LEV_23 + u1_B0(0,3)*rhou0_B0(0,3)*AVG_1_1_LEV_21 +
      u1_B0(0,3)*rhou1_B0(0,3)*AVG_1_1_LEV_22;

    CF_35 = p_B0(0,3)*AVG_1_1_LEV_32 + rhou1_B0(0,3)*AVG_1_1_LEV_30 + p_B0(0,3)*u1_B0(0,3)*AVG_1_1_LEV_33 +
      u1_B0(0,3)*rhoE_B0(0,3)*AVG_1_1_LEV_33 + u1_B0(0,3)*rhou0_B0(0,3)*AVG_1_1_LEV_31 +
      u1_B0(0,3)*rhou1_B0(0,3)*AVG_1_1_LEV_32;

    CS_05 = rho_B0(0,3)*AVG_1_1_LEV_00 + rhoE_B0(0,3)*AVG_1_1_LEV_03 + rhou0_B0(0,3)*AVG_1_1_LEV_01 +
      rhou1_B0(0,3)*AVG_1_1_LEV_02;

   CS_15 = rho_B0(0,3)*AVG_1_1_LEV_10 + rhou0_B0(0,3)*AVG_1_1_LEV_11;

    CS_25 = rho_B0(0,3)*AVG_1_1_LEV_20 + rhoE_B0(0,3)*AVG_1_1_LEV_23 + rhou0_B0(0,3)*AVG_1_1_LEV_21 +
      rhou1_B0(0,3)*AVG_1_1_LEV_22;

    CS_35 = rho_B0(0,3)*AVG_1_1_LEV_30 + rhoE_B0(0,3)*AVG_1_1_LEV_33 + rhou0_B0(0,3)*AVG_1_1_LEV_31 +
      rhou1_B0(0,3)*AVG_1_1_LEV_32;

   max_lambda_00 = shock_filter_control*fmax(fabs(u1_B0(0,0)), fabs(u1_B0(0,1)));

   max_lambda_11 = max_lambda_00;

   max_lambda_22 = shock_filter_control*fmax(fabs(a_B0(0,1) + u1_B0(0,1)), fabs(a_B0(0,0) + u1_B0(0,0)));

   max_lambda_33 = shock_filter_control*fmax(fabs(-u1_B0(0,1) + a_B0(0,1)), fabs(-u1_B0(0,0) + a_B0(0,0)));

    beta_0 = ((1.0/3.0))*((CS_04*max_lambda_00 + CF_04)*(CS_04*max_lambda_00 + CF_04)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_04*max_lambda_00 + CF_04) + ((25.0/6.0))*(CS_03*max_lambda_00 +
      CF_03))*(CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(CS_02*max_lambda_00 +
      CF_02)*(-(31.0/6.0)*(CS_03*max_lambda_00 + CF_03) + ((5.0/3.0))*(CS_02*max_lambda_00 + CF_02) +
      ((11.0/6.0))*(CS_04*max_lambda_00 + CF_04));

    beta_1 = ((1.0/3.0))*((CS_03*max_lambda_00 + CF_03)*(CS_03*max_lambda_00 + CF_03)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_03*max_lambda_00 + CF_03) + ((13.0/6.0))*(CS_02*max_lambda_00 +
      CF_02))*(CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(CS_01*max_lambda_00 +
      CF_01)*(-(13.0/6.0)*(CS_02*max_lambda_00 + CF_02) + ((2.0/3.0))*(CS_01*max_lambda_00 + CF_01) +
      ((5.0/6.0))*(CS_03*max_lambda_00 + CF_03));

    beta_2 = ((5.0/6.0))*((CS_02*max_lambda_00 + CF_02)*(CS_02*max_lambda_00 + CF_02)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_00*max_lambda_00 + CF_00) + ((11.0/6.0))*(CS_02*max_lambda_00 +
      CF_02))*(CS_00*max_lambda_00 + CF_00) + ((1.0/2.0))*(CS_01*max_lambda_00 +
      CF_01)*(-(31.0/6.0)*(CS_02*max_lambda_00 + CF_02) - (19.0/6.0)*(CS_00*max_lambda_00 + CF_00) +
      ((25.0/6.0))*(CS_01*max_lambda_00 + CF_01));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_0 = (-(7.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_00*max_lambda_00 + CF_00) +
      ((11.0/12.0))*(CS_02*max_lambda_00 + CF_02))*omega_2 + (-(1.0/12.0)*(CS_01*max_lambda_00 + CF_01) +
      ((1.0/6.0))*(CS_03*max_lambda_00 + CF_03) + ((5.0/12.0))*(CS_02*max_lambda_00 + CF_02))*omega_1 +
      (-(1.0/12.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(CS_03*max_lambda_00 + CF_03))*omega_0 + Recon_0;

    beta_0 = ((1.0/3.0))*((-CS_05*max_lambda_00 + CF_05)*(-CS_05*max_lambda_00 + CF_05)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_05*max_lambda_00 + CF_05) + ((25.0/6.0))*(-CS_04*max_lambda_00 +
      CF_04))*(-CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(-CS_03*max_lambda_00 +
      CF_03)*(-(31.0/6.0)*(-CS_04*max_lambda_00 + CF_04) + ((5.0/3.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((11.0/6.0))*(-CS_05*max_lambda_00 + CF_05));

    beta_1 = ((1.0/3.0))*((-CS_04*max_lambda_00 + CF_04)*(-CS_04*max_lambda_00 + CF_04)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_04*max_lambda_00 + CF_04) + ((13.0/6.0))*(-CS_03*max_lambda_00 +
      CF_03))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
      CF_02)*(-(13.0/6.0)*(-CS_03*max_lambda_00 + CF_03) + ((2.0/3.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/6.0))*(-CS_04*max_lambda_00 + CF_04));

    beta_2 = ((5.0/6.0))*((-CS_03*max_lambda_00 + CF_03)*(-CS_03*max_lambda_00 + CF_03)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_03*max_lambda_00 + CF_03) + ((25.0/6.0))*(-CS_02*max_lambda_00 +
      CF_02))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01)*(-(19.0/6.0)*(-CS_02*max_lambda_00 + CF_02) + ((2.0/3.0))*(-CS_01*max_lambda_00 + CF_01) +
      ((11.0/6.0))*(-CS_03*max_lambda_00 + CF_03));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_0 = (-(7.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_05*max_lambda_00 + CF_05) +
      ((11.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*omega_0 + (-(1.0/12.0)*(-CS_01*max_lambda_00 + CF_01) +
      ((1.0/6.0))*(-CS_03*max_lambda_00 + CF_03) + ((5.0/12.0))*(-CS_02*max_lambda_00 + CF_02))*omega_2 +
      (-(1.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*omega_1 + Recon_0;

    beta_0 = ((1.0/3.0))*((CS_14*max_lambda_11 + CF_14)*(CS_14*max_lambda_11 + CF_14)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_14*max_lambda_11 + CF_14) + ((25.0/6.0))*(CS_13*max_lambda_11 +
      CF_13))*(CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(CS_12*max_lambda_11 +
      CF_12)*(-(31.0/6.0)*(CS_13*max_lambda_11 + CF_13) + ((5.0/3.0))*(CS_12*max_lambda_11 + CF_12) +
      ((11.0/6.0))*(CS_14*max_lambda_11 + CF_14));

    beta_1 = ((1.0/3.0))*((CS_13*max_lambda_11 + CF_13)*(CS_13*max_lambda_11 + CF_13)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_13*max_lambda_11 + CF_13) + ((13.0/6.0))*(CS_12*max_lambda_11 +
      CF_12))*(CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(CS_11*max_lambda_11 +
      CF_11)*(-(13.0/6.0)*(CS_12*max_lambda_11 + CF_12) + ((2.0/3.0))*(CS_11*max_lambda_11 + CF_11) +
      ((5.0/6.0))*(CS_13*max_lambda_11 + CF_13));

    beta_2 = ((5.0/6.0))*((CS_12*max_lambda_11 + CF_12)*(CS_12*max_lambda_11 + CF_12)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_10*max_lambda_11 + CF_10) + ((11.0/6.0))*(CS_12*max_lambda_11 +
      CF_12))*(CS_10*max_lambda_11 + CF_10) + ((1.0/2.0))*(CS_11*max_lambda_11 +
      CF_11)*(-(31.0/6.0)*(CS_12*max_lambda_11 + CF_12) - (19.0/6.0)*(CS_10*max_lambda_11 + CF_10) +
      ((25.0/6.0))*(CS_11*max_lambda_11 + CF_11));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_1 = (-(7.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_10*max_lambda_11 + CF_10) +
      ((11.0/12.0))*(CS_12*max_lambda_11 + CF_12))*omega_2 + (-(1.0/12.0)*(CS_11*max_lambda_11 + CF_11) +
      ((1.0/6.0))*(CS_13*max_lambda_11 + CF_13) + ((5.0/12.0))*(CS_12*max_lambda_11 + CF_12))*omega_1 +
      (-(1.0/12.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(CS_13*max_lambda_11 + CF_13))*omega_0 + Recon_1;

    beta_0 = ((1.0/3.0))*((-CS_15*max_lambda_11 + CF_15)*(-CS_15*max_lambda_11 + CF_15)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_15*max_lambda_11 + CF_15) + ((25.0/6.0))*(-CS_14*max_lambda_11 +
      CF_14))*(-CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(-CS_13*max_lambda_11 +
      CF_13)*(-(31.0/6.0)*(-CS_14*max_lambda_11 + CF_14) + ((5.0/3.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((11.0/6.0))*(-CS_15*max_lambda_11 + CF_15));

    beta_1 = ((1.0/3.0))*((-CS_14*max_lambda_11 + CF_14)*(-CS_14*max_lambda_11 + CF_14)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_14*max_lambda_11 + CF_14) + ((13.0/6.0))*(-CS_13*max_lambda_11 +
      CF_13))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
      CF_12)*(-(13.0/6.0)*(-CS_13*max_lambda_11 + CF_13) + ((2.0/3.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/6.0))*(-CS_14*max_lambda_11 + CF_14));

    beta_2 = ((5.0/6.0))*((-CS_13*max_lambda_11 + CF_13)*(-CS_13*max_lambda_11 + CF_13)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_13*max_lambda_11 + CF_13) + ((25.0/6.0))*(-CS_12*max_lambda_11 +
      CF_12))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11)*(-(19.0/6.0)*(-CS_12*max_lambda_11 + CF_12) + ((2.0/3.0))*(-CS_11*max_lambda_11 + CF_11) +
      ((11.0/6.0))*(-CS_13*max_lambda_11 + CF_13));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_1 = (-(7.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_15*max_lambda_11 + CF_15) +
      ((11.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*omega_0 + (-(1.0/12.0)*(-CS_11*max_lambda_11 + CF_11) +
      ((1.0/6.0))*(-CS_13*max_lambda_11 + CF_13) + ((5.0/12.0))*(-CS_12*max_lambda_11 + CF_12))*omega_2 +
      (-(1.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*omega_1 + Recon_1;

    beta_0 = ((1.0/3.0))*((CS_24*max_lambda_22 + CF_24)*(CS_24*max_lambda_22 + CF_24)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_24*max_lambda_22 + CF_24) + ((25.0/6.0))*(CS_23*max_lambda_22 +
      CF_23))*(CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(CS_22*max_lambda_22 +
      CF_22)*(-(31.0/6.0)*(CS_23*max_lambda_22 + CF_23) + ((5.0/3.0))*(CS_22*max_lambda_22 + CF_22) +
      ((11.0/6.0))*(CS_24*max_lambda_22 + CF_24));

    beta_1 = ((1.0/3.0))*((CS_23*max_lambda_22 + CF_23)*(CS_23*max_lambda_22 + CF_23)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_23*max_lambda_22 + CF_23) + ((13.0/6.0))*(CS_22*max_lambda_22 +
      CF_22))*(CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(CS_21*max_lambda_22 +
      CF_21)*(-(13.0/6.0)*(CS_22*max_lambda_22 + CF_22) + ((2.0/3.0))*(CS_21*max_lambda_22 + CF_21) +
      ((5.0/6.0))*(CS_23*max_lambda_22 + CF_23));

    beta_2 = ((5.0/6.0))*((CS_22*max_lambda_22 + CF_22)*(CS_22*max_lambda_22 + CF_22)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_20*max_lambda_22 + CF_20) + ((11.0/6.0))*(CS_22*max_lambda_22 +
      CF_22))*(CS_20*max_lambda_22 + CF_20) + ((1.0/2.0))*(CS_21*max_lambda_22 +
      CF_21)*(-(31.0/6.0)*(CS_22*max_lambda_22 + CF_22) - (19.0/6.0)*(CS_20*max_lambda_22 + CF_20) +
      ((25.0/6.0))*(CS_21*max_lambda_22 + CF_21));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_2 = (-(7.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_20*max_lambda_22 + CF_20) +
      ((11.0/12.0))*(CS_22*max_lambda_22 + CF_22))*omega_2 + (-(1.0/12.0)*(CS_21*max_lambda_22 + CF_21) +
      ((1.0/6.0))*(CS_23*max_lambda_22 + CF_23) + ((5.0/12.0))*(CS_22*max_lambda_22 + CF_22))*omega_1 +
      (-(1.0/12.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(CS_23*max_lambda_22 + CF_23))*omega_0 + Recon_2;

    beta_0 = ((1.0/3.0))*((-CS_25*max_lambda_22 + CF_25)*(-CS_25*max_lambda_22 + CF_25)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_25*max_lambda_22 + CF_25) + ((25.0/6.0))*(-CS_24*max_lambda_22 +
      CF_24))*(-CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(-CS_23*max_lambda_22 +
      CF_23)*(-(31.0/6.0)*(-CS_24*max_lambda_22 + CF_24) + ((5.0/3.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((11.0/6.0))*(-CS_25*max_lambda_22 + CF_25));

    beta_1 = ((1.0/3.0))*((-CS_24*max_lambda_22 + CF_24)*(-CS_24*max_lambda_22 + CF_24)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_24*max_lambda_22 + CF_24) + ((13.0/6.0))*(-CS_23*max_lambda_22 +
      CF_23))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
      CF_22)*(-(13.0/6.0)*(-CS_23*max_lambda_22 + CF_23) + ((2.0/3.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/6.0))*(-CS_24*max_lambda_22 + CF_24));

    beta_2 = ((5.0/6.0))*((-CS_23*max_lambda_22 + CF_23)*(-CS_23*max_lambda_22 + CF_23)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_23*max_lambda_22 + CF_23) + ((25.0/6.0))*(-CS_22*max_lambda_22 +
      CF_22))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21)*(-(19.0/6.0)*(-CS_22*max_lambda_22 + CF_22) + ((2.0/3.0))*(-CS_21*max_lambda_22 + CF_21) +
      ((11.0/6.0))*(-CS_23*max_lambda_22 + CF_23));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_2 = (-(7.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_25*max_lambda_22 + CF_25) +
      ((11.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*omega_0 + (-(1.0/12.0)*(-CS_21*max_lambda_22 + CF_21) +
      ((1.0/6.0))*(-CS_23*max_lambda_22 + CF_23) + ((5.0/12.0))*(-CS_22*max_lambda_22 + CF_22))*omega_2 +
      (-(1.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*omega_1 + Recon_2;

    beta_0 = ((1.0/3.0))*((CS_34*max_lambda_33 + CF_34)*(CS_34*max_lambda_33 + CF_34)) +
      ((1.0/2.0))*(-(19.0/6.0)*(CS_34*max_lambda_33 + CF_34) + ((25.0/6.0))*(CS_33*max_lambda_33 +
      CF_33))*(CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(CS_32*max_lambda_33 +
      CF_32)*(-(31.0/6.0)*(CS_33*max_lambda_33 + CF_33) + ((5.0/3.0))*(CS_32*max_lambda_33 + CF_32) +
      ((11.0/6.0))*(CS_34*max_lambda_33 + CF_34));

    beta_1 = ((1.0/3.0))*((CS_33*max_lambda_33 + CF_33)*(CS_33*max_lambda_33 + CF_33)) +
      ((1.0/2.0))*(-(13.0/6.0)*(CS_33*max_lambda_33 + CF_33) + ((13.0/6.0))*(CS_32*max_lambda_33 +
      CF_32))*(CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(CS_31*max_lambda_33 +
      CF_31)*(-(13.0/6.0)*(CS_32*max_lambda_33 + CF_32) + ((2.0/3.0))*(CS_31*max_lambda_33 + CF_31) +
      ((5.0/6.0))*(CS_33*max_lambda_33 + CF_33));

    beta_2 = ((5.0/6.0))*((CS_32*max_lambda_33 + CF_32)*(CS_32*max_lambda_33 + CF_32)) +
      ((1.0/2.0))*(((2.0/3.0))*(CS_30*max_lambda_33 + CF_30) + ((11.0/6.0))*(CS_32*max_lambda_33 +
      CF_32))*(CS_30*max_lambda_33 + CF_30) + ((1.0/2.0))*(CS_31*max_lambda_33 +
      CF_31)*(-(31.0/6.0)*(CS_32*max_lambda_33 + CF_32) - (19.0/6.0)*(CS_30*max_lambda_33 + CF_30) +
      ((25.0/6.0))*(CS_31*max_lambda_33 + CF_31));

    alpha_0 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_3 = (-(7.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_30*max_lambda_33 + CF_30) +
      ((11.0/12.0))*(CS_32*max_lambda_33 + CF_32))*omega_2 + (-(1.0/12.0)*(CS_31*max_lambda_33 + CF_31) +
      ((1.0/6.0))*(CS_33*max_lambda_33 + CF_33) + ((5.0/12.0))*(CS_32*max_lambda_33 + CF_32))*omega_1 +
      (-(1.0/12.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(CS_33*max_lambda_33 + CF_33))*omega_0 + Recon_3;

    beta_0 = ((1.0/3.0))*((-CS_35*max_lambda_33 + CF_35)*(-CS_35*max_lambda_33 + CF_35)) +
      ((1.0/2.0))*(-(19.0/6.0)*(-CS_35*max_lambda_33 + CF_35) + ((25.0/6.0))*(-CS_34*max_lambda_33 +
      CF_34))*(-CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(-CS_33*max_lambda_33 +
      CF_33)*(-(31.0/6.0)*(-CS_34*max_lambda_33 + CF_34) + ((5.0/3.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((11.0/6.0))*(-CS_35*max_lambda_33 + CF_35));

    beta_1 = ((1.0/3.0))*((-CS_34*max_lambda_33 + CF_34)*(-CS_34*max_lambda_33 + CF_34)) +
      ((1.0/2.0))*(-(13.0/6.0)*(-CS_34*max_lambda_33 + CF_34) + ((13.0/6.0))*(-CS_33*max_lambda_33 +
      CF_33))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
      CF_32)*(-(13.0/6.0)*(-CS_33*max_lambda_33 + CF_33) + ((2.0/3.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/6.0))*(-CS_34*max_lambda_33 + CF_34));

    beta_2 = ((5.0/6.0))*((-CS_33*max_lambda_33 + CF_33)*(-CS_33*max_lambda_33 + CF_33)) +
      ((1.0/2.0))*(-(31.0/6.0)*(-CS_33*max_lambda_33 + CF_33) + ((25.0/6.0))*(-CS_32*max_lambda_33 +
      CF_32))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31)*(-(19.0/6.0)*(-CS_32*max_lambda_33 + CF_32) + ((2.0/3.0))*(-CS_31*max_lambda_33 + CF_31) +
      ((11.0/6.0))*(-CS_33*max_lambda_33 + CF_33));

    alpha_0 = 0.1 + ((1.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 +
      beta_0));

   alpha_1 = 0.6 + ((3.0/5.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

    alpha_2 = 0.3 + ((3.0/10.0))*(fabs(-beta_2 + beta_0)*fabs(-beta_2 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 +
      beta_2));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2));

   omega_0 = alpha_0*inv_alpha_sum;

   omega_1 = alpha_1*inv_alpha_sum;

   omega_2 = alpha_2*inv_alpha_sum;

    Recon_3 = (-(7.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_35*max_lambda_33 + CF_35) +
      ((11.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*omega_0 + (-(1.0/12.0)*(-CS_31*max_lambda_33 + CF_31) +
      ((1.0/6.0))*(-CS_33*max_lambda_33 + CF_33) + ((5.0/12.0))*(-CS_32*max_lambda_33 + CF_32))*omega_2 +
      (-(1.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*omega_1 + Recon_3;

    wk4_B0(0,0) = 0.707106781186547*AVG_1_rho*Recon_2*inv_AVG_a + 0.707106781186547*AVG_1_rho*Recon_3*inv_AVG_a +
      Recon_0;

    wk5_B0(0,0) = AVG_1_rho*Recon_1 + AVG_1_u0*Recon_0 + 0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_2*inv_AVG_a +
      0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_3*inv_AVG_a;

    wk6_B0(0,0) = AVG_1_u1*Recon_0 + 0.707106781186547*(-AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a +
      0.707106781186547*(AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_2*inv_AVG_a;

    wk7_B0(0,0) = (((1.0/2.0))*(AVG_1_u0*AVG_1_u0) + ((1.0/2.0))*(AVG_1_u1*AVG_1_u1))*Recon_0 +
      AVG_1_rho*AVG_1_u0*Recon_1 + 0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) +
      (AVG_1_u1*AVG_1_u1))*gamma_m1)*invgamma_m1 + AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_2*inv_AVG_a +
      0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) +
      (AVG_1_u1*AVG_1_u1))*gamma_m1)*invgamma_m1 - AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a;

}

 void opensbliblock00Kernel010(const ACC<double> &wk0_B0, const ACC<double> &wk1_B0, const ACC<double> &wk2_B0, const
ACC<double> &wk3_B0, const ACC<double> &wk4_B0, const ACC<double> &wk5_B0, const ACC<double> &wk6_B0, const ACC<double>
&wk7_B0, ACC<double> &Residual0_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0, ACC<double> &Residual3_B0)
{
   Residual0_B0(0,0) = -(-wk0_B0(-1,0) + wk0_B0(0,0))*invDelta0block0 - (-wk4_B0(0,-1) + wk4_B0(0,0))*invDelta1block0;

   Residual1_B0(0,0) = -(-wk1_B0(-1,0) + wk1_B0(0,0))*invDelta0block0 - (-wk5_B0(0,-1) + wk5_B0(0,0))*invDelta1block0;

   Residual2_B0(0,0) = -(-wk2_B0(-1,0) + wk2_B0(0,0))*invDelta0block0 - (-wk6_B0(0,-1) + wk6_B0(0,0))*invDelta1block0;

   Residual3_B0(0,0) = -(-wk3_B0(-1,0) + wk3_B0(0,0))*invDelta0block0 - (-wk7_B0(0,-1) + wk7_B0(0,0))*invDelta1block0;

}

 void opensbliblock00Kernel018(const ACC<double> &Residual0_B0, const ACC<double> &Residual1_B0, const ACC<double>
&Residual2_B0, const ACC<double> &Residual3_B0, const ACC<double> &rhoE_RKold_B0, const ACC<double> &rho_RKold_B0, const
ACC<double> &rhou0_RKold_B0, const ACC<double> &rhou1_RKold_B0, ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double>
&rhou0_B0, ACC<double> &rhou1_B0, const double *rknew)
{
   rho_B0(0,0) = rknew[0]*dt*Residual0_B0(0,0) + rho_RKold_B0(0,0);

   rhou0_B0(0,0) = rknew[0]*dt*Residual1_B0(0,0) + rhou0_RKold_B0(0,0);

   rhou1_B0(0,0) = rknew[0]*dt*Residual2_B0(0,0) + rhou1_RKold_B0(0,0);

   rhoE_B0(0,0) = rknew[0]*dt*Residual3_B0(0,0) + rhoE_RKold_B0(0,0);

}

 void opensbliblock00Kernel017(const ACC<double> &Residual0_B0, const ACC<double> &Residual1_B0, const ACC<double>
&Residual2_B0, const ACC<double> &Residual3_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double>
&rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, const double *rkold)
{
   rho_RKold_B0(0,0) = rkold[0]*dt*Residual0_B0(0,0) + rho_RKold_B0(0,0);

   rhou0_RKold_B0(0,0) = rkold[0]*dt*Residual1_B0(0,0) + rhou0_RKold_B0(0,0);

   rhou1_RKold_B0(0,0) = rkold[0]*dt*Residual2_B0(0,0) + rhou1_RKold_B0(0,0);

   rhoE_RKold_B0(0,0) = rkold[0]*dt*Residual3_B0(0,0) + rhoE_RKold_B0(0,0);

}

#endif
