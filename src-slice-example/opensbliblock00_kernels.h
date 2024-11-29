#ifndef OPENSBLIBLOCK00_KERNEL_H
#define OPENSBLIBLOCK00_KERNEL_H
 void opensbliblock00Kernel041(ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0,
ACC<double> &rhou2_B0, ACC<double> &x0_B0, ACC<double> &x1_B0, ACC<double> &x2_B0, const int *idx)
{
   double amp = 0.0;
   double b = 0.0;
   double cx0 = 0.0;
   double cx1 = 0.0;
   double cx2 = 0.0;
   double d = 0.0;
   double p = 0.0;
   double sx0 = 0.0;
   double sx1 = 0.0;
   double sx2 = 0.0;
   double u0 = 0.0;
   double u1 = 0.0;
   double u2 = 0.0;
   double ubar = 0.0;
   double vonkar = 0.0;
   x0_B0(0,0,0) = Delta0block0*idx[0];

   x1_B0(0,0,0) = -1.0*tanh((1.0 - 2.0*idx[1]/(-1.0 + block0np1))*stretch)/tanh(stretch);

   x2_B0(0,0,0) = Delta2block0*idx[2];

   b = 5.50000000000000;

   vonkar = 2.50000000000000;

    ubar = (((1 - fabs(x1_B0(0,0,0)))*Re < 10.0) ? (
   (1 - fabs(x1_B0(0,0,0)))*Re
)
: (
   vonkar*log((1 -
      fabs(x1_B0(0,0,0)))*Re) + b
));

   amp = 0.1*b + 0.1*vonkar*log(Re);

   sx0 = sin(4.0*M_PI*invlx0*x0_B0(0,0,0));

   sx1 = sin(M_PI*x1_B0(0,0,0));

   sx2 = sin(2.0*M_PI*invlx2*x2_B0(0,0,0));

   cx0 = cos(4.0*M_PI*invlx0*x0_B0(0,0,0));

   cx1 = 1 + cos(M_PI*x1_B0(0,0,0));

   cx2 = cos(2.0*M_PI*invlx2*x2_B0(0,0,0));

   u0 = 0.5*lx0*amp*cx0*sx1*sx2 + ubar;

   u1 = -amp*cx1*sx0*sx2;

   u2 = -0.5*lx2*amp*cx2*sx0*sx1;

   p = invgama*inv2Minf;

   d = 1;

   rho_B0(0,0,0) = d;

   rhou0_B0(0,0,0) = d*u0;

   rhou1_B0(0,0,0) = d*u1;

   rhou2_B0(0,0,0) = d*u2;

   rhoE_B0(0,0,0) = p/(-1 + gama) + 0.5*((u0*u0) + (u1*u1) + (u2*u2))*d;

}

 void opensbliblock00Kernel043(const ACC<double> &x1_B0, ACC<double> &D11_B0, ACC<double> &detJ_B0, ACC<double> &wk4_B0,
const int *idx)
{
   double d1_x1_dy = 0.0;
    d1_x1_dy = invDelta1block0*((idx[1] == 0) ? (
   3.0*x1_B0(0,1,0) + 0.333333333333333*x1_B0(0,3,0) -
      1.5*x1_B0(0,2,0) - 1.83333333333333*x1_B0(0,0,0)
)
: ((idx[1] == 1) ? (
   0.0394168524399447*x1_B0(0,2,0) +
      0.00571369039775442*x1_B0(0,4,0) + 0.719443173328855*x1_B0(0,1,0) - 0.322484932882161*x1_B0(0,0,0) -
      0.0658051057710389*x1_B0(0,3,0) - 0.376283677513354*x1_B0(0,-1,0)
)
: ((idx[1] == 2) ? (
  
      0.197184333887745*x1_B0(0,0,0) + 0.521455851089587*x1_B0(0,1,0) + 0.113446470384241*x1_B0(0,-2,0) -
      0.00412637789557492*x1_B0(0,3,0) - 0.0367146847001261*x1_B0(0,2,0) - 0.791245592765872*x1_B0(0,-1,0)
)
: ((idx[1]
      == 3) ? (
   0.0451033223343881*x1_B0(0,0,0) + 0.652141084861241*x1_B0(0,1,0) + 0.121937153224065*x1_B0(0,-2,0) -
      0.00932597985049999*x1_B0(0,-3,0) - 0.727822147724592*x1_B0(0,-1,0) - 0.082033432844602*x1_B0(0,2,0)
)
: ((idx[1]
      == -1 + block0np1) ? (
   1.5*x1_B0(0,-2,0) + 1.83333333333333*x1_B0(0,0,0) - 3.0*x1_B0(0,-1,0) -
      0.333333333333333*x1_B0(0,-3,0)
)
: ((idx[1] == -2 + block0np1) ? (
   0.322484932882161*x1_B0(0,0,0) +
      0.0658051057710389*x1_B0(0,-3,0) + 0.376283677513354*x1_B0(0,1,0) - 0.0394168524399447*x1_B0(0,-2,0) -
      0.00571369039775442*x1_B0(0,-4,0) - 0.719443173328855*x1_B0(0,-1,0)
)
: ((idx[1] == -3 + block0np1) ? (
  
      0.00412637789557492*x1_B0(0,-3,0) + 0.0367146847001261*x1_B0(0,-2,0) + 0.791245592765872*x1_B0(0,1,0) -
      0.197184333887745*x1_B0(0,0,0) - 0.521455851089587*x1_B0(0,-1,0) - 0.113446470384241*x1_B0(0,2,0)
)
: ((idx[1] ==
      -4 + block0np1) ? (
   0.00932597985049999*x1_B0(0,3,0) + 0.727822147724592*x1_B0(0,1,0) +
      0.082033432844602*x1_B0(0,-2,0) - 0.0451033223343881*x1_B0(0,0,0) - 0.652141084861241*x1_B0(0,-1,0) -
      0.121937153224065*x1_B0(0,2,0)
)
: (
   -(2.0/3.0)*x1_B0(0,-1,0) - (1.0/12.0)*x1_B0(0,2,0) +
      ((1.0/12.0))*x1_B0(0,-2,0) + ((2.0/3.0))*x1_B0(0,1,0)
)))))))));

   wk4_B0(0,0,0) = d1_x1_dy;

   detJ_B0(0,0,0) = d1_x1_dy;

   D11_B0(0,0,0) = 1.0/(d1_x1_dy);

}

void opensbliblock00Kernel046(ACC<double> &D11_B0, ACC<double> &detJ_B0)
{
   D11_B0(0,-1,0) = D11_B0(0,1,0);

   detJ_B0(0,-1,0) = detJ_B0(0,1,0);

   D11_B0(0,-2,0) = D11_B0(0,2,0);

   detJ_B0(0,-2,0) = detJ_B0(0,2,0);

   D11_B0(0,-3,0) = D11_B0(0,3,0);

   detJ_B0(0,-3,0) = detJ_B0(0,3,0);

   D11_B0(0,-4,0) = D11_B0(0,4,0);

   detJ_B0(0,-4,0) = detJ_B0(0,4,0);

   D11_B0(0,-5,0) = D11_B0(0,5,0);

   detJ_B0(0,-5,0) = detJ_B0(0,5,0);

}

void opensbliblock00Kernel047(ACC<double> &D11_B0, ACC<double> &detJ_B0)
{
   D11_B0(0,1,0) = D11_B0(0,-1,0);

   detJ_B0(0,1,0) = detJ_B0(0,-1,0);

   D11_B0(0,2,0) = D11_B0(0,-2,0);

   detJ_B0(0,2,0) = detJ_B0(0,-2,0);

   D11_B0(0,3,0) = D11_B0(0,-3,0);

   detJ_B0(0,3,0) = detJ_B0(0,-3,0);

   D11_B0(0,4,0) = D11_B0(0,-4,0);

   detJ_B0(0,4,0) = detJ_B0(0,-4,0);

   D11_B0(0,5,0) = D11_B0(0,-5,0);

   detJ_B0(0,5,0) = detJ_B0(0,-5,0);

}

void opensbliblock00Kernel051(const ACC<double> &D11_B0, ACC<double> &SD111_B0, const int *idx)
{
   double d1_D11_dy = 0.0;
    d1_D11_dy = invDelta1block0*((idx[1] == 0) ? (
   3.0*D11_B0(0,1,0) + 0.333333333333333*D11_B0(0,3,0) -
      1.5*D11_B0(0,2,0) - 1.83333333333333*D11_B0(0,0,0)
)
: ((idx[1] == 1) ? (
   0.0394168524399447*D11_B0(0,2,0) +
      0.00571369039775442*D11_B0(0,4,0) + 0.719443173328855*D11_B0(0,1,0) - 0.322484932882161*D11_B0(0,0,0) -
      0.0658051057710389*D11_B0(0,3,0) - 0.376283677513354*D11_B0(0,-1,0)
)
: ((idx[1] == 2) ? (
  
      0.197184333887745*D11_B0(0,0,0) + 0.521455851089587*D11_B0(0,1,0) + 0.113446470384241*D11_B0(0,-2,0) -
      0.00412637789557492*D11_B0(0,3,0) - 0.0367146847001261*D11_B0(0,2,0) - 0.791245592765872*D11_B0(0,-1,0)
)
:
      ((idx[1] == 3) ? (
   0.0451033223343881*D11_B0(0,0,0) + 0.652141084861241*D11_B0(0,1,0) +
      0.121937153224065*D11_B0(0,-2,0) - 0.00932597985049999*D11_B0(0,-3,0) - 0.727822147724592*D11_B0(0,-1,0) -
      0.082033432844602*D11_B0(0,2,0)
)
: ((idx[1] == -1 + block0np1) ? (
   1.5*D11_B0(0,-2,0) +
      1.83333333333333*D11_B0(0,0,0) - 3.0*D11_B0(0,-1,0) - 0.333333333333333*D11_B0(0,-3,0)
)
: ((idx[1] == -2 +
      block0np1) ? (
   0.322484932882161*D11_B0(0,0,0) + 0.0658051057710389*D11_B0(0,-3,0) +
      0.376283677513354*D11_B0(0,1,0) - 0.0394168524399447*D11_B0(0,-2,0) - 0.00571369039775442*D11_B0(0,-4,0) -
      0.719443173328855*D11_B0(0,-1,0)
)
: ((idx[1] == -3 + block0np1) ? (
   0.00412637789557492*D11_B0(0,-3,0) +
      0.0367146847001261*D11_B0(0,-2,0) + 0.791245592765872*D11_B0(0,1,0) - 0.197184333887745*D11_B0(0,0,0) -
      0.521455851089587*D11_B0(0,-1,0) - 0.113446470384241*D11_B0(0,2,0)
)
: ((idx[1] == -4 + block0np1) ? (
  
      0.00932597985049999*D11_B0(0,3,0) + 0.727822147724592*D11_B0(0,1,0) + 0.082033432844602*D11_B0(0,-2,0) -
      0.0451033223343881*D11_B0(0,0,0) - 0.652141084861241*D11_B0(0,-1,0) - 0.121937153224065*D11_B0(0,2,0)
)
: (
  
      -(2.0/3.0)*D11_B0(0,-1,0) - (1.0/12.0)*D11_B0(0,2,0) + ((1.0/12.0))*D11_B0(0,-2,0) +
      ((2.0/3.0))*D11_B0(0,1,0)
)))))))));

   SD111_B0(0,0,0) = d1_D11_dy;

}

 void opensbliblock00Kernel037(ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0,
ACC<double> &rhou2_B0)
{
   double Pwall = 0.0;
   double T1 = 0.0;
   double T2 = 0.0;
   double T3 = 0.0;
   double T4 = 0.0;
   double T5 = 0.0;
   double T_above = 0.0;
   double rho_halo_1 = 0.0;
   double rho_halo_2 = 0.0;
   double rho_halo_3 = 0.0;
   double rho_halo_4 = 0.0;
   double rho_halo_5 = 0.0;
   double u01 = 0.0;
   double u02 = 0.0;
   double u03 = 0.0;
   double u04 = 0.0;
   double u05 = 0.0;
   double u11 = 0.0;
   double u12 = 0.0;
   double u13 = 0.0;
   double u14 = 0.0;
   double u15 = 0.0;
   double u21 = 0.0;
   double u22 = 0.0;
   double u23 = 0.0;
   double u24 = 0.0;
   double u25 = 0.0;
   rhou0_B0(0,0,0) = 0.0;

   rhou1_B0(0,0,0) = 0.0;

   rhou2_B0(0,0,0) = 0.0;

   rhoE_B0(0,0,0) = Twall*invgama*inv2Minf*rho_B0(0,0,0)/(-1 + gama);

    Pwall = (-1 + gama)*(-(((1.0/2.0))*(rhou0_B0(0,0,0)*rhou0_B0(0,0,0)) + ((1.0/2.0))*(rhou1_B0(0,0,0)*rhou1_B0(0,0,0))
      + ((1.0/2.0))*(rhou2_B0(0,0,0)*rhou2_B0(0,0,0)))/rho_B0(0,0,0) + rhoE_B0(0,0,0));

   u01 = rhou0_B0(0,1,0)/rho_B0(0,1,0);

   u02 = rhou0_B0(0,2,0)/rho_B0(0,2,0);

   u03 = rhou0_B0(0,3,0)/rho_B0(0,3,0);

   u04 = rhou0_B0(0,4,0)/rho_B0(0,4,0);

   u05 = rhou0_B0(0,5,0)/rho_B0(0,5,0);

   u11 = rhou1_B0(0,1,0)/rho_B0(0,1,0);

   u12 = rhou1_B0(0,2,0)/rho_B0(0,2,0);

   u13 = rhou1_B0(0,3,0)/rho_B0(0,3,0);

   u14 = rhou1_B0(0,4,0)/rho_B0(0,4,0);

   u15 = rhou1_B0(0,5,0)/rho_B0(0,5,0);

   u21 = rhou2_B0(0,1,0)/rho_B0(0,1,0);

   u22 = rhou2_B0(0,2,0)/rho_B0(0,2,0);

   u23 = rhou2_B0(0,3,0)/rho_B0(0,3,0);

   u24 = rhou2_B0(0,4,0)/rho_B0(0,4,0);

   u25 = rhou2_B0(0,5,0)/rho_B0(0,5,0);

    T_above = (Minf*Minf)*(-1 + gama)*(-(((1.0/2.0))*(rhou0_B0(0,1,0)*rhou0_B0(0,1,0)) +
      ((1.0/2.0))*(rhou1_B0(0,1,0)*rhou1_B0(0,1,0)) + ((1.0/2.0))*(rhou2_B0(0,1,0)*rhou2_B0(0,1,0)))/rho_B0(0,1,0) +
      rhoE_B0(0,1,0))*gama/rho_B0(0,1,0);

   T1 = -T_above + 2*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,-1,0) = rho_halo_1;

   rhou0_B0(0,-1,0) = -rho_halo_1*u01;

   rhou1_B0(0,-1,0) = -rho_halo_1*u11;

   rhou2_B0(0,-1,0) = -rho_halo_1*u21;

   rhoE_B0(0,-1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   T2 = -2*T_above + 3*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,-1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,-2,0) = rho_halo_2;

   rhou0_B0(0,-1,0) = -rho_halo_1*u01;

   rhou1_B0(0,-1,0) = -rho_halo_1*u11;

   rhou2_B0(0,-1,0) = -rho_halo_1*u21;

   rhou0_B0(0,-2,0) = -rho_halo_2*u02;

   rhou1_B0(0,-2,0) = -rho_halo_2*u12;

   rhou2_B0(0,-2,0) = -rho_halo_2*u22;

   rhoE_B0(0,-1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,-2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   T3 = -3*T_above + 4*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,-1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,-2,0) = rho_halo_2;

   rho_halo_3 = (Minf*Minf)*gama*Pwall/T3;

   rho_B0(0,-3,0) = rho_halo_3;

   rhou0_B0(0,-1,0) = -rho_halo_1*u01;

   rhou1_B0(0,-1,0) = -rho_halo_1*u11;

   rhou2_B0(0,-1,0) = -rho_halo_1*u21;

   rhou0_B0(0,-2,0) = -rho_halo_2*u02;

   rhou1_B0(0,-2,0) = -rho_halo_2*u12;

   rhou2_B0(0,-2,0) = -rho_halo_2*u22;

   rhou0_B0(0,-3,0) = -rho_halo_3*u03;

   rhou1_B0(0,-3,0) = -rho_halo_3*u13;

   rhou2_B0(0,-3,0) = -rho_halo_3*u23;

   rhoE_B0(0,-1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,-2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   rhoE_B0(0,-3,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u03*u03) + (u13*u13) + (u23*u23))*rho_halo_3;

   T4 = -4*T_above + 5*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,-1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,-2,0) = rho_halo_2;

   rho_halo_3 = (Minf*Minf)*gama*Pwall/T3;

   rho_B0(0,-3,0) = rho_halo_3;

   rho_halo_4 = (Minf*Minf)*gama*Pwall/T4;

   rho_B0(0,-4,0) = rho_halo_4;

   rhou0_B0(0,-1,0) = -rho_halo_1*u01;

   rhou1_B0(0,-1,0) = -rho_halo_1*u11;

   rhou2_B0(0,-1,0) = -rho_halo_1*u21;

   rhou0_B0(0,-2,0) = -rho_halo_2*u02;

   rhou1_B0(0,-2,0) = -rho_halo_2*u12;

   rhou2_B0(0,-2,0) = -rho_halo_2*u22;

   rhou0_B0(0,-3,0) = -rho_halo_3*u03;

   rhou1_B0(0,-3,0) = -rho_halo_3*u13;

   rhou2_B0(0,-3,0) = -rho_halo_3*u23;

   rhou0_B0(0,-4,0) = -rho_halo_4*u04;

   rhou1_B0(0,-4,0) = -rho_halo_4*u14;

   rhou2_B0(0,-4,0) = -rho_halo_4*u24;

   rhoE_B0(0,-1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,-2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   rhoE_B0(0,-3,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u03*u03) + (u13*u13) + (u23*u23))*rho_halo_3;

   rhoE_B0(0,-4,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u04*u04) + (u14*u14) + (u24*u24))*rho_halo_4;

   T5 = -5*T_above + 6*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,-1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,-2,0) = rho_halo_2;

   rho_halo_3 = (Minf*Minf)*gama*Pwall/T3;

   rho_B0(0,-3,0) = rho_halo_3;

   rho_halo_4 = (Minf*Minf)*gama*Pwall/T4;

   rho_B0(0,-4,0) = rho_halo_4;

   rho_halo_5 = (Minf*Minf)*gama*Pwall/T5;

   rho_B0(0,-5,0) = rho_halo_5;

   rhou0_B0(0,-1,0) = -rho_halo_1*u01;

   rhou1_B0(0,-1,0) = -rho_halo_1*u11;

   rhou2_B0(0,-1,0) = -rho_halo_1*u21;

   rhou0_B0(0,-2,0) = -rho_halo_2*u02;

   rhou1_B0(0,-2,0) = -rho_halo_2*u12;

   rhou2_B0(0,-2,0) = -rho_halo_2*u22;

   rhou0_B0(0,-3,0) = -rho_halo_3*u03;

   rhou1_B0(0,-3,0) = -rho_halo_3*u13;

   rhou2_B0(0,-3,0) = -rho_halo_3*u23;

   rhou0_B0(0,-4,0) = -rho_halo_4*u04;

   rhou1_B0(0,-4,0) = -rho_halo_4*u14;

   rhou2_B0(0,-4,0) = -rho_halo_4*u24;

   rhou0_B0(0,-5,0) = -rho_halo_5*u05;

   rhou1_B0(0,-5,0) = -rho_halo_5*u15;

   rhou2_B0(0,-5,0) = -rho_halo_5*u25;

   rhoE_B0(0,-1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,-2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   rhoE_B0(0,-3,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u03*u03) + (u13*u13) + (u23*u23))*rho_halo_3;

   rhoE_B0(0,-4,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u04*u04) + (u14*u14) + (u24*u24))*rho_halo_4;

   rhoE_B0(0,-5,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u05*u05) + (u15*u15) + (u25*u25))*rho_halo_5;

}

 void opensbliblock00Kernel038(ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0,
ACC<double> &rhou2_B0)
{
   double Pwall = 0.0;
   double T1 = 0.0;
   double T2 = 0.0;
   double T3 = 0.0;
   double T4 = 0.0;
   double T5 = 0.0;
   double T_above = 0.0;
   double rho_halo_1 = 0.0;
   double rho_halo_2 = 0.0;
   double rho_halo_3 = 0.0;
   double rho_halo_4 = 0.0;
   double rho_halo_5 = 0.0;
   double u01 = 0.0;
   double u02 = 0.0;
   double u03 = 0.0;
   double u04 = 0.0;
   double u05 = 0.0;
   double u11 = 0.0;
   double u12 = 0.0;
   double u13 = 0.0;
   double u14 = 0.0;
   double u15 = 0.0;
   double u21 = 0.0;
   double u22 = 0.0;
   double u23 = 0.0;
   double u24 = 0.0;
   double u25 = 0.0;
   rhou0_B0(0,0,0) = 0.0;

   rhou1_B0(0,0,0) = 0.0;

   rhou2_B0(0,0,0) = 0.0;

   rhoE_B0(0,0,0) = Twall*invgama*inv2Minf*rho_B0(0,0,0)/(-1 + gama);

    Pwall = (-1 + gama)*(-(((1.0/2.0))*(rhou0_B0(0,0,0)*rhou0_B0(0,0,0)) + ((1.0/2.0))*(rhou1_B0(0,0,0)*rhou1_B0(0,0,0))
      + ((1.0/2.0))*(rhou2_B0(0,0,0)*rhou2_B0(0,0,0)))/rho_B0(0,0,0) + rhoE_B0(0,0,0));

   u01 = rhou0_B0(0,-1,0)/rho_B0(0,-1,0);

   u02 = rhou0_B0(0,-2,0)/rho_B0(0,-2,0);

   u03 = rhou0_B0(0,-3,0)/rho_B0(0,-3,0);

   u04 = rhou0_B0(0,-4,0)/rho_B0(0,-4,0);

   u05 = rhou0_B0(0,-5,0)/rho_B0(0,-5,0);

   u11 = rhou1_B0(0,-1,0)/rho_B0(0,-1,0);

   u12 = rhou1_B0(0,-2,0)/rho_B0(0,-2,0);

   u13 = rhou1_B0(0,-3,0)/rho_B0(0,-3,0);

   u14 = rhou1_B0(0,-4,0)/rho_B0(0,-4,0);

   u15 = rhou1_B0(0,-5,0)/rho_B0(0,-5,0);

   u21 = rhou2_B0(0,-1,0)/rho_B0(0,-1,0);

   u22 = rhou2_B0(0,-2,0)/rho_B0(0,-2,0);

   u23 = rhou2_B0(0,-3,0)/rho_B0(0,-3,0);

   u24 = rhou2_B0(0,-4,0)/rho_B0(0,-4,0);

   u25 = rhou2_B0(0,-5,0)/rho_B0(0,-5,0);

    T_above = (Minf*Minf)*(-1 + gama)*(-(((1.0/2.0))*(rhou0_B0(0,-1,0)*rhou0_B0(0,-1,0)) +
      ((1.0/2.0))*(rhou1_B0(0,-1,0)*rhou1_B0(0,-1,0)) + ((1.0/2.0))*(rhou2_B0(0,-1,0)*rhou2_B0(0,-1,0)))/rho_B0(0,-1,0)
      + rhoE_B0(0,-1,0))*gama/rho_B0(0,-1,0);

   T1 = -T_above + 2*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,1,0) = rho_halo_1;

   rhou0_B0(0,1,0) = -rho_halo_1*u01;

   rhou1_B0(0,1,0) = -rho_halo_1*u11;

   rhou2_B0(0,1,0) = -rho_halo_1*u21;

   rhoE_B0(0,1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   T2 = -2*T_above + 3*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,2,0) = rho_halo_2;

   rhou0_B0(0,1,0) = -rho_halo_1*u01;

   rhou1_B0(0,1,0) = -rho_halo_1*u11;

   rhou2_B0(0,1,0) = -rho_halo_1*u21;

   rhou0_B0(0,2,0) = -rho_halo_2*u02;

   rhou1_B0(0,2,0) = -rho_halo_2*u12;

   rhou2_B0(0,2,0) = -rho_halo_2*u22;

   rhoE_B0(0,1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   T3 = -3*T_above + 4*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,2,0) = rho_halo_2;

   rho_halo_3 = (Minf*Minf)*gama*Pwall/T3;

   rho_B0(0,3,0) = rho_halo_3;

   rhou0_B0(0,1,0) = -rho_halo_1*u01;

   rhou1_B0(0,1,0) = -rho_halo_1*u11;

   rhou2_B0(0,1,0) = -rho_halo_1*u21;

   rhou0_B0(0,2,0) = -rho_halo_2*u02;

   rhou1_B0(0,2,0) = -rho_halo_2*u12;

   rhou2_B0(0,2,0) = -rho_halo_2*u22;

   rhou0_B0(0,3,0) = -rho_halo_3*u03;

   rhou1_B0(0,3,0) = -rho_halo_3*u13;

   rhou2_B0(0,3,0) = -rho_halo_3*u23;

   rhoE_B0(0,1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   rhoE_B0(0,3,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u03*u03) + (u13*u13) + (u23*u23))*rho_halo_3;

   T4 = -4*T_above + 5*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,2,0) = rho_halo_2;

   rho_halo_3 = (Minf*Minf)*gama*Pwall/T3;

   rho_B0(0,3,0) = rho_halo_3;

   rho_halo_4 = (Minf*Minf)*gama*Pwall/T4;

   rho_B0(0,4,0) = rho_halo_4;

   rhou0_B0(0,1,0) = -rho_halo_1*u01;

   rhou1_B0(0,1,0) = -rho_halo_1*u11;

   rhou2_B0(0,1,0) = -rho_halo_1*u21;

   rhou0_B0(0,2,0) = -rho_halo_2*u02;

   rhou1_B0(0,2,0) = -rho_halo_2*u12;

   rhou2_B0(0,2,0) = -rho_halo_2*u22;

   rhou0_B0(0,3,0) = -rho_halo_3*u03;

   rhou1_B0(0,3,0) = -rho_halo_3*u13;

   rhou2_B0(0,3,0) = -rho_halo_3*u23;

   rhou0_B0(0,4,0) = -rho_halo_4*u04;

   rhou1_B0(0,4,0) = -rho_halo_4*u14;

   rhou2_B0(0,4,0) = -rho_halo_4*u24;

   rhoE_B0(0,1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   rhoE_B0(0,3,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u03*u03) + (u13*u13) + (u23*u23))*rho_halo_3;

   rhoE_B0(0,4,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u04*u04) + (u14*u14) + (u24*u24))*rho_halo_4;

   T5 = -5*T_above + 6*Twall;

   rho_halo_1 = (Minf*Minf)*gama*Pwall/T1;

   rho_B0(0,1,0) = rho_halo_1;

   rho_halo_2 = (Minf*Minf)*gama*Pwall/T2;

   rho_B0(0,2,0) = rho_halo_2;

   rho_halo_3 = (Minf*Minf)*gama*Pwall/T3;

   rho_B0(0,3,0) = rho_halo_3;

   rho_halo_4 = (Minf*Minf)*gama*Pwall/T4;

   rho_B0(0,4,0) = rho_halo_4;

   rho_halo_5 = (Minf*Minf)*gama*Pwall/T5;

   rho_B0(0,5,0) = rho_halo_5;

   rhou0_B0(0,1,0) = -rho_halo_1*u01;

   rhou1_B0(0,1,0) = -rho_halo_1*u11;

   rhou2_B0(0,1,0) = -rho_halo_1*u21;

   rhou0_B0(0,2,0) = -rho_halo_2*u02;

   rhou1_B0(0,2,0) = -rho_halo_2*u12;

   rhou2_B0(0,2,0) = -rho_halo_2*u22;

   rhou0_B0(0,3,0) = -rho_halo_3*u03;

   rhou1_B0(0,3,0) = -rho_halo_3*u13;

   rhou2_B0(0,3,0) = -rho_halo_3*u23;

   rhou0_B0(0,4,0) = -rho_halo_4*u04;

   rhou1_B0(0,4,0) = -rho_halo_4*u14;

   rhou2_B0(0,4,0) = -rho_halo_4*u24;

   rhou0_B0(0,5,0) = -rho_halo_5*u05;

   rhou1_B0(0,5,0) = -rho_halo_5*u15;

   rhou2_B0(0,5,0) = -rho_halo_5*u25;

   rhoE_B0(0,1,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u01*u01) + (u11*u11) + (u21*u21))*rho_halo_1;

   rhoE_B0(0,2,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u02*u02) + (u12*u12) + (u22*u22))*rho_halo_2;

   rhoE_B0(0,3,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u03*u03) + (u13*u13) + (u23*u23))*rho_halo_3;

   rhoE_B0(0,4,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u04*u04) + (u14*u14) + (u24*u24))*rho_halo_4;

   rhoE_B0(0,5,0) = inv_gamma_m1*Pwall + ((1.0/2.0))*((u05*u05) + (u15*u15) + (u25*u25))*rho_halo_5;

}

void opensbliblock00Kernel005(const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, ACC<double> &u0_B0)
{
   u0_B0(0,0,0) = rhou0_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel007(const ACC<double> &rho_B0, const ACC<double> &rhou1_B0, ACC<double> &u1_B0)
{
   u1_B0(0,0,0) = rhou1_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel009(const ACC<double> &rho_B0, const ACC<double> &rhou2_B0, ACC<double> &u2_B0)
{
   u2_B0(0,0,0) = rhou2_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel034(const ACC<double> &x1_B0, ACC<double> &phi_B0)
{
   phi_B0(0,0,0) = tanh(aCF*x1_B0(0,0,0));

}

 void opensbliblock00Kernel023(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double> &p_B0)
{
    p_B0(0,0,0) = (-1 + gama)*(-(1.0/2.0)*(u0_B0(0,0,0)*u0_B0(0,0,0))*rho_B0(0,0,0) -
      (1.0/2.0)*(u1_B0(0,0,0)*u1_B0(0,0,0))*rho_B0(0,0,0) - (1.0/2.0)*(u2_B0(0,0,0)*u2_B0(0,0,0))*rho_B0(0,0,0) +
      rhoE_B0(0,0,0));

}

void opensbliblock00Kernel011(const ACC<double> &p_B0, const ACC<double> &rho_B0, ACC<double> &T_B0)
{
   T_B0(0,0,0) = (Minf*Minf)*gama*p_B0(0,0,0)/rho_B0(0,0,0);

}

 void opensbliblock00Kernel020(const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const ACC<double> &rho_B0,
ACC<double> &H_B0)
{
   H_B0(0,0,0) = (p_B0(0,0,0) + rhoE_B0(0,0,0))/rho_B0(0,0,0);

}

void opensbliblock00Kernel026(const ACC<double> &T_B0, ACC<double> &mu_B0)
{
   mu_B0(0,0,0) = pow(T_B0(0,0,0), 0.7);

}

void opensbliblock00Kernel004(const ACC<double> &u0_B0, ACC<double> &wk0_B0)
{
    wk0_B0(0,0,0) = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel006(const ACC<double> &u1_B0, ACC<double> &wk1_B0)
{
    wk1_B0(0,0,0) = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel008(const ACC<double> &u2_B0, ACC<double> &wk2_B0)
{
    wk2_B0(0,0,0) = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel010(const ACC<double> &T_B0, ACC<double> &wk3_B0)
{
    wk3_B0(0,0,0) = (-(2.0/3.0)*T_B0(-1,0,0) - (1.0/12.0)*T_B0(2,0,0) + ((1.0/12.0))*T_B0(-2,0,0) +
      ((2.0/3.0))*T_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel012(const ACC<double> &u0_B0, ACC<double> &wk4_B0, const int *idx)
{
   if (idx[1] == 0){

       wk4_B0(0,0,0) = (3.0*u0_B0(0,1,0) + 0.333333333333333*u0_B0(0,3,0) - 1.5*u0_B0(0,2,0) -
            1.83333333333333*u0_B0(0,0,0))*invDelta1block0;

   }

   else if (idx[1] == 1){

       wk4_B0(0,0,0) = (0.0394168524399447*u0_B0(0,2,0) + 0.00571369039775442*u0_B0(0,4,0) +
            0.719443173328855*u0_B0(0,1,0) - 0.322484932882161*u0_B0(0,0,0) - 0.0658051057710389*u0_B0(0,3,0) -
            0.376283677513354*u0_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 2){

       wk4_B0(0,0,0) = (0.197184333887745*u0_B0(0,0,0) + 0.521455851089587*u0_B0(0,1,0) +
            0.113446470384241*u0_B0(0,-2,0) - 0.00412637789557492*u0_B0(0,3,0) - 0.0367146847001261*u0_B0(0,2,0) -
            0.791245592765872*u0_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 3){

       wk4_B0(0,0,0) = (0.0451033223343881*u0_B0(0,0,0) + 0.652141084861241*u0_B0(0,1,0) +
            0.121937153224065*u0_B0(0,-2,0) - 0.00932597985049999*u0_B0(0,-3,0) - 0.727822147724592*u0_B0(0,-1,0) -
            0.082033432844602*u0_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       wk4_B0(0,0,0) = (1.5*u0_B0(0,-2,0) + 1.83333333333333*u0_B0(0,0,0) - 3.0*u0_B0(0,-1,0) -
            0.333333333333333*u0_B0(0,-3,0))*invDelta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       wk4_B0(0,0,0) = (0.322484932882161*u0_B0(0,0,0) + 0.0658051057710389*u0_B0(0,-3,0) +
            0.376283677513354*u0_B0(0,1,0) - 0.0394168524399447*u0_B0(0,-2,0) - 0.00571369039775442*u0_B0(0,-4,0) -
            0.719443173328855*u0_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == -3 + block0np1){

       wk4_B0(0,0,0) = (0.00412637789557492*u0_B0(0,-3,0) + 0.0367146847001261*u0_B0(0,-2,0) +
            0.791245592765872*u0_B0(0,1,0) - 0.197184333887745*u0_B0(0,0,0) - 0.521455851089587*u0_B0(0,-1,0) -
            0.113446470384241*u0_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -4 + block0np1){

       wk4_B0(0,0,0) = (0.00932597985049999*u0_B0(0,3,0) + 0.727822147724592*u0_B0(0,1,0) +
            0.082033432844602*u0_B0(0,-2,0) - 0.0451033223343881*u0_B0(0,0,0) - 0.652141084861241*u0_B0(0,-1,0) -
            0.121937153224065*u0_B0(0,2,0))*invDelta1block0;

   }

   else{

       wk4_B0(0,0,0) = (-(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) + ((1.0/12.0))*u0_B0(0,-2,0) +
            ((2.0/3.0))*u0_B0(0,1,0))*invDelta1block0;

   }

}

void opensbliblock00Kernel013(const ACC<double> &u1_B0, ACC<double> &wk5_B0, const int *idx)
{
   if (idx[1] == 0){

       wk5_B0(0,0,0) = (3.0*u1_B0(0,1,0) + 0.333333333333333*u1_B0(0,3,0) - 1.5*u1_B0(0,2,0) -
            1.83333333333333*u1_B0(0,0,0))*invDelta1block0;

   }

   else if (idx[1] == 1){

       wk5_B0(0,0,0) = (0.0394168524399447*u1_B0(0,2,0) + 0.00571369039775442*u1_B0(0,4,0) +
            0.719443173328855*u1_B0(0,1,0) - 0.322484932882161*u1_B0(0,0,0) - 0.0658051057710389*u1_B0(0,3,0) -
            0.376283677513354*u1_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 2){

       wk5_B0(0,0,0) = (0.197184333887745*u1_B0(0,0,0) + 0.521455851089587*u1_B0(0,1,0) +
            0.113446470384241*u1_B0(0,-2,0) - 0.00412637789557492*u1_B0(0,3,0) - 0.0367146847001261*u1_B0(0,2,0) -
            0.791245592765872*u1_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 3){

       wk5_B0(0,0,0) = (0.0451033223343881*u1_B0(0,0,0) + 0.652141084861241*u1_B0(0,1,0) +
            0.121937153224065*u1_B0(0,-2,0) - 0.00932597985049999*u1_B0(0,-3,0) - 0.727822147724592*u1_B0(0,-1,0) -
            0.082033432844602*u1_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       wk5_B0(0,0,0) = (1.5*u1_B0(0,-2,0) + 1.83333333333333*u1_B0(0,0,0) - 3.0*u1_B0(0,-1,0) -
            0.333333333333333*u1_B0(0,-3,0))*invDelta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       wk5_B0(0,0,0) = (0.322484932882161*u1_B0(0,0,0) + 0.0658051057710389*u1_B0(0,-3,0) +
            0.376283677513354*u1_B0(0,1,0) - 0.0394168524399447*u1_B0(0,-2,0) - 0.00571369039775442*u1_B0(0,-4,0) -
            0.719443173328855*u1_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == -3 + block0np1){

       wk5_B0(0,0,0) = (0.00412637789557492*u1_B0(0,-3,0) + 0.0367146847001261*u1_B0(0,-2,0) +
            0.791245592765872*u1_B0(0,1,0) - 0.197184333887745*u1_B0(0,0,0) - 0.521455851089587*u1_B0(0,-1,0) -
            0.113446470384241*u1_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -4 + block0np1){

       wk5_B0(0,0,0) = (0.00932597985049999*u1_B0(0,3,0) + 0.727822147724592*u1_B0(0,1,0) +
            0.082033432844602*u1_B0(0,-2,0) - 0.0451033223343881*u1_B0(0,0,0) - 0.652141084861241*u1_B0(0,-1,0) -
            0.121937153224065*u1_B0(0,2,0))*invDelta1block0;

   }

   else{

       wk5_B0(0,0,0) = (-(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) + ((1.0/12.0))*u1_B0(0,-2,0) +
            ((2.0/3.0))*u1_B0(0,1,0))*invDelta1block0;

   }

}

void opensbliblock00Kernel014(const ACC<double> &u2_B0, ACC<double> &wk6_B0, const int *idx)
{
   if (idx[1] == 0){

       wk6_B0(0,0,0) = (3.0*u2_B0(0,1,0) + 0.333333333333333*u2_B0(0,3,0) - 1.5*u2_B0(0,2,0) -
            1.83333333333333*u2_B0(0,0,0))*invDelta1block0;

   }

   else if (idx[1] == 1){

       wk6_B0(0,0,0) = (0.0394168524399447*u2_B0(0,2,0) + 0.00571369039775442*u2_B0(0,4,0) +
            0.719443173328855*u2_B0(0,1,0) - 0.322484932882161*u2_B0(0,0,0) - 0.0658051057710389*u2_B0(0,3,0) -
            0.376283677513354*u2_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 2){

       wk6_B0(0,0,0) = (0.197184333887745*u2_B0(0,0,0) + 0.521455851089587*u2_B0(0,1,0) +
            0.113446470384241*u2_B0(0,-2,0) - 0.00412637789557492*u2_B0(0,3,0) - 0.0367146847001261*u2_B0(0,2,0) -
            0.791245592765872*u2_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 3){

       wk6_B0(0,0,0) = (0.0451033223343881*u2_B0(0,0,0) + 0.652141084861241*u2_B0(0,1,0) +
            0.121937153224065*u2_B0(0,-2,0) - 0.00932597985049999*u2_B0(0,-3,0) - 0.727822147724592*u2_B0(0,-1,0) -
            0.082033432844602*u2_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       wk6_B0(0,0,0) = (1.5*u2_B0(0,-2,0) + 1.83333333333333*u2_B0(0,0,0) - 3.0*u2_B0(0,-1,0) -
            0.333333333333333*u2_B0(0,-3,0))*invDelta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       wk6_B0(0,0,0) = (0.322484932882161*u2_B0(0,0,0) + 0.0658051057710389*u2_B0(0,-3,0) +
            0.376283677513354*u2_B0(0,1,0) - 0.0394168524399447*u2_B0(0,-2,0) - 0.00571369039775442*u2_B0(0,-4,0) -
            0.719443173328855*u2_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == -3 + block0np1){

       wk6_B0(0,0,0) = (0.00412637789557492*u2_B0(0,-3,0) + 0.0367146847001261*u2_B0(0,-2,0) +
            0.791245592765872*u2_B0(0,1,0) - 0.197184333887745*u2_B0(0,0,0) - 0.521455851089587*u2_B0(0,-1,0) -
            0.113446470384241*u2_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -4 + block0np1){

       wk6_B0(0,0,0) = (0.00932597985049999*u2_B0(0,3,0) + 0.727822147724592*u2_B0(0,1,0) +
            0.082033432844602*u2_B0(0,-2,0) - 0.0451033223343881*u2_B0(0,0,0) - 0.652141084861241*u2_B0(0,-1,0) -
            0.121937153224065*u2_B0(0,2,0))*invDelta1block0;

   }

   else{

       wk6_B0(0,0,0) = (-(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) + ((1.0/12.0))*u2_B0(0,-2,0) +
            ((2.0/3.0))*u2_B0(0,1,0))*invDelta1block0;

   }

}

void opensbliblock00Kernel015(const ACC<double> &T_B0, ACC<double> &wk7_B0, const int *idx)
{
   if (idx[1] == 0){

       wk7_B0(0,0,0) = (3.0*T_B0(0,1,0) + 0.333333333333333*T_B0(0,3,0) - 1.5*T_B0(0,2,0) -
            1.83333333333333*T_B0(0,0,0))*invDelta1block0;

   }

   else if (idx[1] == 1){

       wk7_B0(0,0,0) = (0.0394168524399447*T_B0(0,2,0) + 0.00571369039775442*T_B0(0,4,0) + 0.719443173328855*T_B0(0,1,0)
            - 0.322484932882161*T_B0(0,0,0) - 0.0658051057710389*T_B0(0,3,0) -
            0.376283677513354*T_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 2){

       wk7_B0(0,0,0) = (0.197184333887745*T_B0(0,0,0) + 0.521455851089587*T_B0(0,1,0) + 0.113446470384241*T_B0(0,-2,0) -
            0.00412637789557492*T_B0(0,3,0) - 0.0367146847001261*T_B0(0,2,0) -
            0.791245592765872*T_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 3){

       wk7_B0(0,0,0) = (0.0451033223343881*T_B0(0,0,0) + 0.652141084861241*T_B0(0,1,0) + 0.121937153224065*T_B0(0,-2,0)
            - 0.00932597985049999*T_B0(0,-3,0) - 0.727822147724592*T_B0(0,-1,0) -
            0.082033432844602*T_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       wk7_B0(0,0,0) = (1.5*T_B0(0,-2,0) + 1.83333333333333*T_B0(0,0,0) - 3.0*T_B0(0,-1,0) -
            0.333333333333333*T_B0(0,-3,0))*invDelta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       wk7_B0(0,0,0) = (0.322484932882161*T_B0(0,0,0) + 0.0658051057710389*T_B0(0,-3,0) + 0.376283677513354*T_B0(0,1,0)
            - 0.0394168524399447*T_B0(0,-2,0) - 0.00571369039775442*T_B0(0,-4,0) -
            0.719443173328855*T_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == -3 + block0np1){

       wk7_B0(0,0,0) = (0.00412637789557492*T_B0(0,-3,0) + 0.0367146847001261*T_B0(0,-2,0) +
            0.791245592765872*T_B0(0,1,0) - 0.197184333887745*T_B0(0,0,0) - 0.521455851089587*T_B0(0,-1,0) -
            0.113446470384241*T_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -4 + block0np1){

       wk7_B0(0,0,0) = (0.00932597985049999*T_B0(0,3,0) + 0.727822147724592*T_B0(0,1,0) + 0.082033432844602*T_B0(0,-2,0)
            - 0.0451033223343881*T_B0(0,0,0) - 0.652141084861241*T_B0(0,-1,0) -
            0.121937153224065*T_B0(0,2,0))*invDelta1block0;

   }

   else{

       wk7_B0(0,0,0) = (-(2.0/3.0)*T_B0(0,-1,0) - (1.0/12.0)*T_B0(0,2,0) + ((1.0/12.0))*T_B0(0,-2,0) +
            ((2.0/3.0))*T_B0(0,1,0))*invDelta1block0;

   }

}

void opensbliblock00Kernel016(const ACC<double> &u0_B0, ACC<double> &wk8_B0)
{
    wk8_B0(0,0,0) = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel017(const ACC<double> &u1_B0, ACC<double> &wk9_B0)
{
    wk9_B0(0,0,0) = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel018(const ACC<double> &u2_B0, ACC<double> &wk10_B0)
{
    wk10_B0(0,0,0) = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel019(const ACC<double> &T_B0, ACC<double> &wk11_B0)
{
    wk11_B0(0,0,0) = (-(2.0/3.0)*T_B0(0,0,-1) - (1.0/12.0)*T_B0(0,0,2) + ((1.0/12.0))*T_B0(0,0,-2) +
      ((2.0/3.0))*T_B0(0,0,1))*invDelta2block0;

}

 void opensbliblock00Kernel032(const ACC<double> &D11_B0, const ACC<double> &H_B0, const ACC<double> &p_B0, const
ACC<double> &phi_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, const ACC<double>
&wk0_B0, const ACC<double> &wk10_B0, const ACC<double> &wk1_B0, const ACC<double> &wk2_B0, const ACC<double> &wk4_B0,
const ACC<double> &wk5_B0, const ACC<double> &wk6_B0, const ACC<double> &wk8_B0, const ACC<double> &wk9_B0, ACC<double>
&Residual0_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0, ACC<double> &Residual3_B0, ACC<double>
&Residual4_B0, const int *idx)
{
   double d1_H_dx = 0.0;
   double d1_H_dy = 0.0;
   double d1_H_dz = 0.0;
   double d1_Hrho_dx = 0.0;
   double d1_Hrho_dy = 0.0;
   double d1_Hrho_dz = 0.0;
   double d1_Hrhou0_dx = 0.0;
   double d1_Hrhou1_dy = 0.0;
   double d1_Hrhou2_dz = 0.0;
   double d1_Hu0_dx = 0.0;
   double d1_Hu1_dy = 0.0;
   double d1_Hu2_dz = 0.0;
   double d1_p_dx = 0.0;
   double d1_p_dy = 0.0;
   double d1_p_dz = 0.0;
   double d1_rho_dx = 0.0;
   double d1_rho_dy = 0.0;
   double d1_rho_dz = 0.0;
   double d1_rhou0_dx = 0.0;
   double d1_rhou0_dy = 0.0;
   double d1_rhou0_dz = 0.0;
   double d1_rhou0u0_dx = 0.0;
   double d1_rhou0u1_dx = 0.0;
   double d1_rhou0u2_dx = 0.0;
   double d1_rhou1_dx = 0.0;
   double d1_rhou1_dy = 0.0;
   double d1_rhou1_dz = 0.0;
   double d1_rhou1u0_dy = 0.0;
   double d1_rhou1u1_dy = 0.0;
   double d1_rhou1u2_dy = 0.0;
   double d1_rhou2_dx = 0.0;
   double d1_rhou2_dy = 0.0;
   double d1_rhou2_dz = 0.0;
   double d1_rhou2u0_dz = 0.0;
   double d1_rhou2u1_dz = 0.0;
   double d1_rhou2u2_dz = 0.0;
   double d1_u0u0_dx = 0.0;
   double d1_u0u1_dx = 0.0;
   double d1_u0u1_dy = 0.0;
   double d1_u0u2_dx = 0.0;
   double d1_u0u2_dz = 0.0;
   double d1_u1u1_dy = 0.0;
   double d1_u1u2_dy = 0.0;
   double d1_u1u2_dz = 0.0;
   double d1_u2u2_dz = 0.0;
   if (idx[1] == 0){

       d1_H_dy = (3.0*H_B0(0,1,0) + 0.333333333333333*H_B0(0,3,0) - 1.5*H_B0(0,2,0) -
            1.83333333333333*H_B0(0,0,0))*invDelta1block0;

       d1_Hrho_dy = (3.0*H_B0(0,1,0)*rho_B0(0,1,0) + 0.333333333333333*H_B0(0,3,0)*rho_B0(0,3,0) -
            1.5*H_B0(0,2,0)*rho_B0(0,2,0) - 1.83333333333333*H_B0(0,0,0)*rho_B0(0,0,0))*invDelta1block0;

       d1_Hrhou1_dy = (3.0*H_B0(0,1,0)*rhou1_B0(0,1,0) + 0.333333333333333*H_B0(0,3,0)*rhou1_B0(0,3,0) -
            1.5*H_B0(0,2,0)*rhou1_B0(0,2,0) - 1.83333333333333*H_B0(0,0,0)*rhou1_B0(0,0,0))*invDelta1block0;

       d1_Hu1_dy = (3.0*H_B0(0,1,0)*u1_B0(0,1,0) + 0.333333333333333*H_B0(0,3,0)*u1_B0(0,3,0) -
            1.5*H_B0(0,2,0)*u1_B0(0,2,0) - 1.83333333333333*H_B0(0,0,0)*u1_B0(0,0,0))*invDelta1block0;

       d1_p_dy = (3.0*p_B0(0,1,0) + 0.333333333333333*p_B0(0,3,0) - 1.5*p_B0(0,2,0) -
            1.83333333333333*p_B0(0,0,0))*invDelta1block0;

       d1_rho_dy = (3.0*rho_B0(0,1,0) + 0.333333333333333*rho_B0(0,3,0) - 1.5*rho_B0(0,2,0) -
            1.83333333333333*rho_B0(0,0,0))*invDelta1block0;

       d1_rhou0_dy = (3.0*rhou0_B0(0,1,0) + 0.333333333333333*rhou0_B0(0,3,0) - 1.5*rhou0_B0(0,2,0) -
            1.83333333333333*rhou0_B0(0,0,0))*invDelta1block0;

       d1_rhou1_dy = (3.0*rhou1_B0(0,1,0) + 0.333333333333333*rhou1_B0(0,3,0) - 1.5*rhou1_B0(0,2,0) -
            1.83333333333333*rhou1_B0(0,0,0))*invDelta1block0;

       d1_rhou1u0_dy = (3.0*u0_B0(0,1,0)*rhou1_B0(0,1,0) + 0.333333333333333*u0_B0(0,3,0)*rhou1_B0(0,3,0) -
            1.5*u0_B0(0,2,0)*rhou1_B0(0,2,0) - 1.83333333333333*u0_B0(0,0,0)*rhou1_B0(0,0,0))*invDelta1block0;

       d1_rhou1u1_dy = (3.0*u1_B0(0,1,0)*rhou1_B0(0,1,0) + 0.333333333333333*u1_B0(0,3,0)*rhou1_B0(0,3,0) -
            1.5*u1_B0(0,2,0)*rhou1_B0(0,2,0) - 1.83333333333333*u1_B0(0,0,0)*rhou1_B0(0,0,0))*invDelta1block0;

       d1_rhou1u2_dy = (3.0*u2_B0(0,1,0)*rhou1_B0(0,1,0) + 0.333333333333333*u2_B0(0,3,0)*rhou1_B0(0,3,0) -
            1.5*u2_B0(0,2,0)*rhou1_B0(0,2,0) - 1.83333333333333*u2_B0(0,0,0)*rhou1_B0(0,0,0))*invDelta1block0;

       d1_rhou2_dy = (3.0*rhou2_B0(0,1,0) + 0.333333333333333*rhou2_B0(0,3,0) - 1.5*rhou2_B0(0,2,0) -
            1.83333333333333*rhou2_B0(0,0,0))*invDelta1block0;

       d1_u0u1_dy = (3.0*u0_B0(0,1,0)*u1_B0(0,1,0) + 0.333333333333333*u0_B0(0,3,0)*u1_B0(0,3,0) -
            1.5*u0_B0(0,2,0)*u1_B0(0,2,0) - 1.83333333333333*u0_B0(0,0,0)*u1_B0(0,0,0))*invDelta1block0;

       d1_u1u1_dy = (3.0*(u1_B0(0,1,0)*u1_B0(0,1,0)) + 0.333333333333333*(u1_B0(0,3,0)*u1_B0(0,3,0)) -
            1.5*(u1_B0(0,2,0)*u1_B0(0,2,0)) - 1.83333333333333*(u1_B0(0,0,0)*u1_B0(0,0,0)))*invDelta1block0;

       d1_u1u2_dy = (3.0*u1_B0(0,1,0)*u2_B0(0,1,0) + 0.333333333333333*u1_B0(0,3,0)*u2_B0(0,3,0) -
            1.5*u1_B0(0,2,0)*u2_B0(0,2,0) - 1.83333333333333*u1_B0(0,0,0)*u2_B0(0,0,0))*invDelta1block0;

   }

   else if (idx[1] == 1){

       d1_H_dy = (0.0394168524399447*H_B0(0,2,0) + 0.00571369039775442*H_B0(0,4,0) + 0.719443173328855*H_B0(0,1,0) -
            0.322484932882161*H_B0(0,0,0) - 0.0658051057710389*H_B0(0,3,0) -
            0.376283677513354*H_B0(0,-1,0))*invDelta1block0;

       d1_Hrho_dy = (0.0394168524399447*H_B0(0,2,0)*rho_B0(0,2,0) + 0.00571369039775442*H_B0(0,4,0)*rho_B0(0,4,0) +
            0.719443173328855*H_B0(0,1,0)*rho_B0(0,1,0) - 0.322484932882161*H_B0(0,0,0)*rho_B0(0,0,0) -
            0.0658051057710389*H_B0(0,3,0)*rho_B0(0,3,0) -
            0.376283677513354*H_B0(0,-1,0)*rho_B0(0,-1,0))*invDelta1block0;

       d1_Hrhou1_dy = (0.0394168524399447*H_B0(0,2,0)*rhou1_B0(0,2,0) + 0.00571369039775442*H_B0(0,4,0)*rhou1_B0(0,4,0)
            + 0.719443173328855*H_B0(0,1,0)*rhou1_B0(0,1,0) - 0.322484932882161*H_B0(0,0,0)*rhou1_B0(0,0,0) -
            0.0658051057710389*H_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.376283677513354*H_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_Hu1_dy = (0.0394168524399447*H_B0(0,2,0)*u1_B0(0,2,0) + 0.00571369039775442*H_B0(0,4,0)*u1_B0(0,4,0) +
            0.719443173328855*H_B0(0,1,0)*u1_B0(0,1,0) - 0.322484932882161*H_B0(0,0,0)*u1_B0(0,0,0) -
            0.0658051057710389*H_B0(0,3,0)*u1_B0(0,3,0) -
            0.376283677513354*H_B0(0,-1,0)*u1_B0(0,-1,0))*invDelta1block0;

       d1_p_dy = (0.0394168524399447*p_B0(0,2,0) + 0.00571369039775442*p_B0(0,4,0) + 0.719443173328855*p_B0(0,1,0) -
            0.322484932882161*p_B0(0,0,0) - 0.0658051057710389*p_B0(0,3,0) -
            0.376283677513354*p_B0(0,-1,0))*invDelta1block0;

       d1_rho_dy = (0.0394168524399447*rho_B0(0,2,0) + 0.00571369039775442*rho_B0(0,4,0) +
            0.719443173328855*rho_B0(0,1,0) - 0.322484932882161*rho_B0(0,0,0) - 0.0658051057710389*rho_B0(0,3,0) -
            0.376283677513354*rho_B0(0,-1,0))*invDelta1block0;

       d1_rhou0_dy = (0.0394168524399447*rhou0_B0(0,2,0) + 0.00571369039775442*rhou0_B0(0,4,0) +
            0.719443173328855*rhou0_B0(0,1,0) - 0.322484932882161*rhou0_B0(0,0,0) - 0.0658051057710389*rhou0_B0(0,3,0) -
            0.376283677513354*rhou0_B0(0,-1,0))*invDelta1block0;

       d1_rhou1_dy = (0.0394168524399447*rhou1_B0(0,2,0) + 0.00571369039775442*rhou1_B0(0,4,0) +
            0.719443173328855*rhou1_B0(0,1,0) - 0.322484932882161*rhou1_B0(0,0,0) - 0.0658051057710389*rhou1_B0(0,3,0) -
            0.376283677513354*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u0_dy = (0.0394168524399447*u0_B0(0,2,0)*rhou1_B0(0,2,0) +
            0.00571369039775442*u0_B0(0,4,0)*rhou1_B0(0,4,0) + 0.719443173328855*u0_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.322484932882161*u0_B0(0,0,0)*rhou1_B0(0,0,0) - 0.0658051057710389*u0_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.376283677513354*u0_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u1_dy = (0.0394168524399447*u1_B0(0,2,0)*rhou1_B0(0,2,0) +
            0.00571369039775442*u1_B0(0,4,0)*rhou1_B0(0,4,0) + 0.719443173328855*u1_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.322484932882161*u1_B0(0,0,0)*rhou1_B0(0,0,0) - 0.0658051057710389*u1_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.376283677513354*u1_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u2_dy = (0.0394168524399447*u2_B0(0,2,0)*rhou1_B0(0,2,0) +
            0.00571369039775442*u2_B0(0,4,0)*rhou1_B0(0,4,0) + 0.719443173328855*u2_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.322484932882161*u2_B0(0,0,0)*rhou1_B0(0,0,0) - 0.0658051057710389*u2_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.376283677513354*u2_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou2_dy = (0.0394168524399447*rhou2_B0(0,2,0) + 0.00571369039775442*rhou2_B0(0,4,0) +
            0.719443173328855*rhou2_B0(0,1,0) - 0.322484932882161*rhou2_B0(0,0,0) - 0.0658051057710389*rhou2_B0(0,3,0) -
            0.376283677513354*rhou2_B0(0,-1,0))*invDelta1block0;

       d1_u0u1_dy = (0.0394168524399447*u0_B0(0,2,0)*u1_B0(0,2,0) + 0.00571369039775442*u0_B0(0,4,0)*u1_B0(0,4,0) +
            0.719443173328855*u0_B0(0,1,0)*u1_B0(0,1,0) - 0.322484932882161*u0_B0(0,0,0)*u1_B0(0,0,0) -
            0.0658051057710389*u0_B0(0,3,0)*u1_B0(0,3,0) -
            0.376283677513354*u0_B0(0,-1,0)*u1_B0(0,-1,0))*invDelta1block0;

       d1_u1u1_dy = (0.0394168524399447*(u1_B0(0,2,0)*u1_B0(0,2,0)) + 0.00571369039775442*(u1_B0(0,4,0)*u1_B0(0,4,0)) +
            0.719443173328855*(u1_B0(0,1,0)*u1_B0(0,1,0)) - 0.322484932882161*(u1_B0(0,0,0)*u1_B0(0,0,0)) -
            0.0658051057710389*(u1_B0(0,3,0)*u1_B0(0,3,0)) -
            0.376283677513354*(u1_B0(0,-1,0)*u1_B0(0,-1,0)))*invDelta1block0;

       d1_u1u2_dy = (0.0394168524399447*u1_B0(0,2,0)*u2_B0(0,2,0) + 0.00571369039775442*u1_B0(0,4,0)*u2_B0(0,4,0) +
            0.719443173328855*u1_B0(0,1,0)*u2_B0(0,1,0) - 0.322484932882161*u1_B0(0,0,0)*u2_B0(0,0,0) -
            0.0658051057710389*u1_B0(0,3,0)*u2_B0(0,3,0) -
            0.376283677513354*u1_B0(0,-1,0)*u2_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 2){

       d1_H_dy = (0.197184333887745*H_B0(0,0,0) + 0.521455851089587*H_B0(0,1,0) + 0.113446470384241*H_B0(0,-2,0) -
            0.00412637789557492*H_B0(0,3,0) - 0.0367146847001261*H_B0(0,2,0) -
            0.791245592765872*H_B0(0,-1,0))*invDelta1block0;

       d1_Hrho_dy = (0.197184333887745*H_B0(0,0,0)*rho_B0(0,0,0) + 0.521455851089587*H_B0(0,1,0)*rho_B0(0,1,0) +
            0.113446470384241*H_B0(0,-2,0)*rho_B0(0,-2,0) - 0.00412637789557492*H_B0(0,3,0)*rho_B0(0,3,0) -
            0.0367146847001261*H_B0(0,2,0)*rho_B0(0,2,0) -
            0.791245592765872*H_B0(0,-1,0)*rho_B0(0,-1,0))*invDelta1block0;

       d1_Hrhou1_dy = (0.197184333887745*H_B0(0,0,0)*rhou1_B0(0,0,0) + 0.521455851089587*H_B0(0,1,0)*rhou1_B0(0,1,0) +
            0.113446470384241*H_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00412637789557492*H_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.0367146847001261*H_B0(0,2,0)*rhou1_B0(0,2,0) -
            0.791245592765872*H_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_Hu1_dy = (0.197184333887745*H_B0(0,0,0)*u1_B0(0,0,0) + 0.521455851089587*H_B0(0,1,0)*u1_B0(0,1,0) +
            0.113446470384241*H_B0(0,-2,0)*u1_B0(0,-2,0) - 0.00412637789557492*H_B0(0,3,0)*u1_B0(0,3,0) -
            0.0367146847001261*H_B0(0,2,0)*u1_B0(0,2,0) -
            0.791245592765872*H_B0(0,-1,0)*u1_B0(0,-1,0))*invDelta1block0;

       d1_p_dy = (0.197184333887745*p_B0(0,0,0) + 0.521455851089587*p_B0(0,1,0) + 0.113446470384241*p_B0(0,-2,0) -
            0.00412637789557492*p_B0(0,3,0) - 0.0367146847001261*p_B0(0,2,0) -
            0.791245592765872*p_B0(0,-1,0))*invDelta1block0;

       d1_rho_dy = (0.197184333887745*rho_B0(0,0,0) + 0.521455851089587*rho_B0(0,1,0) + 0.113446470384241*rho_B0(0,-2,0)
            - 0.00412637789557492*rho_B0(0,3,0) - 0.0367146847001261*rho_B0(0,2,0) -
            0.791245592765872*rho_B0(0,-1,0))*invDelta1block0;

       d1_rhou0_dy = (0.197184333887745*rhou0_B0(0,0,0) + 0.521455851089587*rhou0_B0(0,1,0) +
            0.113446470384241*rhou0_B0(0,-2,0) - 0.00412637789557492*rhou0_B0(0,3,0) -
            0.0367146847001261*rhou0_B0(0,2,0) - 0.791245592765872*rhou0_B0(0,-1,0))*invDelta1block0;

       d1_rhou1_dy = (0.197184333887745*rhou1_B0(0,0,0) + 0.521455851089587*rhou1_B0(0,1,0) +
            0.113446470384241*rhou1_B0(0,-2,0) - 0.00412637789557492*rhou1_B0(0,3,0) -
            0.0367146847001261*rhou1_B0(0,2,0) - 0.791245592765872*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u0_dy = (0.197184333887745*u0_B0(0,0,0)*rhou1_B0(0,0,0) + 0.521455851089587*u0_B0(0,1,0)*rhou1_B0(0,1,0)
            + 0.113446470384241*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00412637789557492*u0_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.0367146847001261*u0_B0(0,2,0)*rhou1_B0(0,2,0) -
            0.791245592765872*u0_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u1_dy = (0.197184333887745*u1_B0(0,0,0)*rhou1_B0(0,0,0) + 0.521455851089587*u1_B0(0,1,0)*rhou1_B0(0,1,0)
            + 0.113446470384241*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00412637789557492*u1_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.0367146847001261*u1_B0(0,2,0)*rhou1_B0(0,2,0) -
            0.791245592765872*u1_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u2_dy = (0.197184333887745*u2_B0(0,0,0)*rhou1_B0(0,0,0) + 0.521455851089587*u2_B0(0,1,0)*rhou1_B0(0,1,0)
            + 0.113446470384241*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00412637789557492*u2_B0(0,3,0)*rhou1_B0(0,3,0) -
            0.0367146847001261*u2_B0(0,2,0)*rhou1_B0(0,2,0) -
            0.791245592765872*u2_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou2_dy = (0.197184333887745*rhou2_B0(0,0,0) + 0.521455851089587*rhou2_B0(0,1,0) +
            0.113446470384241*rhou2_B0(0,-2,0) - 0.00412637789557492*rhou2_B0(0,3,0) -
            0.0367146847001261*rhou2_B0(0,2,0) - 0.791245592765872*rhou2_B0(0,-1,0))*invDelta1block0;

       d1_u0u1_dy = (0.197184333887745*u0_B0(0,0,0)*u1_B0(0,0,0) + 0.521455851089587*u0_B0(0,1,0)*u1_B0(0,1,0) +
            0.113446470384241*u0_B0(0,-2,0)*u1_B0(0,-2,0) - 0.00412637789557492*u0_B0(0,3,0)*u1_B0(0,3,0) -
            0.0367146847001261*u0_B0(0,2,0)*u1_B0(0,2,0) -
            0.791245592765872*u0_B0(0,-1,0)*u1_B0(0,-1,0))*invDelta1block0;

       d1_u1u1_dy = (0.197184333887745*(u1_B0(0,0,0)*u1_B0(0,0,0)) + 0.521455851089587*(u1_B0(0,1,0)*u1_B0(0,1,0)) +
            0.113446470384241*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) - 0.00412637789557492*(u1_B0(0,3,0)*u1_B0(0,3,0)) -
            0.0367146847001261*(u1_B0(0,2,0)*u1_B0(0,2,0)) -
            0.791245592765872*(u1_B0(0,-1,0)*u1_B0(0,-1,0)))*invDelta1block0;

       d1_u1u2_dy = (0.197184333887745*u1_B0(0,0,0)*u2_B0(0,0,0) + 0.521455851089587*u1_B0(0,1,0)*u2_B0(0,1,0) +
            0.113446470384241*u1_B0(0,-2,0)*u2_B0(0,-2,0) - 0.00412637789557492*u1_B0(0,3,0)*u2_B0(0,3,0) -
            0.0367146847001261*u1_B0(0,2,0)*u2_B0(0,2,0) -
            0.791245592765872*u1_B0(0,-1,0)*u2_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 3){

       d1_H_dy = (0.0451033223343881*H_B0(0,0,0) + 0.652141084861241*H_B0(0,1,0) + 0.121937153224065*H_B0(0,-2,0) -
            0.00932597985049999*H_B0(0,-3,0) - 0.727822147724592*H_B0(0,-1,0) -
            0.082033432844602*H_B0(0,2,0))*invDelta1block0;

       d1_Hrho_dy = (0.0451033223343881*H_B0(0,0,0)*rho_B0(0,0,0) + 0.652141084861241*H_B0(0,1,0)*rho_B0(0,1,0) +
            0.121937153224065*H_B0(0,-2,0)*rho_B0(0,-2,0) - 0.00932597985049999*H_B0(0,-3,0)*rho_B0(0,-3,0) -
            0.727822147724592*H_B0(0,-1,0)*rho_B0(0,-1,0) -
            0.082033432844602*H_B0(0,2,0)*rho_B0(0,2,0))*invDelta1block0;

       d1_Hrhou1_dy = (0.0451033223343881*H_B0(0,0,0)*rhou1_B0(0,0,0) + 0.652141084861241*H_B0(0,1,0)*rhou1_B0(0,1,0) +
            0.121937153224065*H_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00932597985049999*H_B0(0,-3,0)*rhou1_B0(0,-3,0) -
            0.727822147724592*H_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.082033432844602*H_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_Hu1_dy = (0.0451033223343881*H_B0(0,0,0)*u1_B0(0,0,0) + 0.652141084861241*H_B0(0,1,0)*u1_B0(0,1,0) +
            0.121937153224065*H_B0(0,-2,0)*u1_B0(0,-2,0) - 0.00932597985049999*H_B0(0,-3,0)*u1_B0(0,-3,0) -
            0.727822147724592*H_B0(0,-1,0)*u1_B0(0,-1,0) - 0.082033432844602*H_B0(0,2,0)*u1_B0(0,2,0))*invDelta1block0;

       d1_p_dy = (0.0451033223343881*p_B0(0,0,0) + 0.652141084861241*p_B0(0,1,0) + 0.121937153224065*p_B0(0,-2,0) -
            0.00932597985049999*p_B0(0,-3,0) - 0.727822147724592*p_B0(0,-1,0) -
            0.082033432844602*p_B0(0,2,0))*invDelta1block0;

       d1_rho_dy = (0.0451033223343881*rho_B0(0,0,0) + 0.652141084861241*rho_B0(0,1,0) +
            0.121937153224065*rho_B0(0,-2,0) - 0.00932597985049999*rho_B0(0,-3,0) - 0.727822147724592*rho_B0(0,-1,0) -
            0.082033432844602*rho_B0(0,2,0))*invDelta1block0;

       d1_rhou0_dy = (0.0451033223343881*rhou0_B0(0,0,0) + 0.652141084861241*rhou0_B0(0,1,0) +
            0.121937153224065*rhou0_B0(0,-2,0) - 0.00932597985049999*rhou0_B0(0,-3,0) -
            0.727822147724592*rhou0_B0(0,-1,0) - 0.082033432844602*rhou0_B0(0,2,0))*invDelta1block0;

       d1_rhou1_dy = (0.0451033223343881*rhou1_B0(0,0,0) + 0.652141084861241*rhou1_B0(0,1,0) +
            0.121937153224065*rhou1_B0(0,-2,0) - 0.00932597985049999*rhou1_B0(0,-3,0) -
            0.727822147724592*rhou1_B0(0,-1,0) - 0.082033432844602*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u0_dy = (0.0451033223343881*u0_B0(0,0,0)*rhou1_B0(0,0,0) + 0.652141084861241*u0_B0(0,1,0)*rhou1_B0(0,1,0)
            + 0.121937153224065*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00932597985049999*u0_B0(0,-3,0)*rhou1_B0(0,-3,0) -
            0.727822147724592*u0_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.082033432844602*u0_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u1_dy = (0.0451033223343881*u1_B0(0,0,0)*rhou1_B0(0,0,0) + 0.652141084861241*u1_B0(0,1,0)*rhou1_B0(0,1,0)
            + 0.121937153224065*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00932597985049999*u1_B0(0,-3,0)*rhou1_B0(0,-3,0) -
            0.727822147724592*u1_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.082033432844602*u1_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u2_dy = (0.0451033223343881*u2_B0(0,0,0)*rhou1_B0(0,0,0) + 0.652141084861241*u2_B0(0,1,0)*rhou1_B0(0,1,0)
            + 0.121937153224065*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00932597985049999*u2_B0(0,-3,0)*rhou1_B0(0,-3,0) -
            0.727822147724592*u2_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.082033432844602*u2_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou2_dy = (0.0451033223343881*rhou2_B0(0,0,0) + 0.652141084861241*rhou2_B0(0,1,0) +
            0.121937153224065*rhou2_B0(0,-2,0) - 0.00932597985049999*rhou2_B0(0,-3,0) -
            0.727822147724592*rhou2_B0(0,-1,0) - 0.082033432844602*rhou2_B0(0,2,0))*invDelta1block0;

       d1_u0u1_dy = (0.0451033223343881*u0_B0(0,0,0)*u1_B0(0,0,0) + 0.652141084861241*u0_B0(0,1,0)*u1_B0(0,1,0) +
            0.121937153224065*u0_B0(0,-2,0)*u1_B0(0,-2,0) - 0.00932597985049999*u0_B0(0,-3,0)*u1_B0(0,-3,0) -
            0.727822147724592*u0_B0(0,-1,0)*u1_B0(0,-1,0) -
            0.082033432844602*u0_B0(0,2,0)*u1_B0(0,2,0))*invDelta1block0;

       d1_u1u1_dy = (0.0451033223343881*(u1_B0(0,0,0)*u1_B0(0,0,0)) + 0.652141084861241*(u1_B0(0,1,0)*u1_B0(0,1,0)) +
            0.121937153224065*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) - 0.00932597985049999*(u1_B0(0,-3,0)*u1_B0(0,-3,0)) -
            0.727822147724592*(u1_B0(0,-1,0)*u1_B0(0,-1,0)) -
            0.082033432844602*(u1_B0(0,2,0)*u1_B0(0,2,0)))*invDelta1block0;

       d1_u1u2_dy = (0.0451033223343881*u1_B0(0,0,0)*u2_B0(0,0,0) + 0.652141084861241*u1_B0(0,1,0)*u2_B0(0,1,0) +
            0.121937153224065*u1_B0(0,-2,0)*u2_B0(0,-2,0) - 0.00932597985049999*u1_B0(0,-3,0)*u2_B0(0,-3,0) -
            0.727822147724592*u1_B0(0,-1,0)*u2_B0(0,-1,0) -
            0.082033432844602*u1_B0(0,2,0)*u2_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       d1_H_dy = (1.5*H_B0(0,-2,0) + 1.83333333333333*H_B0(0,0,0) - 3.0*H_B0(0,-1,0) -
            0.333333333333333*H_B0(0,-3,0))*invDelta1block0;

       d1_Hrho_dy = (1.5*H_B0(0,-2,0)*rho_B0(0,-2,0) + 1.83333333333333*H_B0(0,0,0)*rho_B0(0,0,0) -
            3.0*H_B0(0,-1,0)*rho_B0(0,-1,0) - 0.333333333333333*H_B0(0,-3,0)*rho_B0(0,-3,0))*invDelta1block0;

       d1_Hrhou1_dy = (1.5*H_B0(0,-2,0)*rhou1_B0(0,-2,0) + 1.83333333333333*H_B0(0,0,0)*rhou1_B0(0,0,0) -
            3.0*H_B0(0,-1,0)*rhou1_B0(0,-1,0) - 0.333333333333333*H_B0(0,-3,0)*rhou1_B0(0,-3,0))*invDelta1block0;

       d1_Hu1_dy = (1.5*H_B0(0,-2,0)*u1_B0(0,-2,0) + 1.83333333333333*H_B0(0,0,0)*u1_B0(0,0,0) -
            3.0*H_B0(0,-1,0)*u1_B0(0,-1,0) - 0.333333333333333*H_B0(0,-3,0)*u1_B0(0,-3,0))*invDelta1block0;

       d1_p_dy = (1.5*p_B0(0,-2,0) + 1.83333333333333*p_B0(0,0,0) - 3.0*p_B0(0,-1,0) -
            0.333333333333333*p_B0(0,-3,0))*invDelta1block0;

       d1_rho_dy = (1.5*rho_B0(0,-2,0) + 1.83333333333333*rho_B0(0,0,0) - 3.0*rho_B0(0,-1,0) -
            0.333333333333333*rho_B0(0,-3,0))*invDelta1block0;

       d1_rhou0_dy = (1.5*rhou0_B0(0,-2,0) + 1.83333333333333*rhou0_B0(0,0,0) - 3.0*rhou0_B0(0,-1,0) -
            0.333333333333333*rhou0_B0(0,-3,0))*invDelta1block0;

       d1_rhou1_dy = (1.5*rhou1_B0(0,-2,0) + 1.83333333333333*rhou1_B0(0,0,0) - 3.0*rhou1_B0(0,-1,0) -
            0.333333333333333*rhou1_B0(0,-3,0))*invDelta1block0;

       d1_rhou1u0_dy = (1.5*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) + 1.83333333333333*u0_B0(0,0,0)*rhou1_B0(0,0,0) -
            3.0*u0_B0(0,-1,0)*rhou1_B0(0,-1,0) - 0.333333333333333*u0_B0(0,-3,0)*rhou1_B0(0,-3,0))*invDelta1block0;

       d1_rhou1u1_dy = (1.5*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) + 1.83333333333333*u1_B0(0,0,0)*rhou1_B0(0,0,0) -
            3.0*u1_B0(0,-1,0)*rhou1_B0(0,-1,0) - 0.333333333333333*u1_B0(0,-3,0)*rhou1_B0(0,-3,0))*invDelta1block0;

       d1_rhou1u2_dy = (1.5*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) + 1.83333333333333*u2_B0(0,0,0)*rhou1_B0(0,0,0) -
            3.0*u2_B0(0,-1,0)*rhou1_B0(0,-1,0) - 0.333333333333333*u2_B0(0,-3,0)*rhou1_B0(0,-3,0))*invDelta1block0;

       d1_rhou2_dy = (1.5*rhou2_B0(0,-2,0) + 1.83333333333333*rhou2_B0(0,0,0) - 3.0*rhou2_B0(0,-1,0) -
            0.333333333333333*rhou2_B0(0,-3,0))*invDelta1block0;

       d1_u0u1_dy = (1.5*u0_B0(0,-2,0)*u1_B0(0,-2,0) + 1.83333333333333*u0_B0(0,0,0)*u1_B0(0,0,0) -
            3.0*u0_B0(0,-1,0)*u1_B0(0,-1,0) - 0.333333333333333*u0_B0(0,-3,0)*u1_B0(0,-3,0))*invDelta1block0;

       d1_u1u1_dy = (1.5*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) + 1.83333333333333*(u1_B0(0,0,0)*u1_B0(0,0,0)) -
            3.0*(u1_B0(0,-1,0)*u1_B0(0,-1,0)) - 0.333333333333333*(u1_B0(0,-3,0)*u1_B0(0,-3,0)))*invDelta1block0;

       d1_u1u2_dy = (1.5*u1_B0(0,-2,0)*u2_B0(0,-2,0) + 1.83333333333333*u1_B0(0,0,0)*u2_B0(0,0,0) -
            3.0*u1_B0(0,-1,0)*u2_B0(0,-1,0) - 0.333333333333333*u1_B0(0,-3,0)*u2_B0(0,-3,0))*invDelta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       d1_H_dy = (0.322484932882161*H_B0(0,0,0) + 0.0658051057710389*H_B0(0,-3,0) + 0.376283677513354*H_B0(0,1,0) -
            0.0394168524399447*H_B0(0,-2,0) - 0.00571369039775442*H_B0(0,-4,0) -
            0.719443173328855*H_B0(0,-1,0))*invDelta1block0;

       d1_Hrho_dy = (0.322484932882161*H_B0(0,0,0)*rho_B0(0,0,0) + 0.0658051057710389*H_B0(0,-3,0)*rho_B0(0,-3,0) +
            0.376283677513354*H_B0(0,1,0)*rho_B0(0,1,0) - 0.0394168524399447*H_B0(0,-2,0)*rho_B0(0,-2,0) -
            0.00571369039775442*H_B0(0,-4,0)*rho_B0(0,-4,0) -
            0.719443173328855*H_B0(0,-1,0)*rho_B0(0,-1,0))*invDelta1block0;

       d1_Hrhou1_dy = (0.322484932882161*H_B0(0,0,0)*rhou1_B0(0,0,0) + 0.0658051057710389*H_B0(0,-3,0)*rhou1_B0(0,-3,0)
            + 0.376283677513354*H_B0(0,1,0)*rhou1_B0(0,1,0) - 0.0394168524399447*H_B0(0,-2,0)*rhou1_B0(0,-2,0) -
            0.00571369039775442*H_B0(0,-4,0)*rhou1_B0(0,-4,0) -
            0.719443173328855*H_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_Hu1_dy = (0.322484932882161*H_B0(0,0,0)*u1_B0(0,0,0) + 0.0658051057710389*H_B0(0,-3,0)*u1_B0(0,-3,0) +
            0.376283677513354*H_B0(0,1,0)*u1_B0(0,1,0) - 0.0394168524399447*H_B0(0,-2,0)*u1_B0(0,-2,0) -
            0.00571369039775442*H_B0(0,-4,0)*u1_B0(0,-4,0) -
            0.719443173328855*H_B0(0,-1,0)*u1_B0(0,-1,0))*invDelta1block0;

       d1_p_dy = (0.322484932882161*p_B0(0,0,0) + 0.0658051057710389*p_B0(0,-3,0) + 0.376283677513354*p_B0(0,1,0) -
            0.0394168524399447*p_B0(0,-2,0) - 0.00571369039775442*p_B0(0,-4,0) -
            0.719443173328855*p_B0(0,-1,0))*invDelta1block0;

       d1_rho_dy = (0.322484932882161*rho_B0(0,0,0) + 0.0658051057710389*rho_B0(0,-3,0) +
            0.376283677513354*rho_B0(0,1,0) - 0.0394168524399447*rho_B0(0,-2,0) - 0.00571369039775442*rho_B0(0,-4,0) -
            0.719443173328855*rho_B0(0,-1,0))*invDelta1block0;

       d1_rhou0_dy = (0.322484932882161*rhou0_B0(0,0,0) + 0.0658051057710389*rhou0_B0(0,-3,0) +
            0.376283677513354*rhou0_B0(0,1,0) - 0.0394168524399447*rhou0_B0(0,-2,0) -
            0.00571369039775442*rhou0_B0(0,-4,0) - 0.719443173328855*rhou0_B0(0,-1,0))*invDelta1block0;

       d1_rhou1_dy = (0.322484932882161*rhou1_B0(0,0,0) + 0.0658051057710389*rhou1_B0(0,-3,0) +
            0.376283677513354*rhou1_B0(0,1,0) - 0.0394168524399447*rhou1_B0(0,-2,0) -
            0.00571369039775442*rhou1_B0(0,-4,0) - 0.719443173328855*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u0_dy = (0.322484932882161*u0_B0(0,0,0)*rhou1_B0(0,0,0) +
            0.0658051057710389*u0_B0(0,-3,0)*rhou1_B0(0,-3,0) + 0.376283677513354*u0_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.0394168524399447*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00571369039775442*u0_B0(0,-4,0)*rhou1_B0(0,-4,0) -
            0.719443173328855*u0_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u1_dy = (0.322484932882161*u1_B0(0,0,0)*rhou1_B0(0,0,0) +
            0.0658051057710389*u1_B0(0,-3,0)*rhou1_B0(0,-3,0) + 0.376283677513354*u1_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.0394168524399447*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00571369039775442*u1_B0(0,-4,0)*rhou1_B0(0,-4,0) -
            0.719443173328855*u1_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou1u2_dy = (0.322484932882161*u2_B0(0,0,0)*rhou1_B0(0,0,0) +
            0.0658051057710389*u2_B0(0,-3,0)*rhou1_B0(0,-3,0) + 0.376283677513354*u2_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.0394168524399447*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.00571369039775442*u2_B0(0,-4,0)*rhou1_B0(0,-4,0) -
            0.719443173328855*u2_B0(0,-1,0)*rhou1_B0(0,-1,0))*invDelta1block0;

       d1_rhou2_dy = (0.322484932882161*rhou2_B0(0,0,0) + 0.0658051057710389*rhou2_B0(0,-3,0) +
            0.376283677513354*rhou2_B0(0,1,0) - 0.0394168524399447*rhou2_B0(0,-2,0) -
            0.00571369039775442*rhou2_B0(0,-4,0) - 0.719443173328855*rhou2_B0(0,-1,0))*invDelta1block0;

       d1_u0u1_dy = (0.322484932882161*u0_B0(0,0,0)*u1_B0(0,0,0) + 0.0658051057710389*u0_B0(0,-3,0)*u1_B0(0,-3,0) +
            0.376283677513354*u0_B0(0,1,0)*u1_B0(0,1,0) - 0.0394168524399447*u0_B0(0,-2,0)*u1_B0(0,-2,0) -
            0.00571369039775442*u0_B0(0,-4,0)*u1_B0(0,-4,0) -
            0.719443173328855*u0_B0(0,-1,0)*u1_B0(0,-1,0))*invDelta1block0;

       d1_u1u1_dy = (0.322484932882161*(u1_B0(0,0,0)*u1_B0(0,0,0)) + 0.0658051057710389*(u1_B0(0,-3,0)*u1_B0(0,-3,0)) +
            0.376283677513354*(u1_B0(0,1,0)*u1_B0(0,1,0)) - 0.0394168524399447*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) -
            0.00571369039775442*(u1_B0(0,-4,0)*u1_B0(0,-4,0)) -
            0.719443173328855*(u1_B0(0,-1,0)*u1_B0(0,-1,0)))*invDelta1block0;

       d1_u1u2_dy = (0.322484932882161*u1_B0(0,0,0)*u2_B0(0,0,0) + 0.0658051057710389*u1_B0(0,-3,0)*u2_B0(0,-3,0) +
            0.376283677513354*u1_B0(0,1,0)*u2_B0(0,1,0) - 0.0394168524399447*u1_B0(0,-2,0)*u2_B0(0,-2,0) -
            0.00571369039775442*u1_B0(0,-4,0)*u2_B0(0,-4,0) -
            0.719443173328855*u1_B0(0,-1,0)*u2_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == -3 + block0np1){

       d1_H_dy = (0.00412637789557492*H_B0(0,-3,0) + 0.0367146847001261*H_B0(0,-2,0) + 0.791245592765872*H_B0(0,1,0) -
            0.197184333887745*H_B0(0,0,0) - 0.521455851089587*H_B0(0,-1,0) -
            0.113446470384241*H_B0(0,2,0))*invDelta1block0;

       d1_Hrho_dy = (0.00412637789557492*H_B0(0,-3,0)*rho_B0(0,-3,0) + 0.0367146847001261*H_B0(0,-2,0)*rho_B0(0,-2,0) +
            0.791245592765872*H_B0(0,1,0)*rho_B0(0,1,0) - 0.197184333887745*H_B0(0,0,0)*rho_B0(0,0,0) -
            0.521455851089587*H_B0(0,-1,0)*rho_B0(0,-1,0) -
            0.113446470384241*H_B0(0,2,0)*rho_B0(0,2,0))*invDelta1block0;

       d1_Hrhou1_dy = (0.00412637789557492*H_B0(0,-3,0)*rhou1_B0(0,-3,0) +
            0.0367146847001261*H_B0(0,-2,0)*rhou1_B0(0,-2,0) + 0.791245592765872*H_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.197184333887745*H_B0(0,0,0)*rhou1_B0(0,0,0) - 0.521455851089587*H_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.113446470384241*H_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_Hu1_dy = (0.00412637789557492*H_B0(0,-3,0)*u1_B0(0,-3,0) + 0.0367146847001261*H_B0(0,-2,0)*u1_B0(0,-2,0) +
            0.791245592765872*H_B0(0,1,0)*u1_B0(0,1,0) - 0.197184333887745*H_B0(0,0,0)*u1_B0(0,0,0) -
            0.521455851089587*H_B0(0,-1,0)*u1_B0(0,-1,0) - 0.113446470384241*H_B0(0,2,0)*u1_B0(0,2,0))*invDelta1block0;

       d1_p_dy = (0.00412637789557492*p_B0(0,-3,0) + 0.0367146847001261*p_B0(0,-2,0) + 0.791245592765872*p_B0(0,1,0) -
            0.197184333887745*p_B0(0,0,0) - 0.521455851089587*p_B0(0,-1,0) -
            0.113446470384241*p_B0(0,2,0))*invDelta1block0;

       d1_rho_dy = (0.00412637789557492*rho_B0(0,-3,0) + 0.0367146847001261*rho_B0(0,-2,0) +
            0.791245592765872*rho_B0(0,1,0) - 0.197184333887745*rho_B0(0,0,0) - 0.521455851089587*rho_B0(0,-1,0) -
            0.113446470384241*rho_B0(0,2,0))*invDelta1block0;

       d1_rhou0_dy = (0.00412637789557492*rhou0_B0(0,-3,0) + 0.0367146847001261*rhou0_B0(0,-2,0) +
            0.791245592765872*rhou0_B0(0,1,0) - 0.197184333887745*rhou0_B0(0,0,0) - 0.521455851089587*rhou0_B0(0,-1,0) -
            0.113446470384241*rhou0_B0(0,2,0))*invDelta1block0;

       d1_rhou1_dy = (0.00412637789557492*rhou1_B0(0,-3,0) + 0.0367146847001261*rhou1_B0(0,-2,0) +
            0.791245592765872*rhou1_B0(0,1,0) - 0.197184333887745*rhou1_B0(0,0,0) - 0.521455851089587*rhou1_B0(0,-1,0) -
            0.113446470384241*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u0_dy = (0.00412637789557492*u0_B0(0,-3,0)*rhou1_B0(0,-3,0) +
            0.0367146847001261*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) + 0.791245592765872*u0_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.197184333887745*u0_B0(0,0,0)*rhou1_B0(0,0,0) - 0.521455851089587*u0_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.113446470384241*u0_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u1_dy = (0.00412637789557492*u1_B0(0,-3,0)*rhou1_B0(0,-3,0) +
            0.0367146847001261*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) + 0.791245592765872*u1_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.197184333887745*u1_B0(0,0,0)*rhou1_B0(0,0,0) - 0.521455851089587*u1_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.113446470384241*u1_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u2_dy = (0.00412637789557492*u2_B0(0,-3,0)*rhou1_B0(0,-3,0) +
            0.0367146847001261*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) + 0.791245592765872*u2_B0(0,1,0)*rhou1_B0(0,1,0) -
            0.197184333887745*u2_B0(0,0,0)*rhou1_B0(0,0,0) - 0.521455851089587*u2_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.113446470384241*u2_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou2_dy = (0.00412637789557492*rhou2_B0(0,-3,0) + 0.0367146847001261*rhou2_B0(0,-2,0) +
            0.791245592765872*rhou2_B0(0,1,0) - 0.197184333887745*rhou2_B0(0,0,0) - 0.521455851089587*rhou2_B0(0,-1,0) -
            0.113446470384241*rhou2_B0(0,2,0))*invDelta1block0;

       d1_u0u1_dy = (0.00412637789557492*u0_B0(0,-3,0)*u1_B0(0,-3,0) + 0.0367146847001261*u0_B0(0,-2,0)*u1_B0(0,-2,0) +
            0.791245592765872*u0_B0(0,1,0)*u1_B0(0,1,0) - 0.197184333887745*u0_B0(0,0,0)*u1_B0(0,0,0) -
            0.521455851089587*u0_B0(0,-1,0)*u1_B0(0,-1,0) -
            0.113446470384241*u0_B0(0,2,0)*u1_B0(0,2,0))*invDelta1block0;

       d1_u1u1_dy = (0.00412637789557492*(u1_B0(0,-3,0)*u1_B0(0,-3,0)) +
            0.0367146847001261*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) + 0.791245592765872*(u1_B0(0,1,0)*u1_B0(0,1,0)) -
            0.197184333887745*(u1_B0(0,0,0)*u1_B0(0,0,0)) - 0.521455851089587*(u1_B0(0,-1,0)*u1_B0(0,-1,0)) -
            0.113446470384241*(u1_B0(0,2,0)*u1_B0(0,2,0)))*invDelta1block0;

       d1_u1u2_dy = (0.00412637789557492*u1_B0(0,-3,0)*u2_B0(0,-3,0) + 0.0367146847001261*u1_B0(0,-2,0)*u2_B0(0,-2,0) +
            0.791245592765872*u1_B0(0,1,0)*u2_B0(0,1,0) - 0.197184333887745*u1_B0(0,0,0)*u2_B0(0,0,0) -
            0.521455851089587*u1_B0(0,-1,0)*u2_B0(0,-1,0) -
            0.113446470384241*u1_B0(0,2,0)*u2_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -4 + block0np1){

       d1_H_dy = (0.00932597985049999*H_B0(0,3,0) + 0.727822147724592*H_B0(0,1,0) + 0.082033432844602*H_B0(0,-2,0) -
            0.0451033223343881*H_B0(0,0,0) - 0.652141084861241*H_B0(0,-1,0) -
            0.121937153224065*H_B0(0,2,0))*invDelta1block0;

       d1_Hrho_dy = (0.00932597985049999*H_B0(0,3,0)*rho_B0(0,3,0) + 0.727822147724592*H_B0(0,1,0)*rho_B0(0,1,0) +
            0.082033432844602*H_B0(0,-2,0)*rho_B0(0,-2,0) - 0.0451033223343881*H_B0(0,0,0)*rho_B0(0,0,0) -
            0.652141084861241*H_B0(0,-1,0)*rho_B0(0,-1,0) -
            0.121937153224065*H_B0(0,2,0)*rho_B0(0,2,0))*invDelta1block0;

       d1_Hrhou1_dy = (0.00932597985049999*H_B0(0,3,0)*rhou1_B0(0,3,0) + 0.727822147724592*H_B0(0,1,0)*rhou1_B0(0,1,0) +
            0.082033432844602*H_B0(0,-2,0)*rhou1_B0(0,-2,0) - 0.0451033223343881*H_B0(0,0,0)*rhou1_B0(0,0,0) -
            0.652141084861241*H_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.121937153224065*H_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_Hu1_dy = (0.00932597985049999*H_B0(0,3,0)*u1_B0(0,3,0) + 0.727822147724592*H_B0(0,1,0)*u1_B0(0,1,0) +
            0.082033432844602*H_B0(0,-2,0)*u1_B0(0,-2,0) - 0.0451033223343881*H_B0(0,0,0)*u1_B0(0,0,0) -
            0.652141084861241*H_B0(0,-1,0)*u1_B0(0,-1,0) - 0.121937153224065*H_B0(0,2,0)*u1_B0(0,2,0))*invDelta1block0;

       d1_p_dy = (0.00932597985049999*p_B0(0,3,0) + 0.727822147724592*p_B0(0,1,0) + 0.082033432844602*p_B0(0,-2,0) -
            0.0451033223343881*p_B0(0,0,0) - 0.652141084861241*p_B0(0,-1,0) -
            0.121937153224065*p_B0(0,2,0))*invDelta1block0;

       d1_rho_dy = (0.00932597985049999*rho_B0(0,3,0) + 0.727822147724592*rho_B0(0,1,0) +
            0.082033432844602*rho_B0(0,-2,0) - 0.0451033223343881*rho_B0(0,0,0) - 0.652141084861241*rho_B0(0,-1,0) -
            0.121937153224065*rho_B0(0,2,0))*invDelta1block0;

       d1_rhou0_dy = (0.00932597985049999*rhou0_B0(0,3,0) + 0.727822147724592*rhou0_B0(0,1,0) +
            0.082033432844602*rhou0_B0(0,-2,0) - 0.0451033223343881*rhou0_B0(0,0,0) - 0.652141084861241*rhou0_B0(0,-1,0)
            - 0.121937153224065*rhou0_B0(0,2,0))*invDelta1block0;

       d1_rhou1_dy = (0.00932597985049999*rhou1_B0(0,3,0) + 0.727822147724592*rhou1_B0(0,1,0) +
            0.082033432844602*rhou1_B0(0,-2,0) - 0.0451033223343881*rhou1_B0(0,0,0) - 0.652141084861241*rhou1_B0(0,-1,0)
            - 0.121937153224065*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u0_dy = (0.00932597985049999*u0_B0(0,3,0)*rhou1_B0(0,3,0) +
            0.727822147724592*u0_B0(0,1,0)*rhou1_B0(0,1,0) + 0.082033432844602*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) -
            0.0451033223343881*u0_B0(0,0,0)*rhou1_B0(0,0,0) - 0.652141084861241*u0_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.121937153224065*u0_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u1_dy = (0.00932597985049999*u1_B0(0,3,0)*rhou1_B0(0,3,0) +
            0.727822147724592*u1_B0(0,1,0)*rhou1_B0(0,1,0) + 0.082033432844602*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) -
            0.0451033223343881*u1_B0(0,0,0)*rhou1_B0(0,0,0) - 0.652141084861241*u1_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.121937153224065*u1_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou1u2_dy = (0.00932597985049999*u2_B0(0,3,0)*rhou1_B0(0,3,0) +
            0.727822147724592*u2_B0(0,1,0)*rhou1_B0(0,1,0) + 0.082033432844602*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) -
            0.0451033223343881*u2_B0(0,0,0)*rhou1_B0(0,0,0) - 0.652141084861241*u2_B0(0,-1,0)*rhou1_B0(0,-1,0) -
            0.121937153224065*u2_B0(0,2,0)*rhou1_B0(0,2,0))*invDelta1block0;

       d1_rhou2_dy = (0.00932597985049999*rhou2_B0(0,3,0) + 0.727822147724592*rhou2_B0(0,1,0) +
            0.082033432844602*rhou2_B0(0,-2,0) - 0.0451033223343881*rhou2_B0(0,0,0) - 0.652141084861241*rhou2_B0(0,-1,0)
            - 0.121937153224065*rhou2_B0(0,2,0))*invDelta1block0;

       d1_u0u1_dy = (0.00932597985049999*u0_B0(0,3,0)*u1_B0(0,3,0) + 0.727822147724592*u0_B0(0,1,0)*u1_B0(0,1,0) +
            0.082033432844602*u0_B0(0,-2,0)*u1_B0(0,-2,0) - 0.0451033223343881*u0_B0(0,0,0)*u1_B0(0,0,0) -
            0.652141084861241*u0_B0(0,-1,0)*u1_B0(0,-1,0) -
            0.121937153224065*u0_B0(0,2,0)*u1_B0(0,2,0))*invDelta1block0;

       d1_u1u1_dy = (0.00932597985049999*(u1_B0(0,3,0)*u1_B0(0,3,0)) + 0.727822147724592*(u1_B0(0,1,0)*u1_B0(0,1,0)) +
            0.082033432844602*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) - 0.0451033223343881*(u1_B0(0,0,0)*u1_B0(0,0,0)) -
            0.652141084861241*(u1_B0(0,-1,0)*u1_B0(0,-1,0)) -
            0.121937153224065*(u1_B0(0,2,0)*u1_B0(0,2,0)))*invDelta1block0;

       d1_u1u2_dy = (0.00932597985049999*u1_B0(0,3,0)*u2_B0(0,3,0) + 0.727822147724592*u1_B0(0,1,0)*u2_B0(0,1,0) +
            0.082033432844602*u1_B0(0,-2,0)*u2_B0(0,-2,0) - 0.0451033223343881*u1_B0(0,0,0)*u2_B0(0,0,0) -
            0.652141084861241*u1_B0(0,-1,0)*u2_B0(0,-1,0) -
            0.121937153224065*u1_B0(0,2,0)*u2_B0(0,2,0))*invDelta1block0;

   }

   else{

       d1_H_dy = (-(2.0/3.0)*H_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0) + ((1.0/12.0))*H_B0(0,-2,0) +
            ((2.0/3.0))*H_B0(0,1,0))*invDelta1block0;

       d1_Hrho_dy = (-(2.0/3.0)*H_B0(0,-1,0)*rho_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0)*rho_B0(0,2,0) +
            ((1.0/12.0))*H_B0(0,-2,0)*rho_B0(0,-2,0) + ((2.0/3.0))*H_B0(0,1,0)*rho_B0(0,1,0))*invDelta1block0;

       d1_Hrhou1_dy = (-(2.0/3.0)*H_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0)*rhou1_B0(0,2,0) +
            ((1.0/12.0))*H_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*H_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

       d1_Hu1_dy = (-(2.0/3.0)*H_B0(0,-1,0)*u1_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0)*u1_B0(0,2,0) +
            ((1.0/12.0))*H_B0(0,-2,0)*u1_B0(0,-2,0) + ((2.0/3.0))*H_B0(0,1,0)*u1_B0(0,1,0))*invDelta1block0;

       d1_p_dy = (-(2.0/3.0)*p_B0(0,-1,0) - (1.0/12.0)*p_B0(0,2,0) + ((1.0/12.0))*p_B0(0,-2,0) +
            ((2.0/3.0))*p_B0(0,1,0))*invDelta1block0;

       d1_rho_dy = (-(2.0/3.0)*rho_B0(0,-1,0) - (1.0/12.0)*rho_B0(0,2,0) + ((1.0/12.0))*rho_B0(0,-2,0) +
            ((2.0/3.0))*rho_B0(0,1,0))*invDelta1block0;

       d1_rhou0_dy = (-(2.0/3.0)*rhou0_B0(0,-1,0) - (1.0/12.0)*rhou0_B0(0,2,0) + ((1.0/12.0))*rhou0_B0(0,-2,0) +
            ((2.0/3.0))*rhou0_B0(0,1,0))*invDelta1block0;

       d1_rhou1_dy = (-(2.0/3.0)*rhou1_B0(0,-1,0) - (1.0/12.0)*rhou1_B0(0,2,0) + ((1.0/12.0))*rhou1_B0(0,-2,0) +
            ((2.0/3.0))*rhou1_B0(0,1,0))*invDelta1block0;

       d1_rhou1u0_dy = (-(2.0/3.0)*u0_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0)*rhou1_B0(0,2,0) +
            ((1.0/12.0))*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*u0_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

       d1_rhou1u1_dy = (-(2.0/3.0)*u1_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0)*rhou1_B0(0,2,0) +
            ((1.0/12.0))*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*u1_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

       d1_rhou1u2_dy = (-(2.0/3.0)*u2_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0)*rhou1_B0(0,2,0) +
            ((1.0/12.0))*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*u2_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

       d1_rhou2_dy = (-(2.0/3.0)*rhou2_B0(0,-1,0) - (1.0/12.0)*rhou2_B0(0,2,0) + ((1.0/12.0))*rhou2_B0(0,-2,0) +
            ((2.0/3.0))*rhou2_B0(0,1,0))*invDelta1block0;

       d1_u0u1_dy = (-(2.0/3.0)*u0_B0(0,-1,0)*u1_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0)*u1_B0(0,2,0) +
            ((1.0/12.0))*u0_B0(0,-2,0)*u1_B0(0,-2,0) + ((2.0/3.0))*u0_B0(0,1,0)*u1_B0(0,1,0))*invDelta1block0;

       d1_u1u1_dy = (-(2.0/3.0)*(u1_B0(0,-1,0)*u1_B0(0,-1,0)) - (1.0/12.0)*(u1_B0(0,2,0)*u1_B0(0,2,0)) +
            ((1.0/12.0))*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) + ((2.0/3.0))*(u1_B0(0,1,0)*u1_B0(0,1,0)))*invDelta1block0;

       d1_u1u2_dy = (-(2.0/3.0)*u1_B0(0,-1,0)*u2_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0)*u2_B0(0,2,0) +
            ((1.0/12.0))*u1_B0(0,-2,0)*u2_B0(0,-2,0) + ((2.0/3.0))*u1_B0(0,1,0)*u2_B0(0,1,0))*invDelta1block0;

   }

    d1_H_dx = (-(2.0/3.0)*H_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0) + ((1.0/12.0))*H_B0(-2,0,0) +
      ((2.0/3.0))*H_B0(1,0,0))*invDelta0block0;

    d1_Hrho_dx = (-(2.0/3.0)*H_B0(-1,0,0)*rho_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0)*rho_B0(2,0,0) +
      ((1.0/12.0))*H_B0(-2,0,0)*rho_B0(-2,0,0) + ((2.0/3.0))*H_B0(1,0,0)*rho_B0(1,0,0))*invDelta0block0;

    d1_Hrhou0_dx = (-(2.0/3.0)*H_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*H_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*H_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_Hu0_dx = (-(2.0/3.0)*H_B0(-1,0,0)*u0_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0)*u0_B0(2,0,0) +
      ((1.0/12.0))*H_B0(-2,0,0)*u0_B0(-2,0,0) + ((2.0/3.0))*H_B0(1,0,0)*u0_B0(1,0,0))*invDelta0block0;

    d1_p_dx = (-(2.0/3.0)*p_B0(-1,0,0) - (1.0/12.0)*p_B0(2,0,0) + ((1.0/12.0))*p_B0(-2,0,0) +
      ((2.0/3.0))*p_B0(1,0,0))*invDelta0block0;

    d1_rho_dx = (-(2.0/3.0)*rho_B0(-1,0,0) - (1.0/12.0)*rho_B0(2,0,0) + ((1.0/12.0))*rho_B0(-2,0,0) +
      ((2.0/3.0))*rho_B0(1,0,0))*invDelta0block0;

    d1_rhou0_dx = (-(2.0/3.0)*rhou0_B0(-1,0,0) - (1.0/12.0)*rhou0_B0(2,0,0) + ((1.0/12.0))*rhou0_B0(-2,0,0) +
      ((2.0/3.0))*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou0u0_dx = (-(2.0/3.0)*u0_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*u0_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*u0_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou0u1_dx = (-(2.0/3.0)*u1_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*u1_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*u1_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou0u2_dx = (-(2.0/3.0)*u2_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*u2_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*u2_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou1_dx = (-(2.0/3.0)*rhou1_B0(-1,0,0) - (1.0/12.0)*rhou1_B0(2,0,0) + ((1.0/12.0))*rhou1_B0(-2,0,0) +
      ((2.0/3.0))*rhou1_B0(1,0,0))*invDelta0block0;

    d1_rhou2_dx = (-(2.0/3.0)*rhou2_B0(-1,0,0) - (1.0/12.0)*rhou2_B0(2,0,0) + ((1.0/12.0))*rhou2_B0(-2,0,0) +
      ((2.0/3.0))*rhou2_B0(1,0,0))*invDelta0block0;

    d1_u0u0_dx = (-(2.0/3.0)*(u0_B0(-1,0,0)*u0_B0(-1,0,0)) - (1.0/12.0)*(u0_B0(2,0,0)*u0_B0(2,0,0)) +
      ((1.0/12.0))*(u0_B0(-2,0,0)*u0_B0(-2,0,0)) + ((2.0/3.0))*(u0_B0(1,0,0)*u0_B0(1,0,0)))*invDelta0block0;

    d1_u0u1_dx = (-(2.0/3.0)*u0_B0(-1,0,0)*u1_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0)*u1_B0(2,0,0) +
      ((1.0/12.0))*u0_B0(-2,0,0)*u1_B0(-2,0,0) + ((2.0/3.0))*u0_B0(1,0,0)*u1_B0(1,0,0))*invDelta0block0;

    d1_u0u2_dx = (-(2.0/3.0)*u0_B0(-1,0,0)*u2_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0)*u2_B0(2,0,0) +
      ((1.0/12.0))*u0_B0(-2,0,0)*u2_B0(-2,0,0) + ((2.0/3.0))*u0_B0(1,0,0)*u2_B0(1,0,0))*invDelta0block0;

    d1_H_dz = (-(2.0/3.0)*H_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2) + ((1.0/12.0))*H_B0(0,0,-2) +
      ((2.0/3.0))*H_B0(0,0,1))*invDelta2block0;

    d1_Hrho_dz = (-(2.0/3.0)*H_B0(0,0,-1)*rho_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2)*rho_B0(0,0,2) +
      ((1.0/12.0))*H_B0(0,0,-2)*rho_B0(0,0,-2) + ((2.0/3.0))*H_B0(0,0,1)*rho_B0(0,0,1))*invDelta2block0;

    d1_Hrhou2_dz = (-(2.0/3.0)*H_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*H_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*H_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_Hu2_dz = (-(2.0/3.0)*H_B0(0,0,-1)*u2_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2)*u2_B0(0,0,2) +
      ((1.0/12.0))*H_B0(0,0,-2)*u2_B0(0,0,-2) + ((2.0/3.0))*H_B0(0,0,1)*u2_B0(0,0,1))*invDelta2block0;

    d1_p_dz = (-(2.0/3.0)*p_B0(0,0,-1) - (1.0/12.0)*p_B0(0,0,2) + ((1.0/12.0))*p_B0(0,0,-2) +
      ((2.0/3.0))*p_B0(0,0,1))*invDelta2block0;

    d1_rho_dz = (-(2.0/3.0)*rho_B0(0,0,-1) - (1.0/12.0)*rho_B0(0,0,2) + ((1.0/12.0))*rho_B0(0,0,-2) +
      ((2.0/3.0))*rho_B0(0,0,1))*invDelta2block0;

    d1_rhou0_dz = (-(2.0/3.0)*rhou0_B0(0,0,-1) - (1.0/12.0)*rhou0_B0(0,0,2) + ((1.0/12.0))*rhou0_B0(0,0,-2) +
      ((2.0/3.0))*rhou0_B0(0,0,1))*invDelta2block0;

    d1_rhou1_dz = (-(2.0/3.0)*rhou1_B0(0,0,-1) - (1.0/12.0)*rhou1_B0(0,0,2) + ((1.0/12.0))*rhou1_B0(0,0,-2) +
      ((2.0/3.0))*rhou1_B0(0,0,1))*invDelta2block0;

    d1_rhou2_dz = (-(2.0/3.0)*rhou2_B0(0,0,-1) - (1.0/12.0)*rhou2_B0(0,0,2) + ((1.0/12.0))*rhou2_B0(0,0,-2) +
      ((2.0/3.0))*rhou2_B0(0,0,1))*invDelta2block0;

    d1_rhou2u0_dz = (-(2.0/3.0)*u0_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*u0_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*u0_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_rhou2u1_dz = (-(2.0/3.0)*u1_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*u1_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*u1_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_rhou2u2_dz = (-(2.0/3.0)*u2_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*u2_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*u2_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_u0u2_dz = (-(2.0/3.0)*u0_B0(0,0,-1)*u2_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2)*u2_B0(0,0,2) +
      ((1.0/12.0))*u0_B0(0,0,-2)*u2_B0(0,0,-2) + ((2.0/3.0))*u0_B0(0,0,1)*u2_B0(0,0,1))*invDelta2block0;

    d1_u1u2_dz = (-(2.0/3.0)*u1_B0(0,0,-1)*u2_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2)*u2_B0(0,0,2) +
      ((1.0/12.0))*u1_B0(0,0,-2)*u2_B0(0,0,-2) + ((2.0/3.0))*u1_B0(0,0,1)*u2_B0(0,0,1))*invDelta2block0;

    d1_u2u2_dz = (-(2.0/3.0)*(u2_B0(0,0,-1)*u2_B0(0,0,-1)) - (1.0/12.0)*(u2_B0(0,0,2)*u2_B0(0,0,2)) +
      ((1.0/12.0))*(u2_B0(0,0,-2)*u2_B0(0,0,-2)) + ((2.0/3.0))*(u2_B0(0,0,1)*u2_B0(0,0,1)))*invDelta2block0;

    Residual0_B0(0,0,0) = -(1.0/2.0)*d1_rhou0_dx - (1.0/2.0)*d1_rhou2_dz - (1.0/2.0)*(D11_B0(0,0,0)*wk5_B0(0,0,0) +
      wk0_B0(0,0,0) + wk10_B0(0,0,0))*rho_B0(0,0,0) - (1.0/2.0)*u0_B0(0,0,0)*d1_rho_dx -
      (1.0/2.0)*u2_B0(0,0,0)*d1_rho_dz - (1.0/2.0)*D11_B0(0,0,0)*d1_rhou1_dy -
      (1.0/2.0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_rho_dy;

    Residual1_B0(0,0,0) = -d1_p_dx - (1.0/4.0)*d1_rhou0u0_dx - (1.0/4.0)*d1_rhou2u0_dz - c0*phi_B0(0,0,0) -
      (1.0/2.0)*u0_B0(0,0,0)*d1_rhou0_dx - (1.0/2.0)*wk0_B0(0,0,0)*rhou0_B0(0,0,0) -
      (1.0/4.0)*(u0_B0(0,0,0)*u0_B0(0,0,0))*d1_rho_dx - (1.0/4.0)*u0_B0(0,0,0)*d1_rhou2_dz -
      (1.0/4.0)*u2_B0(0,0,0)*d1_rhou0_dz - (1.0/4.0)*D11_B0(0,0,0)*d1_rhou1u0_dy - (1.0/4.0)*rho_B0(0,0,0)*d1_u0u0_dx -
      (1.0/4.0)*rho_B0(0,0,0)*d1_u0u2_dz - (1.0/4.0)*wk8_B0(0,0,0)*rhou2_B0(0,0,0) -
      (1.0/4.0)*wk10_B0(0,0,0)*rhou0_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dz -
      (1.0/4.0)*u0_B0(0,0,0)*D11_B0(0,0,0)*d1_rhou1_dy - (1.0/4.0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_rhou0_dy -
      (1.0/4.0)*D11_B0(0,0,0)*rho_B0(0,0,0)*d1_u0u1_dy - (1.0/4.0)*D11_B0(0,0,0)*wk4_B0(0,0,0)*rhou1_B0(0,0,0) -
      (1.0/4.0)*D11_B0(0,0,0)*wk5_B0(0,0,0)*rhou0_B0(0,0,0) -
      (1.0/4.0)*u0_B0(0,0,0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_rho_dy;

    Residual2_B0(0,0,0) = -(1.0/4.0)*d1_rhou0u1_dx - (1.0/4.0)*d1_rhou2u1_dz - c1*phi_B0(0,0,0) - D11_B0(0,0,0)*d1_p_dy
      - (1.0/4.0)*(D11_B0(0,0,0)*wk5_B0(0,0,0) + wk0_B0(0,0,0) + wk10_B0(0,0,0))*rhou1_B0(0,0,0) -
      (1.0/4.0)*(D11_B0(0,0,0)*d1_rhou1_dy + d1_rhou0_dx + d1_rhou2_dz)*u1_B0(0,0,0) -
      (1.0/4.0)*(D11_B0(0,0,0)*d1_u1u1_dy + d1_u0u1_dx + d1_u1u2_dz)*rho_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*d1_rhou1_dx
      - (1.0/4.0)*u2_B0(0,0,0)*d1_rhou1_dz - (1.0/4.0)*D11_B0(0,0,0)*d1_rhou1u1_dy -
      (1.0/4.0)*wk1_B0(0,0,0)*rhou0_B0(0,0,0) - (1.0/4.0)*wk9_B0(0,0,0)*rhou2_B0(0,0,0) -
      (1.0/4.0)*(u1_B0(0,0,0)*u1_B0(0,0,0))*D11_B0(0,0,0)*d1_rho_dy - (1.0/4.0)*u0_B0(0,0,0)*u1_B0(0,0,0)*d1_rho_dx -
      (1.0/4.0)*u1_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dz - (1.0/4.0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_rhou1_dy -
      (1.0/4.0)*D11_B0(0,0,0)*wk5_B0(0,0,0)*rhou1_B0(0,0,0);

    Residual3_B0(0,0,0) = -d1_p_dz - (1.0/4.0)*d1_rhou0u2_dx - (1.0/4.0)*d1_rhou2u2_dz - c2*phi_B0(0,0,0) -
      (1.0/4.0)*(u2_B0(0,0,0)*u2_B0(0,0,0))*d1_rho_dz - (1.0/4.0)*(D11_B0(0,0,0)*wk5_B0(0,0,0) + wk0_B0(0,0,0) +
      wk10_B0(0,0,0))*rhou2_B0(0,0,0) - (1.0/4.0)*(D11_B0(0,0,0)*d1_rhou1_dy + d1_rhou0_dx + d1_rhou2_dz)*u2_B0(0,0,0) -
      (1.0/4.0)*(D11_B0(0,0,0)*d1_u1u2_dy + d1_u0u2_dx + d1_u2u2_dz)*rho_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*d1_rhou2_dx
      - (1.0/4.0)*u2_B0(0,0,0)*d1_rhou2_dz - (1.0/4.0)*D11_B0(0,0,0)*d1_rhou1u2_dy -
      (1.0/4.0)*wk2_B0(0,0,0)*rhou0_B0(0,0,0) - (1.0/4.0)*wk10_B0(0,0,0)*rhou2_B0(0,0,0) -
      (1.0/4.0)*u0_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dx - (1.0/4.0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_rhou2_dy -
      (1.0/4.0)*D11_B0(0,0,0)*wk6_B0(0,0,0)*rhou1_B0(0,0,0) -
      (1.0/4.0)*u1_B0(0,0,0)*u2_B0(0,0,0)*D11_B0(0,0,0)*d1_rho_dy;

    Residual4_B0(0,0,0) = -(1.0/4.0)*d1_Hrhou0_dx - (1.0/4.0)*d1_Hrhou2_dz - (1.0/4.0)*(D11_B0(0,0,0)*d1_Hu1_dy +
      d1_Hu0_dx + d1_Hu2_dz)*rho_B0(0,0,0) - (1.0/4.0)*(D11_B0(0,0,0)*d1_rhou1_dy + d1_rhou0_dx +
      d1_rhou2_dz)*H_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*d1_Hrho_dx - (1.0/4.0)*u2_B0(0,0,0)*d1_Hrho_dz -
      (1.0/4.0)*D11_B0(0,0,0)*d1_Hrhou1_dy - (1.0/4.0)*rhou0_B0(0,0,0)*d1_H_dx - (1.0/4.0)*rhou2_B0(0,0,0)*d1_H_dz -
      c0*u0_B0(0,0,0)*phi_B0(0,0,0) - c1*u1_B0(0,0,0)*phi_B0(0,0,0) - c2*u2_B0(0,0,0)*phi_B0(0,0,0) -
      (1.0/4.0)*(D11_B0(0,0,0)*wk5_B0(0,0,0) + wk0_B0(0,0,0) + wk10_B0(0,0,0))*H_B0(0,0,0)*rho_B0(0,0,0) -
      (1.0/4.0)*H_B0(0,0,0)*u0_B0(0,0,0)*d1_rho_dx - (1.0/4.0)*H_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dz -
      (1.0/4.0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_Hrho_dy - (1.0/4.0)*D11_B0(0,0,0)*rhou1_B0(0,0,0)*d1_H_dy -
      (1.0/4.0)*H_B0(0,0,0)*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_rho_dy;

}

 void opensbliblock00Kernel033(const ACC<double> &D11_B0, const ACC<double> &SD111_B0, const ACC<double> &T_B0, const
ACC<double> &mu_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, const ACC<double>
&wk0_B0, const ACC<double> &wk10_B0, const ACC<double> &wk11_B0, const ACC<double> &wk1_B0, const ACC<double> &wk2_B0,
const ACC<double> &wk3_B0, const ACC<double> &wk4_B0, const ACC<double> &wk5_B0, const ACC<double> &wk6_B0, const
ACC<double> &wk7_B0, const ACC<double> &wk8_B0, const ACC<double> &wk9_B0, ACC<double> &Residual1_B0, ACC<double>
&Residual2_B0, ACC<double> &Residual3_B0, ACC<double> &Residual4_B0, const int *idx)
{
   double d1_mu_dx = 0.0;
   double d1_mu_dy = 0.0;
   double d1_mu_dz = 0.0;
   double d1_wk0_dy = 0.0;
   double d1_wk0_dz = 0.0;
   double d1_wk1_dy = 0.0;
   double d1_wk2_dz = 0.0;
   double d1_wk5_dz = 0.0;
   double d1_wk6_dz = 0.0;
   double d2_T_dx = 0.0;
   double d2_T_dy = 0.0;
   double d2_T_dz = 0.0;
   double d2_u0_dx = 0.0;
   double d2_u0_dy = 0.0;
   double d2_u0_dz = 0.0;
   double d2_u1_dx = 0.0;
   double d2_u1_dy = 0.0;
   double d2_u1_dz = 0.0;
   double d2_u2_dx = 0.0;
   double d2_u2_dy = 0.0;
   double d2_u2_dz = 0.0;
   if (idx[1] == 0){

       d1_mu_dy = (3.0*mu_B0(0,1,0) + 0.333333333333333*mu_B0(0,3,0) - 1.5*mu_B0(0,2,0) -
            1.83333333333333*mu_B0(0,0,0))*invDelta1block0;

       d1_wk0_dy = (3.0*wk0_B0(0,1,0) + 0.333333333333333*wk0_B0(0,3,0) - 1.5*wk0_B0(0,2,0) -
            1.83333333333333*wk0_B0(0,0,0))*invDelta1block0;

       d1_wk1_dy = (3.0*wk1_B0(0,1,0) + 0.333333333333333*wk1_B0(0,3,0) - 1.5*wk1_B0(0,2,0) -
            1.83333333333333*wk1_B0(0,0,0))*invDelta1block0;

   }

   else if (idx[1] == 1){

       d1_mu_dy = (0.0394168524399447*mu_B0(0,2,0) + 0.00571369039775442*mu_B0(0,4,0) + 0.719443173328855*mu_B0(0,1,0) -
            0.322484932882161*mu_B0(0,0,0) - 0.0658051057710389*mu_B0(0,3,0) -
            0.376283677513354*mu_B0(0,-1,0))*invDelta1block0;

       d1_wk0_dy = (0.0394168524399447*wk0_B0(0,2,0) + 0.00571369039775442*wk0_B0(0,4,0) +
            0.719443173328855*wk0_B0(0,1,0) - 0.322484932882161*wk0_B0(0,0,0) - 0.0658051057710389*wk0_B0(0,3,0) -
            0.376283677513354*wk0_B0(0,-1,0))*invDelta1block0;

       d1_wk1_dy = (0.0394168524399447*wk1_B0(0,2,0) + 0.00571369039775442*wk1_B0(0,4,0) +
            0.719443173328855*wk1_B0(0,1,0) - 0.322484932882161*wk1_B0(0,0,0) - 0.0658051057710389*wk1_B0(0,3,0) -
            0.376283677513354*wk1_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 2){

       d1_mu_dy = (0.197184333887745*mu_B0(0,0,0) + 0.521455851089587*mu_B0(0,1,0) + 0.113446470384241*mu_B0(0,-2,0) -
            0.00412637789557492*mu_B0(0,3,0) - 0.0367146847001261*mu_B0(0,2,0) -
            0.791245592765872*mu_B0(0,-1,0))*invDelta1block0;

       d1_wk0_dy = (0.197184333887745*wk0_B0(0,0,0) + 0.521455851089587*wk0_B0(0,1,0) + 0.113446470384241*wk0_B0(0,-2,0)
            - 0.00412637789557492*wk0_B0(0,3,0) - 0.0367146847001261*wk0_B0(0,2,0) -
            0.791245592765872*wk0_B0(0,-1,0))*invDelta1block0;

       d1_wk1_dy = (0.197184333887745*wk1_B0(0,0,0) + 0.521455851089587*wk1_B0(0,1,0) + 0.113446470384241*wk1_B0(0,-2,0)
            - 0.00412637789557492*wk1_B0(0,3,0) - 0.0367146847001261*wk1_B0(0,2,0) -
            0.791245592765872*wk1_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == 3){

       d1_mu_dy = (0.0451033223343881*mu_B0(0,0,0) + 0.652141084861241*mu_B0(0,1,0) + 0.121937153224065*mu_B0(0,-2,0) -
            0.00932597985049999*mu_B0(0,-3,0) - 0.727822147724592*mu_B0(0,-1,0) -
            0.082033432844602*mu_B0(0,2,0))*invDelta1block0;

       d1_wk0_dy = (0.0451033223343881*wk0_B0(0,0,0) + 0.652141084861241*wk0_B0(0,1,0) +
            0.121937153224065*wk0_B0(0,-2,0) - 0.00932597985049999*wk0_B0(0,-3,0) - 0.727822147724592*wk0_B0(0,-1,0) -
            0.082033432844602*wk0_B0(0,2,0))*invDelta1block0;

       d1_wk1_dy = (0.0451033223343881*wk1_B0(0,0,0) + 0.652141084861241*wk1_B0(0,1,0) +
            0.121937153224065*wk1_B0(0,-2,0) - 0.00932597985049999*wk1_B0(0,-3,0) - 0.727822147724592*wk1_B0(0,-1,0) -
            0.082033432844602*wk1_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       d1_mu_dy = (1.5*mu_B0(0,-2,0) + 1.83333333333333*mu_B0(0,0,0) - 3.0*mu_B0(0,-1,0) -
            0.333333333333333*mu_B0(0,-3,0))*invDelta1block0;

       d1_wk0_dy = (1.5*wk0_B0(0,-2,0) + 1.83333333333333*wk0_B0(0,0,0) - 3.0*wk0_B0(0,-1,0) -
            0.333333333333333*wk0_B0(0,-3,0))*invDelta1block0;

       d1_wk1_dy = (1.5*wk1_B0(0,-2,0) + 1.83333333333333*wk1_B0(0,0,0) - 3.0*wk1_B0(0,-1,0) -
            0.333333333333333*wk1_B0(0,-3,0))*invDelta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       d1_mu_dy = (0.322484932882161*mu_B0(0,0,0) + 0.0658051057710389*mu_B0(0,-3,0) + 0.376283677513354*mu_B0(0,1,0) -
            0.0394168524399447*mu_B0(0,-2,0) - 0.00571369039775442*mu_B0(0,-4,0) -
            0.719443173328855*mu_B0(0,-1,0))*invDelta1block0;

       d1_wk0_dy = (0.322484932882161*wk0_B0(0,0,0) + 0.0658051057710389*wk0_B0(0,-3,0) +
            0.376283677513354*wk0_B0(0,1,0) - 0.0394168524399447*wk0_B0(0,-2,0) - 0.00571369039775442*wk0_B0(0,-4,0) -
            0.719443173328855*wk0_B0(0,-1,0))*invDelta1block0;

       d1_wk1_dy = (0.322484932882161*wk1_B0(0,0,0) + 0.0658051057710389*wk1_B0(0,-3,0) +
            0.376283677513354*wk1_B0(0,1,0) - 0.0394168524399447*wk1_B0(0,-2,0) - 0.00571369039775442*wk1_B0(0,-4,0) -
            0.719443173328855*wk1_B0(0,-1,0))*invDelta1block0;

   }

   else if (idx[1] == -3 + block0np1){

       d1_mu_dy = (0.00412637789557492*mu_B0(0,-3,0) + 0.0367146847001261*mu_B0(0,-2,0) + 0.791245592765872*mu_B0(0,1,0)
            - 0.197184333887745*mu_B0(0,0,0) - 0.521455851089587*mu_B0(0,-1,0) -
            0.113446470384241*mu_B0(0,2,0))*invDelta1block0;

       d1_wk0_dy = (0.00412637789557492*wk0_B0(0,-3,0) + 0.0367146847001261*wk0_B0(0,-2,0) +
            0.791245592765872*wk0_B0(0,1,0) - 0.197184333887745*wk0_B0(0,0,0) - 0.521455851089587*wk0_B0(0,-1,0) -
            0.113446470384241*wk0_B0(0,2,0))*invDelta1block0;

       d1_wk1_dy = (0.00412637789557492*wk1_B0(0,-3,0) + 0.0367146847001261*wk1_B0(0,-2,0) +
            0.791245592765872*wk1_B0(0,1,0) - 0.197184333887745*wk1_B0(0,0,0) - 0.521455851089587*wk1_B0(0,-1,0) -
            0.113446470384241*wk1_B0(0,2,0))*invDelta1block0;

   }

   else if (idx[1] == -4 + block0np1){

       d1_mu_dy = (0.00932597985049999*mu_B0(0,3,0) + 0.727822147724592*mu_B0(0,1,0) + 0.082033432844602*mu_B0(0,-2,0) -
            0.0451033223343881*mu_B0(0,0,0) - 0.652141084861241*mu_B0(0,-1,0) -
            0.121937153224065*mu_B0(0,2,0))*invDelta1block0;

       d1_wk0_dy = (0.00932597985049999*wk0_B0(0,3,0) + 0.727822147724592*wk0_B0(0,1,0) +
            0.082033432844602*wk0_B0(0,-2,0) - 0.0451033223343881*wk0_B0(0,0,0) - 0.652141084861241*wk0_B0(0,-1,0) -
            0.121937153224065*wk0_B0(0,2,0))*invDelta1block0;

       d1_wk1_dy = (0.00932597985049999*wk1_B0(0,3,0) + 0.727822147724592*wk1_B0(0,1,0) +
            0.082033432844602*wk1_B0(0,-2,0) - 0.0451033223343881*wk1_B0(0,0,0) - 0.652141084861241*wk1_B0(0,-1,0) -
            0.121937153224065*wk1_B0(0,2,0))*invDelta1block0;

   }

   else{

       d1_mu_dy = (-(2.0/3.0)*mu_B0(0,-1,0) - (1.0/12.0)*mu_B0(0,2,0) + ((1.0/12.0))*mu_B0(0,-2,0) +
            ((2.0/3.0))*mu_B0(0,1,0))*invDelta1block0;

       d1_wk0_dy = (-(2.0/3.0)*wk0_B0(0,-1,0) - (1.0/12.0)*wk0_B0(0,2,0) + ((1.0/12.0))*wk0_B0(0,-2,0) +
            ((2.0/3.0))*wk0_B0(0,1,0))*invDelta1block0;

       d1_wk1_dy = (-(2.0/3.0)*wk1_B0(0,-1,0) - (1.0/12.0)*wk1_B0(0,2,0) + ((1.0/12.0))*wk1_B0(0,-2,0) +
            ((2.0/3.0))*wk1_B0(0,1,0))*invDelta1block0;

   }

   if (idx[1] == 0){

       d2_T_dy = (-(26.0/3.0)*T_B0(0,1,0) - (14.0/3.0)*T_B0(0,3,0) + ((11.0/12.0))*T_B0(0,4,0) +
            ((19.0/2.0))*T_B0(0,2,0) + ((35.0/12.0))*T_B0(0,0,0))*inv2Delta1block0;

       d2_u0_dy = (-(26.0/3.0)*u0_B0(0,1,0) - (14.0/3.0)*u0_B0(0,3,0) + ((11.0/12.0))*u0_B0(0,4,0) +
            ((19.0/2.0))*u0_B0(0,2,0) + ((35.0/12.0))*u0_B0(0,0,0))*inv2Delta1block0;

       d2_u1_dy = (-(26.0/3.0)*u1_B0(0,1,0) - (14.0/3.0)*u1_B0(0,3,0) + ((11.0/12.0))*u1_B0(0,4,0) +
            ((19.0/2.0))*u1_B0(0,2,0) + ((35.0/12.0))*u1_B0(0,0,0))*inv2Delta1block0;

       d2_u2_dy = (-(26.0/3.0)*u2_B0(0,1,0) - (14.0/3.0)*u2_B0(0,3,0) + ((11.0/12.0))*u2_B0(0,4,0) +
            ((19.0/2.0))*u2_B0(0,2,0) + ((35.0/12.0))*u2_B0(0,0,0))*inv2Delta1block0;

   }

   else if (idx[1] == 1){

       d2_T_dy = (((1.0/2.0))*T_B0(0,1,0) - (5.0/3.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(0,3,0) + ((1.0/3.0))*T_B0(0,2,0) +
            ((11.0/12.0))*T_B0(0,-1,0))*inv2Delta1block0;

       d2_u0_dy = (((1.0/2.0))*u0_B0(0,1,0) - (5.0/3.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(0,3,0) +
            ((1.0/3.0))*u0_B0(0,2,0) + ((11.0/12.0))*u0_B0(0,-1,0))*inv2Delta1block0;

       d2_u1_dy = (((1.0/2.0))*u1_B0(0,1,0) - (5.0/3.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(0,3,0) +
            ((1.0/3.0))*u1_B0(0,2,0) + ((11.0/12.0))*u1_B0(0,-1,0))*inv2Delta1block0;

       d2_u2_dy = (((1.0/2.0))*u2_B0(0,1,0) - (5.0/3.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(0,3,0) +
            ((1.0/3.0))*u2_B0(0,2,0) + ((11.0/12.0))*u2_B0(0,-1,0))*inv2Delta1block0;

   }

   else if (idx[1] == -1 + block0np1){

       d2_T_dy = (-(26.0/3.0)*T_B0(0,-1,0) - (14.0/3.0)*T_B0(0,-3,0) + ((11.0/12.0))*T_B0(0,-4,0) +
            ((19.0/2.0))*T_B0(0,-2,0) + ((35.0/12.0))*T_B0(0,0,0))*inv2Delta1block0;

       d2_u0_dy = (-(26.0/3.0)*u0_B0(0,-1,0) - (14.0/3.0)*u0_B0(0,-3,0) + ((11.0/12.0))*u0_B0(0,-4,0) +
            ((19.0/2.0))*u0_B0(0,-2,0) + ((35.0/12.0))*u0_B0(0,0,0))*inv2Delta1block0;

       d2_u1_dy = (-(26.0/3.0)*u1_B0(0,-1,0) - (14.0/3.0)*u1_B0(0,-3,0) + ((11.0/12.0))*u1_B0(0,-4,0) +
            ((19.0/2.0))*u1_B0(0,-2,0) + ((35.0/12.0))*u1_B0(0,0,0))*inv2Delta1block0;

       d2_u2_dy = (-(26.0/3.0)*u2_B0(0,-1,0) - (14.0/3.0)*u2_B0(0,-3,0) + ((11.0/12.0))*u2_B0(0,-4,0) +
            ((19.0/2.0))*u2_B0(0,-2,0) + ((35.0/12.0))*u2_B0(0,0,0))*inv2Delta1block0;

   }

   else if (idx[1] == -2 + block0np1){

       d2_T_dy = (((1.0/2.0))*T_B0(0,-1,0) - (5.0/3.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(0,-3,0) + ((1.0/3.0))*T_B0(0,-2,0)
            + ((11.0/12.0))*T_B0(0,1,0))*inv2Delta1block0;

       d2_u0_dy = (((1.0/2.0))*u0_B0(0,-1,0) - (5.0/3.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(0,-3,0) +
            ((1.0/3.0))*u0_B0(0,-2,0) + ((11.0/12.0))*u0_B0(0,1,0))*inv2Delta1block0;

       d2_u1_dy = (((1.0/2.0))*u1_B0(0,-1,0) - (5.0/3.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(0,-3,0) +
            ((1.0/3.0))*u1_B0(0,-2,0) + ((11.0/12.0))*u1_B0(0,1,0))*inv2Delta1block0;

       d2_u2_dy = (((1.0/2.0))*u2_B0(0,-1,0) - (5.0/3.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(0,-3,0) +
            ((1.0/3.0))*u2_B0(0,-2,0) + ((11.0/12.0))*u2_B0(0,1,0))*inv2Delta1block0;

   }

   else{

       d2_T_dy = (-(5.0/2.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(0,-2,0) - (1.0/12.0)*T_B0(0,2,0) + ((4.0/3.0))*T_B0(0,1,0) +
            ((4.0/3.0))*T_B0(0,-1,0))*inv2Delta1block0;

       d2_u0_dy = (-(5.0/2.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(0,-2,0) - (1.0/12.0)*u0_B0(0,2,0) +
            ((4.0/3.0))*u0_B0(0,1,0) + ((4.0/3.0))*u0_B0(0,-1,0))*inv2Delta1block0;

       d2_u1_dy = (-(5.0/2.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(0,-2,0) - (1.0/12.0)*u1_B0(0,2,0) +
            ((4.0/3.0))*u1_B0(0,1,0) + ((4.0/3.0))*u1_B0(0,-1,0))*inv2Delta1block0;

       d2_u2_dy = (-(5.0/2.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(0,-2,0) - (1.0/12.0)*u2_B0(0,2,0) +
            ((4.0/3.0))*u2_B0(0,1,0) + ((4.0/3.0))*u2_B0(0,-1,0))*inv2Delta1block0;

   }

    d2_T_dx = (-(5.0/2.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(-2,0,0) - (1.0/12.0)*T_B0(2,0,0) + ((4.0/3.0))*T_B0(1,0,0) +
      ((4.0/3.0))*T_B0(-1,0,0))*inv2Delta0block0;

    d1_mu_dx = (-(2.0/3.0)*mu_B0(-1,0,0) - (1.0/12.0)*mu_B0(2,0,0) + ((1.0/12.0))*mu_B0(-2,0,0) +
      ((2.0/3.0))*mu_B0(1,0,0))*invDelta0block0;

    d2_u0_dx = (-(5.0/2.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(-2,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((4.0/3.0))*u0_B0(1,0,0)
      + ((4.0/3.0))*u0_B0(-1,0,0))*inv2Delta0block0;

    d2_u1_dx = (-(5.0/2.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(-2,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((4.0/3.0))*u1_B0(1,0,0)
      + ((4.0/3.0))*u1_B0(-1,0,0))*inv2Delta0block0;

    d2_u2_dx = (-(5.0/2.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(-2,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((4.0/3.0))*u2_B0(1,0,0)
      + ((4.0/3.0))*u2_B0(-1,0,0))*inv2Delta0block0;

    d2_T_dz = (-(5.0/2.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(0,0,-2) - (1.0/12.0)*T_B0(0,0,2) + ((4.0/3.0))*T_B0(0,0,1) +
      ((4.0/3.0))*T_B0(0,0,-1))*inv2Delta2block0;

    d1_mu_dz = (-(2.0/3.0)*mu_B0(0,0,-1) - (1.0/12.0)*mu_B0(0,0,2) + ((1.0/12.0))*mu_B0(0,0,-2) +
      ((2.0/3.0))*mu_B0(0,0,1))*invDelta2block0;

    d2_u0_dz = (-(5.0/2.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(0,0,-2) - (1.0/12.0)*u0_B0(0,0,2) + ((4.0/3.0))*u0_B0(0,0,1)
      + ((4.0/3.0))*u0_B0(0,0,-1))*inv2Delta2block0;

    d2_u1_dz = (-(5.0/2.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(0,0,-2) - (1.0/12.0)*u1_B0(0,0,2) + ((4.0/3.0))*u1_B0(0,0,1)
      + ((4.0/3.0))*u1_B0(0,0,-1))*inv2Delta2block0;

    d2_u2_dz = (-(5.0/2.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(0,0,-2) - (1.0/12.0)*u2_B0(0,0,2) + ((4.0/3.0))*u2_B0(0,0,1)
      + ((4.0/3.0))*u2_B0(0,0,-1))*inv2Delta2block0;

    d1_wk0_dz = (-(2.0/3.0)*wk0_B0(0,0,-1) - (1.0/12.0)*wk0_B0(0,0,2) + ((1.0/12.0))*wk0_B0(0,0,-2) +
      ((2.0/3.0))*wk0_B0(0,0,1))*invDelta2block0;

    d1_wk2_dz = (-(2.0/3.0)*wk2_B0(0,0,-1) - (1.0/12.0)*wk2_B0(0,0,2) + ((1.0/12.0))*wk2_B0(0,0,-2) +
      ((2.0/3.0))*wk2_B0(0,0,1))*invDelta2block0;

    d1_wk5_dz = (-(2.0/3.0)*wk5_B0(0,0,-1) - (1.0/12.0)*wk5_B0(0,0,2) + ((1.0/12.0))*wk5_B0(0,0,-2) +
      ((2.0/3.0))*wk5_B0(0,0,1))*invDelta2block0;

    d1_wk6_dz = (-(2.0/3.0)*wk6_B0(0,0,-1) - (1.0/12.0)*wk6_B0(0,0,2) + ((1.0/12.0))*wk6_B0(0,0,-2) +
      ((2.0/3.0))*wk6_B0(0,0,1))*invDelta2block0;

    Residual1_B0(0,0,0) = (wk2_B0(0,0,0) + wk8_B0(0,0,0))*invRe*d1_mu_dz + (-(2.0/3.0)*wk10_B0(0,0,0) +
      ((4.0/3.0))*wk0_B0(0,0,0) - (2.0/3.0)*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*d1_mu_dx + (((1.0/3.0))*d1_wk2_dz +
      ((4.0/3.0))*d2_u0_dx + (D11_B0(0,0,0)*D11_B0(0,0,0))*d2_u0_dy + ((1.0/3.0))*D11_B0(0,0,0)*d1_wk1_dy +
      D11_B0(0,0,0)*wk4_B0(0,0,0)*SD111_B0(0,0,0) + d2_u0_dz)*invRe*mu_B0(0,0,0) + (D11_B0(0,0,0)*wk4_B0(0,0,0) +
      wk1_B0(0,0,0))*invRe*D11_B0(0,0,0)*d1_mu_dy + Residual1_B0(0,0,0);

    Residual2_B0(0,0,0) = (D11_B0(0,0,0)*wk4_B0(0,0,0) + wk1_B0(0,0,0))*invRe*d1_mu_dx + (D11_B0(0,0,0)*wk6_B0(0,0,0) +
      wk9_B0(0,0,0))*invRe*d1_mu_dz + (((1.0/3.0))*D11_B0(0,0,0)*d1_wk0_dy + ((1.0/3.0))*D11_B0(0,0,0)*d1_wk6_dz +
      ((4.0/3.0))*(D11_B0(0,0,0)*D11_B0(0,0,0))*d2_u1_dy + ((4.0/3.0))*D11_B0(0,0,0)*wk5_B0(0,0,0)*SD111_B0(0,0,0) +
      d2_u1_dx + d2_u1_dz)*invRe*mu_B0(0,0,0) + (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk10_B0(0,0,0) +
      ((4.0/3.0))*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*D11_B0(0,0,0)*d1_mu_dy + Residual2_B0(0,0,0);

    Residual3_B0(0,0,0) = (wk2_B0(0,0,0) + wk8_B0(0,0,0))*invRe*d1_mu_dx + (-(2.0/3.0)*wk0_B0(0,0,0) +
      ((4.0/3.0))*wk10_B0(0,0,0) - (2.0/3.0)*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*d1_mu_dz + (((1.0/3.0))*d1_wk0_dz +
      ((4.0/3.0))*d2_u2_dz + (D11_B0(0,0,0)*D11_B0(0,0,0))*d2_u2_dy + ((1.0/3.0))*D11_B0(0,0,0)*d1_wk5_dz +
      D11_B0(0,0,0)*wk6_B0(0,0,0)*SD111_B0(0,0,0) + d2_u2_dx)*invRe*mu_B0(0,0,0) + (D11_B0(0,0,0)*wk6_B0(0,0,0) +
      wk9_B0(0,0,0))*invRe*D11_B0(0,0,0)*d1_mu_dy + Residual3_B0(0,0,0);

    Residual4_B0(0,0,0) = (D11_B0(0,0,0)*wk4_B0(0,0,0) + wk1_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk1_B0(0,0,0) +
      (D11_B0(0,0,0)*wk4_B0(0,0,0) + wk1_B0(0,0,0))*invRe*u1_B0(0,0,0)*d1_mu_dx + (D11_B0(0,0,0)*wk6_B0(0,0,0) +
      wk9_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk9_B0(0,0,0) + (D11_B0(0,0,0)*wk6_B0(0,0,0) +
      wk9_B0(0,0,0))*invRe*u1_B0(0,0,0)*d1_mu_dz + (wk2_B0(0,0,0) + wk8_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk2_B0(0,0,0) +
      (wk2_B0(0,0,0) + wk8_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk8_B0(0,0,0) + (wk2_B0(0,0,0) +
      wk8_B0(0,0,0))*invRe*u0_B0(0,0,0)*d1_mu_dz + (wk2_B0(0,0,0) + wk8_B0(0,0,0))*invRe*u2_B0(0,0,0)*d1_mu_dx +
      (-(2.0/3.0)*wk0_B0(0,0,0) + ((4.0/3.0))*wk10_B0(0,0,0) -
      (2.0/3.0)*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk10_B0(0,0,0) + (-(2.0/3.0)*wk0_B0(0,0,0) +
      ((4.0/3.0))*wk10_B0(0,0,0) - (2.0/3.0)*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*u2_B0(0,0,0)*d1_mu_dz +
      (-(2.0/3.0)*wk10_B0(0,0,0) + ((4.0/3.0))*wk0_B0(0,0,0) -
      (2.0/3.0)*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk0_B0(0,0,0) + (-(2.0/3.0)*wk10_B0(0,0,0) +
      ((4.0/3.0))*wk0_B0(0,0,0) - (2.0/3.0)*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*u0_B0(0,0,0)*d1_mu_dx +
      (((1.0/3.0))*d1_wk0_dz + ((4.0/3.0))*d2_u2_dz + (D11_B0(0,0,0)*D11_B0(0,0,0))*d2_u2_dy +
      ((1.0/3.0))*D11_B0(0,0,0)*d1_wk5_dz + D11_B0(0,0,0)*wk6_B0(0,0,0)*SD111_B0(0,0,0) +
      d2_u2_dx)*invRe*mu_B0(0,0,0)*u2_B0(0,0,0) + (((1.0/3.0))*d1_wk2_dz + ((4.0/3.0))*d2_u0_dx +
      (D11_B0(0,0,0)*D11_B0(0,0,0))*d2_u0_dy + ((1.0/3.0))*D11_B0(0,0,0)*d1_wk1_dy +
      D11_B0(0,0,0)*wk4_B0(0,0,0)*SD111_B0(0,0,0) + d2_u0_dz)*invRe*mu_B0(0,0,0)*u0_B0(0,0,0) +
      (((1.0/3.0))*D11_B0(0,0,0)*d1_wk0_dy + ((1.0/3.0))*D11_B0(0,0,0)*d1_wk6_dz +
      ((4.0/3.0))*(D11_B0(0,0,0)*D11_B0(0,0,0))*d2_u1_dy + ((4.0/3.0))*D11_B0(0,0,0)*wk5_B0(0,0,0)*SD111_B0(0,0,0) +
      d2_u1_dx + d2_u1_dz)*invRe*mu_B0(0,0,0)*u1_B0(0,0,0) + (D11_B0(0,0,0)*wk4_B0(0,0,0) +
      wk1_B0(0,0,0))*invRe*mu_B0(0,0,0)*D11_B0(0,0,0)*wk4_B0(0,0,0) + (D11_B0(0,0,0)*wk4_B0(0,0,0) +
      wk1_B0(0,0,0))*invRe*u0_B0(0,0,0)*D11_B0(0,0,0)*d1_mu_dy + (D11_B0(0,0,0)*wk6_B0(0,0,0) +
      wk9_B0(0,0,0))*invRe*mu_B0(0,0,0)*D11_B0(0,0,0)*wk6_B0(0,0,0) + (D11_B0(0,0,0)*wk6_B0(0,0,0) +
      wk9_B0(0,0,0))*invRe*u2_B0(0,0,0)*D11_B0(0,0,0)*d1_mu_dy + (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk10_B0(0,0,0) +
      ((4.0/3.0))*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*mu_B0(0,0,0)*D11_B0(0,0,0)*wk5_B0(0,0,0) +
      (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk10_B0(0,0,0) +
      ((4.0/3.0))*D11_B0(0,0,0)*wk5_B0(0,0,0))*invRe*u1_B0(0,0,0)*D11_B0(0,0,0)*d1_mu_dy +
      ((D11_B0(0,0,0)*D11_B0(0,0,0))*d2_T_dy + D11_B0(0,0,0)*wk7_B0(0,0,0)*SD111_B0(0,0,0) + d2_T_dx +
      d2_T_dz)*invPr*invRe*inv2Minf*inv_gamma_m1*mu_B0(0,0,0) + invPr*invRe*inv2Minf*inv_gamma_m1*wk3_B0(0,0,0)*d1_mu_dx
      + invPr*invRe*inv2Minf*inv_gamma_m1*wk11_B0(0,0,0)*d1_mu_dz +
      (D11_B0(0,0,0)*D11_B0(0,0,0))*invPr*invRe*inv2Minf*inv_gamma_m1*wk7_B0(0,0,0)*d1_mu_dy + Residual4_B0(0,0,0);

}

 void opensbliblock00Kernel071(const ACC<double> &Residual0_B0, const ACC<double> &Residual1_B0, const ACC<double>
&Residual2_B0, const ACC<double> &Residual3_B0, const ACC<double> &Residual4_B0, ACC<double> &rhoE_B0, ACC<double>
&rhoE_RKold_B0, ACC<double> &rho_B0, ACC<double> &rho_RKold_B0, ACC<double> &rhou0_B0, ACC<double> &rhou0_RKold_B0,
ACC<double> &rhou1_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_B0, ACC<double> &rhou2_RKold_B0, const double
*rkA, const double *rkB)
{
   rho_RKold_B0(0,0,0) = rkA[0]*rho_RKold_B0(0,0,0) + dt*Residual0_B0(0,0,0);

   rho_B0(0,0,0) = rkB[0]*rho_RKold_B0(0,0,0) + rho_B0(0,0,0);

   rhou0_RKold_B0(0,0,0) = rkA[0]*rhou0_RKold_B0(0,0,0) + dt*Residual1_B0(0,0,0);

   rhou0_B0(0,0,0) = rkB[0]*rhou0_RKold_B0(0,0,0) + rhou0_B0(0,0,0);

   rhou1_RKold_B0(0,0,0) = rkA[0]*rhou1_RKold_B0(0,0,0) + dt*Residual2_B0(0,0,0);

   rhou1_B0(0,0,0) = rkB[0]*rhou1_RKold_B0(0,0,0) + rhou1_B0(0,0,0);

   rhou2_RKold_B0(0,0,0) = rkA[0]*rhou2_RKold_B0(0,0,0) + dt*Residual3_B0(0,0,0);

   rhou2_B0(0,0,0) = rkB[0]*rhou2_RKold_B0(0,0,0) + rhou2_B0(0,0,0);

   rhoE_RKold_B0(0,0,0) = rkA[0]*rhoE_RKold_B0(0,0,0) + dt*Residual4_B0(0,0,0);

   rhoE_B0(0,0,0) = rkB[0]*rhoE_RKold_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel062(ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double> &rhou0_RKold_B0,
ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_RKold_B0)
{
   rho_RKold_B0(0,0,0) = 0.0;

   rhou0_RKold_B0(0,0,0) = 0.0;

   rhou1_RKold_B0(0,0,0) = 0.0;

   rhou2_RKold_B0(0,0,0) = 0.0;

   rhoE_RKold_B0(0,0,0) = 0.0;

}

 void opensbliblock00Kernel069(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, ACC<double> &E_mean_B0, ACC<double> &M_mean_B0, ACC<double>
&TT_mean_B0, ACC<double> &T_mean_B0, ACC<double> &a_mean_B0, ACC<double> &mu_mean_B0, ACC<double> &p_mean_B0,
ACC<double> &pp_mean_B0, ACC<double> &rhomean_B0, ACC<double> &rhou0mean_B0, ACC<double> &rhou0u0mean_B0, ACC<double>
&rhou1mean_B0, ACC<double> &rhou1u0mean_B0, ACC<double> &rhou1u1mean_B0, ACC<double> &rhou2mean_B0, ACC<double>
&rhou2u0mean_B0, ACC<double> &rhou2u1mean_B0, ACC<double> &rhou2u2mean_B0, ACC<double> &u0mean_B0, ACC<double>
&u0u0mean_B0, ACC<double> &u1mean_B0, ACC<double> &u1u0mean_B0, ACC<double> &u1u1mean_B0, ACC<double> &u2mean_B0,
ACC<double> &u2u0mean_B0, ACC<double> &u2u1mean_B0, ACC<double> &u2u2mean_B0)
{
   double rhoinv = 0.0;
   double u0 = 0.0;
   double u1 = 0.0;
   double u2 = 0.0;
   rhomean_B0(0,0,0) = rho_B0(0,0,0) + rhomean_B0(0,0,0);

   rhou0mean_B0(0,0,0) = rhou0_B0(0,0,0) + rhou0mean_B0(0,0,0);

   rhou1mean_B0(0,0,0) = rhou1_B0(0,0,0) + rhou1mean_B0(0,0,0);

   rhou2mean_B0(0,0,0) = rhou2_B0(0,0,0) + rhou2mean_B0(0,0,0);

   E_mean_B0(0,0,0) = rhoE_B0(0,0,0)/rho_B0(0,0,0) + E_mean_B0(0,0,0);

   u0 = rhou0_B0(0,0,0)/rho_B0(0,0,0);

   u1 = rhou1_B0(0,0,0)/rho_B0(0,0,0);

   u2 = rhou2_B0(0,0,0)/rho_B0(0,0,0);

   u0mean_B0(0,0,0) = u0mean_B0(0,0,0) + u0;

   u1mean_B0(0,0,0) = u1mean_B0(0,0,0) + u1;

   u2mean_B0(0,0,0) = u2mean_B0(0,0,0) + u2;

   u0u0mean_B0(0,0,0) = (u0*u0) + u0u0mean_B0(0,0,0);

   u1u0mean_B0(0,0,0) = u0*u1 + u1u0mean_B0(0,0,0);

   u1u1mean_B0(0,0,0) = (u1*u1) + u1u1mean_B0(0,0,0);

   u2u0mean_B0(0,0,0) = u0*u2 + u2u0mean_B0(0,0,0);

   u2u1mean_B0(0,0,0) = u1*u2 + u2u1mean_B0(0,0,0);

   u2u2mean_B0(0,0,0) = (u2*u2) + u2u2mean_B0(0,0,0);

    p_mean_B0(0,0,0) = 0.4*rhoE_B0(0,0,0) - 0.2*(u0*u0)*rho_B0(0,0,0) - 0.2*(u1*u1)*rho_B0(0,0,0) -
      0.2*(u2*u2)*rho_B0(0,0,0) + p_mean_B0(0,0,0);

    pp_mean_B0(0,0,0) = ((0.4*rhoE_B0(0,0,0) - 0.2*(u0*u0)*rho_B0(0,0,0) - 0.2*(u1*u1)*rho_B0(0,0,0) -
      0.2*(u2*u2)*rho_B0(0,0,0))*(0.4*rhoE_B0(0,0,0) - 0.2*(u0*u0)*rho_B0(0,0,0) - 0.2*(u1*u1)*rho_B0(0,0,0) -
      0.2*(u2*u2)*rho_B0(0,0,0))) + pp_mean_B0(0,0,0);

    a_mean_B0(0,0,0) = sqrt(((0.56*rhoE_B0(0,0,0) - 0.28*(u0*u0)*rho_B0(0,0,0) - 0.28*(u1*u1)*rho_B0(0,0,0) -
      0.28*(u2*u2)*rho_B0(0,0,0))/rho_B0(0,0,0))) + a_mean_B0(0,0,0);

    T_mean_B0(0,0,0) = (0.00510734*rhoE_B0(0,0,0) - 0.00255367*(u0*u0)*rho_B0(0,0,0) - 0.00255367*(u1*u1)*rho_B0(0,0,0)
      - 0.00255367*(u2*u2)*rho_B0(0,0,0))/rho_B0(0,0,0) + T_mean_B0(0,0,0);

    TT_mean_B0(0,0,0) = ((0.00510734*rhoE_B0(0,0,0) - 0.00255367*(u0*u0)*rho_B0(0,0,0) -
      0.00255367*(u1*u1)*rho_B0(0,0,0) - 0.00255367*(u2*u2)*rho_B0(0,0,0))*(0.00510734*rhoE_B0(0,0,0) -
      0.00255367*(u0*u0)*rho_B0(0,0,0) - 0.00255367*(u1*u1)*rho_B0(0,0,0) -
      0.00255367*(u2*u2)*rho_B0(0,0,0)))/(rho_B0(0,0,0)*rho_B0(0,0,0)) + TT_mean_B0(0,0,0);

    mu_mean_B0(0,0,0) = pow((0.00510734*rhoE_B0(0,0,0) - 0.00255367*(u0*u0)*rho_B0(0,0,0) -
      0.00255367*(u1*u1)*rho_B0(0,0,0) - 0.00255367*(u2*u2)*rho_B0(0,0,0))/rho_B0(0,0,0), 0.7) + mu_mean_B0(0,0,0);

    M_mean_B0(0,0,0) = pow((0.56*rhoE_B0(0,0,0) - 0.28*(u0*u0)*rho_B0(0,0,0) - 0.28*(u1*u1)*rho_B0(0,0,0) -
      0.28*(u2*u2)*rho_B0(0,0,0))/rho_B0(0,0,0), -0.5)*sqrt(((u0*u0) + (u1*u1) + (u2*u2))) + M_mean_B0(0,0,0);

   rhoinv = 1.0/rho_B0(0,0,0);

   rhou0u0mean_B0(0,0,0) = (rhou0_B0(0,0,0)*rhou0_B0(0,0,0))*rhoinv + rhou0u0mean_B0(0,0,0);

   rhou1u0mean_B0(0,0,0) = rhou0_B0(0,0,0)*rhou1_B0(0,0,0)*rhoinv + rhou1u0mean_B0(0,0,0);

   rhou1u1mean_B0(0,0,0) = (rhou1_B0(0,0,0)*rhou1_B0(0,0,0))*rhoinv + rhou1u1mean_B0(0,0,0);

   rhou2u0mean_B0(0,0,0) = rhou0_B0(0,0,0)*rhou2_B0(0,0,0)*rhoinv + rhou2u0mean_B0(0,0,0);

   rhou2u1mean_B0(0,0,0) = rhou1_B0(0,0,0)*rhou2_B0(0,0,0)*rhoinv + rhou2u1mean_B0(0,0,0);

   rhou2u2mean_B0(0,0,0) = (rhou2_B0(0,0,0)*rhou2_B0(0,0,0))*rhoinv + rhou2u2mean_B0(0,0,0);

}

 void opensbliblock00Kernel063(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double>
&rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_RKold_B0)
{
    rho_RKold_B0(0,0,0) = -(7.0/32.0)*rho_B0(1,0,0) - (7.0/32.0)*rho_B0(-1,0,0) - (1.0/32.0)*rho_B0(-3,0,0) -
      (1.0/32.0)*rho_B0(3,0,0) + ((1.0/256.0))*rho_B0(-4,0,0) + ((1.0/256.0))*rho_B0(4,0,0) +
      ((7.0/64.0))*rho_B0(-2,0,0) + ((7.0/64.0))*rho_B0(2,0,0) + ((35.0/128.0))*rho_B0(0,0,0);

    rhou0_RKold_B0(0,0,0) = -(7.0/32.0)*rhou0_B0(1,0,0) - (7.0/32.0)*rhou0_B0(-1,0,0) - (1.0/32.0)*rhou0_B0(-3,0,0) -
      (1.0/32.0)*rhou0_B0(3,0,0) + ((1.0/256.0))*rhou0_B0(-4,0,0) + ((1.0/256.0))*rhou0_B0(4,0,0) +
      ((7.0/64.0))*rhou0_B0(-2,0,0) + ((7.0/64.0))*rhou0_B0(2,0,0) + ((35.0/128.0))*rhou0_B0(0,0,0);

    rhou1_RKold_B0(0,0,0) = -(7.0/32.0)*rhou1_B0(1,0,0) - (7.0/32.0)*rhou1_B0(-1,0,0) - (1.0/32.0)*rhou1_B0(-3,0,0) -
      (1.0/32.0)*rhou1_B0(3,0,0) + ((1.0/256.0))*rhou1_B0(-4,0,0) + ((1.0/256.0))*rhou1_B0(4,0,0) +
      ((7.0/64.0))*rhou1_B0(-2,0,0) + ((7.0/64.0))*rhou1_B0(2,0,0) + ((35.0/128.0))*rhou1_B0(0,0,0);

    rhou2_RKold_B0(0,0,0) = -(7.0/32.0)*rhou2_B0(1,0,0) - (7.0/32.0)*rhou2_B0(-1,0,0) - (1.0/32.0)*rhou2_B0(-3,0,0) -
      (1.0/32.0)*rhou2_B0(3,0,0) + ((1.0/256.0))*rhou2_B0(-4,0,0) + ((1.0/256.0))*rhou2_B0(4,0,0) +
      ((7.0/64.0))*rhou2_B0(-2,0,0) + ((7.0/64.0))*rhou2_B0(2,0,0) + ((35.0/128.0))*rhou2_B0(0,0,0);

    rhoE_RKold_B0(0,0,0) = -(7.0/32.0)*rhoE_B0(1,0,0) - (7.0/32.0)*rhoE_B0(-1,0,0) - (1.0/32.0)*rhoE_B0(-3,0,0) -
      (1.0/32.0)*rhoE_B0(3,0,0) + ((1.0/256.0))*rhoE_B0(-4,0,0) + ((1.0/256.0))*rhoE_B0(4,0,0) +
      ((7.0/64.0))*rhoE_B0(-2,0,0) + ((7.0/64.0))*rhoE_B0(2,0,0) + ((35.0/128.0))*rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel064(const ACC<double> &rhoE_RKold_B0, const ACC<double> &rho_RKold_B0, const ACC<double>
&rhou0_RKold_B0, const ACC<double> &rhou1_RKold_B0, const ACC<double> &rhou2_RKold_B0, ACC<double> &rhoE_B0, ACC<double>
&rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0, ACC<double> &rhou2_B0)
{
   rho_B0(0,0,0) = -DRP_filt*rho_RKold_B0(0,0,0) + rho_B0(0,0,0);

   rhou0_B0(0,0,0) = -DRP_filt*rhou0_RKold_B0(0,0,0) + rhou0_B0(0,0,0);

   rhou1_B0(0,0,0) = -DRP_filt*rhou1_RKold_B0(0,0,0) + rhou1_B0(0,0,0);

   rhou2_B0(0,0,0) = -DRP_filt*rhou2_RKold_B0(0,0,0) + rhou2_B0(0,0,0);

   rhoE_B0(0,0,0) = -DRP_filt*rhoE_RKold_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel065(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double>
&rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_RKold_B0, const int *idx)
{
   if (idx[1] > -5 + block0np1 || idx[1] < 4){

      rho_RKold_B0(0,0,0) = 0.0;

      rhou0_RKold_B0(0,0,0) = 0.0;

      rhou1_RKold_B0(0,0,0) = 0.0;

      rhou2_RKold_B0(0,0,0) = 0.0;

      rhoE_RKold_B0(0,0,0) = 0.0;

   }

   else{

       rho_RKold_B0(0,0,0) = -(7.0/32.0)*rho_B0(0,1,0) - (7.0/32.0)*rho_B0(0,-1,0) - (1.0/32.0)*rho_B0(0,-3,0) -
            (1.0/32.0)*rho_B0(0,3,0) + ((1.0/256.0))*rho_B0(0,-4,0) + ((1.0/256.0))*rho_B0(0,4,0) +
            ((7.0/64.0))*rho_B0(0,-2,0) + ((7.0/64.0))*rho_B0(0,2,0) + ((35.0/128.0))*rho_B0(0,0,0);

       rhou0_RKold_B0(0,0,0) = -(7.0/32.0)*rhou0_B0(0,1,0) - (7.0/32.0)*rhou0_B0(0,-1,0) - (1.0/32.0)*rhou0_B0(0,-3,0) -
            (1.0/32.0)*rhou0_B0(0,3,0) + ((1.0/256.0))*rhou0_B0(0,-4,0) + ((1.0/256.0))*rhou0_B0(0,4,0) +
            ((7.0/64.0))*rhou0_B0(0,-2,0) + ((7.0/64.0))*rhou0_B0(0,2,0) + ((35.0/128.0))*rhou0_B0(0,0,0);

       rhou1_RKold_B0(0,0,0) = -(7.0/32.0)*rhou1_B0(0,1,0) - (7.0/32.0)*rhou1_B0(0,-1,0) - (1.0/32.0)*rhou1_B0(0,-3,0) -
            (1.0/32.0)*rhou1_B0(0,3,0) + ((1.0/256.0))*rhou1_B0(0,-4,0) + ((1.0/256.0))*rhou1_B0(0,4,0) +
            ((7.0/64.0))*rhou1_B0(0,-2,0) + ((7.0/64.0))*rhou1_B0(0,2,0) + ((35.0/128.0))*rhou1_B0(0,0,0);

       rhou2_RKold_B0(0,0,0) = -(7.0/32.0)*rhou2_B0(0,1,0) - (7.0/32.0)*rhou2_B0(0,-1,0) - (1.0/32.0)*rhou2_B0(0,-3,0) -
            (1.0/32.0)*rhou2_B0(0,3,0) + ((1.0/256.0))*rhou2_B0(0,-4,0) + ((1.0/256.0))*rhou2_B0(0,4,0) +
            ((7.0/64.0))*rhou2_B0(0,-2,0) + ((7.0/64.0))*rhou2_B0(0,2,0) + ((35.0/128.0))*rhou2_B0(0,0,0);

       rhoE_RKold_B0(0,0,0) = -(7.0/32.0)*rhoE_B0(0,1,0) - (7.0/32.0)*rhoE_B0(0,-1,0) - (1.0/32.0)*rhoE_B0(0,-3,0) -
            (1.0/32.0)*rhoE_B0(0,3,0) + ((1.0/256.0))*rhoE_B0(0,-4,0) + ((1.0/256.0))*rhoE_B0(0,4,0) +
            ((7.0/64.0))*rhoE_B0(0,-2,0) + ((7.0/64.0))*rhoE_B0(0,2,0) + ((35.0/128.0))*rhoE_B0(0,0,0);

   }

}

 void opensbliblock00Kernel066(const ACC<double> &rhoE_RKold_B0, const ACC<double> &rho_RKold_B0, const ACC<double>
&rhou0_RKold_B0, const ACC<double> &rhou1_RKold_B0, const ACC<double> &rhou2_RKold_B0, ACC<double> &rhoE_B0, ACC<double>
&rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0, ACC<double> &rhou2_B0)
{
   rho_B0(0,0,0) = -DRP_filt*rho_RKold_B0(0,0,0) + rho_B0(0,0,0);

   rhou0_B0(0,0,0) = -DRP_filt*rhou0_RKold_B0(0,0,0) + rhou0_B0(0,0,0);

   rhou1_B0(0,0,0) = -DRP_filt*rhou1_RKold_B0(0,0,0) + rhou1_B0(0,0,0);

   rhou2_B0(0,0,0) = -DRP_filt*rhou2_RKold_B0(0,0,0) + rhou2_B0(0,0,0);

   rhoE_B0(0,0,0) = -DRP_filt*rhoE_RKold_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel067(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double>
&rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_RKold_B0)
{
    rho_RKold_B0(0,0,0) = -(7.0/32.0)*rho_B0(0,0,1) - (7.0/32.0)*rho_B0(0,0,-1) - (1.0/32.0)*rho_B0(0,0,-3) -
      (1.0/32.0)*rho_B0(0,0,3) + ((1.0/256.0))*rho_B0(0,0,-4) + ((1.0/256.0))*rho_B0(0,0,4) +
      ((7.0/64.0))*rho_B0(0,0,-2) + ((7.0/64.0))*rho_B0(0,0,2) + ((35.0/128.0))*rho_B0(0,0,0);

    rhou0_RKold_B0(0,0,0) = -(7.0/32.0)*rhou0_B0(0,0,1) - (7.0/32.0)*rhou0_B0(0,0,-1) - (1.0/32.0)*rhou0_B0(0,0,-3) -
      (1.0/32.0)*rhou0_B0(0,0,3) + ((1.0/256.0))*rhou0_B0(0,0,-4) + ((1.0/256.0))*rhou0_B0(0,0,4) +
      ((7.0/64.0))*rhou0_B0(0,0,-2) + ((7.0/64.0))*rhou0_B0(0,0,2) + ((35.0/128.0))*rhou0_B0(0,0,0);

    rhou1_RKold_B0(0,0,0) = -(7.0/32.0)*rhou1_B0(0,0,1) - (7.0/32.0)*rhou1_B0(0,0,-1) - (1.0/32.0)*rhou1_B0(0,0,-3) -
      (1.0/32.0)*rhou1_B0(0,0,3) + ((1.0/256.0))*rhou1_B0(0,0,-4) + ((1.0/256.0))*rhou1_B0(0,0,4) +
      ((7.0/64.0))*rhou1_B0(0,0,-2) + ((7.0/64.0))*rhou1_B0(0,0,2) + ((35.0/128.0))*rhou1_B0(0,0,0);

    rhou2_RKold_B0(0,0,0) = -(7.0/32.0)*rhou2_B0(0,0,1) - (7.0/32.0)*rhou2_B0(0,0,-1) - (1.0/32.0)*rhou2_B0(0,0,-3) -
      (1.0/32.0)*rhou2_B0(0,0,3) + ((1.0/256.0))*rhou2_B0(0,0,-4) + ((1.0/256.0))*rhou2_B0(0,0,4) +
      ((7.0/64.0))*rhou2_B0(0,0,-2) + ((7.0/64.0))*rhou2_B0(0,0,2) + ((35.0/128.0))*rhou2_B0(0,0,0);

    rhoE_RKold_B0(0,0,0) = -(7.0/32.0)*rhoE_B0(0,0,1) - (7.0/32.0)*rhoE_B0(0,0,-1) - (1.0/32.0)*rhoE_B0(0,0,-3) -
      (1.0/32.0)*rhoE_B0(0,0,3) + ((1.0/256.0))*rhoE_B0(0,0,-4) + ((1.0/256.0))*rhoE_B0(0,0,4) +
      ((7.0/64.0))*rhoE_B0(0,0,-2) + ((7.0/64.0))*rhoE_B0(0,0,2) + ((35.0/128.0))*rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel068(const ACC<double> &rhoE_RKold_B0, const ACC<double> &rho_RKold_B0, const ACC<double>
&rhou0_RKold_B0, const ACC<double> &rhou1_RKold_B0, const ACC<double> &rhou2_RKold_B0, ACC<double> &rhoE_B0, ACC<double>
&rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0, ACC<double> &rhou2_B0)
{
   rho_B0(0,0,0) = -DRP_filt*rho_RKold_B0(0,0,0) + rho_B0(0,0,0);

   rhou0_B0(0,0,0) = -DRP_filt*rhou0_RKold_B0(0,0,0) + rhou0_B0(0,0,0);

   rhou1_B0(0,0,0) = -DRP_filt*rhou1_RKold_B0(0,0,0) + rhou1_B0(0,0,0);

   rhou2_B0(0,0,0) = -DRP_filt*rhou2_RKold_B0(0,0,0) + rhou2_B0(0,0,0);

   rhoE_B0(0,0,0) = -DRP_filt*rhoE_RKold_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel052(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, ACC<double> &a_B0, ACC<double> &u0_B0, ACC<double> &u1_B0,
ACC<double> &u2_B0, ACC<double> &p_B0)
{
   double inv_rho = 0.0;
   inv_rho = 1.0/rho_B0(0,0,0);

   u0_B0(0,0,0) = rhou0_B0(0,0,0)*inv_rho;

   u1_B0(0,0,0) = rhou1_B0(0,0,0)*inv_rho;

   u2_B0(0,0,0) = rhou2_B0(0,0,0)*inv_rho;

    p_B0(0,0,0) = (-1 + gama)*(-(0.5*(rhou0_B0(0,0,0)*rhou0_B0(0,0,0)) + 0.5*(rhou1_B0(0,0,0)*rhou1_B0(0,0,0)) +
      0.5*(rhou2_B0(0,0,0)*rhou2_B0(0,0,0)))*inv_rho + rhoE_B0(0,0,0));

   a_B0(0,0,0) = sqrt((gama*p_B0(0,0,0)*inv_rho));

}

 void opensbliblock00Kernel056(const ACC<double> &D11_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const
ACC<double> &u2_B0, ACC<double> &kappa_B0, const int *idx)
{
   double d1_u0_dx = 0.0;
   double d1_u0_dy = 0.0;
   double d1_u0_dz = 0.0;
   double d1_u1_dx = 0.0;
   double d1_u1_dy = 0.0;
   double d1_u1_dz = 0.0;
   double d1_u2_dx = 0.0;
   double d1_u2_dy = 0.0;
   double d1_u2_dz = 0.0;
    d1_u0_dz = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

    d1_u0_dy = invDelta1block0*((idx[1] == 0) ? (
   3.0*u0_B0(0,1,0) + 0.333333333333333*u0_B0(0,3,0) -
      1.5*u0_B0(0,2,0) - 1.83333333333333*u0_B0(0,0,0)
)
: ((idx[1] == 1) ? (
   0.0394168524399447*u0_B0(0,2,0) +
      0.00571369039775442*u0_B0(0,4,0) + 0.719443173328855*u0_B0(0,1,0) - 0.322484932882161*u0_B0(0,0,0) -
      0.0658051057710389*u0_B0(0,3,0) - 0.376283677513354*u0_B0(0,-1,0)
)
: ((idx[1] == 2) ? (
  
      0.197184333887745*u0_B0(0,0,0) + 0.521455851089587*u0_B0(0,1,0) + 0.113446470384241*u0_B0(0,-2,0) -
      0.00412637789557492*u0_B0(0,3,0) - 0.0367146847001261*u0_B0(0,2,0) - 0.791245592765872*u0_B0(0,-1,0)
)
: ((idx[1]
      == 3) ? (
   0.0451033223343881*u0_B0(0,0,0) + 0.652141084861241*u0_B0(0,1,0) + 0.121937153224065*u0_B0(0,-2,0) -
      0.00932597985049999*u0_B0(0,-3,0) - 0.727822147724592*u0_B0(0,-1,0) - 0.082033432844602*u0_B0(0,2,0)
)
: ((idx[1]
      == -1 + block0np1) ? (
   1.5*u0_B0(0,-2,0) + 1.83333333333333*u0_B0(0,0,0) - 3.0*u0_B0(0,-1,0) -
      0.333333333333333*u0_B0(0,-3,0)
)
: ((idx[1] == -2 + block0np1) ? (
   0.322484932882161*u0_B0(0,0,0) +
      0.0658051057710389*u0_B0(0,-3,0) + 0.376283677513354*u0_B0(0,1,0) - 0.0394168524399447*u0_B0(0,-2,0) -
      0.00571369039775442*u0_B0(0,-4,0) - 0.719443173328855*u0_B0(0,-1,0)
)
: ((idx[1] == -3 + block0np1) ? (
  
      0.00412637789557492*u0_B0(0,-3,0) + 0.0367146847001261*u0_B0(0,-2,0) + 0.791245592765872*u0_B0(0,1,0) -
      0.197184333887745*u0_B0(0,0,0) - 0.521455851089587*u0_B0(0,-1,0) - 0.113446470384241*u0_B0(0,2,0)
)
: ((idx[1] ==
      -4 + block0np1) ? (
   0.00932597985049999*u0_B0(0,3,0) + 0.727822147724592*u0_B0(0,1,0) +
      0.082033432844602*u0_B0(0,-2,0) - 0.0451033223343881*u0_B0(0,0,0) - 0.652141084861241*u0_B0(0,-1,0) -
      0.121937153224065*u0_B0(0,2,0)
)
: (
   -(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) +
      ((1.0/12.0))*u0_B0(0,-2,0) + ((2.0/3.0))*u0_B0(0,1,0)
)))))))));

    d1_u1_dy = invDelta1block0*((idx[1] == 0) ? (
   3.0*u1_B0(0,1,0) + 0.333333333333333*u1_B0(0,3,0) -
      1.5*u1_B0(0,2,0) - 1.83333333333333*u1_B0(0,0,0)
)
: ((idx[1] == 1) ? (
   0.0394168524399447*u1_B0(0,2,0) +
      0.00571369039775442*u1_B0(0,4,0) + 0.719443173328855*u1_B0(0,1,0) - 0.322484932882161*u1_B0(0,0,0) -
      0.0658051057710389*u1_B0(0,3,0) - 0.376283677513354*u1_B0(0,-1,0)
)
: ((idx[1] == 2) ? (
  
      0.197184333887745*u1_B0(0,0,0) + 0.521455851089587*u1_B0(0,1,0) + 0.113446470384241*u1_B0(0,-2,0) -
      0.00412637789557492*u1_B0(0,3,0) - 0.0367146847001261*u1_B0(0,2,0) - 0.791245592765872*u1_B0(0,-1,0)
)
: ((idx[1]
      == 3) ? (
   0.0451033223343881*u1_B0(0,0,0) + 0.652141084861241*u1_B0(0,1,0) + 0.121937153224065*u1_B0(0,-2,0) -
      0.00932597985049999*u1_B0(0,-3,0) - 0.727822147724592*u1_B0(0,-1,0) - 0.082033432844602*u1_B0(0,2,0)
)
: ((idx[1]
      == -1 + block0np1) ? (
   1.5*u1_B0(0,-2,0) + 1.83333333333333*u1_B0(0,0,0) - 3.0*u1_B0(0,-1,0) -
      0.333333333333333*u1_B0(0,-3,0)
)
: ((idx[1] == -2 + block0np1) ? (
   0.322484932882161*u1_B0(0,0,0) +
      0.0658051057710389*u1_B0(0,-3,0) + 0.376283677513354*u1_B0(0,1,0) - 0.0394168524399447*u1_B0(0,-2,0) -
      0.00571369039775442*u1_B0(0,-4,0) - 0.719443173328855*u1_B0(0,-1,0)
)
: ((idx[1] == -3 + block0np1) ? (
  
      0.00412637789557492*u1_B0(0,-3,0) + 0.0367146847001261*u1_B0(0,-2,0) + 0.791245592765872*u1_B0(0,1,0) -
      0.197184333887745*u1_B0(0,0,0) - 0.521455851089587*u1_B0(0,-1,0) - 0.113446470384241*u1_B0(0,2,0)
)
: ((idx[1] ==
      -4 + block0np1) ? (
   0.00932597985049999*u1_B0(0,3,0) + 0.727822147724592*u1_B0(0,1,0) +
      0.082033432844602*u1_B0(0,-2,0) - 0.0451033223343881*u1_B0(0,0,0) - 0.652141084861241*u1_B0(0,-1,0) -
      0.121937153224065*u1_B0(0,2,0)
)
: (
   -(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) +
      ((1.0/12.0))*u1_B0(0,-2,0) + ((2.0/3.0))*u1_B0(0,1,0)
)))))))));

    d1_u1_dz = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

    d1_u0_dx = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

    d1_u2_dy = invDelta1block0*((idx[1] == 0) ? (
   3.0*u2_B0(0,1,0) + 0.333333333333333*u2_B0(0,3,0) -
      1.5*u2_B0(0,2,0) - 1.83333333333333*u2_B0(0,0,0)
)
: ((idx[1] == 1) ? (
   0.0394168524399447*u2_B0(0,2,0) +
      0.00571369039775442*u2_B0(0,4,0) + 0.719443173328855*u2_B0(0,1,0) - 0.322484932882161*u2_B0(0,0,0) -
      0.0658051057710389*u2_B0(0,3,0) - 0.376283677513354*u2_B0(0,-1,0)
)
: ((idx[1] == 2) ? (
  
      0.197184333887745*u2_B0(0,0,0) + 0.521455851089587*u2_B0(0,1,0) + 0.113446470384241*u2_B0(0,-2,0) -
      0.00412637789557492*u2_B0(0,3,0) - 0.0367146847001261*u2_B0(0,2,0) - 0.791245592765872*u2_B0(0,-1,0)
)
: ((idx[1]
      == 3) ? (
   0.0451033223343881*u2_B0(0,0,0) + 0.652141084861241*u2_B0(0,1,0) + 0.121937153224065*u2_B0(0,-2,0) -
      0.00932597985049999*u2_B0(0,-3,0) - 0.727822147724592*u2_B0(0,-1,0) - 0.082033432844602*u2_B0(0,2,0)
)
: ((idx[1]
      == -1 + block0np1) ? (
   1.5*u2_B0(0,-2,0) + 1.83333333333333*u2_B0(0,0,0) - 3.0*u2_B0(0,-1,0) -
      0.333333333333333*u2_B0(0,-3,0)
)
: ((idx[1] == -2 + block0np1) ? (
   0.322484932882161*u2_B0(0,0,0) +
      0.0658051057710389*u2_B0(0,-3,0) + 0.376283677513354*u2_B0(0,1,0) - 0.0394168524399447*u2_B0(0,-2,0) -
      0.00571369039775442*u2_B0(0,-4,0) - 0.719443173328855*u2_B0(0,-1,0)
)
: ((idx[1] == -3 + block0np1) ? (
  
      0.00412637789557492*u2_B0(0,-3,0) + 0.0367146847001261*u2_B0(0,-2,0) + 0.791245592765872*u2_B0(0,1,0) -
      0.197184333887745*u2_B0(0,0,0) - 0.521455851089587*u2_B0(0,-1,0) - 0.113446470384241*u2_B0(0,2,0)
)
: ((idx[1] ==
      -4 + block0np1) ? (
   0.00932597985049999*u2_B0(0,3,0) + 0.727822147724592*u2_B0(0,1,0) +
      0.082033432844602*u2_B0(0,-2,0) - 0.0451033223343881*u2_B0(0,0,0) - 0.652141084861241*u2_B0(0,-1,0) -
      0.121937153224065*u2_B0(0,2,0)
)
: (
   -(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) +
      ((1.0/12.0))*u2_B0(0,-2,0) + ((2.0/3.0))*u2_B0(0,1,0)
)))))))));

    d1_u2_dz = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

    d1_u2_dx = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

    d1_u1_dx = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

    kappa_B0(0,0,0) = ((D11_B0(0,0,0)*d1_u1_dy + d1_u0_dx + d1_u2_dz)*(D11_B0(0,0,0)*d1_u1_dy + d1_u0_dx +
      d1_u2_dz))*(0.5 - 0.5*tanh(2.5 + 500.0*(D11_B0(0,0,0)*d1_u1_dy + d1_u0_dx +
      d1_u2_dz)/sqrt(((Delta0block0*Delta0block0) + (Delta1block0*Delta1block0) +
      (Delta2block0*Delta2block0)))))/(1.0e-40 + ((-d1_u1_dz + D11_B0(0,0,0)*d1_u2_dy)*(-d1_u1_dz +
      D11_B0(0,0,0)*d1_u2_dy)) + ((-d1_u2_dx + d1_u0_dz)*(-d1_u2_dx + d1_u0_dz)) + ((-D11_B0(0,0,0)*d1_u0_dy +
      d1_u1_dx)*(-D11_B0(0,0,0)*d1_u0_dy + d1_u1_dx)) + ((D11_B0(0,0,0)*d1_u1_dy + d1_u0_dx +
      d1_u2_dz)*(D11_B0(0,0,0)*d1_u1_dy + d1_u0_dx + d1_u2_dz)));

}

 void opensbliblock00Kernel057(ACC<double> &Residual0_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0,
ACC<double> &Residual3_B0, ACC<double> &Residual4_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double>
&rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_RKold_B0, ACC<double> &wk0_B0, ACC<double> &wk1_B0,
ACC<double> &wk2_B0, ACC<double> &wk3_B0, ACC<double> &wk4_B0)
{
   wk0_B0(0,0,0) = 0.0;

   wk1_B0(0,0,0) = 0.0;

   wk2_B0(0,0,0) = 0.0;

   wk3_B0(0,0,0) = 0.0;

   wk4_B0(0,0,0) = 0.0;

   Residual0_B0(0,0,0) = 0.0;

   Residual1_B0(0,0,0) = 0.0;

   Residual2_B0(0,0,0) = 0.0;

   Residual3_B0(0,0,0) = 0.0;

   Residual4_B0(0,0,0) = 0.0;

   rho_RKold_B0(0,0,0) = 0.0;

   rhou0_RKold_B0(0,0,0) = 0.0;

   rhou1_RKold_B0(0,0,0) = 0.0;

   rhou2_RKold_B0(0,0,0) = 0.0;

   rhoE_RKold_B0(0,0,0) = 0.0;

}

 void opensbliblock00Kernel058(const ACC<double> &a_B0, const ACC<double> &kappa_B0, const ACC<double> &p_B0, const
ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&wk0_B0, ACC<double> &wk1_B0, ACC<double> &wk2_B0, ACC<double> &wk3_B0, ACC<double> &wk4_B0)
{
   double AVG_0_0_LEV_00 = 0.0;
   double AVG_0_0_LEV_01 = 0.0;
   double AVG_0_0_LEV_02 = 0.0;
   double AVG_0_0_LEV_03 = 0.0;
   double AVG_0_0_LEV_04 = 0.0;
   double AVG_0_0_LEV_10 = 0.0;
   double AVG_0_0_LEV_13 = 0.0;
   double AVG_0_0_LEV_20 = 0.0;
   double AVG_0_0_LEV_22 = 0.0;
   double AVG_0_0_LEV_30 = 0.0;
   double AVG_0_0_LEV_31 = 0.0;
   double AVG_0_0_LEV_32 = 0.0;
   double AVG_0_0_LEV_33 = 0.0;
   double AVG_0_0_LEV_34 = 0.0;
   double AVG_0_0_LEV_40 = 0.0;
   double AVG_0_0_LEV_41 = 0.0;
   double AVG_0_0_LEV_42 = 0.0;
   double AVG_0_0_LEV_43 = 0.0;
   double AVG_0_0_LEV_44 = 0.0;
   double AVG_0_a = 0.0;
   double AVG_0_inv_rho = 0.0;
   double AVG_0_rho = 0.0;
   double AVG_0_u0 = 0.0;
   double AVG_0_u1 = 0.0;
   double AVG_0_u2 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_06 = 0.0;
   double CF_07 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_16 = 0.0;
   double CF_17 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_26 = 0.0;
   double CF_27 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CF_36 = 0.0;
   double CF_37 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CF_44 = 0.0;
   double CF_45 = 0.0;
   double CF_46 = 0.0;
   double CF_47 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_06 = 0.0;
   double CS_07 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_16 = 0.0;
   double CS_17 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_26 = 0.0;
   double CS_27 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double CS_36 = 0.0;
   double CS_37 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double CS_44 = 0.0;
   double CS_45 = 0.0;
   double CS_46 = 0.0;
   double CS_47 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double alpha_3 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double beta_3 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double omega_2 = 0.0;
   double omega_3 = 0.0;
   double rj0 = 0.0;
   double rj1 = 0.0;
   double rj2 = 0.0;
   double rj3 = 0.0;
   double rj4 = 0.0;
    if (fmax(kappa_B0(2,0,0), fmax(kappa_B0(-2,0,0), fmax(kappa_B0(-1,0,0), fmax(kappa_B0(1,0,0), fmax(kappa_B0(-3,0,0),
      kappa_B0(0,0,0)))))) > Ducros_check){

      AVG_0_rho = sqrt((rho_B0(0,0,0)*rho_B0(1,0,0)));

      AVG_0_inv_rho = 1.0/((sqrt(rho_B0(0,0,0)) + sqrt(rho_B0(1,0,0))));

      AVG_0_u0 = (sqrt(rho_B0(0,0,0))*u0_B0(0,0,0) + sqrt(rho_B0(1,0,0))*u0_B0(1,0,0))*AVG_0_inv_rho;

      AVG_0_u1 = (sqrt(rho_B0(0,0,0))*u1_B0(0,0,0) + sqrt(rho_B0(1,0,0))*u1_B0(1,0,0))*AVG_0_inv_rho;

      AVG_0_u2 = (sqrt(rho_B0(0,0,0))*u2_B0(0,0,0) + sqrt(rho_B0(1,0,0))*u2_B0(1,0,0))*AVG_0_inv_rho;

       AVG_0_a = sqrt(((-(1.0/2.0)*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) + (AVG_0_u2*AVG_0_u2)) + ((p_B0(0,0,0) +
            rhoE_B0(0,0,0))/sqrt(rho_B0(0,0,0)) + (p_B0(1,0,0) +
            rhoE_B0(1,0,0))/sqrt(rho_B0(1,0,0)))*AVG_0_inv_rho)*gamma_m1));

      inv_AVG_a = 1.0/(AVG_0_a);

      inv_AVG_rho = 1.0/(AVG_0_rho);

       AVG_0_0_LEV_00 = -(1.0/2.0)*(-2 - (AVG_0_u0*AVG_0_u0)*(inv_AVG_a*inv_AVG_a) -
            (AVG_0_u1*AVG_0_u1)*(inv_AVG_a*inv_AVG_a) - (AVG_0_u2*AVG_0_u2)*(inv_AVG_a*inv_AVG_a) +
            (AVG_0_u0*AVG_0_u0)*(inv_AVG_a*inv_AVG_a)*gama + (AVG_0_u1*AVG_0_u1)*(inv_AVG_a*inv_AVG_a)*gama +
            (AVG_0_u2*AVG_0_u2)*(inv_AVG_a*inv_AVG_a)*gama);

      AVG_0_0_LEV_01 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u0;

      AVG_0_0_LEV_02 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u1;

      AVG_0_0_LEV_03 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u2;

      AVG_0_0_LEV_04 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

      AVG_0_0_LEV_10 = -AVG_0_u2*inv_AVG_rho;

      AVG_0_0_LEV_13 = inv_AVG_rho;

      AVG_0_0_LEV_20 = AVG_0_u1*inv_AVG_rho;

      AVG_0_0_LEV_22 = -inv_AVG_rho;

       AVG_0_0_LEV_30 = -0.353553390593274*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) + (AVG_0_u2*AVG_0_u2) -
            (AVG_0_u0*AVG_0_u0)*gama - (AVG_0_u1*AVG_0_u1)*gama - (AVG_0_u2*AVG_0_u2)*gama +
            2*AVG_0_a*AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_31 = 0.707106781186547*(-gama*AVG_0_u0 + AVG_0_a + AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_32 = -0.707106781186547*gamma_m1*AVG_0_u1*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_33 = -0.707106781186547*gamma_m1*AVG_0_u2*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_34 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       AVG_0_0_LEV_40 = 0.353553390593274*(-(AVG_0_u0*AVG_0_u0) - (AVG_0_u1*AVG_0_u1) - (AVG_0_u2*AVG_0_u2) +
            (AVG_0_u0*AVG_0_u0)*gama + (AVG_0_u1*AVG_0_u1)*gama + (AVG_0_u2*AVG_0_u2)*gama +
            2*AVG_0_a*AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_41 = -0.707106781186547*(-AVG_0_u0 + gama*AVG_0_u0 + AVG_0_a)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_42 = -0.707106781186547*gamma_m1*AVG_0_u1*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_43 = -0.707106781186547*gamma_m1*AVG_0_u2*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_44 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       CF_00 = p_B0(-3,0,0)*AVG_0_0_LEV_01 + p_B0(-3,0,0)*u0_B0(-3,0,0)*AVG_0_0_LEV_04 +
            u0_B0(-3,0,0)*rho_B0(-3,0,0)*AVG_0_0_LEV_00 + u0_B0(-3,0,0)*rhoE_B0(-3,0,0)*AVG_0_0_LEV_04 +
            u0_B0(-3,0,0)*rhou0_B0(-3,0,0)*AVG_0_0_LEV_01 + u0_B0(-3,0,0)*rhou1_B0(-3,0,0)*AVG_0_0_LEV_02 +
            u0_B0(-3,0,0)*rhou2_B0(-3,0,0)*AVG_0_0_LEV_03;

      CF_10 = (rho_B0(-3,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-3,0,0)*AVG_0_0_LEV_13)*u0_B0(-3,0,0);

      CF_20 = (rho_B0(-3,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-3,0,0)*AVG_0_0_LEV_22)*u0_B0(-3,0,0);

       CF_30 = p_B0(-3,0,0)*AVG_0_0_LEV_31 + p_B0(-3,0,0)*u0_B0(-3,0,0)*AVG_0_0_LEV_34 +
            u0_B0(-3,0,0)*rho_B0(-3,0,0)*AVG_0_0_LEV_30 + u0_B0(-3,0,0)*rhoE_B0(-3,0,0)*AVG_0_0_LEV_34 +
            u0_B0(-3,0,0)*rhou0_B0(-3,0,0)*AVG_0_0_LEV_31 + u0_B0(-3,0,0)*rhou1_B0(-3,0,0)*AVG_0_0_LEV_32 +
            u0_B0(-3,0,0)*rhou2_B0(-3,0,0)*AVG_0_0_LEV_33;

       CF_40 = p_B0(-3,0,0)*AVG_0_0_LEV_41 + p_B0(-3,0,0)*u0_B0(-3,0,0)*AVG_0_0_LEV_44 +
            u0_B0(-3,0,0)*rho_B0(-3,0,0)*AVG_0_0_LEV_40 + u0_B0(-3,0,0)*rhoE_B0(-3,0,0)*AVG_0_0_LEV_44 +
            u0_B0(-3,0,0)*rhou0_B0(-3,0,0)*AVG_0_0_LEV_41 + u0_B0(-3,0,0)*rhou1_B0(-3,0,0)*AVG_0_0_LEV_42 +
            u0_B0(-3,0,0)*rhou2_B0(-3,0,0)*AVG_0_0_LEV_43;

       CS_00 = rho_B0(-3,0,0)*AVG_0_0_LEV_00 + rhoE_B0(-3,0,0)*AVG_0_0_LEV_04 + rhou0_B0(-3,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(-3,0,0)*AVG_0_0_LEV_02 + rhou2_B0(-3,0,0)*AVG_0_0_LEV_03;

      CS_10 = rho_B0(-3,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-3,0,0)*AVG_0_0_LEV_13;

      CS_20 = rho_B0(-3,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-3,0,0)*AVG_0_0_LEV_22;

       CS_30 = rho_B0(-3,0,0)*AVG_0_0_LEV_30 + rhoE_B0(-3,0,0)*AVG_0_0_LEV_34 + rhou0_B0(-3,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(-3,0,0)*AVG_0_0_LEV_32 + rhou2_B0(-3,0,0)*AVG_0_0_LEV_33;

       CS_40 = rho_B0(-3,0,0)*AVG_0_0_LEV_40 + rhoE_B0(-3,0,0)*AVG_0_0_LEV_44 + rhou0_B0(-3,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(-3,0,0)*AVG_0_0_LEV_42 + rhou2_B0(-3,0,0)*AVG_0_0_LEV_43;

       CF_01 = p_B0(-2,0,0)*AVG_0_0_LEV_01 + p_B0(-2,0,0)*u0_B0(-2,0,0)*AVG_0_0_LEV_04 +
            u0_B0(-2,0,0)*rho_B0(-2,0,0)*AVG_0_0_LEV_00 + u0_B0(-2,0,0)*rhoE_B0(-2,0,0)*AVG_0_0_LEV_04 +
            u0_B0(-2,0,0)*rhou0_B0(-2,0,0)*AVG_0_0_LEV_01 + u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_02 +
            u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_03;

      CF_11 = (rho_B0(-2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_13)*u0_B0(-2,0,0);

      CF_21 = (rho_B0(-2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-2,0,0)*AVG_0_0_LEV_22)*u0_B0(-2,0,0);

       CF_31 = p_B0(-2,0,0)*AVG_0_0_LEV_31 + p_B0(-2,0,0)*u0_B0(-2,0,0)*AVG_0_0_LEV_34 +
            u0_B0(-2,0,0)*rho_B0(-2,0,0)*AVG_0_0_LEV_30 + u0_B0(-2,0,0)*rhoE_B0(-2,0,0)*AVG_0_0_LEV_34 +
            u0_B0(-2,0,0)*rhou0_B0(-2,0,0)*AVG_0_0_LEV_31 + u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_32 +
            u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_33;

       CF_41 = p_B0(-2,0,0)*AVG_0_0_LEV_41 + p_B0(-2,0,0)*u0_B0(-2,0,0)*AVG_0_0_LEV_44 +
            u0_B0(-2,0,0)*rho_B0(-2,0,0)*AVG_0_0_LEV_40 + u0_B0(-2,0,0)*rhoE_B0(-2,0,0)*AVG_0_0_LEV_44 +
            u0_B0(-2,0,0)*rhou0_B0(-2,0,0)*AVG_0_0_LEV_41 + u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_42 +
            u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_43;

       CS_01 = rho_B0(-2,0,0)*AVG_0_0_LEV_00 + rhoE_B0(-2,0,0)*AVG_0_0_LEV_04 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(-2,0,0)*AVG_0_0_LEV_02 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_03;

      CS_11 = rho_B0(-2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_13;

      CS_21 = rho_B0(-2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-2,0,0)*AVG_0_0_LEV_22;

       CS_31 = rho_B0(-2,0,0)*AVG_0_0_LEV_30 + rhoE_B0(-2,0,0)*AVG_0_0_LEV_34 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(-2,0,0)*AVG_0_0_LEV_32 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_33;

       CS_41 = rho_B0(-2,0,0)*AVG_0_0_LEV_40 + rhoE_B0(-2,0,0)*AVG_0_0_LEV_44 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(-2,0,0)*AVG_0_0_LEV_42 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_43;

       CF_02 = p_B0(-1,0,0)*AVG_0_0_LEV_01 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_04 +
            u0_B0(-1,0,0)*rho_B0(-1,0,0)*AVG_0_0_LEV_00 + u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_04 +
            u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_01 + u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_02 +
            u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_03;

      CF_12 = (rho_B0(-1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_13)*u0_B0(-1,0,0);

      CF_22 = (rho_B0(-1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-1,0,0)*AVG_0_0_LEV_22)*u0_B0(-1,0,0);

       CF_32 = p_B0(-1,0,0)*AVG_0_0_LEV_31 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_34 +
            u0_B0(-1,0,0)*rho_B0(-1,0,0)*AVG_0_0_LEV_30 + u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_34 +
            u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_31 + u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_32 +
            u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_33;

       CF_42 = p_B0(-1,0,0)*AVG_0_0_LEV_41 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_44 +
            u0_B0(-1,0,0)*rho_B0(-1,0,0)*AVG_0_0_LEV_40 + u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_44 +
            u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_41 + u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_42 +
            u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_43;

       CS_02 = rho_B0(-1,0,0)*AVG_0_0_LEV_00 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_04 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(-1,0,0)*AVG_0_0_LEV_02 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_03;

      CS_12 = rho_B0(-1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_13;

      CS_22 = rho_B0(-1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-1,0,0)*AVG_0_0_LEV_22;

       CS_32 = rho_B0(-1,0,0)*AVG_0_0_LEV_30 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_34 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(-1,0,0)*AVG_0_0_LEV_32 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_33;

       CS_42 = rho_B0(-1,0,0)*AVG_0_0_LEV_40 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_44 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(-1,0,0)*AVG_0_0_LEV_42 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_43;

       CF_03 = p_B0(0,0,0)*AVG_0_0_LEV_01 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_04 +
            u0_B0(0,0,0)*rho_B0(0,0,0)*AVG_0_0_LEV_00 + u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_04 +
            u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_01 + u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_02 +
            u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_03;

      CF_13 = (rho_B0(0,0,0)*AVG_0_0_LEV_10 + rhou2_B0(0,0,0)*AVG_0_0_LEV_13)*u0_B0(0,0,0);

      CF_23 = (rho_B0(0,0,0)*AVG_0_0_LEV_20 + rhou1_B0(0,0,0)*AVG_0_0_LEV_22)*u0_B0(0,0,0);

       CF_33 = p_B0(0,0,0)*AVG_0_0_LEV_31 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_34 +
            u0_B0(0,0,0)*rho_B0(0,0,0)*AVG_0_0_LEV_30 + u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_34 +
            u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_31 + u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_32 +
            u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_33;

       CF_43 = p_B0(0,0,0)*AVG_0_0_LEV_41 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_44 +
            u0_B0(0,0,0)*rho_B0(0,0,0)*AVG_0_0_LEV_40 + u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_44 +
            u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_41 + u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_42 +
            u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_43;

       CS_03 = rho_B0(0,0,0)*AVG_0_0_LEV_00 + rhoE_B0(0,0,0)*AVG_0_0_LEV_04 + rhou0_B0(0,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(0,0,0)*AVG_0_0_LEV_02 + rhou2_B0(0,0,0)*AVG_0_0_LEV_03;

      CS_13 = rho_B0(0,0,0)*AVG_0_0_LEV_10 + rhou2_B0(0,0,0)*AVG_0_0_LEV_13;

      CS_23 = rho_B0(0,0,0)*AVG_0_0_LEV_20 + rhou1_B0(0,0,0)*AVG_0_0_LEV_22;

       CS_33 = rho_B0(0,0,0)*AVG_0_0_LEV_30 + rhoE_B0(0,0,0)*AVG_0_0_LEV_34 + rhou0_B0(0,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(0,0,0)*AVG_0_0_LEV_32 + rhou2_B0(0,0,0)*AVG_0_0_LEV_33;

       CS_43 = rho_B0(0,0,0)*AVG_0_0_LEV_40 + rhoE_B0(0,0,0)*AVG_0_0_LEV_44 + rhou0_B0(0,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(0,0,0)*AVG_0_0_LEV_42 + rhou2_B0(0,0,0)*AVG_0_0_LEV_43;

       CF_04 = p_B0(1,0,0)*AVG_0_0_LEV_01 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_04 +
            u0_B0(1,0,0)*rho_B0(1,0,0)*AVG_0_0_LEV_00 + u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_04 +
            u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_01 + u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_02 +
            u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_03;

      CF_14 = (rho_B0(1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(1,0,0)*AVG_0_0_LEV_13)*u0_B0(1,0,0);

      CF_24 = (rho_B0(1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(1,0,0)*AVG_0_0_LEV_22)*u0_B0(1,0,0);

       CF_34 = p_B0(1,0,0)*AVG_0_0_LEV_31 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_34 +
            u0_B0(1,0,0)*rho_B0(1,0,0)*AVG_0_0_LEV_30 + u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_34 +
            u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_31 + u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_32 +
            u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_33;

       CF_44 = p_B0(1,0,0)*AVG_0_0_LEV_41 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_44 +
            u0_B0(1,0,0)*rho_B0(1,0,0)*AVG_0_0_LEV_40 + u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_44 +
            u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_41 + u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_42 +
            u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_43;

       CS_04 = rho_B0(1,0,0)*AVG_0_0_LEV_00 + rhoE_B0(1,0,0)*AVG_0_0_LEV_04 + rhou0_B0(1,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(1,0,0)*AVG_0_0_LEV_02 + rhou2_B0(1,0,0)*AVG_0_0_LEV_03;

      CS_14 = rho_B0(1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(1,0,0)*AVG_0_0_LEV_13;

      CS_24 = rho_B0(1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(1,0,0)*AVG_0_0_LEV_22;

       CS_34 = rho_B0(1,0,0)*AVG_0_0_LEV_30 + rhoE_B0(1,0,0)*AVG_0_0_LEV_34 + rhou0_B0(1,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(1,0,0)*AVG_0_0_LEV_32 + rhou2_B0(1,0,0)*AVG_0_0_LEV_33;

       CS_44 = rho_B0(1,0,0)*AVG_0_0_LEV_40 + rhoE_B0(1,0,0)*AVG_0_0_LEV_44 + rhou0_B0(1,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(1,0,0)*AVG_0_0_LEV_42 + rhou2_B0(1,0,0)*AVG_0_0_LEV_43;

       CF_05 = p_B0(2,0,0)*AVG_0_0_LEV_01 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_04 +
            u0_B0(2,0,0)*rho_B0(2,0,0)*AVG_0_0_LEV_00 + u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_04 +
            u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_01 + u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_02 +
            u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_03;

      CF_15 = (rho_B0(2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(2,0,0)*AVG_0_0_LEV_13)*u0_B0(2,0,0);

      CF_25 = (rho_B0(2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(2,0,0)*AVG_0_0_LEV_22)*u0_B0(2,0,0);

       CF_35 = p_B0(2,0,0)*AVG_0_0_LEV_31 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_34 +
            u0_B0(2,0,0)*rho_B0(2,0,0)*AVG_0_0_LEV_30 + u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_34 +
            u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_31 + u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_32 +
            u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_33;

       CF_45 = p_B0(2,0,0)*AVG_0_0_LEV_41 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_44 +
            u0_B0(2,0,0)*rho_B0(2,0,0)*AVG_0_0_LEV_40 + u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_44 +
            u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_41 + u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_42 +
            u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_43;

       CS_05 = rho_B0(2,0,0)*AVG_0_0_LEV_00 + rhoE_B0(2,0,0)*AVG_0_0_LEV_04 + rhou0_B0(2,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(2,0,0)*AVG_0_0_LEV_02 + rhou2_B0(2,0,0)*AVG_0_0_LEV_03;

      CS_15 = rho_B0(2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(2,0,0)*AVG_0_0_LEV_13;

      CS_25 = rho_B0(2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(2,0,0)*AVG_0_0_LEV_22;

       CS_35 = rho_B0(2,0,0)*AVG_0_0_LEV_30 + rhoE_B0(2,0,0)*AVG_0_0_LEV_34 + rhou0_B0(2,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(2,0,0)*AVG_0_0_LEV_32 + rhou2_B0(2,0,0)*AVG_0_0_LEV_33;

       CS_45 = rho_B0(2,0,0)*AVG_0_0_LEV_40 + rhoE_B0(2,0,0)*AVG_0_0_LEV_44 + rhou0_B0(2,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(2,0,0)*AVG_0_0_LEV_42 + rhou2_B0(2,0,0)*AVG_0_0_LEV_43;

       CF_06 = p_B0(3,0,0)*AVG_0_0_LEV_01 + p_B0(3,0,0)*u0_B0(3,0,0)*AVG_0_0_LEV_04 +
            u0_B0(3,0,0)*rho_B0(3,0,0)*AVG_0_0_LEV_00 + u0_B0(3,0,0)*rhoE_B0(3,0,0)*AVG_0_0_LEV_04 +
            u0_B0(3,0,0)*rhou0_B0(3,0,0)*AVG_0_0_LEV_01 + u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_02 +
            u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_03;

      CF_16 = (rho_B0(3,0,0)*AVG_0_0_LEV_10 + rhou2_B0(3,0,0)*AVG_0_0_LEV_13)*u0_B0(3,0,0);

      CF_26 = (rho_B0(3,0,0)*AVG_0_0_LEV_20 + rhou1_B0(3,0,0)*AVG_0_0_LEV_22)*u0_B0(3,0,0);

       CF_36 = p_B0(3,0,0)*AVG_0_0_LEV_31 + p_B0(3,0,0)*u0_B0(3,0,0)*AVG_0_0_LEV_34 +
            u0_B0(3,0,0)*rho_B0(3,0,0)*AVG_0_0_LEV_30 + u0_B0(3,0,0)*rhoE_B0(3,0,0)*AVG_0_0_LEV_34 +
            u0_B0(3,0,0)*rhou0_B0(3,0,0)*AVG_0_0_LEV_31 + u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_32 +
            u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_33;

       CF_46 = p_B0(3,0,0)*AVG_0_0_LEV_41 + p_B0(3,0,0)*u0_B0(3,0,0)*AVG_0_0_LEV_44 +
            u0_B0(3,0,0)*rho_B0(3,0,0)*AVG_0_0_LEV_40 + u0_B0(3,0,0)*rhoE_B0(3,0,0)*AVG_0_0_LEV_44 +
            u0_B0(3,0,0)*rhou0_B0(3,0,0)*AVG_0_0_LEV_41 + u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_42 +
            u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_43;

       CS_06 = rho_B0(3,0,0)*AVG_0_0_LEV_00 + rhoE_B0(3,0,0)*AVG_0_0_LEV_04 + rhou0_B0(3,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(3,0,0)*AVG_0_0_LEV_02 + rhou2_B0(3,0,0)*AVG_0_0_LEV_03;

      CS_16 = rho_B0(3,0,0)*AVG_0_0_LEV_10 + rhou2_B0(3,0,0)*AVG_0_0_LEV_13;

      CS_26 = rho_B0(3,0,0)*AVG_0_0_LEV_20 + rhou1_B0(3,0,0)*AVG_0_0_LEV_22;

       CS_36 = rho_B0(3,0,0)*AVG_0_0_LEV_30 + rhoE_B0(3,0,0)*AVG_0_0_LEV_34 + rhou0_B0(3,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(3,0,0)*AVG_0_0_LEV_32 + rhou2_B0(3,0,0)*AVG_0_0_LEV_33;

       CS_46 = rho_B0(3,0,0)*AVG_0_0_LEV_40 + rhoE_B0(3,0,0)*AVG_0_0_LEV_44 + rhou0_B0(3,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(3,0,0)*AVG_0_0_LEV_42 + rhou2_B0(3,0,0)*AVG_0_0_LEV_43;

       CF_07 = p_B0(4,0,0)*AVG_0_0_LEV_01 + p_B0(4,0,0)*u0_B0(4,0,0)*AVG_0_0_LEV_04 +
            u0_B0(4,0,0)*rho_B0(4,0,0)*AVG_0_0_LEV_00 + u0_B0(4,0,0)*rhoE_B0(4,0,0)*AVG_0_0_LEV_04 +
            u0_B0(4,0,0)*rhou0_B0(4,0,0)*AVG_0_0_LEV_01 + u0_B0(4,0,0)*rhou1_B0(4,0,0)*AVG_0_0_LEV_02 +
            u0_B0(4,0,0)*rhou2_B0(4,0,0)*AVG_0_0_LEV_03;

      CF_17 = (rho_B0(4,0,0)*AVG_0_0_LEV_10 + rhou2_B0(4,0,0)*AVG_0_0_LEV_13)*u0_B0(4,0,0);

      CF_27 = (rho_B0(4,0,0)*AVG_0_0_LEV_20 + rhou1_B0(4,0,0)*AVG_0_0_LEV_22)*u0_B0(4,0,0);

       CF_37 = p_B0(4,0,0)*AVG_0_0_LEV_31 + p_B0(4,0,0)*u0_B0(4,0,0)*AVG_0_0_LEV_34 +
            u0_B0(4,0,0)*rho_B0(4,0,0)*AVG_0_0_LEV_30 + u0_B0(4,0,0)*rhoE_B0(4,0,0)*AVG_0_0_LEV_34 +
            u0_B0(4,0,0)*rhou0_B0(4,0,0)*AVG_0_0_LEV_31 + u0_B0(4,0,0)*rhou1_B0(4,0,0)*AVG_0_0_LEV_32 +
            u0_B0(4,0,0)*rhou2_B0(4,0,0)*AVG_0_0_LEV_33;

       CF_47 = p_B0(4,0,0)*AVG_0_0_LEV_41 + p_B0(4,0,0)*u0_B0(4,0,0)*AVG_0_0_LEV_44 +
            u0_B0(4,0,0)*rho_B0(4,0,0)*AVG_0_0_LEV_40 + u0_B0(4,0,0)*rhoE_B0(4,0,0)*AVG_0_0_LEV_44 +
            u0_B0(4,0,0)*rhou0_B0(4,0,0)*AVG_0_0_LEV_41 + u0_B0(4,0,0)*rhou1_B0(4,0,0)*AVG_0_0_LEV_42 +
            u0_B0(4,0,0)*rhou2_B0(4,0,0)*AVG_0_0_LEV_43;

       CS_07 = rho_B0(4,0,0)*AVG_0_0_LEV_00 + rhoE_B0(4,0,0)*AVG_0_0_LEV_04 + rhou0_B0(4,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(4,0,0)*AVG_0_0_LEV_02 + rhou2_B0(4,0,0)*AVG_0_0_LEV_03;

      CS_17 = rho_B0(4,0,0)*AVG_0_0_LEV_10 + rhou2_B0(4,0,0)*AVG_0_0_LEV_13;

      CS_27 = rho_B0(4,0,0)*AVG_0_0_LEV_20 + rhou1_B0(4,0,0)*AVG_0_0_LEV_22;

       CS_37 = rho_B0(4,0,0)*AVG_0_0_LEV_30 + rhoE_B0(4,0,0)*AVG_0_0_LEV_34 + rhou0_B0(4,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(4,0,0)*AVG_0_0_LEV_32 + rhou2_B0(4,0,0)*AVG_0_0_LEV_33;

       CS_47 = rho_B0(4,0,0)*AVG_0_0_LEV_40 + rhoE_B0(4,0,0)*AVG_0_0_LEV_44 + rhou0_B0(4,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(4,0,0)*AVG_0_0_LEV_42 + rhou2_B0(4,0,0)*AVG_0_0_LEV_43;

      max_lambda_00 = shock_filter_control*fmax(fabs(u0_B0(1,0,0)), fabs(u0_B0(0,0,0)));

      max_lambda_11 = max_lambda_00;

      max_lambda_22 = max_lambda_00;

      max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(1,0,0) + u0_B0(1,0,0)), fabs(a_B0(0,0,0) + u0_B0(0,0,0)));

      max_lambda_44 = shock_filter_control*fmax(fabs(-u0_B0(0,0,0) + a_B0(0,0,0)), fabs(-u0_B0(1,0,0) + a_B0(1,0,0)));

       beta_0 = ((547.0/960.0))*((CS_06*max_lambda_00 + CF_06)*(CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_06*max_lambda_00 + CF_06) + ((7043.0/480.0))*(CS_05*max_lambda_00 +
            CF_05))*(CS_05*max_lambda_00 + CF_05) + ((1.0/2.0))*(CS_03*max_lambda_00 +
            CF_03)*(-(1567.0/80.0)*(CS_04*max_lambda_00 + CF_04) - (309.0/80.0)*(CS_06*max_lambda_00 + CF_06) +
            ((2107.0/480.0))*(CS_03*max_lambda_00 + CF_03) + ((3521.0/240.0))*(CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04)*(-(8623.0/240.0)*(CS_05*max_lambda_00 + CF_05) +
            ((2321.0/240.0))*(CS_06*max_lambda_00 + CF_06) + ((11003.0/480.0))*(CS_04*max_lambda_00 + CF_04));

       beta_1 = ((89.0/320.0))*((CS_05*max_lambda_00 + CF_05)*(CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_05*max_lambda_00 + CF_05) + ((2843.0/480.0))*(CS_04*max_lambda_00 +
            CF_04))*(CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(1261.0/240.0)*(CS_03*max_lambda_00 + CF_03) - (247.0/240.0)*(CS_05*max_lambda_00 + CF_05) +
            ((547.0/480.0))*(CS_02*max_lambda_00 + CF_02) + ((961.0/240.0))*(CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(CS_03*max_lambda_00 + CF_03)*(-(2983.0/240.0)*(CS_04*max_lambda_00 + CF_04) +
            ((267.0/80.0))*(CS_05*max_lambda_00 + CF_05) + ((3443.0/480.0))*(CS_03*max_lambda_00 + CF_03));

       beta_2 = ((547.0/960.0))*((CS_04*max_lambda_00 + CF_04)*(CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_04*max_lambda_00 + CF_04) + ((3443.0/480.0))*(CS_03*max_lambda_00 +
            CF_03))*(CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(CS_01*max_lambda_00 +
            CF_01)*(-(247.0/240.0)*(CS_04*max_lambda_00 + CF_04) + ((89.0/160.0))*(CS_01*max_lambda_00 + CF_01) +
            ((267.0/80.0))*(CS_03*max_lambda_00 + CF_03)) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(2983.0/240.0)*(CS_03*max_lambda_00 + CF_03) - (821.0/240.0)*(CS_01*max_lambda_00 + CF_01) +
            ((961.0/240.0))*(CS_04*max_lambda_00 + CF_04) + ((2843.0/480.0))*(CS_02*max_lambda_00 + CF_02));

       beta_3 = ((2107.0/960.0))*((CS_03*max_lambda_00 + CF_03)*(CS_03*max_lambda_00 + CF_03)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_03*max_lambda_00 + CF_03) + ((547.0/480.0))*(CS_00*max_lambda_00 +
            CF_00))*(CS_00*max_lambda_00 + CF_00) + ((1.0/2.0))*(CS_01*max_lambda_00 +
            CF_01)*(-(647.0/80.0)*(CS_00*max_lambda_00 + CF_00) + ((3521.0/240.0))*(CS_03*max_lambda_00 + CF_03) +
            ((7043.0/480.0))*(CS_01*max_lambda_00 + CF_01)) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(8623.0/240.0)*(CS_01*max_lambda_00 + CF_01) - (1567.0/80.0)*(CS_03*max_lambda_00 + CF_03) +
            ((2321.0/240.0))*(CS_00*max_lambda_00 + CF_00) + ((11003.0/480.0))*(CS_02*max_lambda_00 + CF_02));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj0 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_0 = (-(23.0/24.0)*(CS_02*max_lambda_00 + CF_02) - (1.0/8.0)*(CS_00*max_lambda_00 + CF_00) +
            ((13.0/24.0))*(CS_01*max_lambda_00 + CF_01) + ((25.0/24.0))*(CS_03*max_lambda_00 + CF_03))*omega_3 +
            (-(5.0/24.0)*(CS_02*max_lambda_00 + CF_02) + ((1.0/8.0))*(CS_04*max_lambda_00 + CF_04) +
            ((1.0/24.0))*(CS_01*max_lambda_00 + CF_01) + ((13.0/24.0))*(CS_03*max_lambda_00 + CF_03))*omega_2 +
            (-(5.0/24.0)*(CS_05*max_lambda_00 + CF_05) + ((1.0/8.0))*(CS_03*max_lambda_00 + CF_03) +
            ((1.0/24.0))*(CS_06*max_lambda_00 + CF_06) + ((13.0/24.0))*(CS_04*max_lambda_00 + CF_04))*omega_0 +
            (-(1.0/24.0)*(CS_02*max_lambda_00 + CF_02) - (1.0/24.0)*(CS_05*max_lambda_00 + CF_05) +
            ((7.0/24.0))*(CS_03*max_lambda_00 + CF_03) + ((7.0/24.0))*(CS_04*max_lambda_00 + CF_04))*omega_1 + Recon_0;

       beta_0 = ((547.0/960.0))*((-CS_07*max_lambda_00 + CF_07)*(-CS_07*max_lambda_00 + CF_07)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_07*max_lambda_00 + CF_07) + ((7043.0/480.0))*(-CS_06*max_lambda_00 +
            CF_06))*(-CS_06*max_lambda_00 + CF_06) + ((1.0/2.0))*(-CS_04*max_lambda_00 +
            CF_04)*(-(1567.0/80.0)*(-CS_05*max_lambda_00 + CF_05) - (309.0/80.0)*(-CS_07*max_lambda_00 + CF_07) +
            ((2107.0/480.0))*(-CS_04*max_lambda_00 + CF_04) + ((3521.0/240.0))*(-CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05)*(-(8623.0/240.0)*(-CS_06*max_lambda_00 + CF_06) +
            ((2321.0/240.0))*(-CS_07*max_lambda_00 + CF_07) + ((11003.0/480.0))*(-CS_05*max_lambda_00 + CF_05));

       beta_1 = ((89.0/320.0))*((-CS_06*max_lambda_00 + CF_06)*(-CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_06*max_lambda_00 + CF_06) + ((2843.0/480.0))*(-CS_05*max_lambda_00 +
            CF_05))*(-CS_05*max_lambda_00 + CF_05) + ((1.0/2.0))*(-CS_03*max_lambda_00 +
            CF_03)*(-(1261.0/240.0)*(-CS_04*max_lambda_00 + CF_04) - (247.0/240.0)*(-CS_06*max_lambda_00 + CF_06) +
            ((547.0/480.0))*(-CS_03*max_lambda_00 + CF_03) + ((961.0/240.0))*(-CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-CS_04*max_lambda_00 + CF_04)*(-(2983.0/240.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((267.0/80.0))*(-CS_06*max_lambda_00 + CF_06) + ((3443.0/480.0))*(-CS_04*max_lambda_00 + CF_04));

       beta_2 = ((547.0/960.0))*((-CS_05*max_lambda_00 + CF_05)*(-CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_05*max_lambda_00 + CF_05) + ((3443.0/480.0))*(-CS_04*max_lambda_00 +
            CF_04))*(-CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
            CF_02)*(-(821.0/240.0)*(-CS_03*max_lambda_00 + CF_03) - (247.0/240.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((89.0/160.0))*(-CS_02*max_lambda_00 + CF_02) + ((267.0/80.0))*(-CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03)*(-(2983.0/240.0)*(-CS_04*max_lambda_00 + CF_04) +
            ((961.0/240.0))*(-CS_05*max_lambda_00 + CF_05) + ((2843.0/480.0))*(-CS_03*max_lambda_00 + CF_03));

       beta_3 = ((2107.0/960.0))*((-CS_04*max_lambda_00 + CF_04)*(-CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_04*max_lambda_00 + CF_04) + ((11003.0/480.0))*(-CS_03*max_lambda_00 +
            CF_03))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_01*max_lambda_00 +
            CF_01)*(-(309.0/80.0)*(-CS_04*max_lambda_00 + CF_04) + ((547.0/480.0))*(-CS_01*max_lambda_00 + CF_01) +
            ((2321.0/240.0))*(-CS_03*max_lambda_00 + CF_03)) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
            CF_02)*(-(8623.0/240.0)*(-CS_03*max_lambda_00 + CF_03) - (647.0/80.0)*(-CS_01*max_lambda_00 + CF_01) +
            ((3521.0/240.0))*(-CS_04*max_lambda_00 + CF_04) + ((7043.0/480.0))*(-CS_02*max_lambda_00 + CF_02));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj0 = fmax(rj0, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_0 = (-(23.0/24.0)*(-CS_05*max_lambda_00 + CF_05) - (1.0/8.0)*(-CS_07*max_lambda_00 + CF_07) +
            ((13.0/24.0))*(-CS_06*max_lambda_00 + CF_06) + ((25.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_0 +
            (-(5.0/24.0)*(-CS_02*max_lambda_00 + CF_02) + ((1.0/8.0))*(-CS_04*max_lambda_00 + CF_04) +
            ((1.0/24.0))*(-CS_01*max_lambda_00 + CF_01) + ((13.0/24.0))*(-CS_03*max_lambda_00 + CF_03))*omega_3 +
            (-(5.0/24.0)*(-CS_05*max_lambda_00 + CF_05) + ((1.0/8.0))*(-CS_03*max_lambda_00 + CF_03) +
            ((1.0/24.0))*(-CS_06*max_lambda_00 + CF_06) + ((13.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_1 +
            (-(1.0/24.0)*(-CS_02*max_lambda_00 + CF_02) - (1.0/24.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((7.0/24.0))*(-CS_03*max_lambda_00 + CF_03) + ((7.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_2 +
            Recon_0;

       beta_0 = ((547.0/960.0))*((CS_16*max_lambda_11 + CF_16)*(CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_16*max_lambda_11 + CF_16) + ((7043.0/480.0))*(CS_15*max_lambda_11 +
            CF_15))*(CS_15*max_lambda_11 + CF_15) + ((1.0/2.0))*(CS_13*max_lambda_11 +
            CF_13)*(-(1567.0/80.0)*(CS_14*max_lambda_11 + CF_14) - (309.0/80.0)*(CS_16*max_lambda_11 + CF_16) +
            ((2107.0/480.0))*(CS_13*max_lambda_11 + CF_13) + ((3521.0/240.0))*(CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14)*(-(8623.0/240.0)*(CS_15*max_lambda_11 + CF_15) +
            ((2321.0/240.0))*(CS_16*max_lambda_11 + CF_16) + ((11003.0/480.0))*(CS_14*max_lambda_11 + CF_14));

       beta_1 = ((89.0/320.0))*((CS_15*max_lambda_11 + CF_15)*(CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_15*max_lambda_11 + CF_15) + ((2843.0/480.0))*(CS_14*max_lambda_11 +
            CF_14))*(CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(1261.0/240.0)*(CS_13*max_lambda_11 + CF_13) - (247.0/240.0)*(CS_15*max_lambda_11 + CF_15) +
            ((547.0/480.0))*(CS_12*max_lambda_11 + CF_12) + ((961.0/240.0))*(CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(CS_13*max_lambda_11 + CF_13)*(-(2983.0/240.0)*(CS_14*max_lambda_11 + CF_14) +
            ((267.0/80.0))*(CS_15*max_lambda_11 + CF_15) + ((3443.0/480.0))*(CS_13*max_lambda_11 + CF_13));

       beta_2 = ((547.0/960.0))*((CS_14*max_lambda_11 + CF_14)*(CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_14*max_lambda_11 + CF_14) + ((3443.0/480.0))*(CS_13*max_lambda_11 +
            CF_13))*(CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(CS_11*max_lambda_11 +
            CF_11)*(-(247.0/240.0)*(CS_14*max_lambda_11 + CF_14) + ((89.0/160.0))*(CS_11*max_lambda_11 + CF_11) +
            ((267.0/80.0))*(CS_13*max_lambda_11 + CF_13)) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(2983.0/240.0)*(CS_13*max_lambda_11 + CF_13) - (821.0/240.0)*(CS_11*max_lambda_11 + CF_11) +
            ((961.0/240.0))*(CS_14*max_lambda_11 + CF_14) + ((2843.0/480.0))*(CS_12*max_lambda_11 + CF_12));

       beta_3 = ((2107.0/960.0))*((CS_13*max_lambda_11 + CF_13)*(CS_13*max_lambda_11 + CF_13)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_13*max_lambda_11 + CF_13) + ((547.0/480.0))*(CS_10*max_lambda_11 +
            CF_10))*(CS_10*max_lambda_11 + CF_10) + ((1.0/2.0))*(CS_11*max_lambda_11 +
            CF_11)*(-(647.0/80.0)*(CS_10*max_lambda_11 + CF_10) + ((3521.0/240.0))*(CS_13*max_lambda_11 + CF_13) +
            ((7043.0/480.0))*(CS_11*max_lambda_11 + CF_11)) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(8623.0/240.0)*(CS_11*max_lambda_11 + CF_11) - (1567.0/80.0)*(CS_13*max_lambda_11 + CF_13) +
            ((2321.0/240.0))*(CS_10*max_lambda_11 + CF_10) + ((11003.0/480.0))*(CS_12*max_lambda_11 + CF_12));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj1 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_1 = (-(23.0/24.0)*(CS_12*max_lambda_11 + CF_12) - (1.0/8.0)*(CS_10*max_lambda_11 + CF_10) +
            ((13.0/24.0))*(CS_11*max_lambda_11 + CF_11) + ((25.0/24.0))*(CS_13*max_lambda_11 + CF_13))*omega_3 +
            (-(5.0/24.0)*(CS_12*max_lambda_11 + CF_12) + ((1.0/8.0))*(CS_14*max_lambda_11 + CF_14) +
            ((1.0/24.0))*(CS_11*max_lambda_11 + CF_11) + ((13.0/24.0))*(CS_13*max_lambda_11 + CF_13))*omega_2 +
            (-(5.0/24.0)*(CS_15*max_lambda_11 + CF_15) + ((1.0/8.0))*(CS_13*max_lambda_11 + CF_13) +
            ((1.0/24.0))*(CS_16*max_lambda_11 + CF_16) + ((13.0/24.0))*(CS_14*max_lambda_11 + CF_14))*omega_0 +
            (-(1.0/24.0)*(CS_12*max_lambda_11 + CF_12) - (1.0/24.0)*(CS_15*max_lambda_11 + CF_15) +
            ((7.0/24.0))*(CS_13*max_lambda_11 + CF_13) + ((7.0/24.0))*(CS_14*max_lambda_11 + CF_14))*omega_1 + Recon_1;

       beta_0 = ((547.0/960.0))*((-CS_17*max_lambda_11 + CF_17)*(-CS_17*max_lambda_11 + CF_17)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_17*max_lambda_11 + CF_17) + ((7043.0/480.0))*(-CS_16*max_lambda_11 +
            CF_16))*(-CS_16*max_lambda_11 + CF_16) + ((1.0/2.0))*(-CS_14*max_lambda_11 +
            CF_14)*(-(1567.0/80.0)*(-CS_15*max_lambda_11 + CF_15) - (309.0/80.0)*(-CS_17*max_lambda_11 + CF_17) +
            ((2107.0/480.0))*(-CS_14*max_lambda_11 + CF_14) + ((3521.0/240.0))*(-CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15)*(-(8623.0/240.0)*(-CS_16*max_lambda_11 + CF_16) +
            ((2321.0/240.0))*(-CS_17*max_lambda_11 + CF_17) + ((11003.0/480.0))*(-CS_15*max_lambda_11 + CF_15));

       beta_1 = ((89.0/320.0))*((-CS_16*max_lambda_11 + CF_16)*(-CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_16*max_lambda_11 + CF_16) + ((2843.0/480.0))*(-CS_15*max_lambda_11 +
            CF_15))*(-CS_15*max_lambda_11 + CF_15) + ((1.0/2.0))*(-CS_13*max_lambda_11 +
            CF_13)*(-(1261.0/240.0)*(-CS_14*max_lambda_11 + CF_14) - (247.0/240.0)*(-CS_16*max_lambda_11 + CF_16) +
            ((547.0/480.0))*(-CS_13*max_lambda_11 + CF_13) + ((961.0/240.0))*(-CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-CS_14*max_lambda_11 + CF_14)*(-(2983.0/240.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((267.0/80.0))*(-CS_16*max_lambda_11 + CF_16) + ((3443.0/480.0))*(-CS_14*max_lambda_11 + CF_14));

       beta_2 = ((547.0/960.0))*((-CS_15*max_lambda_11 + CF_15)*(-CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_15*max_lambda_11 + CF_15) + ((3443.0/480.0))*(-CS_14*max_lambda_11 +
            CF_14))*(-CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
            CF_12)*(-(821.0/240.0)*(-CS_13*max_lambda_11 + CF_13) - (247.0/240.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((89.0/160.0))*(-CS_12*max_lambda_11 + CF_12) + ((267.0/80.0))*(-CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13)*(-(2983.0/240.0)*(-CS_14*max_lambda_11 + CF_14) +
            ((961.0/240.0))*(-CS_15*max_lambda_11 + CF_15) + ((2843.0/480.0))*(-CS_13*max_lambda_11 + CF_13));

       beta_3 = ((2107.0/960.0))*((-CS_14*max_lambda_11 + CF_14)*(-CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_14*max_lambda_11 + CF_14) + ((11003.0/480.0))*(-CS_13*max_lambda_11 +
            CF_13))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_11*max_lambda_11 +
            CF_11)*(-(309.0/80.0)*(-CS_14*max_lambda_11 + CF_14) + ((547.0/480.0))*(-CS_11*max_lambda_11 + CF_11) +
            ((2321.0/240.0))*(-CS_13*max_lambda_11 + CF_13)) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
            CF_12)*(-(8623.0/240.0)*(-CS_13*max_lambda_11 + CF_13) - (647.0/80.0)*(-CS_11*max_lambda_11 + CF_11) +
            ((3521.0/240.0))*(-CS_14*max_lambda_11 + CF_14) + ((7043.0/480.0))*(-CS_12*max_lambda_11 + CF_12));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj1 = fmax(rj1, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_1 = (-(23.0/24.0)*(-CS_15*max_lambda_11 + CF_15) - (1.0/8.0)*(-CS_17*max_lambda_11 + CF_17) +
            ((13.0/24.0))*(-CS_16*max_lambda_11 + CF_16) + ((25.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_0 +
            (-(5.0/24.0)*(-CS_12*max_lambda_11 + CF_12) + ((1.0/8.0))*(-CS_14*max_lambda_11 + CF_14) +
            ((1.0/24.0))*(-CS_11*max_lambda_11 + CF_11) + ((13.0/24.0))*(-CS_13*max_lambda_11 + CF_13))*omega_3 +
            (-(5.0/24.0)*(-CS_15*max_lambda_11 + CF_15) + ((1.0/8.0))*(-CS_13*max_lambda_11 + CF_13) +
            ((1.0/24.0))*(-CS_16*max_lambda_11 + CF_16) + ((13.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_1 +
            (-(1.0/24.0)*(-CS_12*max_lambda_11 + CF_12) - (1.0/24.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((7.0/24.0))*(-CS_13*max_lambda_11 + CF_13) + ((7.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_2 +
            Recon_1;

       beta_0 = ((547.0/960.0))*((CS_26*max_lambda_22 + CF_26)*(CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_26*max_lambda_22 + CF_26) + ((7043.0/480.0))*(CS_25*max_lambda_22 +
            CF_25))*(CS_25*max_lambda_22 + CF_25) + ((1.0/2.0))*(CS_23*max_lambda_22 +
            CF_23)*(-(1567.0/80.0)*(CS_24*max_lambda_22 + CF_24) - (309.0/80.0)*(CS_26*max_lambda_22 + CF_26) +
            ((2107.0/480.0))*(CS_23*max_lambda_22 + CF_23) + ((3521.0/240.0))*(CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24)*(-(8623.0/240.0)*(CS_25*max_lambda_22 + CF_25) +
            ((2321.0/240.0))*(CS_26*max_lambda_22 + CF_26) + ((11003.0/480.0))*(CS_24*max_lambda_22 + CF_24));

       beta_1 = ((89.0/320.0))*((CS_25*max_lambda_22 + CF_25)*(CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_25*max_lambda_22 + CF_25) + ((2843.0/480.0))*(CS_24*max_lambda_22 +
            CF_24))*(CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(1261.0/240.0)*(CS_23*max_lambda_22 + CF_23) - (247.0/240.0)*(CS_25*max_lambda_22 + CF_25) +
            ((547.0/480.0))*(CS_22*max_lambda_22 + CF_22) + ((961.0/240.0))*(CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(CS_23*max_lambda_22 + CF_23)*(-(2983.0/240.0)*(CS_24*max_lambda_22 + CF_24) +
            ((267.0/80.0))*(CS_25*max_lambda_22 + CF_25) + ((3443.0/480.0))*(CS_23*max_lambda_22 + CF_23));

       beta_2 = ((547.0/960.0))*((CS_24*max_lambda_22 + CF_24)*(CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_24*max_lambda_22 + CF_24) + ((3443.0/480.0))*(CS_23*max_lambda_22 +
            CF_23))*(CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(CS_21*max_lambda_22 +
            CF_21)*(-(247.0/240.0)*(CS_24*max_lambda_22 + CF_24) + ((89.0/160.0))*(CS_21*max_lambda_22 + CF_21) +
            ((267.0/80.0))*(CS_23*max_lambda_22 + CF_23)) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(2983.0/240.0)*(CS_23*max_lambda_22 + CF_23) - (821.0/240.0)*(CS_21*max_lambda_22 + CF_21) +
            ((961.0/240.0))*(CS_24*max_lambda_22 + CF_24) + ((2843.0/480.0))*(CS_22*max_lambda_22 + CF_22));

       beta_3 = ((2107.0/960.0))*((CS_23*max_lambda_22 + CF_23)*(CS_23*max_lambda_22 + CF_23)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_23*max_lambda_22 + CF_23) + ((547.0/480.0))*(CS_20*max_lambda_22 +
            CF_20))*(CS_20*max_lambda_22 + CF_20) + ((1.0/2.0))*(CS_21*max_lambda_22 +
            CF_21)*(-(647.0/80.0)*(CS_20*max_lambda_22 + CF_20) + ((3521.0/240.0))*(CS_23*max_lambda_22 + CF_23) +
            ((7043.0/480.0))*(CS_21*max_lambda_22 + CF_21)) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(8623.0/240.0)*(CS_21*max_lambda_22 + CF_21) - (1567.0/80.0)*(CS_23*max_lambda_22 + CF_23) +
            ((2321.0/240.0))*(CS_20*max_lambda_22 + CF_20) + ((11003.0/480.0))*(CS_22*max_lambda_22 + CF_22));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj2 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_2 = (-(23.0/24.0)*(CS_22*max_lambda_22 + CF_22) - (1.0/8.0)*(CS_20*max_lambda_22 + CF_20) +
            ((13.0/24.0))*(CS_21*max_lambda_22 + CF_21) + ((25.0/24.0))*(CS_23*max_lambda_22 + CF_23))*omega_3 +
            (-(5.0/24.0)*(CS_22*max_lambda_22 + CF_22) + ((1.0/8.0))*(CS_24*max_lambda_22 + CF_24) +
            ((1.0/24.0))*(CS_21*max_lambda_22 + CF_21) + ((13.0/24.0))*(CS_23*max_lambda_22 + CF_23))*omega_2 +
            (-(5.0/24.0)*(CS_25*max_lambda_22 + CF_25) + ((1.0/8.0))*(CS_23*max_lambda_22 + CF_23) +
            ((1.0/24.0))*(CS_26*max_lambda_22 + CF_26) + ((13.0/24.0))*(CS_24*max_lambda_22 + CF_24))*omega_0 +
            (-(1.0/24.0)*(CS_22*max_lambda_22 + CF_22) - (1.0/24.0)*(CS_25*max_lambda_22 + CF_25) +
            ((7.0/24.0))*(CS_23*max_lambda_22 + CF_23) + ((7.0/24.0))*(CS_24*max_lambda_22 + CF_24))*omega_1 + Recon_2;

       beta_0 = ((547.0/960.0))*((-CS_27*max_lambda_22 + CF_27)*(-CS_27*max_lambda_22 + CF_27)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_27*max_lambda_22 + CF_27) + ((7043.0/480.0))*(-CS_26*max_lambda_22 +
            CF_26))*(-CS_26*max_lambda_22 + CF_26) + ((1.0/2.0))*(-CS_24*max_lambda_22 +
            CF_24)*(-(1567.0/80.0)*(-CS_25*max_lambda_22 + CF_25) - (309.0/80.0)*(-CS_27*max_lambda_22 + CF_27) +
            ((2107.0/480.0))*(-CS_24*max_lambda_22 + CF_24) + ((3521.0/240.0))*(-CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25)*(-(8623.0/240.0)*(-CS_26*max_lambda_22 + CF_26) +
            ((2321.0/240.0))*(-CS_27*max_lambda_22 + CF_27) + ((11003.0/480.0))*(-CS_25*max_lambda_22 + CF_25));

       beta_1 = ((89.0/320.0))*((-CS_26*max_lambda_22 + CF_26)*(-CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_26*max_lambda_22 + CF_26) + ((2843.0/480.0))*(-CS_25*max_lambda_22 +
            CF_25))*(-CS_25*max_lambda_22 + CF_25) + ((1.0/2.0))*(-CS_23*max_lambda_22 +
            CF_23)*(-(1261.0/240.0)*(-CS_24*max_lambda_22 + CF_24) - (247.0/240.0)*(-CS_26*max_lambda_22 + CF_26) +
            ((547.0/480.0))*(-CS_23*max_lambda_22 + CF_23) + ((961.0/240.0))*(-CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-CS_24*max_lambda_22 + CF_24)*(-(2983.0/240.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((267.0/80.0))*(-CS_26*max_lambda_22 + CF_26) + ((3443.0/480.0))*(-CS_24*max_lambda_22 + CF_24));

       beta_2 = ((547.0/960.0))*((-CS_25*max_lambda_22 + CF_25)*(-CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_25*max_lambda_22 + CF_25) + ((3443.0/480.0))*(-CS_24*max_lambda_22 +
            CF_24))*(-CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
            CF_22)*(-(821.0/240.0)*(-CS_23*max_lambda_22 + CF_23) - (247.0/240.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((89.0/160.0))*(-CS_22*max_lambda_22 + CF_22) + ((267.0/80.0))*(-CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23)*(-(2983.0/240.0)*(-CS_24*max_lambda_22 + CF_24) +
            ((961.0/240.0))*(-CS_25*max_lambda_22 + CF_25) + ((2843.0/480.0))*(-CS_23*max_lambda_22 + CF_23));

       beta_3 = ((2107.0/960.0))*((-CS_24*max_lambda_22 + CF_24)*(-CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_24*max_lambda_22 + CF_24) + ((11003.0/480.0))*(-CS_23*max_lambda_22 +
            CF_23))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_21*max_lambda_22 +
            CF_21)*(-(309.0/80.0)*(-CS_24*max_lambda_22 + CF_24) + ((547.0/480.0))*(-CS_21*max_lambda_22 + CF_21) +
            ((2321.0/240.0))*(-CS_23*max_lambda_22 + CF_23)) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
            CF_22)*(-(8623.0/240.0)*(-CS_23*max_lambda_22 + CF_23) - (647.0/80.0)*(-CS_21*max_lambda_22 + CF_21) +
            ((3521.0/240.0))*(-CS_24*max_lambda_22 + CF_24) + ((7043.0/480.0))*(-CS_22*max_lambda_22 + CF_22));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj2 = fmax(rj2, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_2 = (-(23.0/24.0)*(-CS_25*max_lambda_22 + CF_25) - (1.0/8.0)*(-CS_27*max_lambda_22 + CF_27) +
            ((13.0/24.0))*(-CS_26*max_lambda_22 + CF_26) + ((25.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_0 +
            (-(5.0/24.0)*(-CS_22*max_lambda_22 + CF_22) + ((1.0/8.0))*(-CS_24*max_lambda_22 + CF_24) +
            ((1.0/24.0))*(-CS_21*max_lambda_22 + CF_21) + ((13.0/24.0))*(-CS_23*max_lambda_22 + CF_23))*omega_3 +
            (-(5.0/24.0)*(-CS_25*max_lambda_22 + CF_25) + ((1.0/8.0))*(-CS_23*max_lambda_22 + CF_23) +
            ((1.0/24.0))*(-CS_26*max_lambda_22 + CF_26) + ((13.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_1 +
            (-(1.0/24.0)*(-CS_22*max_lambda_22 + CF_22) - (1.0/24.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((7.0/24.0))*(-CS_23*max_lambda_22 + CF_23) + ((7.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_2 +
            Recon_2;

       beta_0 = ((547.0/960.0))*((CS_36*max_lambda_33 + CF_36)*(CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_36*max_lambda_33 + CF_36) + ((7043.0/480.0))*(CS_35*max_lambda_33 +
            CF_35))*(CS_35*max_lambda_33 + CF_35) + ((1.0/2.0))*(CS_33*max_lambda_33 +
            CF_33)*(-(1567.0/80.0)*(CS_34*max_lambda_33 + CF_34) - (309.0/80.0)*(CS_36*max_lambda_33 + CF_36) +
            ((2107.0/480.0))*(CS_33*max_lambda_33 + CF_33) + ((3521.0/240.0))*(CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34)*(-(8623.0/240.0)*(CS_35*max_lambda_33 + CF_35) +
            ((2321.0/240.0))*(CS_36*max_lambda_33 + CF_36) + ((11003.0/480.0))*(CS_34*max_lambda_33 + CF_34));

       beta_1 = ((89.0/320.0))*((CS_35*max_lambda_33 + CF_35)*(CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_35*max_lambda_33 + CF_35) + ((2843.0/480.0))*(CS_34*max_lambda_33 +
            CF_34))*(CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(1261.0/240.0)*(CS_33*max_lambda_33 + CF_33) - (247.0/240.0)*(CS_35*max_lambda_33 + CF_35) +
            ((547.0/480.0))*(CS_32*max_lambda_33 + CF_32) + ((961.0/240.0))*(CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(CS_33*max_lambda_33 + CF_33)*(-(2983.0/240.0)*(CS_34*max_lambda_33 + CF_34) +
            ((267.0/80.0))*(CS_35*max_lambda_33 + CF_35) + ((3443.0/480.0))*(CS_33*max_lambda_33 + CF_33));

       beta_2 = ((547.0/960.0))*((CS_34*max_lambda_33 + CF_34)*(CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_34*max_lambda_33 + CF_34) + ((3443.0/480.0))*(CS_33*max_lambda_33 +
            CF_33))*(CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(CS_31*max_lambda_33 +
            CF_31)*(-(247.0/240.0)*(CS_34*max_lambda_33 + CF_34) + ((89.0/160.0))*(CS_31*max_lambda_33 + CF_31) +
            ((267.0/80.0))*(CS_33*max_lambda_33 + CF_33)) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(2983.0/240.0)*(CS_33*max_lambda_33 + CF_33) - (821.0/240.0)*(CS_31*max_lambda_33 + CF_31) +
            ((961.0/240.0))*(CS_34*max_lambda_33 + CF_34) + ((2843.0/480.0))*(CS_32*max_lambda_33 + CF_32));

       beta_3 = ((2107.0/960.0))*((CS_33*max_lambda_33 + CF_33)*(CS_33*max_lambda_33 + CF_33)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_33*max_lambda_33 + CF_33) + ((547.0/480.0))*(CS_30*max_lambda_33 +
            CF_30))*(CS_30*max_lambda_33 + CF_30) + ((1.0/2.0))*(CS_31*max_lambda_33 +
            CF_31)*(-(647.0/80.0)*(CS_30*max_lambda_33 + CF_30) + ((3521.0/240.0))*(CS_33*max_lambda_33 + CF_33) +
            ((7043.0/480.0))*(CS_31*max_lambda_33 + CF_31)) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(8623.0/240.0)*(CS_31*max_lambda_33 + CF_31) - (1567.0/80.0)*(CS_33*max_lambda_33 + CF_33) +
            ((2321.0/240.0))*(CS_30*max_lambda_33 + CF_30) + ((11003.0/480.0))*(CS_32*max_lambda_33 + CF_32));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj3 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_3 = (-(23.0/24.0)*(CS_32*max_lambda_33 + CF_32) - (1.0/8.0)*(CS_30*max_lambda_33 + CF_30) +
            ((13.0/24.0))*(CS_31*max_lambda_33 + CF_31) + ((25.0/24.0))*(CS_33*max_lambda_33 + CF_33))*omega_3 +
            (-(5.0/24.0)*(CS_32*max_lambda_33 + CF_32) + ((1.0/8.0))*(CS_34*max_lambda_33 + CF_34) +
            ((1.0/24.0))*(CS_31*max_lambda_33 + CF_31) + ((13.0/24.0))*(CS_33*max_lambda_33 + CF_33))*omega_2 +
            (-(5.0/24.0)*(CS_35*max_lambda_33 + CF_35) + ((1.0/8.0))*(CS_33*max_lambda_33 + CF_33) +
            ((1.0/24.0))*(CS_36*max_lambda_33 + CF_36) + ((13.0/24.0))*(CS_34*max_lambda_33 + CF_34))*omega_0 +
            (-(1.0/24.0)*(CS_32*max_lambda_33 + CF_32) - (1.0/24.0)*(CS_35*max_lambda_33 + CF_35) +
            ((7.0/24.0))*(CS_33*max_lambda_33 + CF_33) + ((7.0/24.0))*(CS_34*max_lambda_33 + CF_34))*omega_1 + Recon_3;

       beta_0 = ((547.0/960.0))*((-CS_37*max_lambda_33 + CF_37)*(-CS_37*max_lambda_33 + CF_37)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_37*max_lambda_33 + CF_37) + ((7043.0/480.0))*(-CS_36*max_lambda_33 +
            CF_36))*(-CS_36*max_lambda_33 + CF_36) + ((1.0/2.0))*(-CS_34*max_lambda_33 +
            CF_34)*(-(1567.0/80.0)*(-CS_35*max_lambda_33 + CF_35) - (309.0/80.0)*(-CS_37*max_lambda_33 + CF_37) +
            ((2107.0/480.0))*(-CS_34*max_lambda_33 + CF_34) + ((3521.0/240.0))*(-CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35)*(-(8623.0/240.0)*(-CS_36*max_lambda_33 + CF_36) +
            ((2321.0/240.0))*(-CS_37*max_lambda_33 + CF_37) + ((11003.0/480.0))*(-CS_35*max_lambda_33 + CF_35));

       beta_1 = ((89.0/320.0))*((-CS_36*max_lambda_33 + CF_36)*(-CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_36*max_lambda_33 + CF_36) + ((2843.0/480.0))*(-CS_35*max_lambda_33 +
            CF_35))*(-CS_35*max_lambda_33 + CF_35) + ((1.0/2.0))*(-CS_33*max_lambda_33 +
            CF_33)*(-(1261.0/240.0)*(-CS_34*max_lambda_33 + CF_34) - (247.0/240.0)*(-CS_36*max_lambda_33 + CF_36) +
            ((547.0/480.0))*(-CS_33*max_lambda_33 + CF_33) + ((961.0/240.0))*(-CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-CS_34*max_lambda_33 + CF_34)*(-(2983.0/240.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((267.0/80.0))*(-CS_36*max_lambda_33 + CF_36) + ((3443.0/480.0))*(-CS_34*max_lambda_33 + CF_34));

       beta_2 = ((547.0/960.0))*((-CS_35*max_lambda_33 + CF_35)*(-CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_35*max_lambda_33 + CF_35) + ((3443.0/480.0))*(-CS_34*max_lambda_33 +
            CF_34))*(-CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
            CF_32)*(-(821.0/240.0)*(-CS_33*max_lambda_33 + CF_33) - (247.0/240.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((89.0/160.0))*(-CS_32*max_lambda_33 + CF_32) + ((267.0/80.0))*(-CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33)*(-(2983.0/240.0)*(-CS_34*max_lambda_33 + CF_34) +
            ((961.0/240.0))*(-CS_35*max_lambda_33 + CF_35) + ((2843.0/480.0))*(-CS_33*max_lambda_33 + CF_33));

       beta_3 = ((2107.0/960.0))*((-CS_34*max_lambda_33 + CF_34)*(-CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_34*max_lambda_33 + CF_34) + ((11003.0/480.0))*(-CS_33*max_lambda_33 +
            CF_33))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_31*max_lambda_33 +
            CF_31)*(-(309.0/80.0)*(-CS_34*max_lambda_33 + CF_34) + ((547.0/480.0))*(-CS_31*max_lambda_33 + CF_31) +
            ((2321.0/240.0))*(-CS_33*max_lambda_33 + CF_33)) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
            CF_32)*(-(8623.0/240.0)*(-CS_33*max_lambda_33 + CF_33) - (647.0/80.0)*(-CS_31*max_lambda_33 + CF_31) +
            ((3521.0/240.0))*(-CS_34*max_lambda_33 + CF_34) + ((7043.0/480.0))*(-CS_32*max_lambda_33 + CF_32));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj3 = fmax(rj3, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_3 = (-(23.0/24.0)*(-CS_35*max_lambda_33 + CF_35) - (1.0/8.0)*(-CS_37*max_lambda_33 + CF_37) +
            ((13.0/24.0))*(-CS_36*max_lambda_33 + CF_36) + ((25.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_0 +
            (-(5.0/24.0)*(-CS_32*max_lambda_33 + CF_32) + ((1.0/8.0))*(-CS_34*max_lambda_33 + CF_34) +
            ((1.0/24.0))*(-CS_31*max_lambda_33 + CF_31) + ((13.0/24.0))*(-CS_33*max_lambda_33 + CF_33))*omega_3 +
            (-(5.0/24.0)*(-CS_35*max_lambda_33 + CF_35) + ((1.0/8.0))*(-CS_33*max_lambda_33 + CF_33) +
            ((1.0/24.0))*(-CS_36*max_lambda_33 + CF_36) + ((13.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_1 +
            (-(1.0/24.0)*(-CS_32*max_lambda_33 + CF_32) - (1.0/24.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((7.0/24.0))*(-CS_33*max_lambda_33 + CF_33) + ((7.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_2 +
            Recon_3;

       beta_0 = ((547.0/960.0))*((CS_46*max_lambda_44 + CF_46)*(CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_46*max_lambda_44 + CF_46) + ((7043.0/480.0))*(CS_45*max_lambda_44 +
            CF_45))*(CS_45*max_lambda_44 + CF_45) + ((1.0/2.0))*(CS_43*max_lambda_44 +
            CF_43)*(-(1567.0/80.0)*(CS_44*max_lambda_44 + CF_44) - (309.0/80.0)*(CS_46*max_lambda_44 + CF_46) +
            ((2107.0/480.0))*(CS_43*max_lambda_44 + CF_43) + ((3521.0/240.0))*(CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44)*(-(8623.0/240.0)*(CS_45*max_lambda_44 + CF_45) +
            ((2321.0/240.0))*(CS_46*max_lambda_44 + CF_46) + ((11003.0/480.0))*(CS_44*max_lambda_44 + CF_44));

       beta_1 = ((89.0/320.0))*((CS_45*max_lambda_44 + CF_45)*(CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_45*max_lambda_44 + CF_45) + ((2843.0/480.0))*(CS_44*max_lambda_44 +
            CF_44))*(CS_44*max_lambda_44 + CF_44) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(1261.0/240.0)*(CS_43*max_lambda_44 + CF_43) - (247.0/240.0)*(CS_45*max_lambda_44 + CF_45) +
            ((547.0/480.0))*(CS_42*max_lambda_44 + CF_42) + ((961.0/240.0))*(CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(CS_43*max_lambda_44 + CF_43)*(-(2983.0/240.0)*(CS_44*max_lambda_44 + CF_44) +
            ((267.0/80.0))*(CS_45*max_lambda_44 + CF_45) + ((3443.0/480.0))*(CS_43*max_lambda_44 + CF_43));

       beta_2 = ((547.0/960.0))*((CS_44*max_lambda_44 + CF_44)*(CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_44*max_lambda_44 + CF_44) + ((3443.0/480.0))*(CS_43*max_lambda_44 +
            CF_43))*(CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(CS_41*max_lambda_44 +
            CF_41)*(-(247.0/240.0)*(CS_44*max_lambda_44 + CF_44) + ((89.0/160.0))*(CS_41*max_lambda_44 + CF_41) +
            ((267.0/80.0))*(CS_43*max_lambda_44 + CF_43)) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(2983.0/240.0)*(CS_43*max_lambda_44 + CF_43) - (821.0/240.0)*(CS_41*max_lambda_44 + CF_41) +
            ((961.0/240.0))*(CS_44*max_lambda_44 + CF_44) + ((2843.0/480.0))*(CS_42*max_lambda_44 + CF_42));

       beta_3 = ((2107.0/960.0))*((CS_43*max_lambda_44 + CF_43)*(CS_43*max_lambda_44 + CF_43)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_43*max_lambda_44 + CF_43) + ((547.0/480.0))*(CS_40*max_lambda_44 +
            CF_40))*(CS_40*max_lambda_44 + CF_40) + ((1.0/2.0))*(CS_41*max_lambda_44 +
            CF_41)*(-(647.0/80.0)*(CS_40*max_lambda_44 + CF_40) + ((3521.0/240.0))*(CS_43*max_lambda_44 + CF_43) +
            ((7043.0/480.0))*(CS_41*max_lambda_44 + CF_41)) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(8623.0/240.0)*(CS_41*max_lambda_44 + CF_41) - (1567.0/80.0)*(CS_43*max_lambda_44 + CF_43) +
            ((2321.0/240.0))*(CS_40*max_lambda_44 + CF_40) + ((11003.0/480.0))*(CS_42*max_lambda_44 + CF_42));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj4 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_4 = (-(23.0/24.0)*(CS_42*max_lambda_44 + CF_42) - (1.0/8.0)*(CS_40*max_lambda_44 + CF_40) +
            ((13.0/24.0))*(CS_41*max_lambda_44 + CF_41) + ((25.0/24.0))*(CS_43*max_lambda_44 + CF_43))*omega_3 +
            (-(5.0/24.0)*(CS_42*max_lambda_44 + CF_42) + ((1.0/8.0))*(CS_44*max_lambda_44 + CF_44) +
            ((1.0/24.0))*(CS_41*max_lambda_44 + CF_41) + ((13.0/24.0))*(CS_43*max_lambda_44 + CF_43))*omega_2 +
            (-(5.0/24.0)*(CS_45*max_lambda_44 + CF_45) + ((1.0/8.0))*(CS_43*max_lambda_44 + CF_43) +
            ((1.0/24.0))*(CS_46*max_lambda_44 + CF_46) + ((13.0/24.0))*(CS_44*max_lambda_44 + CF_44))*omega_0 +
            (-(1.0/24.0)*(CS_42*max_lambda_44 + CF_42) - (1.0/24.0)*(CS_45*max_lambda_44 + CF_45) +
            ((7.0/24.0))*(CS_43*max_lambda_44 + CF_43) + ((7.0/24.0))*(CS_44*max_lambda_44 + CF_44))*omega_1 + Recon_4;

       beta_0 = ((547.0/960.0))*((-CS_47*max_lambda_44 + CF_47)*(-CS_47*max_lambda_44 + CF_47)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_47*max_lambda_44 + CF_47) + ((7043.0/480.0))*(-CS_46*max_lambda_44 +
            CF_46))*(-CS_46*max_lambda_44 + CF_46) + ((1.0/2.0))*(-CS_44*max_lambda_44 +
            CF_44)*(-(1567.0/80.0)*(-CS_45*max_lambda_44 + CF_45) - (309.0/80.0)*(-CS_47*max_lambda_44 + CF_47) +
            ((2107.0/480.0))*(-CS_44*max_lambda_44 + CF_44) + ((3521.0/240.0))*(-CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45)*(-(8623.0/240.0)*(-CS_46*max_lambda_44 + CF_46) +
            ((2321.0/240.0))*(-CS_47*max_lambda_44 + CF_47) + ((11003.0/480.0))*(-CS_45*max_lambda_44 + CF_45));

       beta_1 = ((89.0/320.0))*((-CS_46*max_lambda_44 + CF_46)*(-CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_46*max_lambda_44 + CF_46) + ((2843.0/480.0))*(-CS_45*max_lambda_44 +
            CF_45))*(-CS_45*max_lambda_44 + CF_45) + ((1.0/2.0))*(-CS_43*max_lambda_44 +
            CF_43)*(-(1261.0/240.0)*(-CS_44*max_lambda_44 + CF_44) - (247.0/240.0)*(-CS_46*max_lambda_44 + CF_46) +
            ((547.0/480.0))*(-CS_43*max_lambda_44 + CF_43) + ((961.0/240.0))*(-CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-CS_44*max_lambda_44 + CF_44)*(-(2983.0/240.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((267.0/80.0))*(-CS_46*max_lambda_44 + CF_46) + ((3443.0/480.0))*(-CS_44*max_lambda_44 + CF_44));

       beta_2 = ((547.0/960.0))*((-CS_45*max_lambda_44 + CF_45)*(-CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_45*max_lambda_44 + CF_45) + ((3443.0/480.0))*(-CS_44*max_lambda_44 +
            CF_44))*(-CS_44*max_lambda_44 + CF_44) + ((1.0/2.0))*(-CS_42*max_lambda_44 +
            CF_42)*(-(821.0/240.0)*(-CS_43*max_lambda_44 + CF_43) - (247.0/240.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((89.0/160.0))*(-CS_42*max_lambda_44 + CF_42) + ((267.0/80.0))*(-CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43)*(-(2983.0/240.0)*(-CS_44*max_lambda_44 + CF_44) +
            ((961.0/240.0))*(-CS_45*max_lambda_44 + CF_45) + ((2843.0/480.0))*(-CS_43*max_lambda_44 + CF_43));

       beta_3 = ((2107.0/960.0))*((-CS_44*max_lambda_44 + CF_44)*(-CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_44*max_lambda_44 + CF_44) + ((11003.0/480.0))*(-CS_43*max_lambda_44 +
            CF_43))*(-CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(-CS_41*max_lambda_44 +
            CF_41)*(-(309.0/80.0)*(-CS_44*max_lambda_44 + CF_44) + ((547.0/480.0))*(-CS_41*max_lambda_44 + CF_41) +
            ((2321.0/240.0))*(-CS_43*max_lambda_44 + CF_43)) + ((1.0/2.0))*(-CS_42*max_lambda_44 +
            CF_42)*(-(8623.0/240.0)*(-CS_43*max_lambda_44 + CF_43) - (647.0/80.0)*(-CS_41*max_lambda_44 + CF_41) +
            ((3521.0/240.0))*(-CS_44*max_lambda_44 + CF_44) + ((7043.0/480.0))*(-CS_42*max_lambda_44 + CF_42));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj4 = fmax(rj4, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_4 = (-(23.0/24.0)*(-CS_45*max_lambda_44 + CF_45) - (1.0/8.0)*(-CS_47*max_lambda_44 + CF_47) +
            ((13.0/24.0))*(-CS_46*max_lambda_44 + CF_46) + ((25.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_0 +
            (-(5.0/24.0)*(-CS_42*max_lambda_44 + CF_42) + ((1.0/8.0))*(-CS_44*max_lambda_44 + CF_44) +
            ((1.0/24.0))*(-CS_41*max_lambda_44 + CF_41) + ((13.0/24.0))*(-CS_43*max_lambda_44 + CF_43))*omega_3 +
            (-(5.0/24.0)*(-CS_45*max_lambda_44 + CF_45) + ((1.0/8.0))*(-CS_43*max_lambda_44 + CF_43) +
            ((1.0/24.0))*(-CS_46*max_lambda_44 + CF_46) + ((13.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_1 +
            (-(1.0/24.0)*(-CS_42*max_lambda_44 + CF_42) - (1.0/24.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((7.0/24.0))*(-CS_43*max_lambda_44 + CF_43) + ((7.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_2 +
            Recon_4;

       Recon_0 = (((1.0/840.0))*(-533*CF_03 - 533*CF_04 - 29*CF_01 - 29*CF_06 + 3*CF_00 + 3*CF_07 + 139*CF_02 +
            139*CF_05) + Recon_0)*rj0;

       Recon_1 = (((1.0/840.0))*(-533*CF_13 - 533*CF_14 - 29*CF_11 - 29*CF_16 + 3*CF_10 + 3*CF_17 + 139*CF_12 +
            139*CF_15) + Recon_1)*rj1;

       Recon_2 = (((1.0/840.0))*(-533*CF_23 - 533*CF_24 - 29*CF_21 - 29*CF_26 + 3*CF_20 + 3*CF_27 + 139*CF_22 +
            139*CF_25) + Recon_2)*rj2;

       Recon_3 = (((1.0/840.0))*(-533*CF_33 - 533*CF_34 - 29*CF_31 - 29*CF_36 + 3*CF_30 + 3*CF_37 + 139*CF_32 +
            139*CF_35) + Recon_3)*rj3;

       Recon_4 = (((1.0/840.0))*(-533*CF_43 - 533*CF_44 - 29*CF_41 - 29*CF_46 + 3*CF_40 + 3*CF_47 + 139*CF_42 +
            139*CF_45) + Recon_4)*rj4;

       wk0_B0(0,0,0) = 0.707106781186547*AVG_0_rho*Recon_3*inv_AVG_a + 0.707106781186547*AVG_0_rho*Recon_4*inv_AVG_a +
            Recon_0;

       wk1_B0(0,0,0) = AVG_0_u0*Recon_0 + 0.707106781186547*(-AVG_0_a + AVG_0_u0)*AVG_0_rho*Recon_4*inv_AVG_a +
            0.707106781186547*(AVG_0_a + AVG_0_u0)*AVG_0_rho*Recon_3*inv_AVG_a;

       wk2_B0(0,0,0) = AVG_0_u1*Recon_0 - AVG_0_rho*Recon_2 + 0.707106781186547*AVG_0_rho*AVG_0_u1*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_0_rho*AVG_0_u1*Recon_4*inv_AVG_a;

       wk3_B0(0,0,0) = AVG_0_rho*Recon_1 + AVG_0_u2*Recon_0 + 0.707106781186547*AVG_0_rho*AVG_0_u2*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_0_rho*AVG_0_u2*Recon_4*inv_AVG_a;

       wk4_B0(0,0,0) = (((1.0/2.0))*(AVG_0_u0*AVG_0_u0) + ((1.0/2.0))*(AVG_0_u1*AVG_0_u1) +
            ((1.0/2.0))*(AVG_0_u2*AVG_0_u2))*Recon_0 + AVG_0_rho*AVG_0_u2*Recon_1 - AVG_0_rho*AVG_0_u1*Recon_2 +
            0.707106781186547*(((AVG_0_a*AVG_0_a) + ((1.0/2.0))*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) +
            (AVG_0_u2*AVG_0_u2))*gamma_m1)*invgamma_m1 + AVG_0_a*AVG_0_u0)*AVG_0_rho*Recon_3*inv_AVG_a +
            0.707106781186547*(((AVG_0_a*AVG_0_a) + ((1.0/2.0))*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) +
            (AVG_0_u2*AVG_0_u2))*gamma_m1)*invgamma_m1 - AVG_0_a*AVG_0_u0)*AVG_0_rho*Recon_4*inv_AVG_a;

   }

   else{

      wk0_B0(0,0,0) = 0.0;

      wk1_B0(0,0,0) = 0.0;

      wk2_B0(0,0,0) = 0.0;

      wk3_B0(0,0,0) = 0.0;

      wk4_B0(0,0,0) = 0.0;

   }

}

 void opensbliblock00Kernel059(const ACC<double> &a_B0, const ACC<double> &kappa_B0, const ACC<double> &p_B0, const
ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&Residual0_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0, ACC<double> &Residual3_B0, ACC<double>
&Residual4_B0)
{
   double AVG_1_1_LEV_00 = 0.0;
   double AVG_1_1_LEV_03 = 0.0;
   double AVG_1_1_LEV_10 = 0.0;
   double AVG_1_1_LEV_11 = 0.0;
   double AVG_1_1_LEV_12 = 0.0;
   double AVG_1_1_LEV_13 = 0.0;
   double AVG_1_1_LEV_14 = 0.0;
   double AVG_1_1_LEV_20 = 0.0;
   double AVG_1_1_LEV_21 = 0.0;
   double AVG_1_1_LEV_30 = 0.0;
   double AVG_1_1_LEV_31 = 0.0;
   double AVG_1_1_LEV_32 = 0.0;
   double AVG_1_1_LEV_33 = 0.0;
   double AVG_1_1_LEV_34 = 0.0;
   double AVG_1_1_LEV_40 = 0.0;
   double AVG_1_1_LEV_41 = 0.0;
   double AVG_1_1_LEV_42 = 0.0;
   double AVG_1_1_LEV_43 = 0.0;
   double AVG_1_1_LEV_44 = 0.0;
   double AVG_1_a = 0.0;
   double AVG_1_inv_rho = 0.0;
   double AVG_1_rho = 0.0;
   double AVG_1_u0 = 0.0;
   double AVG_1_u1 = 0.0;
   double AVG_1_u2 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_06 = 0.0;
   double CF_07 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_16 = 0.0;
   double CF_17 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_26 = 0.0;
   double CF_27 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CF_36 = 0.0;
   double CF_37 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CF_44 = 0.0;
   double CF_45 = 0.0;
   double CF_46 = 0.0;
   double CF_47 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_06 = 0.0;
   double CS_07 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_16 = 0.0;
   double CS_17 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_26 = 0.0;
   double CS_27 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double CS_36 = 0.0;
   double CS_37 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double CS_44 = 0.0;
   double CS_45 = 0.0;
   double CS_46 = 0.0;
   double CS_47 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double alpha_3 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double beta_3 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double omega_2 = 0.0;
   double omega_3 = 0.0;
   double rj0 = 0.0;
   double rj1 = 0.0;
   double rj2 = 0.0;
   double rj3 = 0.0;
   double rj4 = 0.0;
    if (fmax(kappa_B0(0,2,0), fmax(kappa_B0(0,-2,0), fmax(kappa_B0(0,-1,0), fmax(kappa_B0(0,1,0), fmax(kappa_B0(0,-3,0),
      kappa_B0(0,0,0)))))) > Ducros_check){

      AVG_1_rho = sqrt((rho_B0(0,0,0)*rho_B0(0,1,0)));

      AVG_1_inv_rho = 1.0/((sqrt(rho_B0(0,0,0)) + sqrt(rho_B0(0,1,0))));

      AVG_1_u0 = (sqrt(rho_B0(0,0,0))*u0_B0(0,0,0) + sqrt(rho_B0(0,1,0))*u0_B0(0,1,0))*AVG_1_inv_rho;

      AVG_1_u1 = (sqrt(rho_B0(0,0,0))*u1_B0(0,0,0) + sqrt(rho_B0(0,1,0))*u1_B0(0,1,0))*AVG_1_inv_rho;

      AVG_1_u2 = (sqrt(rho_B0(0,0,0))*u2_B0(0,0,0) + sqrt(rho_B0(0,1,0))*u2_B0(0,1,0))*AVG_1_inv_rho;

       AVG_1_a = sqrt(((-(1.0/2.0)*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) + (AVG_1_u2*AVG_1_u2)) + ((p_B0(0,0,0) +
            rhoE_B0(0,0,0))/sqrt(rho_B0(0,0,0)) + (p_B0(0,1,0) +
            rhoE_B0(0,1,0))/sqrt(rho_B0(0,1,0)))*AVG_1_inv_rho)*gamma_m1));

      inv_AVG_a = 1.0/(AVG_1_a);

      inv_AVG_rho = 1.0/(AVG_1_rho);

      AVG_1_1_LEV_00 = AVG_1_u2*inv_AVG_rho;

      AVG_1_1_LEV_03 = -inv_AVG_rho;

       AVG_1_1_LEV_10 = -(1.0/2.0)*(-2 - (AVG_1_u0*AVG_1_u0)*(inv_AVG_a*inv_AVG_a) -
            (AVG_1_u1*AVG_1_u1)*(inv_AVG_a*inv_AVG_a) - (AVG_1_u2*AVG_1_u2)*(inv_AVG_a*inv_AVG_a) +
            (AVG_1_u0*AVG_1_u0)*(inv_AVG_a*inv_AVG_a)*gama + (AVG_1_u1*AVG_1_u1)*(inv_AVG_a*inv_AVG_a)*gama +
            (AVG_1_u2*AVG_1_u2)*(inv_AVG_a*inv_AVG_a)*gama);

      AVG_1_1_LEV_11 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u0;

      AVG_1_1_LEV_12 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u1;

      AVG_1_1_LEV_13 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u2;

      AVG_1_1_LEV_14 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

      AVG_1_1_LEV_20 = -AVG_1_u0*inv_AVG_rho;

      AVG_1_1_LEV_21 = inv_AVG_rho;

       AVG_1_1_LEV_30 = -0.353553390593274*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) + (AVG_1_u2*AVG_1_u2) -
            (AVG_1_u0*AVG_1_u0)*gama - (AVG_1_u1*AVG_1_u1)*gama - (AVG_1_u2*AVG_1_u2)*gama +
            2*AVG_1_a*AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_31 = -0.707106781186547*gamma_m1*AVG_1_u0*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_32 = 0.707106781186547*(-gama*AVG_1_u1 + AVG_1_a + AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_33 = -0.707106781186547*gamma_m1*AVG_1_u2*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_34 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       AVG_1_1_LEV_40 = 0.353553390593274*(-(AVG_1_u0*AVG_1_u0) - (AVG_1_u1*AVG_1_u1) - (AVG_1_u2*AVG_1_u2) +
            (AVG_1_u0*AVG_1_u0)*gama + (AVG_1_u1*AVG_1_u1)*gama + (AVG_1_u2*AVG_1_u2)*gama +
            2*AVG_1_a*AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_41 = -0.707106781186547*gamma_m1*AVG_1_u0*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_42 = -0.707106781186547*(-AVG_1_u1 + gama*AVG_1_u1 + AVG_1_a)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_43 = -0.707106781186547*gamma_m1*AVG_1_u2*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_44 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

      CF_00 = (rho_B0(0,-3,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-3,0)*AVG_1_1_LEV_03)*u1_B0(0,-3,0);

       CF_10 = p_B0(0,-3,0)*AVG_1_1_LEV_12 + p_B0(0,-3,0)*u1_B0(0,-3,0)*AVG_1_1_LEV_14 +
            u1_B0(0,-3,0)*rho_B0(0,-3,0)*AVG_1_1_LEV_10 + u1_B0(0,-3,0)*rhoE_B0(0,-3,0)*AVG_1_1_LEV_14 +
            u1_B0(0,-3,0)*rhou0_B0(0,-3,0)*AVG_1_1_LEV_11 + u1_B0(0,-3,0)*rhou1_B0(0,-3,0)*AVG_1_1_LEV_12 +
            u1_B0(0,-3,0)*rhou2_B0(0,-3,0)*AVG_1_1_LEV_13;

      CF_20 = (rho_B0(0,-3,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-3,0)*AVG_1_1_LEV_21)*u1_B0(0,-3,0);

       CF_30 = p_B0(0,-3,0)*AVG_1_1_LEV_32 + p_B0(0,-3,0)*u1_B0(0,-3,0)*AVG_1_1_LEV_34 +
            u1_B0(0,-3,0)*rho_B0(0,-3,0)*AVG_1_1_LEV_30 + u1_B0(0,-3,0)*rhoE_B0(0,-3,0)*AVG_1_1_LEV_34 +
            u1_B0(0,-3,0)*rhou0_B0(0,-3,0)*AVG_1_1_LEV_31 + u1_B0(0,-3,0)*rhou1_B0(0,-3,0)*AVG_1_1_LEV_32 +
            u1_B0(0,-3,0)*rhou2_B0(0,-3,0)*AVG_1_1_LEV_33;

       CF_40 = p_B0(0,-3,0)*AVG_1_1_LEV_42 + p_B0(0,-3,0)*u1_B0(0,-3,0)*AVG_1_1_LEV_44 +
            u1_B0(0,-3,0)*rho_B0(0,-3,0)*AVG_1_1_LEV_40 + u1_B0(0,-3,0)*rhoE_B0(0,-3,0)*AVG_1_1_LEV_44 +
            u1_B0(0,-3,0)*rhou0_B0(0,-3,0)*AVG_1_1_LEV_41 + u1_B0(0,-3,0)*rhou1_B0(0,-3,0)*AVG_1_1_LEV_42 +
            u1_B0(0,-3,0)*rhou2_B0(0,-3,0)*AVG_1_1_LEV_43;

      CS_00 = rho_B0(0,-3,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-3,0)*AVG_1_1_LEV_03;

       CS_10 = rho_B0(0,-3,0)*AVG_1_1_LEV_10 + rhoE_B0(0,-3,0)*AVG_1_1_LEV_14 + rhou0_B0(0,-3,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,-3,0)*AVG_1_1_LEV_12 + rhou2_B0(0,-3,0)*AVG_1_1_LEV_13;

      CS_20 = rho_B0(0,-3,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-3,0)*AVG_1_1_LEV_21;

       CS_30 = rho_B0(0,-3,0)*AVG_1_1_LEV_30 + rhoE_B0(0,-3,0)*AVG_1_1_LEV_34 + rhou0_B0(0,-3,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,-3,0)*AVG_1_1_LEV_32 + rhou2_B0(0,-3,0)*AVG_1_1_LEV_33;

       CS_40 = rho_B0(0,-3,0)*AVG_1_1_LEV_40 + rhoE_B0(0,-3,0)*AVG_1_1_LEV_44 + rhou0_B0(0,-3,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,-3,0)*AVG_1_1_LEV_42 + rhou2_B0(0,-3,0)*AVG_1_1_LEV_43;

      CF_01 = (rho_B0(0,-2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_03)*u1_B0(0,-2,0);

       CF_11 = p_B0(0,-2,0)*AVG_1_1_LEV_12 + p_B0(0,-2,0)*u1_B0(0,-2,0)*AVG_1_1_LEV_14 +
            u1_B0(0,-2,0)*rho_B0(0,-2,0)*AVG_1_1_LEV_10 + u1_B0(0,-2,0)*rhoE_B0(0,-2,0)*AVG_1_1_LEV_14 +
            u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_11 + u1_B0(0,-2,0)*rhou1_B0(0,-2,0)*AVG_1_1_LEV_12 +
            u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_13;

      CF_21 = (rho_B0(0,-2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_21)*u1_B0(0,-2,0);

       CF_31 = p_B0(0,-2,0)*AVG_1_1_LEV_32 + p_B0(0,-2,0)*u1_B0(0,-2,0)*AVG_1_1_LEV_34 +
            u1_B0(0,-2,0)*rho_B0(0,-2,0)*AVG_1_1_LEV_30 + u1_B0(0,-2,0)*rhoE_B0(0,-2,0)*AVG_1_1_LEV_34 +
            u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_31 + u1_B0(0,-2,0)*rhou1_B0(0,-2,0)*AVG_1_1_LEV_32 +
            u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_33;

       CF_41 = p_B0(0,-2,0)*AVG_1_1_LEV_42 + p_B0(0,-2,0)*u1_B0(0,-2,0)*AVG_1_1_LEV_44 +
            u1_B0(0,-2,0)*rho_B0(0,-2,0)*AVG_1_1_LEV_40 + u1_B0(0,-2,0)*rhoE_B0(0,-2,0)*AVG_1_1_LEV_44 +
            u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_41 + u1_B0(0,-2,0)*rhou1_B0(0,-2,0)*AVG_1_1_LEV_42 +
            u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_43;

      CS_01 = rho_B0(0,-2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_03;

       CS_11 = rho_B0(0,-2,0)*AVG_1_1_LEV_10 + rhoE_B0(0,-2,0)*AVG_1_1_LEV_14 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,-2,0)*AVG_1_1_LEV_12 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_13;

      CS_21 = rho_B0(0,-2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_21;

       CS_31 = rho_B0(0,-2,0)*AVG_1_1_LEV_30 + rhoE_B0(0,-2,0)*AVG_1_1_LEV_34 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,-2,0)*AVG_1_1_LEV_32 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_33;

       CS_41 = rho_B0(0,-2,0)*AVG_1_1_LEV_40 + rhoE_B0(0,-2,0)*AVG_1_1_LEV_44 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,-2,0)*AVG_1_1_LEV_42 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_43;

      CF_02 = (rho_B0(0,-1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_03)*u1_B0(0,-1,0);

       CF_12 = p_B0(0,-1,0)*AVG_1_1_LEV_12 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_14 +
            u1_B0(0,-1,0)*rho_B0(0,-1,0)*AVG_1_1_LEV_10 + u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_14 +
            u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_11 + u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_12 +
            u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_13;

      CF_22 = (rho_B0(0,-1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_21)*u1_B0(0,-1,0);

       CF_32 = p_B0(0,-1,0)*AVG_1_1_LEV_32 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_34 +
            u1_B0(0,-1,0)*rho_B0(0,-1,0)*AVG_1_1_LEV_30 + u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_34 +
            u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_31 + u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_32 +
            u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_33;

       CF_42 = p_B0(0,-1,0)*AVG_1_1_LEV_42 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_44 +
            u1_B0(0,-1,0)*rho_B0(0,-1,0)*AVG_1_1_LEV_40 + u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_44 +
            u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_41 + u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_42 +
            u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_43;

      CS_02 = rho_B0(0,-1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_03;

       CS_12 = rho_B0(0,-1,0)*AVG_1_1_LEV_10 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_14 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,-1,0)*AVG_1_1_LEV_12 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_13;

      CS_22 = rho_B0(0,-1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_21;

       CS_32 = rho_B0(0,-1,0)*AVG_1_1_LEV_30 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_34 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,-1,0)*AVG_1_1_LEV_32 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_33;

       CS_42 = rho_B0(0,-1,0)*AVG_1_1_LEV_40 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_44 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,-1,0)*AVG_1_1_LEV_42 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_43;

      CF_03 = (rho_B0(0,0,0)*AVG_1_1_LEV_00 + rhou2_B0(0,0,0)*AVG_1_1_LEV_03)*u1_B0(0,0,0);

       CF_13 = p_B0(0,0,0)*AVG_1_1_LEV_12 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_14 +
            u1_B0(0,0,0)*rho_B0(0,0,0)*AVG_1_1_LEV_10 + u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_14 +
            u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_11 + u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_12 +
            u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_13;

      CF_23 = (rho_B0(0,0,0)*AVG_1_1_LEV_20 + rhou0_B0(0,0,0)*AVG_1_1_LEV_21)*u1_B0(0,0,0);

       CF_33 = p_B0(0,0,0)*AVG_1_1_LEV_32 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_34 +
            u1_B0(0,0,0)*rho_B0(0,0,0)*AVG_1_1_LEV_30 + u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_34 +
            u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_31 + u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_32 +
            u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_33;

       CF_43 = p_B0(0,0,0)*AVG_1_1_LEV_42 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_44 +
            u1_B0(0,0,0)*rho_B0(0,0,0)*AVG_1_1_LEV_40 + u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_44 +
            u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_41 + u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_42 +
            u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_43;

      CS_03 = rho_B0(0,0,0)*AVG_1_1_LEV_00 + rhou2_B0(0,0,0)*AVG_1_1_LEV_03;

       CS_13 = rho_B0(0,0,0)*AVG_1_1_LEV_10 + rhoE_B0(0,0,0)*AVG_1_1_LEV_14 + rhou0_B0(0,0,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,0,0)*AVG_1_1_LEV_12 + rhou2_B0(0,0,0)*AVG_1_1_LEV_13;

      CS_23 = rho_B0(0,0,0)*AVG_1_1_LEV_20 + rhou0_B0(0,0,0)*AVG_1_1_LEV_21;

       CS_33 = rho_B0(0,0,0)*AVG_1_1_LEV_30 + rhoE_B0(0,0,0)*AVG_1_1_LEV_34 + rhou0_B0(0,0,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,0,0)*AVG_1_1_LEV_32 + rhou2_B0(0,0,0)*AVG_1_1_LEV_33;

       CS_43 = rho_B0(0,0,0)*AVG_1_1_LEV_40 + rhoE_B0(0,0,0)*AVG_1_1_LEV_44 + rhou0_B0(0,0,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,0,0)*AVG_1_1_LEV_42 + rhou2_B0(0,0,0)*AVG_1_1_LEV_43;

      CF_04 = (rho_B0(0,1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,1,0)*AVG_1_1_LEV_03)*u1_B0(0,1,0);

       CF_14 = p_B0(0,1,0)*AVG_1_1_LEV_12 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_14 +
            u1_B0(0,1,0)*rho_B0(0,1,0)*AVG_1_1_LEV_10 + u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_14 +
            u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_11 + u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_12 +
            u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_13;

      CF_24 = (rho_B0(0,1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,1,0)*AVG_1_1_LEV_21)*u1_B0(0,1,0);

       CF_34 = p_B0(0,1,0)*AVG_1_1_LEV_32 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_34 +
            u1_B0(0,1,0)*rho_B0(0,1,0)*AVG_1_1_LEV_30 + u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_34 +
            u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_31 + u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_32 +
            u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_33;

       CF_44 = p_B0(0,1,0)*AVG_1_1_LEV_42 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_44 +
            u1_B0(0,1,0)*rho_B0(0,1,0)*AVG_1_1_LEV_40 + u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_44 +
            u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_41 + u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_42 +
            u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_43;

      CS_04 = rho_B0(0,1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,1,0)*AVG_1_1_LEV_03;

       CS_14 = rho_B0(0,1,0)*AVG_1_1_LEV_10 + rhoE_B0(0,1,0)*AVG_1_1_LEV_14 + rhou0_B0(0,1,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,1,0)*AVG_1_1_LEV_12 + rhou2_B0(0,1,0)*AVG_1_1_LEV_13;

      CS_24 = rho_B0(0,1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,1,0)*AVG_1_1_LEV_21;

       CS_34 = rho_B0(0,1,0)*AVG_1_1_LEV_30 + rhoE_B0(0,1,0)*AVG_1_1_LEV_34 + rhou0_B0(0,1,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,1,0)*AVG_1_1_LEV_32 + rhou2_B0(0,1,0)*AVG_1_1_LEV_33;

       CS_44 = rho_B0(0,1,0)*AVG_1_1_LEV_40 + rhoE_B0(0,1,0)*AVG_1_1_LEV_44 + rhou0_B0(0,1,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,1,0)*AVG_1_1_LEV_42 + rhou2_B0(0,1,0)*AVG_1_1_LEV_43;

      CF_05 = (rho_B0(0,2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,2,0)*AVG_1_1_LEV_03)*u1_B0(0,2,0);

       CF_15 = p_B0(0,2,0)*AVG_1_1_LEV_12 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_14 +
            u1_B0(0,2,0)*rho_B0(0,2,0)*AVG_1_1_LEV_10 + u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_14 +
            u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_11 + u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_12 +
            u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_13;

      CF_25 = (rho_B0(0,2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,2,0)*AVG_1_1_LEV_21)*u1_B0(0,2,0);

       CF_35 = p_B0(0,2,0)*AVG_1_1_LEV_32 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_34 +
            u1_B0(0,2,0)*rho_B0(0,2,0)*AVG_1_1_LEV_30 + u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_34 +
            u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_31 + u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_32 +
            u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_33;

       CF_45 = p_B0(0,2,0)*AVG_1_1_LEV_42 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_44 +
            u1_B0(0,2,0)*rho_B0(0,2,0)*AVG_1_1_LEV_40 + u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_44 +
            u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_41 + u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_42 +
            u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_43;

      CS_05 = rho_B0(0,2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,2,0)*AVG_1_1_LEV_03;

       CS_15 = rho_B0(0,2,0)*AVG_1_1_LEV_10 + rhoE_B0(0,2,0)*AVG_1_1_LEV_14 + rhou0_B0(0,2,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,2,0)*AVG_1_1_LEV_12 + rhou2_B0(0,2,0)*AVG_1_1_LEV_13;

      CS_25 = rho_B0(0,2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,2,0)*AVG_1_1_LEV_21;

       CS_35 = rho_B0(0,2,0)*AVG_1_1_LEV_30 + rhoE_B0(0,2,0)*AVG_1_1_LEV_34 + rhou0_B0(0,2,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,2,0)*AVG_1_1_LEV_32 + rhou2_B0(0,2,0)*AVG_1_1_LEV_33;

       CS_45 = rho_B0(0,2,0)*AVG_1_1_LEV_40 + rhoE_B0(0,2,0)*AVG_1_1_LEV_44 + rhou0_B0(0,2,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,2,0)*AVG_1_1_LEV_42 + rhou2_B0(0,2,0)*AVG_1_1_LEV_43;

      CF_06 = (rho_B0(0,3,0)*AVG_1_1_LEV_00 + rhou2_B0(0,3,0)*AVG_1_1_LEV_03)*u1_B0(0,3,0);

       CF_16 = p_B0(0,3,0)*AVG_1_1_LEV_12 + p_B0(0,3,0)*u1_B0(0,3,0)*AVG_1_1_LEV_14 +
            u1_B0(0,3,0)*rho_B0(0,3,0)*AVG_1_1_LEV_10 + u1_B0(0,3,0)*rhoE_B0(0,3,0)*AVG_1_1_LEV_14 +
            u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_11 + u1_B0(0,3,0)*rhou1_B0(0,3,0)*AVG_1_1_LEV_12 +
            u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_13;

      CF_26 = (rho_B0(0,3,0)*AVG_1_1_LEV_20 + rhou0_B0(0,3,0)*AVG_1_1_LEV_21)*u1_B0(0,3,0);

       CF_36 = p_B0(0,3,0)*AVG_1_1_LEV_32 + p_B0(0,3,0)*u1_B0(0,3,0)*AVG_1_1_LEV_34 +
            u1_B0(0,3,0)*rho_B0(0,3,0)*AVG_1_1_LEV_30 + u1_B0(0,3,0)*rhoE_B0(0,3,0)*AVG_1_1_LEV_34 +
            u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_31 + u1_B0(0,3,0)*rhou1_B0(0,3,0)*AVG_1_1_LEV_32 +
            u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_33;

       CF_46 = p_B0(0,3,0)*AVG_1_1_LEV_42 + p_B0(0,3,0)*u1_B0(0,3,0)*AVG_1_1_LEV_44 +
            u1_B0(0,3,0)*rho_B0(0,3,0)*AVG_1_1_LEV_40 + u1_B0(0,3,0)*rhoE_B0(0,3,0)*AVG_1_1_LEV_44 +
            u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_41 + u1_B0(0,3,0)*rhou1_B0(0,3,0)*AVG_1_1_LEV_42 +
            u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_43;

      CS_06 = rho_B0(0,3,0)*AVG_1_1_LEV_00 + rhou2_B0(0,3,0)*AVG_1_1_LEV_03;

       CS_16 = rho_B0(0,3,0)*AVG_1_1_LEV_10 + rhoE_B0(0,3,0)*AVG_1_1_LEV_14 + rhou0_B0(0,3,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,3,0)*AVG_1_1_LEV_12 + rhou2_B0(0,3,0)*AVG_1_1_LEV_13;

      CS_26 = rho_B0(0,3,0)*AVG_1_1_LEV_20 + rhou0_B0(0,3,0)*AVG_1_1_LEV_21;

       CS_36 = rho_B0(0,3,0)*AVG_1_1_LEV_30 + rhoE_B0(0,3,0)*AVG_1_1_LEV_34 + rhou0_B0(0,3,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,3,0)*AVG_1_1_LEV_32 + rhou2_B0(0,3,0)*AVG_1_1_LEV_33;

       CS_46 = rho_B0(0,3,0)*AVG_1_1_LEV_40 + rhoE_B0(0,3,0)*AVG_1_1_LEV_44 + rhou0_B0(0,3,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,3,0)*AVG_1_1_LEV_42 + rhou2_B0(0,3,0)*AVG_1_1_LEV_43;

      CF_07 = (rho_B0(0,4,0)*AVG_1_1_LEV_00 + rhou2_B0(0,4,0)*AVG_1_1_LEV_03)*u1_B0(0,4,0);

       CF_17 = p_B0(0,4,0)*AVG_1_1_LEV_12 + p_B0(0,4,0)*u1_B0(0,4,0)*AVG_1_1_LEV_14 +
            u1_B0(0,4,0)*rho_B0(0,4,0)*AVG_1_1_LEV_10 + u1_B0(0,4,0)*rhoE_B0(0,4,0)*AVG_1_1_LEV_14 +
            u1_B0(0,4,0)*rhou0_B0(0,4,0)*AVG_1_1_LEV_11 + u1_B0(0,4,0)*rhou1_B0(0,4,0)*AVG_1_1_LEV_12 +
            u1_B0(0,4,0)*rhou2_B0(0,4,0)*AVG_1_1_LEV_13;

      CF_27 = (rho_B0(0,4,0)*AVG_1_1_LEV_20 + rhou0_B0(0,4,0)*AVG_1_1_LEV_21)*u1_B0(0,4,0);

       CF_37 = p_B0(0,4,0)*AVG_1_1_LEV_32 + p_B0(0,4,0)*u1_B0(0,4,0)*AVG_1_1_LEV_34 +
            u1_B0(0,4,0)*rho_B0(0,4,0)*AVG_1_1_LEV_30 + u1_B0(0,4,0)*rhoE_B0(0,4,0)*AVG_1_1_LEV_34 +
            u1_B0(0,4,0)*rhou0_B0(0,4,0)*AVG_1_1_LEV_31 + u1_B0(0,4,0)*rhou1_B0(0,4,0)*AVG_1_1_LEV_32 +
            u1_B0(0,4,0)*rhou2_B0(0,4,0)*AVG_1_1_LEV_33;

       CF_47 = p_B0(0,4,0)*AVG_1_1_LEV_42 + p_B0(0,4,0)*u1_B0(0,4,0)*AVG_1_1_LEV_44 +
            u1_B0(0,4,0)*rho_B0(0,4,0)*AVG_1_1_LEV_40 + u1_B0(0,4,0)*rhoE_B0(0,4,0)*AVG_1_1_LEV_44 +
            u1_B0(0,4,0)*rhou0_B0(0,4,0)*AVG_1_1_LEV_41 + u1_B0(0,4,0)*rhou1_B0(0,4,0)*AVG_1_1_LEV_42 +
            u1_B0(0,4,0)*rhou2_B0(0,4,0)*AVG_1_1_LEV_43;

      CS_07 = rho_B0(0,4,0)*AVG_1_1_LEV_00 + rhou2_B0(0,4,0)*AVG_1_1_LEV_03;

       CS_17 = rho_B0(0,4,0)*AVG_1_1_LEV_10 + rhoE_B0(0,4,0)*AVG_1_1_LEV_14 + rhou0_B0(0,4,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,4,0)*AVG_1_1_LEV_12 + rhou2_B0(0,4,0)*AVG_1_1_LEV_13;

      CS_27 = rho_B0(0,4,0)*AVG_1_1_LEV_20 + rhou0_B0(0,4,0)*AVG_1_1_LEV_21;

       CS_37 = rho_B0(0,4,0)*AVG_1_1_LEV_30 + rhoE_B0(0,4,0)*AVG_1_1_LEV_34 + rhou0_B0(0,4,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,4,0)*AVG_1_1_LEV_32 + rhou2_B0(0,4,0)*AVG_1_1_LEV_33;

       CS_47 = rho_B0(0,4,0)*AVG_1_1_LEV_40 + rhoE_B0(0,4,0)*AVG_1_1_LEV_44 + rhou0_B0(0,4,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,4,0)*AVG_1_1_LEV_42 + rhou2_B0(0,4,0)*AVG_1_1_LEV_43;

      max_lambda_00 = shock_filter_control*fmax(fabs(u1_B0(0,1,0)), fabs(u1_B0(0,0,0)));

      max_lambda_11 = max_lambda_00;

      max_lambda_22 = max_lambda_00;

      max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(0,1,0) + u1_B0(0,1,0)), fabs(a_B0(0,0,0) + u1_B0(0,0,0)));

      max_lambda_44 = shock_filter_control*fmax(fabs(-u1_B0(0,1,0) + a_B0(0,1,0)), fabs(-u1_B0(0,0,0) + a_B0(0,0,0)));

       beta_0 = ((547.0/960.0))*((CS_06*max_lambda_00 + CF_06)*(CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_06*max_lambda_00 + CF_06) + ((7043.0/480.0))*(CS_05*max_lambda_00 +
            CF_05))*(CS_05*max_lambda_00 + CF_05) + ((1.0/2.0))*(CS_03*max_lambda_00 +
            CF_03)*(-(1567.0/80.0)*(CS_04*max_lambda_00 + CF_04) - (309.0/80.0)*(CS_06*max_lambda_00 + CF_06) +
            ((2107.0/480.0))*(CS_03*max_lambda_00 + CF_03) + ((3521.0/240.0))*(CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04)*(-(8623.0/240.0)*(CS_05*max_lambda_00 + CF_05) +
            ((2321.0/240.0))*(CS_06*max_lambda_00 + CF_06) + ((11003.0/480.0))*(CS_04*max_lambda_00 + CF_04));

       beta_1 = ((89.0/320.0))*((CS_05*max_lambda_00 + CF_05)*(CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_05*max_lambda_00 + CF_05) + ((2843.0/480.0))*(CS_04*max_lambda_00 +
            CF_04))*(CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(1261.0/240.0)*(CS_03*max_lambda_00 + CF_03) - (247.0/240.0)*(CS_05*max_lambda_00 + CF_05) +
            ((547.0/480.0))*(CS_02*max_lambda_00 + CF_02) + ((961.0/240.0))*(CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(CS_03*max_lambda_00 + CF_03)*(-(2983.0/240.0)*(CS_04*max_lambda_00 + CF_04) +
            ((267.0/80.0))*(CS_05*max_lambda_00 + CF_05) + ((3443.0/480.0))*(CS_03*max_lambda_00 + CF_03));

       beta_2 = ((547.0/960.0))*((CS_04*max_lambda_00 + CF_04)*(CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_04*max_lambda_00 + CF_04) + ((3443.0/480.0))*(CS_03*max_lambda_00 +
            CF_03))*(CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(CS_01*max_lambda_00 +
            CF_01)*(-(247.0/240.0)*(CS_04*max_lambda_00 + CF_04) + ((89.0/160.0))*(CS_01*max_lambda_00 + CF_01) +
            ((267.0/80.0))*(CS_03*max_lambda_00 + CF_03)) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(2983.0/240.0)*(CS_03*max_lambda_00 + CF_03) - (821.0/240.0)*(CS_01*max_lambda_00 + CF_01) +
            ((961.0/240.0))*(CS_04*max_lambda_00 + CF_04) + ((2843.0/480.0))*(CS_02*max_lambda_00 + CF_02));

       beta_3 = ((2107.0/960.0))*((CS_03*max_lambda_00 + CF_03)*(CS_03*max_lambda_00 + CF_03)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_03*max_lambda_00 + CF_03) + ((547.0/480.0))*(CS_00*max_lambda_00 +
            CF_00))*(CS_00*max_lambda_00 + CF_00) + ((1.0/2.0))*(CS_01*max_lambda_00 +
            CF_01)*(-(647.0/80.0)*(CS_00*max_lambda_00 + CF_00) + ((3521.0/240.0))*(CS_03*max_lambda_00 + CF_03) +
            ((7043.0/480.0))*(CS_01*max_lambda_00 + CF_01)) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(8623.0/240.0)*(CS_01*max_lambda_00 + CF_01) - (1567.0/80.0)*(CS_03*max_lambda_00 + CF_03) +
            ((2321.0/240.0))*(CS_00*max_lambda_00 + CF_00) + ((11003.0/480.0))*(CS_02*max_lambda_00 + CF_02));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj0 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_0 = (-(23.0/24.0)*(CS_02*max_lambda_00 + CF_02) - (1.0/8.0)*(CS_00*max_lambda_00 + CF_00) +
            ((13.0/24.0))*(CS_01*max_lambda_00 + CF_01) + ((25.0/24.0))*(CS_03*max_lambda_00 + CF_03))*omega_3 +
            (-(5.0/24.0)*(CS_02*max_lambda_00 + CF_02) + ((1.0/8.0))*(CS_04*max_lambda_00 + CF_04) +
            ((1.0/24.0))*(CS_01*max_lambda_00 + CF_01) + ((13.0/24.0))*(CS_03*max_lambda_00 + CF_03))*omega_2 +
            (-(5.0/24.0)*(CS_05*max_lambda_00 + CF_05) + ((1.0/8.0))*(CS_03*max_lambda_00 + CF_03) +
            ((1.0/24.0))*(CS_06*max_lambda_00 + CF_06) + ((13.0/24.0))*(CS_04*max_lambda_00 + CF_04))*omega_0 +
            (-(1.0/24.0)*(CS_02*max_lambda_00 + CF_02) - (1.0/24.0)*(CS_05*max_lambda_00 + CF_05) +
            ((7.0/24.0))*(CS_03*max_lambda_00 + CF_03) + ((7.0/24.0))*(CS_04*max_lambda_00 + CF_04))*omega_1 + Recon_0;

       beta_0 = ((547.0/960.0))*((-CS_07*max_lambda_00 + CF_07)*(-CS_07*max_lambda_00 + CF_07)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_07*max_lambda_00 + CF_07) + ((7043.0/480.0))*(-CS_06*max_lambda_00 +
            CF_06))*(-CS_06*max_lambda_00 + CF_06) + ((1.0/2.0))*(-CS_04*max_lambda_00 +
            CF_04)*(-(1567.0/80.0)*(-CS_05*max_lambda_00 + CF_05) - (309.0/80.0)*(-CS_07*max_lambda_00 + CF_07) +
            ((2107.0/480.0))*(-CS_04*max_lambda_00 + CF_04) + ((3521.0/240.0))*(-CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05)*(-(8623.0/240.0)*(-CS_06*max_lambda_00 + CF_06) +
            ((2321.0/240.0))*(-CS_07*max_lambda_00 + CF_07) + ((11003.0/480.0))*(-CS_05*max_lambda_00 + CF_05));

       beta_1 = ((89.0/320.0))*((-CS_06*max_lambda_00 + CF_06)*(-CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_06*max_lambda_00 + CF_06) + ((2843.0/480.0))*(-CS_05*max_lambda_00 +
            CF_05))*(-CS_05*max_lambda_00 + CF_05) + ((1.0/2.0))*(-CS_03*max_lambda_00 +
            CF_03)*(-(1261.0/240.0)*(-CS_04*max_lambda_00 + CF_04) - (247.0/240.0)*(-CS_06*max_lambda_00 + CF_06) +
            ((547.0/480.0))*(-CS_03*max_lambda_00 + CF_03) + ((961.0/240.0))*(-CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-CS_04*max_lambda_00 + CF_04)*(-(2983.0/240.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((267.0/80.0))*(-CS_06*max_lambda_00 + CF_06) + ((3443.0/480.0))*(-CS_04*max_lambda_00 + CF_04));

       beta_2 = ((547.0/960.0))*((-CS_05*max_lambda_00 + CF_05)*(-CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_05*max_lambda_00 + CF_05) + ((3443.0/480.0))*(-CS_04*max_lambda_00 +
            CF_04))*(-CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
            CF_02)*(-(821.0/240.0)*(-CS_03*max_lambda_00 + CF_03) - (247.0/240.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((89.0/160.0))*(-CS_02*max_lambda_00 + CF_02) + ((267.0/80.0))*(-CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03)*(-(2983.0/240.0)*(-CS_04*max_lambda_00 + CF_04) +
            ((961.0/240.0))*(-CS_05*max_lambda_00 + CF_05) + ((2843.0/480.0))*(-CS_03*max_lambda_00 + CF_03));

       beta_3 = ((2107.0/960.0))*((-CS_04*max_lambda_00 + CF_04)*(-CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_04*max_lambda_00 + CF_04) + ((11003.0/480.0))*(-CS_03*max_lambda_00 +
            CF_03))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_01*max_lambda_00 +
            CF_01)*(-(309.0/80.0)*(-CS_04*max_lambda_00 + CF_04) + ((547.0/480.0))*(-CS_01*max_lambda_00 + CF_01) +
            ((2321.0/240.0))*(-CS_03*max_lambda_00 + CF_03)) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
            CF_02)*(-(8623.0/240.0)*(-CS_03*max_lambda_00 + CF_03) - (647.0/80.0)*(-CS_01*max_lambda_00 + CF_01) +
            ((3521.0/240.0))*(-CS_04*max_lambda_00 + CF_04) + ((7043.0/480.0))*(-CS_02*max_lambda_00 + CF_02));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj0 = fmax(rj0, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_0 = (-(23.0/24.0)*(-CS_05*max_lambda_00 + CF_05) - (1.0/8.0)*(-CS_07*max_lambda_00 + CF_07) +
            ((13.0/24.0))*(-CS_06*max_lambda_00 + CF_06) + ((25.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_0 +
            (-(5.0/24.0)*(-CS_02*max_lambda_00 + CF_02) + ((1.0/8.0))*(-CS_04*max_lambda_00 + CF_04) +
            ((1.0/24.0))*(-CS_01*max_lambda_00 + CF_01) + ((13.0/24.0))*(-CS_03*max_lambda_00 + CF_03))*omega_3 +
            (-(5.0/24.0)*(-CS_05*max_lambda_00 + CF_05) + ((1.0/8.0))*(-CS_03*max_lambda_00 + CF_03) +
            ((1.0/24.0))*(-CS_06*max_lambda_00 + CF_06) + ((13.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_1 +
            (-(1.0/24.0)*(-CS_02*max_lambda_00 + CF_02) - (1.0/24.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((7.0/24.0))*(-CS_03*max_lambda_00 + CF_03) + ((7.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_2 +
            Recon_0;

       beta_0 = ((547.0/960.0))*((CS_16*max_lambda_11 + CF_16)*(CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_16*max_lambda_11 + CF_16) + ((7043.0/480.0))*(CS_15*max_lambda_11 +
            CF_15))*(CS_15*max_lambda_11 + CF_15) + ((1.0/2.0))*(CS_13*max_lambda_11 +
            CF_13)*(-(1567.0/80.0)*(CS_14*max_lambda_11 + CF_14) - (309.0/80.0)*(CS_16*max_lambda_11 + CF_16) +
            ((2107.0/480.0))*(CS_13*max_lambda_11 + CF_13) + ((3521.0/240.0))*(CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14)*(-(8623.0/240.0)*(CS_15*max_lambda_11 + CF_15) +
            ((2321.0/240.0))*(CS_16*max_lambda_11 + CF_16) + ((11003.0/480.0))*(CS_14*max_lambda_11 + CF_14));

       beta_1 = ((89.0/320.0))*((CS_15*max_lambda_11 + CF_15)*(CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_15*max_lambda_11 + CF_15) + ((2843.0/480.0))*(CS_14*max_lambda_11 +
            CF_14))*(CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(1261.0/240.0)*(CS_13*max_lambda_11 + CF_13) - (247.0/240.0)*(CS_15*max_lambda_11 + CF_15) +
            ((547.0/480.0))*(CS_12*max_lambda_11 + CF_12) + ((961.0/240.0))*(CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(CS_13*max_lambda_11 + CF_13)*(-(2983.0/240.0)*(CS_14*max_lambda_11 + CF_14) +
            ((267.0/80.0))*(CS_15*max_lambda_11 + CF_15) + ((3443.0/480.0))*(CS_13*max_lambda_11 + CF_13));

       beta_2 = ((547.0/960.0))*((CS_14*max_lambda_11 + CF_14)*(CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_14*max_lambda_11 + CF_14) + ((3443.0/480.0))*(CS_13*max_lambda_11 +
            CF_13))*(CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(CS_11*max_lambda_11 +
            CF_11)*(-(247.0/240.0)*(CS_14*max_lambda_11 + CF_14) + ((89.0/160.0))*(CS_11*max_lambda_11 + CF_11) +
            ((267.0/80.0))*(CS_13*max_lambda_11 + CF_13)) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(2983.0/240.0)*(CS_13*max_lambda_11 + CF_13) - (821.0/240.0)*(CS_11*max_lambda_11 + CF_11) +
            ((961.0/240.0))*(CS_14*max_lambda_11 + CF_14) + ((2843.0/480.0))*(CS_12*max_lambda_11 + CF_12));

       beta_3 = ((2107.0/960.0))*((CS_13*max_lambda_11 + CF_13)*(CS_13*max_lambda_11 + CF_13)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_13*max_lambda_11 + CF_13) + ((547.0/480.0))*(CS_10*max_lambda_11 +
            CF_10))*(CS_10*max_lambda_11 + CF_10) + ((1.0/2.0))*(CS_11*max_lambda_11 +
            CF_11)*(-(647.0/80.0)*(CS_10*max_lambda_11 + CF_10) + ((3521.0/240.0))*(CS_13*max_lambda_11 + CF_13) +
            ((7043.0/480.0))*(CS_11*max_lambda_11 + CF_11)) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(8623.0/240.0)*(CS_11*max_lambda_11 + CF_11) - (1567.0/80.0)*(CS_13*max_lambda_11 + CF_13) +
            ((2321.0/240.0))*(CS_10*max_lambda_11 + CF_10) + ((11003.0/480.0))*(CS_12*max_lambda_11 + CF_12));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj1 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_1 = (-(23.0/24.0)*(CS_12*max_lambda_11 + CF_12) - (1.0/8.0)*(CS_10*max_lambda_11 + CF_10) +
            ((13.0/24.0))*(CS_11*max_lambda_11 + CF_11) + ((25.0/24.0))*(CS_13*max_lambda_11 + CF_13))*omega_3 +
            (-(5.0/24.0)*(CS_12*max_lambda_11 + CF_12) + ((1.0/8.0))*(CS_14*max_lambda_11 + CF_14) +
            ((1.0/24.0))*(CS_11*max_lambda_11 + CF_11) + ((13.0/24.0))*(CS_13*max_lambda_11 + CF_13))*omega_2 +
            (-(5.0/24.0)*(CS_15*max_lambda_11 + CF_15) + ((1.0/8.0))*(CS_13*max_lambda_11 + CF_13) +
            ((1.0/24.0))*(CS_16*max_lambda_11 + CF_16) + ((13.0/24.0))*(CS_14*max_lambda_11 + CF_14))*omega_0 +
            (-(1.0/24.0)*(CS_12*max_lambda_11 + CF_12) - (1.0/24.0)*(CS_15*max_lambda_11 + CF_15) +
            ((7.0/24.0))*(CS_13*max_lambda_11 + CF_13) + ((7.0/24.0))*(CS_14*max_lambda_11 + CF_14))*omega_1 + Recon_1;

       beta_0 = ((547.0/960.0))*((-CS_17*max_lambda_11 + CF_17)*(-CS_17*max_lambda_11 + CF_17)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_17*max_lambda_11 + CF_17) + ((7043.0/480.0))*(-CS_16*max_lambda_11 +
            CF_16))*(-CS_16*max_lambda_11 + CF_16) + ((1.0/2.0))*(-CS_14*max_lambda_11 +
            CF_14)*(-(1567.0/80.0)*(-CS_15*max_lambda_11 + CF_15) - (309.0/80.0)*(-CS_17*max_lambda_11 + CF_17) +
            ((2107.0/480.0))*(-CS_14*max_lambda_11 + CF_14) + ((3521.0/240.0))*(-CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15)*(-(8623.0/240.0)*(-CS_16*max_lambda_11 + CF_16) +
            ((2321.0/240.0))*(-CS_17*max_lambda_11 + CF_17) + ((11003.0/480.0))*(-CS_15*max_lambda_11 + CF_15));

       beta_1 = ((89.0/320.0))*((-CS_16*max_lambda_11 + CF_16)*(-CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_16*max_lambda_11 + CF_16) + ((2843.0/480.0))*(-CS_15*max_lambda_11 +
            CF_15))*(-CS_15*max_lambda_11 + CF_15) + ((1.0/2.0))*(-CS_13*max_lambda_11 +
            CF_13)*(-(1261.0/240.0)*(-CS_14*max_lambda_11 + CF_14) - (247.0/240.0)*(-CS_16*max_lambda_11 + CF_16) +
            ((547.0/480.0))*(-CS_13*max_lambda_11 + CF_13) + ((961.0/240.0))*(-CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-CS_14*max_lambda_11 + CF_14)*(-(2983.0/240.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((267.0/80.0))*(-CS_16*max_lambda_11 + CF_16) + ((3443.0/480.0))*(-CS_14*max_lambda_11 + CF_14));

       beta_2 = ((547.0/960.0))*((-CS_15*max_lambda_11 + CF_15)*(-CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_15*max_lambda_11 + CF_15) + ((3443.0/480.0))*(-CS_14*max_lambda_11 +
            CF_14))*(-CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
            CF_12)*(-(821.0/240.0)*(-CS_13*max_lambda_11 + CF_13) - (247.0/240.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((89.0/160.0))*(-CS_12*max_lambda_11 + CF_12) + ((267.0/80.0))*(-CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13)*(-(2983.0/240.0)*(-CS_14*max_lambda_11 + CF_14) +
            ((961.0/240.0))*(-CS_15*max_lambda_11 + CF_15) + ((2843.0/480.0))*(-CS_13*max_lambda_11 + CF_13));

       beta_3 = ((2107.0/960.0))*((-CS_14*max_lambda_11 + CF_14)*(-CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_14*max_lambda_11 + CF_14) + ((11003.0/480.0))*(-CS_13*max_lambda_11 +
            CF_13))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_11*max_lambda_11 +
            CF_11)*(-(309.0/80.0)*(-CS_14*max_lambda_11 + CF_14) + ((547.0/480.0))*(-CS_11*max_lambda_11 + CF_11) +
            ((2321.0/240.0))*(-CS_13*max_lambda_11 + CF_13)) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
            CF_12)*(-(8623.0/240.0)*(-CS_13*max_lambda_11 + CF_13) - (647.0/80.0)*(-CS_11*max_lambda_11 + CF_11) +
            ((3521.0/240.0))*(-CS_14*max_lambda_11 + CF_14) + ((7043.0/480.0))*(-CS_12*max_lambda_11 + CF_12));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj1 = fmax(rj1, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_1 = (-(23.0/24.0)*(-CS_15*max_lambda_11 + CF_15) - (1.0/8.0)*(-CS_17*max_lambda_11 + CF_17) +
            ((13.0/24.0))*(-CS_16*max_lambda_11 + CF_16) + ((25.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_0 +
            (-(5.0/24.0)*(-CS_12*max_lambda_11 + CF_12) + ((1.0/8.0))*(-CS_14*max_lambda_11 + CF_14) +
            ((1.0/24.0))*(-CS_11*max_lambda_11 + CF_11) + ((13.0/24.0))*(-CS_13*max_lambda_11 + CF_13))*omega_3 +
            (-(5.0/24.0)*(-CS_15*max_lambda_11 + CF_15) + ((1.0/8.0))*(-CS_13*max_lambda_11 + CF_13) +
            ((1.0/24.0))*(-CS_16*max_lambda_11 + CF_16) + ((13.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_1 +
            (-(1.0/24.0)*(-CS_12*max_lambda_11 + CF_12) - (1.0/24.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((7.0/24.0))*(-CS_13*max_lambda_11 + CF_13) + ((7.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_2 +
            Recon_1;

       beta_0 = ((547.0/960.0))*((CS_26*max_lambda_22 + CF_26)*(CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_26*max_lambda_22 + CF_26) + ((7043.0/480.0))*(CS_25*max_lambda_22 +
            CF_25))*(CS_25*max_lambda_22 + CF_25) + ((1.0/2.0))*(CS_23*max_lambda_22 +
            CF_23)*(-(1567.0/80.0)*(CS_24*max_lambda_22 + CF_24) - (309.0/80.0)*(CS_26*max_lambda_22 + CF_26) +
            ((2107.0/480.0))*(CS_23*max_lambda_22 + CF_23) + ((3521.0/240.0))*(CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24)*(-(8623.0/240.0)*(CS_25*max_lambda_22 + CF_25) +
            ((2321.0/240.0))*(CS_26*max_lambda_22 + CF_26) + ((11003.0/480.0))*(CS_24*max_lambda_22 + CF_24));

       beta_1 = ((89.0/320.0))*((CS_25*max_lambda_22 + CF_25)*(CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_25*max_lambda_22 + CF_25) + ((2843.0/480.0))*(CS_24*max_lambda_22 +
            CF_24))*(CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(1261.0/240.0)*(CS_23*max_lambda_22 + CF_23) - (247.0/240.0)*(CS_25*max_lambda_22 + CF_25) +
            ((547.0/480.0))*(CS_22*max_lambda_22 + CF_22) + ((961.0/240.0))*(CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(CS_23*max_lambda_22 + CF_23)*(-(2983.0/240.0)*(CS_24*max_lambda_22 + CF_24) +
            ((267.0/80.0))*(CS_25*max_lambda_22 + CF_25) + ((3443.0/480.0))*(CS_23*max_lambda_22 + CF_23));

       beta_2 = ((547.0/960.0))*((CS_24*max_lambda_22 + CF_24)*(CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_24*max_lambda_22 + CF_24) + ((3443.0/480.0))*(CS_23*max_lambda_22 +
            CF_23))*(CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(CS_21*max_lambda_22 +
            CF_21)*(-(247.0/240.0)*(CS_24*max_lambda_22 + CF_24) + ((89.0/160.0))*(CS_21*max_lambda_22 + CF_21) +
            ((267.0/80.0))*(CS_23*max_lambda_22 + CF_23)) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(2983.0/240.0)*(CS_23*max_lambda_22 + CF_23) - (821.0/240.0)*(CS_21*max_lambda_22 + CF_21) +
            ((961.0/240.0))*(CS_24*max_lambda_22 + CF_24) + ((2843.0/480.0))*(CS_22*max_lambda_22 + CF_22));

       beta_3 = ((2107.0/960.0))*((CS_23*max_lambda_22 + CF_23)*(CS_23*max_lambda_22 + CF_23)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_23*max_lambda_22 + CF_23) + ((547.0/480.0))*(CS_20*max_lambda_22 +
            CF_20))*(CS_20*max_lambda_22 + CF_20) + ((1.0/2.0))*(CS_21*max_lambda_22 +
            CF_21)*(-(647.0/80.0)*(CS_20*max_lambda_22 + CF_20) + ((3521.0/240.0))*(CS_23*max_lambda_22 + CF_23) +
            ((7043.0/480.0))*(CS_21*max_lambda_22 + CF_21)) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(8623.0/240.0)*(CS_21*max_lambda_22 + CF_21) - (1567.0/80.0)*(CS_23*max_lambda_22 + CF_23) +
            ((2321.0/240.0))*(CS_20*max_lambda_22 + CF_20) + ((11003.0/480.0))*(CS_22*max_lambda_22 + CF_22));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj2 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_2 = (-(23.0/24.0)*(CS_22*max_lambda_22 + CF_22) - (1.0/8.0)*(CS_20*max_lambda_22 + CF_20) +
            ((13.0/24.0))*(CS_21*max_lambda_22 + CF_21) + ((25.0/24.0))*(CS_23*max_lambda_22 + CF_23))*omega_3 +
            (-(5.0/24.0)*(CS_22*max_lambda_22 + CF_22) + ((1.0/8.0))*(CS_24*max_lambda_22 + CF_24) +
            ((1.0/24.0))*(CS_21*max_lambda_22 + CF_21) + ((13.0/24.0))*(CS_23*max_lambda_22 + CF_23))*omega_2 +
            (-(5.0/24.0)*(CS_25*max_lambda_22 + CF_25) + ((1.0/8.0))*(CS_23*max_lambda_22 + CF_23) +
            ((1.0/24.0))*(CS_26*max_lambda_22 + CF_26) + ((13.0/24.0))*(CS_24*max_lambda_22 + CF_24))*omega_0 +
            (-(1.0/24.0)*(CS_22*max_lambda_22 + CF_22) - (1.0/24.0)*(CS_25*max_lambda_22 + CF_25) +
            ((7.0/24.0))*(CS_23*max_lambda_22 + CF_23) + ((7.0/24.0))*(CS_24*max_lambda_22 + CF_24))*omega_1 + Recon_2;

       beta_0 = ((547.0/960.0))*((-CS_27*max_lambda_22 + CF_27)*(-CS_27*max_lambda_22 + CF_27)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_27*max_lambda_22 + CF_27) + ((7043.0/480.0))*(-CS_26*max_lambda_22 +
            CF_26))*(-CS_26*max_lambda_22 + CF_26) + ((1.0/2.0))*(-CS_24*max_lambda_22 +
            CF_24)*(-(1567.0/80.0)*(-CS_25*max_lambda_22 + CF_25) - (309.0/80.0)*(-CS_27*max_lambda_22 + CF_27) +
            ((2107.0/480.0))*(-CS_24*max_lambda_22 + CF_24) + ((3521.0/240.0))*(-CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25)*(-(8623.0/240.0)*(-CS_26*max_lambda_22 + CF_26) +
            ((2321.0/240.0))*(-CS_27*max_lambda_22 + CF_27) + ((11003.0/480.0))*(-CS_25*max_lambda_22 + CF_25));

       beta_1 = ((89.0/320.0))*((-CS_26*max_lambda_22 + CF_26)*(-CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_26*max_lambda_22 + CF_26) + ((2843.0/480.0))*(-CS_25*max_lambda_22 +
            CF_25))*(-CS_25*max_lambda_22 + CF_25) + ((1.0/2.0))*(-CS_23*max_lambda_22 +
            CF_23)*(-(1261.0/240.0)*(-CS_24*max_lambda_22 + CF_24) - (247.0/240.0)*(-CS_26*max_lambda_22 + CF_26) +
            ((547.0/480.0))*(-CS_23*max_lambda_22 + CF_23) + ((961.0/240.0))*(-CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-CS_24*max_lambda_22 + CF_24)*(-(2983.0/240.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((267.0/80.0))*(-CS_26*max_lambda_22 + CF_26) + ((3443.0/480.0))*(-CS_24*max_lambda_22 + CF_24));

       beta_2 = ((547.0/960.0))*((-CS_25*max_lambda_22 + CF_25)*(-CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_25*max_lambda_22 + CF_25) + ((3443.0/480.0))*(-CS_24*max_lambda_22 +
            CF_24))*(-CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
            CF_22)*(-(821.0/240.0)*(-CS_23*max_lambda_22 + CF_23) - (247.0/240.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((89.0/160.0))*(-CS_22*max_lambda_22 + CF_22) + ((267.0/80.0))*(-CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23)*(-(2983.0/240.0)*(-CS_24*max_lambda_22 + CF_24) +
            ((961.0/240.0))*(-CS_25*max_lambda_22 + CF_25) + ((2843.0/480.0))*(-CS_23*max_lambda_22 + CF_23));

       beta_3 = ((2107.0/960.0))*((-CS_24*max_lambda_22 + CF_24)*(-CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_24*max_lambda_22 + CF_24) + ((11003.0/480.0))*(-CS_23*max_lambda_22 +
            CF_23))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_21*max_lambda_22 +
            CF_21)*(-(309.0/80.0)*(-CS_24*max_lambda_22 + CF_24) + ((547.0/480.0))*(-CS_21*max_lambda_22 + CF_21) +
            ((2321.0/240.0))*(-CS_23*max_lambda_22 + CF_23)) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
            CF_22)*(-(8623.0/240.0)*(-CS_23*max_lambda_22 + CF_23) - (647.0/80.0)*(-CS_21*max_lambda_22 + CF_21) +
            ((3521.0/240.0))*(-CS_24*max_lambda_22 + CF_24) + ((7043.0/480.0))*(-CS_22*max_lambda_22 + CF_22));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj2 = fmax(rj2, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_2 = (-(23.0/24.0)*(-CS_25*max_lambda_22 + CF_25) - (1.0/8.0)*(-CS_27*max_lambda_22 + CF_27) +
            ((13.0/24.0))*(-CS_26*max_lambda_22 + CF_26) + ((25.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_0 +
            (-(5.0/24.0)*(-CS_22*max_lambda_22 + CF_22) + ((1.0/8.0))*(-CS_24*max_lambda_22 + CF_24) +
            ((1.0/24.0))*(-CS_21*max_lambda_22 + CF_21) + ((13.0/24.0))*(-CS_23*max_lambda_22 + CF_23))*omega_3 +
            (-(5.0/24.0)*(-CS_25*max_lambda_22 + CF_25) + ((1.0/8.0))*(-CS_23*max_lambda_22 + CF_23) +
            ((1.0/24.0))*(-CS_26*max_lambda_22 + CF_26) + ((13.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_1 +
            (-(1.0/24.0)*(-CS_22*max_lambda_22 + CF_22) - (1.0/24.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((7.0/24.0))*(-CS_23*max_lambda_22 + CF_23) + ((7.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_2 +
            Recon_2;

       beta_0 = ((547.0/960.0))*((CS_36*max_lambda_33 + CF_36)*(CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_36*max_lambda_33 + CF_36) + ((7043.0/480.0))*(CS_35*max_lambda_33 +
            CF_35))*(CS_35*max_lambda_33 + CF_35) + ((1.0/2.0))*(CS_33*max_lambda_33 +
            CF_33)*(-(1567.0/80.0)*(CS_34*max_lambda_33 + CF_34) - (309.0/80.0)*(CS_36*max_lambda_33 + CF_36) +
            ((2107.0/480.0))*(CS_33*max_lambda_33 + CF_33) + ((3521.0/240.0))*(CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34)*(-(8623.0/240.0)*(CS_35*max_lambda_33 + CF_35) +
            ((2321.0/240.0))*(CS_36*max_lambda_33 + CF_36) + ((11003.0/480.0))*(CS_34*max_lambda_33 + CF_34));

       beta_1 = ((89.0/320.0))*((CS_35*max_lambda_33 + CF_35)*(CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_35*max_lambda_33 + CF_35) + ((2843.0/480.0))*(CS_34*max_lambda_33 +
            CF_34))*(CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(1261.0/240.0)*(CS_33*max_lambda_33 + CF_33) - (247.0/240.0)*(CS_35*max_lambda_33 + CF_35) +
            ((547.0/480.0))*(CS_32*max_lambda_33 + CF_32) + ((961.0/240.0))*(CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(CS_33*max_lambda_33 + CF_33)*(-(2983.0/240.0)*(CS_34*max_lambda_33 + CF_34) +
            ((267.0/80.0))*(CS_35*max_lambda_33 + CF_35) + ((3443.0/480.0))*(CS_33*max_lambda_33 + CF_33));

       beta_2 = ((547.0/960.0))*((CS_34*max_lambda_33 + CF_34)*(CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_34*max_lambda_33 + CF_34) + ((3443.0/480.0))*(CS_33*max_lambda_33 +
            CF_33))*(CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(CS_31*max_lambda_33 +
            CF_31)*(-(247.0/240.0)*(CS_34*max_lambda_33 + CF_34) + ((89.0/160.0))*(CS_31*max_lambda_33 + CF_31) +
            ((267.0/80.0))*(CS_33*max_lambda_33 + CF_33)) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(2983.0/240.0)*(CS_33*max_lambda_33 + CF_33) - (821.0/240.0)*(CS_31*max_lambda_33 + CF_31) +
            ((961.0/240.0))*(CS_34*max_lambda_33 + CF_34) + ((2843.0/480.0))*(CS_32*max_lambda_33 + CF_32));

       beta_3 = ((2107.0/960.0))*((CS_33*max_lambda_33 + CF_33)*(CS_33*max_lambda_33 + CF_33)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_33*max_lambda_33 + CF_33) + ((547.0/480.0))*(CS_30*max_lambda_33 +
            CF_30))*(CS_30*max_lambda_33 + CF_30) + ((1.0/2.0))*(CS_31*max_lambda_33 +
            CF_31)*(-(647.0/80.0)*(CS_30*max_lambda_33 + CF_30) + ((3521.0/240.0))*(CS_33*max_lambda_33 + CF_33) +
            ((7043.0/480.0))*(CS_31*max_lambda_33 + CF_31)) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(8623.0/240.0)*(CS_31*max_lambda_33 + CF_31) - (1567.0/80.0)*(CS_33*max_lambda_33 + CF_33) +
            ((2321.0/240.0))*(CS_30*max_lambda_33 + CF_30) + ((11003.0/480.0))*(CS_32*max_lambda_33 + CF_32));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj3 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_3 = (-(23.0/24.0)*(CS_32*max_lambda_33 + CF_32) - (1.0/8.0)*(CS_30*max_lambda_33 + CF_30) +
            ((13.0/24.0))*(CS_31*max_lambda_33 + CF_31) + ((25.0/24.0))*(CS_33*max_lambda_33 + CF_33))*omega_3 +
            (-(5.0/24.0)*(CS_32*max_lambda_33 + CF_32) + ((1.0/8.0))*(CS_34*max_lambda_33 + CF_34) +
            ((1.0/24.0))*(CS_31*max_lambda_33 + CF_31) + ((13.0/24.0))*(CS_33*max_lambda_33 + CF_33))*omega_2 +
            (-(5.0/24.0)*(CS_35*max_lambda_33 + CF_35) + ((1.0/8.0))*(CS_33*max_lambda_33 + CF_33) +
            ((1.0/24.0))*(CS_36*max_lambda_33 + CF_36) + ((13.0/24.0))*(CS_34*max_lambda_33 + CF_34))*omega_0 +
            (-(1.0/24.0)*(CS_32*max_lambda_33 + CF_32) - (1.0/24.0)*(CS_35*max_lambda_33 + CF_35) +
            ((7.0/24.0))*(CS_33*max_lambda_33 + CF_33) + ((7.0/24.0))*(CS_34*max_lambda_33 + CF_34))*omega_1 + Recon_3;

       beta_0 = ((547.0/960.0))*((-CS_37*max_lambda_33 + CF_37)*(-CS_37*max_lambda_33 + CF_37)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_37*max_lambda_33 + CF_37) + ((7043.0/480.0))*(-CS_36*max_lambda_33 +
            CF_36))*(-CS_36*max_lambda_33 + CF_36) + ((1.0/2.0))*(-CS_34*max_lambda_33 +
            CF_34)*(-(1567.0/80.0)*(-CS_35*max_lambda_33 + CF_35) - (309.0/80.0)*(-CS_37*max_lambda_33 + CF_37) +
            ((2107.0/480.0))*(-CS_34*max_lambda_33 + CF_34) + ((3521.0/240.0))*(-CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35)*(-(8623.0/240.0)*(-CS_36*max_lambda_33 + CF_36) +
            ((2321.0/240.0))*(-CS_37*max_lambda_33 + CF_37) + ((11003.0/480.0))*(-CS_35*max_lambda_33 + CF_35));

       beta_1 = ((89.0/320.0))*((-CS_36*max_lambda_33 + CF_36)*(-CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_36*max_lambda_33 + CF_36) + ((2843.0/480.0))*(-CS_35*max_lambda_33 +
            CF_35))*(-CS_35*max_lambda_33 + CF_35) + ((1.0/2.0))*(-CS_33*max_lambda_33 +
            CF_33)*(-(1261.0/240.0)*(-CS_34*max_lambda_33 + CF_34) - (247.0/240.0)*(-CS_36*max_lambda_33 + CF_36) +
            ((547.0/480.0))*(-CS_33*max_lambda_33 + CF_33) + ((961.0/240.0))*(-CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-CS_34*max_lambda_33 + CF_34)*(-(2983.0/240.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((267.0/80.0))*(-CS_36*max_lambda_33 + CF_36) + ((3443.0/480.0))*(-CS_34*max_lambda_33 + CF_34));

       beta_2 = ((547.0/960.0))*((-CS_35*max_lambda_33 + CF_35)*(-CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_35*max_lambda_33 + CF_35) + ((3443.0/480.0))*(-CS_34*max_lambda_33 +
            CF_34))*(-CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
            CF_32)*(-(821.0/240.0)*(-CS_33*max_lambda_33 + CF_33) - (247.0/240.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((89.0/160.0))*(-CS_32*max_lambda_33 + CF_32) + ((267.0/80.0))*(-CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33)*(-(2983.0/240.0)*(-CS_34*max_lambda_33 + CF_34) +
            ((961.0/240.0))*(-CS_35*max_lambda_33 + CF_35) + ((2843.0/480.0))*(-CS_33*max_lambda_33 + CF_33));

       beta_3 = ((2107.0/960.0))*((-CS_34*max_lambda_33 + CF_34)*(-CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_34*max_lambda_33 + CF_34) + ((11003.0/480.0))*(-CS_33*max_lambda_33 +
            CF_33))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_31*max_lambda_33 +
            CF_31)*(-(309.0/80.0)*(-CS_34*max_lambda_33 + CF_34) + ((547.0/480.0))*(-CS_31*max_lambda_33 + CF_31) +
            ((2321.0/240.0))*(-CS_33*max_lambda_33 + CF_33)) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
            CF_32)*(-(8623.0/240.0)*(-CS_33*max_lambda_33 + CF_33) - (647.0/80.0)*(-CS_31*max_lambda_33 + CF_31) +
            ((3521.0/240.0))*(-CS_34*max_lambda_33 + CF_34) + ((7043.0/480.0))*(-CS_32*max_lambda_33 + CF_32));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj3 = fmax(rj3, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_3 = (-(23.0/24.0)*(-CS_35*max_lambda_33 + CF_35) - (1.0/8.0)*(-CS_37*max_lambda_33 + CF_37) +
            ((13.0/24.0))*(-CS_36*max_lambda_33 + CF_36) + ((25.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_0 +
            (-(5.0/24.0)*(-CS_32*max_lambda_33 + CF_32) + ((1.0/8.0))*(-CS_34*max_lambda_33 + CF_34) +
            ((1.0/24.0))*(-CS_31*max_lambda_33 + CF_31) + ((13.0/24.0))*(-CS_33*max_lambda_33 + CF_33))*omega_3 +
            (-(5.0/24.0)*(-CS_35*max_lambda_33 + CF_35) + ((1.0/8.0))*(-CS_33*max_lambda_33 + CF_33) +
            ((1.0/24.0))*(-CS_36*max_lambda_33 + CF_36) + ((13.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_1 +
            (-(1.0/24.0)*(-CS_32*max_lambda_33 + CF_32) - (1.0/24.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((7.0/24.0))*(-CS_33*max_lambda_33 + CF_33) + ((7.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_2 +
            Recon_3;

       beta_0 = ((547.0/960.0))*((CS_46*max_lambda_44 + CF_46)*(CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_46*max_lambda_44 + CF_46) + ((7043.0/480.0))*(CS_45*max_lambda_44 +
            CF_45))*(CS_45*max_lambda_44 + CF_45) + ((1.0/2.0))*(CS_43*max_lambda_44 +
            CF_43)*(-(1567.0/80.0)*(CS_44*max_lambda_44 + CF_44) - (309.0/80.0)*(CS_46*max_lambda_44 + CF_46) +
            ((2107.0/480.0))*(CS_43*max_lambda_44 + CF_43) + ((3521.0/240.0))*(CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44)*(-(8623.0/240.0)*(CS_45*max_lambda_44 + CF_45) +
            ((2321.0/240.0))*(CS_46*max_lambda_44 + CF_46) + ((11003.0/480.0))*(CS_44*max_lambda_44 + CF_44));

       beta_1 = ((89.0/320.0))*((CS_45*max_lambda_44 + CF_45)*(CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_45*max_lambda_44 + CF_45) + ((2843.0/480.0))*(CS_44*max_lambda_44 +
            CF_44))*(CS_44*max_lambda_44 + CF_44) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(1261.0/240.0)*(CS_43*max_lambda_44 + CF_43) - (247.0/240.0)*(CS_45*max_lambda_44 + CF_45) +
            ((547.0/480.0))*(CS_42*max_lambda_44 + CF_42) + ((961.0/240.0))*(CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(CS_43*max_lambda_44 + CF_43)*(-(2983.0/240.0)*(CS_44*max_lambda_44 + CF_44) +
            ((267.0/80.0))*(CS_45*max_lambda_44 + CF_45) + ((3443.0/480.0))*(CS_43*max_lambda_44 + CF_43));

       beta_2 = ((547.0/960.0))*((CS_44*max_lambda_44 + CF_44)*(CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_44*max_lambda_44 + CF_44) + ((3443.0/480.0))*(CS_43*max_lambda_44 +
            CF_43))*(CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(CS_41*max_lambda_44 +
            CF_41)*(-(247.0/240.0)*(CS_44*max_lambda_44 + CF_44) + ((89.0/160.0))*(CS_41*max_lambda_44 + CF_41) +
            ((267.0/80.0))*(CS_43*max_lambda_44 + CF_43)) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(2983.0/240.0)*(CS_43*max_lambda_44 + CF_43) - (821.0/240.0)*(CS_41*max_lambda_44 + CF_41) +
            ((961.0/240.0))*(CS_44*max_lambda_44 + CF_44) + ((2843.0/480.0))*(CS_42*max_lambda_44 + CF_42));

       beta_3 = ((2107.0/960.0))*((CS_43*max_lambda_44 + CF_43)*(CS_43*max_lambda_44 + CF_43)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_43*max_lambda_44 + CF_43) + ((547.0/480.0))*(CS_40*max_lambda_44 +
            CF_40))*(CS_40*max_lambda_44 + CF_40) + ((1.0/2.0))*(CS_41*max_lambda_44 +
            CF_41)*(-(647.0/80.0)*(CS_40*max_lambda_44 + CF_40) + ((3521.0/240.0))*(CS_43*max_lambda_44 + CF_43) +
            ((7043.0/480.0))*(CS_41*max_lambda_44 + CF_41)) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(8623.0/240.0)*(CS_41*max_lambda_44 + CF_41) - (1567.0/80.0)*(CS_43*max_lambda_44 + CF_43) +
            ((2321.0/240.0))*(CS_40*max_lambda_44 + CF_40) + ((11003.0/480.0))*(CS_42*max_lambda_44 + CF_42));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj4 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_4 = (-(23.0/24.0)*(CS_42*max_lambda_44 + CF_42) - (1.0/8.0)*(CS_40*max_lambda_44 + CF_40) +
            ((13.0/24.0))*(CS_41*max_lambda_44 + CF_41) + ((25.0/24.0))*(CS_43*max_lambda_44 + CF_43))*omega_3 +
            (-(5.0/24.0)*(CS_42*max_lambda_44 + CF_42) + ((1.0/8.0))*(CS_44*max_lambda_44 + CF_44) +
            ((1.0/24.0))*(CS_41*max_lambda_44 + CF_41) + ((13.0/24.0))*(CS_43*max_lambda_44 + CF_43))*omega_2 +
            (-(5.0/24.0)*(CS_45*max_lambda_44 + CF_45) + ((1.0/8.0))*(CS_43*max_lambda_44 + CF_43) +
            ((1.0/24.0))*(CS_46*max_lambda_44 + CF_46) + ((13.0/24.0))*(CS_44*max_lambda_44 + CF_44))*omega_0 +
            (-(1.0/24.0)*(CS_42*max_lambda_44 + CF_42) - (1.0/24.0)*(CS_45*max_lambda_44 + CF_45) +
            ((7.0/24.0))*(CS_43*max_lambda_44 + CF_43) + ((7.0/24.0))*(CS_44*max_lambda_44 + CF_44))*omega_1 + Recon_4;

       beta_0 = ((547.0/960.0))*((-CS_47*max_lambda_44 + CF_47)*(-CS_47*max_lambda_44 + CF_47)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_47*max_lambda_44 + CF_47) + ((7043.0/480.0))*(-CS_46*max_lambda_44 +
            CF_46))*(-CS_46*max_lambda_44 + CF_46) + ((1.0/2.0))*(-CS_44*max_lambda_44 +
            CF_44)*(-(1567.0/80.0)*(-CS_45*max_lambda_44 + CF_45) - (309.0/80.0)*(-CS_47*max_lambda_44 + CF_47) +
            ((2107.0/480.0))*(-CS_44*max_lambda_44 + CF_44) + ((3521.0/240.0))*(-CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45)*(-(8623.0/240.0)*(-CS_46*max_lambda_44 + CF_46) +
            ((2321.0/240.0))*(-CS_47*max_lambda_44 + CF_47) + ((11003.0/480.0))*(-CS_45*max_lambda_44 + CF_45));

       beta_1 = ((89.0/320.0))*((-CS_46*max_lambda_44 + CF_46)*(-CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_46*max_lambda_44 + CF_46) + ((2843.0/480.0))*(-CS_45*max_lambda_44 +
            CF_45))*(-CS_45*max_lambda_44 + CF_45) + ((1.0/2.0))*(-CS_43*max_lambda_44 +
            CF_43)*(-(1261.0/240.0)*(-CS_44*max_lambda_44 + CF_44) - (247.0/240.0)*(-CS_46*max_lambda_44 + CF_46) +
            ((547.0/480.0))*(-CS_43*max_lambda_44 + CF_43) + ((961.0/240.0))*(-CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-CS_44*max_lambda_44 + CF_44)*(-(2983.0/240.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((267.0/80.0))*(-CS_46*max_lambda_44 + CF_46) + ((3443.0/480.0))*(-CS_44*max_lambda_44 + CF_44));

       beta_2 = ((547.0/960.0))*((-CS_45*max_lambda_44 + CF_45)*(-CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_45*max_lambda_44 + CF_45) + ((3443.0/480.0))*(-CS_44*max_lambda_44 +
            CF_44))*(-CS_44*max_lambda_44 + CF_44) + ((1.0/2.0))*(-CS_42*max_lambda_44 +
            CF_42)*(-(821.0/240.0)*(-CS_43*max_lambda_44 + CF_43) - (247.0/240.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((89.0/160.0))*(-CS_42*max_lambda_44 + CF_42) + ((267.0/80.0))*(-CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43)*(-(2983.0/240.0)*(-CS_44*max_lambda_44 + CF_44) +
            ((961.0/240.0))*(-CS_45*max_lambda_44 + CF_45) + ((2843.0/480.0))*(-CS_43*max_lambda_44 + CF_43));

       beta_3 = ((2107.0/960.0))*((-CS_44*max_lambda_44 + CF_44)*(-CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_44*max_lambda_44 + CF_44) + ((11003.0/480.0))*(-CS_43*max_lambda_44 +
            CF_43))*(-CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(-CS_41*max_lambda_44 +
            CF_41)*(-(309.0/80.0)*(-CS_44*max_lambda_44 + CF_44) + ((547.0/480.0))*(-CS_41*max_lambda_44 + CF_41) +
            ((2321.0/240.0))*(-CS_43*max_lambda_44 + CF_43)) + ((1.0/2.0))*(-CS_42*max_lambda_44 +
            CF_42)*(-(8623.0/240.0)*(-CS_43*max_lambda_44 + CF_43) - (647.0/80.0)*(-CS_41*max_lambda_44 + CF_41) +
            ((3521.0/240.0))*(-CS_44*max_lambda_44 + CF_44) + ((7043.0/480.0))*(-CS_42*max_lambda_44 + CF_42));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj4 = fmax(rj4, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_4 = (-(23.0/24.0)*(-CS_45*max_lambda_44 + CF_45) - (1.0/8.0)*(-CS_47*max_lambda_44 + CF_47) +
            ((13.0/24.0))*(-CS_46*max_lambda_44 + CF_46) + ((25.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_0 +
            (-(5.0/24.0)*(-CS_42*max_lambda_44 + CF_42) + ((1.0/8.0))*(-CS_44*max_lambda_44 + CF_44) +
            ((1.0/24.0))*(-CS_41*max_lambda_44 + CF_41) + ((13.0/24.0))*(-CS_43*max_lambda_44 + CF_43))*omega_3 +
            (-(5.0/24.0)*(-CS_45*max_lambda_44 + CF_45) + ((1.0/8.0))*(-CS_43*max_lambda_44 + CF_43) +
            ((1.0/24.0))*(-CS_46*max_lambda_44 + CF_46) + ((13.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_1 +
            (-(1.0/24.0)*(-CS_42*max_lambda_44 + CF_42) - (1.0/24.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((7.0/24.0))*(-CS_43*max_lambda_44 + CF_43) + ((7.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_2 +
            Recon_4;

       Recon_0 = (((1.0/840.0))*(-533*CF_03 - 533*CF_04 - 29*CF_01 - 29*CF_06 + 3*CF_00 + 3*CF_07 + 139*CF_02 +
            139*CF_05) + Recon_0)*rj0;

       Recon_1 = (((1.0/840.0))*(-533*CF_13 - 533*CF_14 - 29*CF_11 - 29*CF_16 + 3*CF_10 + 3*CF_17 + 139*CF_12 +
            139*CF_15) + Recon_1)*rj1;

       Recon_2 = (((1.0/840.0))*(-533*CF_23 - 533*CF_24 - 29*CF_21 - 29*CF_26 + 3*CF_20 + 3*CF_27 + 139*CF_22 +
            139*CF_25) + Recon_2)*rj2;

       Recon_3 = (((1.0/840.0))*(-533*CF_33 - 533*CF_34 - 29*CF_31 - 29*CF_36 + 3*CF_30 + 3*CF_37 + 139*CF_32 +
            139*CF_35) + Recon_3)*rj3;

       Recon_4 = (((1.0/840.0))*(-533*CF_43 - 533*CF_44 - 29*CF_41 - 29*CF_46 + 3*CF_40 + 3*CF_47 + 139*CF_42 +
            139*CF_45) + Recon_4)*rj4;

       Residual0_B0(0,0,0) = 0.707106781186547*AVG_1_rho*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_1_rho*Recon_4*inv_AVG_a + Recon_1;

       Residual1_B0(0,0,0) = AVG_1_rho*Recon_2 + AVG_1_u0*Recon_1 +
            0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_4*inv_AVG_a;

       Residual2_B0(0,0,0) = AVG_1_u1*Recon_1 + 0.707106781186547*(-AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_4*inv_AVG_a +
            0.707106781186547*(AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a;

       Residual3_B0(0,0,0) = AVG_1_u2*Recon_1 - AVG_1_rho*Recon_0 +
            0.707106781186547*AVG_1_rho*AVG_1_u2*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_1_rho*AVG_1_u2*Recon_4*inv_AVG_a;

       Residual4_B0(0,0,0) = (((1.0/2.0))*(AVG_1_u0*AVG_1_u0) + ((1.0/2.0))*(AVG_1_u1*AVG_1_u1) +
            ((1.0/2.0))*(AVG_1_u2*AVG_1_u2))*Recon_1 + AVG_1_rho*AVG_1_u0*Recon_2 - AVG_1_rho*AVG_1_u2*Recon_0 +
            0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) +
            (AVG_1_u2*AVG_1_u2))*gamma_m1)*invgamma_m1 + AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a +
            0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) +
            (AVG_1_u2*AVG_1_u2))*gamma_m1)*invgamma_m1 - AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_4*inv_AVG_a;

   }

   else{

      Residual0_B0(0,0,0) = 0.0;

      Residual1_B0(0,0,0) = 0.0;

      Residual2_B0(0,0,0) = 0.0;

      Residual3_B0(0,0,0) = 0.0;

      Residual4_B0(0,0,0) = 0.0;

   }

}

 void opensbliblock00Kernel060(const ACC<double> &a_B0, const ACC<double> &kappa_B0, const ACC<double> &p_B0, const
ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double> &rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double>
&rhou2_RKold_B0)
{
   double AVG_2_2_LEV_00 = 0.0;
   double AVG_2_2_LEV_02 = 0.0;
   double AVG_2_2_LEV_10 = 0.0;
   double AVG_2_2_LEV_11 = 0.0;
   double AVG_2_2_LEV_20 = 0.0;
   double AVG_2_2_LEV_21 = 0.0;
   double AVG_2_2_LEV_22 = 0.0;
   double AVG_2_2_LEV_23 = 0.0;
   double AVG_2_2_LEV_24 = 0.0;
   double AVG_2_2_LEV_30 = 0.0;
   double AVG_2_2_LEV_31 = 0.0;
   double AVG_2_2_LEV_32 = 0.0;
   double AVG_2_2_LEV_33 = 0.0;
   double AVG_2_2_LEV_34 = 0.0;
   double AVG_2_2_LEV_40 = 0.0;
   double AVG_2_2_LEV_41 = 0.0;
   double AVG_2_2_LEV_42 = 0.0;
   double AVG_2_2_LEV_43 = 0.0;
   double AVG_2_2_LEV_44 = 0.0;
   double AVG_2_a = 0.0;
   double AVG_2_inv_rho = 0.0;
   double AVG_2_rho = 0.0;
   double AVG_2_u0 = 0.0;
   double AVG_2_u1 = 0.0;
   double AVG_2_u2 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_06 = 0.0;
   double CF_07 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_16 = 0.0;
   double CF_17 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_26 = 0.0;
   double CF_27 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CF_36 = 0.0;
   double CF_37 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CF_44 = 0.0;
   double CF_45 = 0.0;
   double CF_46 = 0.0;
   double CF_47 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_06 = 0.0;
   double CS_07 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_16 = 0.0;
   double CS_17 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_26 = 0.0;
   double CS_27 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double CS_36 = 0.0;
   double CS_37 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double CS_44 = 0.0;
   double CS_45 = 0.0;
   double CS_46 = 0.0;
   double CS_47 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double alpha_3 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double beta_3 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double omega_2 = 0.0;
   double omega_3 = 0.0;
   double rj0 = 0.0;
   double rj1 = 0.0;
   double rj2 = 0.0;
   double rj3 = 0.0;
   double rj4 = 0.0;
    if (fmax(kappa_B0(0,0,1), fmax(kappa_B0(0,0,-3), fmax(kappa_B0(0,0,-1), fmax(kappa_B0(0,0,2), fmax(kappa_B0(0,0,0),
      kappa_B0(0,0,-2)))))) > Ducros_check){

      AVG_2_rho = sqrt((rho_B0(0,0,0)*rho_B0(0,0,1)));

      AVG_2_inv_rho = 1.0/((sqrt(rho_B0(0,0,0)) + sqrt(rho_B0(0,0,1))));

      AVG_2_u0 = (sqrt(rho_B0(0,0,0))*u0_B0(0,0,0) + sqrt(rho_B0(0,0,1))*u0_B0(0,0,1))*AVG_2_inv_rho;

      AVG_2_u1 = (sqrt(rho_B0(0,0,0))*u1_B0(0,0,0) + sqrt(rho_B0(0,0,1))*u1_B0(0,0,1))*AVG_2_inv_rho;

      AVG_2_u2 = (sqrt(rho_B0(0,0,0))*u2_B0(0,0,0) + sqrt(rho_B0(0,0,1))*u2_B0(0,0,1))*AVG_2_inv_rho;

       AVG_2_a = sqrt(((-(1.0/2.0)*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) + (AVG_2_u2*AVG_2_u2)) + ((p_B0(0,0,0) +
            rhoE_B0(0,0,0))/sqrt(rho_B0(0,0,0)) + (p_B0(0,0,1) +
            rhoE_B0(0,0,1))/sqrt(rho_B0(0,0,1)))*AVG_2_inv_rho)*gamma_m1));

      inv_AVG_a = 1.0/(AVG_2_a);

      inv_AVG_rho = 1.0/(AVG_2_rho);

      AVG_2_2_LEV_00 = -AVG_2_u1*inv_AVG_rho;

      AVG_2_2_LEV_02 = inv_AVG_rho;

      AVG_2_2_LEV_10 = AVG_2_u0*inv_AVG_rho;

      AVG_2_2_LEV_11 = -inv_AVG_rho;

       AVG_2_2_LEV_20 = -(1.0/2.0)*(-2 - (AVG_2_u0*AVG_2_u0)*(inv_AVG_a*inv_AVG_a) -
            (AVG_2_u1*AVG_2_u1)*(inv_AVG_a*inv_AVG_a) - (AVG_2_u2*AVG_2_u2)*(inv_AVG_a*inv_AVG_a) +
            (AVG_2_u0*AVG_2_u0)*(inv_AVG_a*inv_AVG_a)*gama + (AVG_2_u1*AVG_2_u1)*(inv_AVG_a*inv_AVG_a)*gama +
            (AVG_2_u2*AVG_2_u2)*(inv_AVG_a*inv_AVG_a)*gama);

      AVG_2_2_LEV_21 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_2_u0;

      AVG_2_2_LEV_22 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_2_u1;

      AVG_2_2_LEV_23 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_2_u2;

      AVG_2_2_LEV_24 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

       AVG_2_2_LEV_30 = -0.353553390593274*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) + (AVG_2_u2*AVG_2_u2) -
            (AVG_2_u0*AVG_2_u0)*gama - (AVG_2_u1*AVG_2_u1)*gama - (AVG_2_u2*AVG_2_u2)*gama +
            2*AVG_2_a*AVG_2_u2)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_31 = -0.707106781186547*gamma_m1*AVG_2_u0*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_32 = -0.707106781186547*gamma_m1*AVG_2_u1*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_33 = 0.707106781186547*(-gama*AVG_2_u2 + AVG_2_a + AVG_2_u2)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_34 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       AVG_2_2_LEV_40 = 0.353553390593274*(-(AVG_2_u0*AVG_2_u0) - (AVG_2_u1*AVG_2_u1) - (AVG_2_u2*AVG_2_u2) +
            (AVG_2_u0*AVG_2_u0)*gama + (AVG_2_u1*AVG_2_u1)*gama + (AVG_2_u2*AVG_2_u2)*gama +
            2*AVG_2_a*AVG_2_u2)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_41 = -0.707106781186547*gamma_m1*AVG_2_u0*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_42 = -0.707106781186547*gamma_m1*AVG_2_u1*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_43 = -0.707106781186547*(-AVG_2_u2 + gama*AVG_2_u2 + AVG_2_a)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_44 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

      CF_00 = (rho_B0(0,0,-3)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-3)*AVG_2_2_LEV_02)*u2_B0(0,0,-3);

      CF_10 = (rho_B0(0,0,-3)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-3)*AVG_2_2_LEV_11)*u2_B0(0,0,-3);

       CF_20 = p_B0(0,0,-3)*AVG_2_2_LEV_23 + p_B0(0,0,-3)*u2_B0(0,0,-3)*AVG_2_2_LEV_24 +
            u2_B0(0,0,-3)*rho_B0(0,0,-3)*AVG_2_2_LEV_20 + u2_B0(0,0,-3)*rhoE_B0(0,0,-3)*AVG_2_2_LEV_24 +
            u2_B0(0,0,-3)*rhou0_B0(0,0,-3)*AVG_2_2_LEV_21 + u2_B0(0,0,-3)*rhou1_B0(0,0,-3)*AVG_2_2_LEV_22 +
            u2_B0(0,0,-3)*rhou2_B0(0,0,-3)*AVG_2_2_LEV_23;

       CF_30 = p_B0(0,0,-3)*AVG_2_2_LEV_33 + p_B0(0,0,-3)*u2_B0(0,0,-3)*AVG_2_2_LEV_34 +
            u2_B0(0,0,-3)*rho_B0(0,0,-3)*AVG_2_2_LEV_30 + u2_B0(0,0,-3)*rhoE_B0(0,0,-3)*AVG_2_2_LEV_34 +
            u2_B0(0,0,-3)*rhou0_B0(0,0,-3)*AVG_2_2_LEV_31 + u2_B0(0,0,-3)*rhou1_B0(0,0,-3)*AVG_2_2_LEV_32 +
            u2_B0(0,0,-3)*rhou2_B0(0,0,-3)*AVG_2_2_LEV_33;

       CF_40 = p_B0(0,0,-3)*AVG_2_2_LEV_43 + p_B0(0,0,-3)*u2_B0(0,0,-3)*AVG_2_2_LEV_44 +
            u2_B0(0,0,-3)*rho_B0(0,0,-3)*AVG_2_2_LEV_40 + u2_B0(0,0,-3)*rhoE_B0(0,0,-3)*AVG_2_2_LEV_44 +
            u2_B0(0,0,-3)*rhou0_B0(0,0,-3)*AVG_2_2_LEV_41 + u2_B0(0,0,-3)*rhou1_B0(0,0,-3)*AVG_2_2_LEV_42 +
            u2_B0(0,0,-3)*rhou2_B0(0,0,-3)*AVG_2_2_LEV_43;

      CS_00 = rho_B0(0,0,-3)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-3)*AVG_2_2_LEV_02;

      CS_10 = rho_B0(0,0,-3)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-3)*AVG_2_2_LEV_11;

       CS_20 = rho_B0(0,0,-3)*AVG_2_2_LEV_20 + rhoE_B0(0,0,-3)*AVG_2_2_LEV_24 + rhou0_B0(0,0,-3)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,-3)*AVG_2_2_LEV_22 + rhou2_B0(0,0,-3)*AVG_2_2_LEV_23;

       CS_30 = rho_B0(0,0,-3)*AVG_2_2_LEV_30 + rhoE_B0(0,0,-3)*AVG_2_2_LEV_34 + rhou0_B0(0,0,-3)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,-3)*AVG_2_2_LEV_32 + rhou2_B0(0,0,-3)*AVG_2_2_LEV_33;

       CS_40 = rho_B0(0,0,-3)*AVG_2_2_LEV_40 + rhoE_B0(0,0,-3)*AVG_2_2_LEV_44 + rhou0_B0(0,0,-3)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,-3)*AVG_2_2_LEV_42 + rhou2_B0(0,0,-3)*AVG_2_2_LEV_43;

      CF_01 = (rho_B0(0,0,-2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-2)*AVG_2_2_LEV_02)*u2_B0(0,0,-2);

      CF_11 = (rho_B0(0,0,-2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_11)*u2_B0(0,0,-2);

       CF_21 = p_B0(0,0,-2)*AVG_2_2_LEV_23 + p_B0(0,0,-2)*u2_B0(0,0,-2)*AVG_2_2_LEV_24 +
            u2_B0(0,0,-2)*rho_B0(0,0,-2)*AVG_2_2_LEV_20 + u2_B0(0,0,-2)*rhoE_B0(0,0,-2)*AVG_2_2_LEV_24 +
            u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_21 + u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_22 +
            u2_B0(0,0,-2)*rhou2_B0(0,0,-2)*AVG_2_2_LEV_23;

       CF_31 = p_B0(0,0,-2)*AVG_2_2_LEV_33 + p_B0(0,0,-2)*u2_B0(0,0,-2)*AVG_2_2_LEV_34 +
            u2_B0(0,0,-2)*rho_B0(0,0,-2)*AVG_2_2_LEV_30 + u2_B0(0,0,-2)*rhoE_B0(0,0,-2)*AVG_2_2_LEV_34 +
            u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_31 + u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_32 +
            u2_B0(0,0,-2)*rhou2_B0(0,0,-2)*AVG_2_2_LEV_33;

       CF_41 = p_B0(0,0,-2)*AVG_2_2_LEV_43 + p_B0(0,0,-2)*u2_B0(0,0,-2)*AVG_2_2_LEV_44 +
            u2_B0(0,0,-2)*rho_B0(0,0,-2)*AVG_2_2_LEV_40 + u2_B0(0,0,-2)*rhoE_B0(0,0,-2)*AVG_2_2_LEV_44 +
            u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_41 + u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_42 +
            u2_B0(0,0,-2)*rhou2_B0(0,0,-2)*AVG_2_2_LEV_43;

      CS_01 = rho_B0(0,0,-2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-2)*AVG_2_2_LEV_02;

      CS_11 = rho_B0(0,0,-2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_11;

       CS_21 = rho_B0(0,0,-2)*AVG_2_2_LEV_20 + rhoE_B0(0,0,-2)*AVG_2_2_LEV_24 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,-2)*AVG_2_2_LEV_22 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_23;

       CS_31 = rho_B0(0,0,-2)*AVG_2_2_LEV_30 + rhoE_B0(0,0,-2)*AVG_2_2_LEV_34 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,-2)*AVG_2_2_LEV_32 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_33;

       CS_41 = rho_B0(0,0,-2)*AVG_2_2_LEV_40 + rhoE_B0(0,0,-2)*AVG_2_2_LEV_44 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,-2)*AVG_2_2_LEV_42 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_43;

      CF_02 = (rho_B0(0,0,-1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-1)*AVG_2_2_LEV_02)*u2_B0(0,0,-1);

      CF_12 = (rho_B0(0,0,-1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_11)*u2_B0(0,0,-1);

       CF_22 = p_B0(0,0,-1)*AVG_2_2_LEV_23 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_24 +
            u2_B0(0,0,-1)*rho_B0(0,0,-1)*AVG_2_2_LEV_20 + u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_24 +
            u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_21 + u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_22 +
            u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_23;

       CF_32 = p_B0(0,0,-1)*AVG_2_2_LEV_33 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_34 +
            u2_B0(0,0,-1)*rho_B0(0,0,-1)*AVG_2_2_LEV_30 + u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_34 +
            u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_31 + u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_32 +
            u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_33;

       CF_42 = p_B0(0,0,-1)*AVG_2_2_LEV_43 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_44 +
            u2_B0(0,0,-1)*rho_B0(0,0,-1)*AVG_2_2_LEV_40 + u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_44 +
            u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_41 + u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_42 +
            u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_43;

      CS_02 = rho_B0(0,0,-1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-1)*AVG_2_2_LEV_02;

      CS_12 = rho_B0(0,0,-1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_11;

       CS_22 = rho_B0(0,0,-1)*AVG_2_2_LEV_20 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_24 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,-1)*AVG_2_2_LEV_22 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_23;

       CS_32 = rho_B0(0,0,-1)*AVG_2_2_LEV_30 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_34 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,-1)*AVG_2_2_LEV_32 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_33;

       CS_42 = rho_B0(0,0,-1)*AVG_2_2_LEV_40 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_44 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,-1)*AVG_2_2_LEV_42 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_43;

      CF_03 = (rho_B0(0,0,0)*AVG_2_2_LEV_00 + rhou1_B0(0,0,0)*AVG_2_2_LEV_02)*u2_B0(0,0,0);

      CF_13 = (rho_B0(0,0,0)*AVG_2_2_LEV_10 + rhou0_B0(0,0,0)*AVG_2_2_LEV_11)*u2_B0(0,0,0);

       CF_23 = p_B0(0,0,0)*AVG_2_2_LEV_23 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_24 +
            u2_B0(0,0,0)*rho_B0(0,0,0)*AVG_2_2_LEV_20 + u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_24 +
            u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_21 + u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_22 +
            u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_23;

       CF_33 = p_B0(0,0,0)*AVG_2_2_LEV_33 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_34 +
            u2_B0(0,0,0)*rho_B0(0,0,0)*AVG_2_2_LEV_30 + u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_34 +
            u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_31 + u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_32 +
            u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_33;

       CF_43 = p_B0(0,0,0)*AVG_2_2_LEV_43 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_44 +
            u2_B0(0,0,0)*rho_B0(0,0,0)*AVG_2_2_LEV_40 + u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_44 +
            u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_41 + u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_42 +
            u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_43;

      CS_03 = rho_B0(0,0,0)*AVG_2_2_LEV_00 + rhou1_B0(0,0,0)*AVG_2_2_LEV_02;

      CS_13 = rho_B0(0,0,0)*AVG_2_2_LEV_10 + rhou0_B0(0,0,0)*AVG_2_2_LEV_11;

       CS_23 = rho_B0(0,0,0)*AVG_2_2_LEV_20 + rhoE_B0(0,0,0)*AVG_2_2_LEV_24 + rhou0_B0(0,0,0)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,0)*AVG_2_2_LEV_22 + rhou2_B0(0,0,0)*AVG_2_2_LEV_23;

       CS_33 = rho_B0(0,0,0)*AVG_2_2_LEV_30 + rhoE_B0(0,0,0)*AVG_2_2_LEV_34 + rhou0_B0(0,0,0)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,0)*AVG_2_2_LEV_32 + rhou2_B0(0,0,0)*AVG_2_2_LEV_33;

       CS_43 = rho_B0(0,0,0)*AVG_2_2_LEV_40 + rhoE_B0(0,0,0)*AVG_2_2_LEV_44 + rhou0_B0(0,0,0)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,0)*AVG_2_2_LEV_42 + rhou2_B0(0,0,0)*AVG_2_2_LEV_43;

      CF_04 = (rho_B0(0,0,1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,1)*AVG_2_2_LEV_02)*u2_B0(0,0,1);

      CF_14 = (rho_B0(0,0,1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,1)*AVG_2_2_LEV_11)*u2_B0(0,0,1);

       CF_24 = p_B0(0,0,1)*AVG_2_2_LEV_23 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_24 +
            u2_B0(0,0,1)*rho_B0(0,0,1)*AVG_2_2_LEV_20 + u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_24 +
            u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_21 + u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_22 +
            u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_23;

       CF_34 = p_B0(0,0,1)*AVG_2_2_LEV_33 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_34 +
            u2_B0(0,0,1)*rho_B0(0,0,1)*AVG_2_2_LEV_30 + u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_34 +
            u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_31 + u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_32 +
            u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_33;

       CF_44 = p_B0(0,0,1)*AVG_2_2_LEV_43 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_44 +
            u2_B0(0,0,1)*rho_B0(0,0,1)*AVG_2_2_LEV_40 + u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_44 +
            u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_41 + u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_42 +
            u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_43;

      CS_04 = rho_B0(0,0,1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,1)*AVG_2_2_LEV_02;

      CS_14 = rho_B0(0,0,1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,1)*AVG_2_2_LEV_11;

       CS_24 = rho_B0(0,0,1)*AVG_2_2_LEV_20 + rhoE_B0(0,0,1)*AVG_2_2_LEV_24 + rhou0_B0(0,0,1)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,1)*AVG_2_2_LEV_22 + rhou2_B0(0,0,1)*AVG_2_2_LEV_23;

       CS_34 = rho_B0(0,0,1)*AVG_2_2_LEV_30 + rhoE_B0(0,0,1)*AVG_2_2_LEV_34 + rhou0_B0(0,0,1)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,1)*AVG_2_2_LEV_32 + rhou2_B0(0,0,1)*AVG_2_2_LEV_33;

       CS_44 = rho_B0(0,0,1)*AVG_2_2_LEV_40 + rhoE_B0(0,0,1)*AVG_2_2_LEV_44 + rhou0_B0(0,0,1)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,1)*AVG_2_2_LEV_42 + rhou2_B0(0,0,1)*AVG_2_2_LEV_43;

      CF_05 = (rho_B0(0,0,2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,2)*AVG_2_2_LEV_02)*u2_B0(0,0,2);

      CF_15 = (rho_B0(0,0,2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,2)*AVG_2_2_LEV_11)*u2_B0(0,0,2);

       CF_25 = p_B0(0,0,2)*AVG_2_2_LEV_23 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_24 +
            u2_B0(0,0,2)*rho_B0(0,0,2)*AVG_2_2_LEV_20 + u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_24 +
            u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_21 + u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_22 +
            u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_23;

       CF_35 = p_B0(0,0,2)*AVG_2_2_LEV_33 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_34 +
            u2_B0(0,0,2)*rho_B0(0,0,2)*AVG_2_2_LEV_30 + u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_34 +
            u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_31 + u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_32 +
            u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_33;

       CF_45 = p_B0(0,0,2)*AVG_2_2_LEV_43 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_44 +
            u2_B0(0,0,2)*rho_B0(0,0,2)*AVG_2_2_LEV_40 + u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_44 +
            u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_41 + u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_42 +
            u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_43;

      CS_05 = rho_B0(0,0,2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,2)*AVG_2_2_LEV_02;

      CS_15 = rho_B0(0,0,2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,2)*AVG_2_2_LEV_11;

       CS_25 = rho_B0(0,0,2)*AVG_2_2_LEV_20 + rhoE_B0(0,0,2)*AVG_2_2_LEV_24 + rhou0_B0(0,0,2)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,2)*AVG_2_2_LEV_22 + rhou2_B0(0,0,2)*AVG_2_2_LEV_23;

       CS_35 = rho_B0(0,0,2)*AVG_2_2_LEV_30 + rhoE_B0(0,0,2)*AVG_2_2_LEV_34 + rhou0_B0(0,0,2)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,2)*AVG_2_2_LEV_32 + rhou2_B0(0,0,2)*AVG_2_2_LEV_33;

       CS_45 = rho_B0(0,0,2)*AVG_2_2_LEV_40 + rhoE_B0(0,0,2)*AVG_2_2_LEV_44 + rhou0_B0(0,0,2)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,2)*AVG_2_2_LEV_42 + rhou2_B0(0,0,2)*AVG_2_2_LEV_43;

      CF_06 = (rho_B0(0,0,3)*AVG_2_2_LEV_00 + rhou1_B0(0,0,3)*AVG_2_2_LEV_02)*u2_B0(0,0,3);

      CF_16 = (rho_B0(0,0,3)*AVG_2_2_LEV_10 + rhou0_B0(0,0,3)*AVG_2_2_LEV_11)*u2_B0(0,0,3);

       CF_26 = p_B0(0,0,3)*AVG_2_2_LEV_23 + p_B0(0,0,3)*u2_B0(0,0,3)*AVG_2_2_LEV_24 +
            u2_B0(0,0,3)*rho_B0(0,0,3)*AVG_2_2_LEV_20 + u2_B0(0,0,3)*rhoE_B0(0,0,3)*AVG_2_2_LEV_24 +
            u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_21 + u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_22 +
            u2_B0(0,0,3)*rhou2_B0(0,0,3)*AVG_2_2_LEV_23;

       CF_36 = p_B0(0,0,3)*AVG_2_2_LEV_33 + p_B0(0,0,3)*u2_B0(0,0,3)*AVG_2_2_LEV_34 +
            u2_B0(0,0,3)*rho_B0(0,0,3)*AVG_2_2_LEV_30 + u2_B0(0,0,3)*rhoE_B0(0,0,3)*AVG_2_2_LEV_34 +
            u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_31 + u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_32 +
            u2_B0(0,0,3)*rhou2_B0(0,0,3)*AVG_2_2_LEV_33;

       CF_46 = p_B0(0,0,3)*AVG_2_2_LEV_43 + p_B0(0,0,3)*u2_B0(0,0,3)*AVG_2_2_LEV_44 +
            u2_B0(0,0,3)*rho_B0(0,0,3)*AVG_2_2_LEV_40 + u2_B0(0,0,3)*rhoE_B0(0,0,3)*AVG_2_2_LEV_44 +
            u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_41 + u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_42 +
            u2_B0(0,0,3)*rhou2_B0(0,0,3)*AVG_2_2_LEV_43;

      CS_06 = rho_B0(0,0,3)*AVG_2_2_LEV_00 + rhou1_B0(0,0,3)*AVG_2_2_LEV_02;

      CS_16 = rho_B0(0,0,3)*AVG_2_2_LEV_10 + rhou0_B0(0,0,3)*AVG_2_2_LEV_11;

       CS_26 = rho_B0(0,0,3)*AVG_2_2_LEV_20 + rhoE_B0(0,0,3)*AVG_2_2_LEV_24 + rhou0_B0(0,0,3)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,3)*AVG_2_2_LEV_22 + rhou2_B0(0,0,3)*AVG_2_2_LEV_23;

       CS_36 = rho_B0(0,0,3)*AVG_2_2_LEV_30 + rhoE_B0(0,0,3)*AVG_2_2_LEV_34 + rhou0_B0(0,0,3)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,3)*AVG_2_2_LEV_32 + rhou2_B0(0,0,3)*AVG_2_2_LEV_33;

       CS_46 = rho_B0(0,0,3)*AVG_2_2_LEV_40 + rhoE_B0(0,0,3)*AVG_2_2_LEV_44 + rhou0_B0(0,0,3)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,3)*AVG_2_2_LEV_42 + rhou2_B0(0,0,3)*AVG_2_2_LEV_43;

      CF_07 = (rho_B0(0,0,4)*AVG_2_2_LEV_00 + rhou1_B0(0,0,4)*AVG_2_2_LEV_02)*u2_B0(0,0,4);

      CF_17 = (rho_B0(0,0,4)*AVG_2_2_LEV_10 + rhou0_B0(0,0,4)*AVG_2_2_LEV_11)*u2_B0(0,0,4);

       CF_27 = p_B0(0,0,4)*AVG_2_2_LEV_23 + p_B0(0,0,4)*u2_B0(0,0,4)*AVG_2_2_LEV_24 +
            u2_B0(0,0,4)*rho_B0(0,0,4)*AVG_2_2_LEV_20 + u2_B0(0,0,4)*rhoE_B0(0,0,4)*AVG_2_2_LEV_24 +
            u2_B0(0,0,4)*rhou0_B0(0,0,4)*AVG_2_2_LEV_21 + u2_B0(0,0,4)*rhou1_B0(0,0,4)*AVG_2_2_LEV_22 +
            u2_B0(0,0,4)*rhou2_B0(0,0,4)*AVG_2_2_LEV_23;

       CF_37 = p_B0(0,0,4)*AVG_2_2_LEV_33 + p_B0(0,0,4)*u2_B0(0,0,4)*AVG_2_2_LEV_34 +
            u2_B0(0,0,4)*rho_B0(0,0,4)*AVG_2_2_LEV_30 + u2_B0(0,0,4)*rhoE_B0(0,0,4)*AVG_2_2_LEV_34 +
            u2_B0(0,0,4)*rhou0_B0(0,0,4)*AVG_2_2_LEV_31 + u2_B0(0,0,4)*rhou1_B0(0,0,4)*AVG_2_2_LEV_32 +
            u2_B0(0,0,4)*rhou2_B0(0,0,4)*AVG_2_2_LEV_33;

       CF_47 = p_B0(0,0,4)*AVG_2_2_LEV_43 + p_B0(0,0,4)*u2_B0(0,0,4)*AVG_2_2_LEV_44 +
            u2_B0(0,0,4)*rho_B0(0,0,4)*AVG_2_2_LEV_40 + u2_B0(0,0,4)*rhoE_B0(0,0,4)*AVG_2_2_LEV_44 +
            u2_B0(0,0,4)*rhou0_B0(0,0,4)*AVG_2_2_LEV_41 + u2_B0(0,0,4)*rhou1_B0(0,0,4)*AVG_2_2_LEV_42 +
            u2_B0(0,0,4)*rhou2_B0(0,0,4)*AVG_2_2_LEV_43;

      CS_07 = rho_B0(0,0,4)*AVG_2_2_LEV_00 + rhou1_B0(0,0,4)*AVG_2_2_LEV_02;

      CS_17 = rho_B0(0,0,4)*AVG_2_2_LEV_10 + rhou0_B0(0,0,4)*AVG_2_2_LEV_11;

       CS_27 = rho_B0(0,0,4)*AVG_2_2_LEV_20 + rhoE_B0(0,0,4)*AVG_2_2_LEV_24 + rhou0_B0(0,0,4)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,4)*AVG_2_2_LEV_22 + rhou2_B0(0,0,4)*AVG_2_2_LEV_23;

       CS_37 = rho_B0(0,0,4)*AVG_2_2_LEV_30 + rhoE_B0(0,0,4)*AVG_2_2_LEV_34 + rhou0_B0(0,0,4)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,4)*AVG_2_2_LEV_32 + rhou2_B0(0,0,4)*AVG_2_2_LEV_33;

       CS_47 = rho_B0(0,0,4)*AVG_2_2_LEV_40 + rhoE_B0(0,0,4)*AVG_2_2_LEV_44 + rhou0_B0(0,0,4)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,4)*AVG_2_2_LEV_42 + rhou2_B0(0,0,4)*AVG_2_2_LEV_43;

      max_lambda_00 = shock_filter_control*fmax(fabs(u2_B0(0,0,0)), fabs(u2_B0(0,0,1)));

      max_lambda_11 = max_lambda_00;

      max_lambda_22 = max_lambda_00;

      max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(0,0,0) + u2_B0(0,0,0)), fabs(a_B0(0,0,1) + u2_B0(0,0,1)));

      max_lambda_44 = shock_filter_control*fmax(fabs(-u2_B0(0,0,1) + a_B0(0,0,1)), fabs(-u2_B0(0,0,0) + a_B0(0,0,0)));

       beta_0 = ((547.0/960.0))*((CS_06*max_lambda_00 + CF_06)*(CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_06*max_lambda_00 + CF_06) + ((7043.0/480.0))*(CS_05*max_lambda_00 +
            CF_05))*(CS_05*max_lambda_00 + CF_05) + ((1.0/2.0))*(CS_03*max_lambda_00 +
            CF_03)*(-(1567.0/80.0)*(CS_04*max_lambda_00 + CF_04) - (309.0/80.0)*(CS_06*max_lambda_00 + CF_06) +
            ((2107.0/480.0))*(CS_03*max_lambda_00 + CF_03) + ((3521.0/240.0))*(CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04)*(-(8623.0/240.0)*(CS_05*max_lambda_00 + CF_05) +
            ((2321.0/240.0))*(CS_06*max_lambda_00 + CF_06) + ((11003.0/480.0))*(CS_04*max_lambda_00 + CF_04));

       beta_1 = ((89.0/320.0))*((CS_05*max_lambda_00 + CF_05)*(CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_05*max_lambda_00 + CF_05) + ((2843.0/480.0))*(CS_04*max_lambda_00 +
            CF_04))*(CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(1261.0/240.0)*(CS_03*max_lambda_00 + CF_03) - (247.0/240.0)*(CS_05*max_lambda_00 + CF_05) +
            ((547.0/480.0))*(CS_02*max_lambda_00 + CF_02) + ((961.0/240.0))*(CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(CS_03*max_lambda_00 + CF_03)*(-(2983.0/240.0)*(CS_04*max_lambda_00 + CF_04) +
            ((267.0/80.0))*(CS_05*max_lambda_00 + CF_05) + ((3443.0/480.0))*(CS_03*max_lambda_00 + CF_03));

       beta_2 = ((547.0/960.0))*((CS_04*max_lambda_00 + CF_04)*(CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_04*max_lambda_00 + CF_04) + ((3443.0/480.0))*(CS_03*max_lambda_00 +
            CF_03))*(CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(CS_01*max_lambda_00 +
            CF_01)*(-(247.0/240.0)*(CS_04*max_lambda_00 + CF_04) + ((89.0/160.0))*(CS_01*max_lambda_00 + CF_01) +
            ((267.0/80.0))*(CS_03*max_lambda_00 + CF_03)) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(2983.0/240.0)*(CS_03*max_lambda_00 + CF_03) - (821.0/240.0)*(CS_01*max_lambda_00 + CF_01) +
            ((961.0/240.0))*(CS_04*max_lambda_00 + CF_04) + ((2843.0/480.0))*(CS_02*max_lambda_00 + CF_02));

       beta_3 = ((2107.0/960.0))*((CS_03*max_lambda_00 + CF_03)*(CS_03*max_lambda_00 + CF_03)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_03*max_lambda_00 + CF_03) + ((547.0/480.0))*(CS_00*max_lambda_00 +
            CF_00))*(CS_00*max_lambda_00 + CF_00) + ((1.0/2.0))*(CS_01*max_lambda_00 +
            CF_01)*(-(647.0/80.0)*(CS_00*max_lambda_00 + CF_00) + ((3521.0/240.0))*(CS_03*max_lambda_00 + CF_03) +
            ((7043.0/480.0))*(CS_01*max_lambda_00 + CF_01)) + ((1.0/2.0))*(CS_02*max_lambda_00 +
            CF_02)*(-(8623.0/240.0)*(CS_01*max_lambda_00 + CF_01) - (1567.0/80.0)*(CS_03*max_lambda_00 + CF_03) +
            ((2321.0/240.0))*(CS_00*max_lambda_00 + CF_00) + ((11003.0/480.0))*(CS_02*max_lambda_00 + CF_02));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj0 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_0 = (-(23.0/24.0)*(CS_02*max_lambda_00 + CF_02) - (1.0/8.0)*(CS_00*max_lambda_00 + CF_00) +
            ((13.0/24.0))*(CS_01*max_lambda_00 + CF_01) + ((25.0/24.0))*(CS_03*max_lambda_00 + CF_03))*omega_3 +
            (-(5.0/24.0)*(CS_02*max_lambda_00 + CF_02) + ((1.0/8.0))*(CS_04*max_lambda_00 + CF_04) +
            ((1.0/24.0))*(CS_01*max_lambda_00 + CF_01) + ((13.0/24.0))*(CS_03*max_lambda_00 + CF_03))*omega_2 +
            (-(5.0/24.0)*(CS_05*max_lambda_00 + CF_05) + ((1.0/8.0))*(CS_03*max_lambda_00 + CF_03) +
            ((1.0/24.0))*(CS_06*max_lambda_00 + CF_06) + ((13.0/24.0))*(CS_04*max_lambda_00 + CF_04))*omega_0 +
            (-(1.0/24.0)*(CS_02*max_lambda_00 + CF_02) - (1.0/24.0)*(CS_05*max_lambda_00 + CF_05) +
            ((7.0/24.0))*(CS_03*max_lambda_00 + CF_03) + ((7.0/24.0))*(CS_04*max_lambda_00 + CF_04))*omega_1 + Recon_0;

       beta_0 = ((547.0/960.0))*((-CS_07*max_lambda_00 + CF_07)*(-CS_07*max_lambda_00 + CF_07)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_07*max_lambda_00 + CF_07) + ((7043.0/480.0))*(-CS_06*max_lambda_00 +
            CF_06))*(-CS_06*max_lambda_00 + CF_06) + ((1.0/2.0))*(-CS_04*max_lambda_00 +
            CF_04)*(-(1567.0/80.0)*(-CS_05*max_lambda_00 + CF_05) - (309.0/80.0)*(-CS_07*max_lambda_00 + CF_07) +
            ((2107.0/480.0))*(-CS_04*max_lambda_00 + CF_04) + ((3521.0/240.0))*(-CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05)*(-(8623.0/240.0)*(-CS_06*max_lambda_00 + CF_06) +
            ((2321.0/240.0))*(-CS_07*max_lambda_00 + CF_07) + ((11003.0/480.0))*(-CS_05*max_lambda_00 + CF_05));

       beta_1 = ((89.0/320.0))*((-CS_06*max_lambda_00 + CF_06)*(-CS_06*max_lambda_00 + CF_06)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_06*max_lambda_00 + CF_06) + ((2843.0/480.0))*(-CS_05*max_lambda_00 +
            CF_05))*(-CS_05*max_lambda_00 + CF_05) + ((1.0/2.0))*(-CS_03*max_lambda_00 +
            CF_03)*(-(1261.0/240.0)*(-CS_04*max_lambda_00 + CF_04) - (247.0/240.0)*(-CS_06*max_lambda_00 + CF_06) +
            ((547.0/480.0))*(-CS_03*max_lambda_00 + CF_03) + ((961.0/240.0))*(-CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-CS_04*max_lambda_00 + CF_04)*(-(2983.0/240.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((267.0/80.0))*(-CS_06*max_lambda_00 + CF_06) + ((3443.0/480.0))*(-CS_04*max_lambda_00 + CF_04));

       beta_2 = ((547.0/960.0))*((-CS_05*max_lambda_00 + CF_05)*(-CS_05*max_lambda_00 + CF_05)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_05*max_lambda_00 + CF_05) + ((3443.0/480.0))*(-CS_04*max_lambda_00 +
            CF_04))*(-CS_04*max_lambda_00 + CF_04) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
            CF_02)*(-(821.0/240.0)*(-CS_03*max_lambda_00 + CF_03) - (247.0/240.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((89.0/160.0))*(-CS_02*max_lambda_00 + CF_02) + ((267.0/80.0))*(-CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03)*(-(2983.0/240.0)*(-CS_04*max_lambda_00 + CF_04) +
            ((961.0/240.0))*(-CS_05*max_lambda_00 + CF_05) + ((2843.0/480.0))*(-CS_03*max_lambda_00 + CF_03));

       beta_3 = ((2107.0/960.0))*((-CS_04*max_lambda_00 + CF_04)*(-CS_04*max_lambda_00 + CF_04)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_04*max_lambda_00 + CF_04) + ((11003.0/480.0))*(-CS_03*max_lambda_00 +
            CF_03))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_01*max_lambda_00 +
            CF_01)*(-(309.0/80.0)*(-CS_04*max_lambda_00 + CF_04) + ((547.0/480.0))*(-CS_01*max_lambda_00 + CF_01) +
            ((2321.0/240.0))*(-CS_03*max_lambda_00 + CF_03)) + ((1.0/2.0))*(-CS_02*max_lambda_00 +
            CF_02)*(-(8623.0/240.0)*(-CS_03*max_lambda_00 + CF_03) - (647.0/80.0)*(-CS_01*max_lambda_00 + CF_01) +
            ((3521.0/240.0))*(-CS_04*max_lambda_00 + CF_04) + ((7043.0/480.0))*(-CS_02*max_lambda_00 + CF_02));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj0 = fmax(rj0, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_0 = (-(23.0/24.0)*(-CS_05*max_lambda_00 + CF_05) - (1.0/8.0)*(-CS_07*max_lambda_00 + CF_07) +
            ((13.0/24.0))*(-CS_06*max_lambda_00 + CF_06) + ((25.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_0 +
            (-(5.0/24.0)*(-CS_02*max_lambda_00 + CF_02) + ((1.0/8.0))*(-CS_04*max_lambda_00 + CF_04) +
            ((1.0/24.0))*(-CS_01*max_lambda_00 + CF_01) + ((13.0/24.0))*(-CS_03*max_lambda_00 + CF_03))*omega_3 +
            (-(5.0/24.0)*(-CS_05*max_lambda_00 + CF_05) + ((1.0/8.0))*(-CS_03*max_lambda_00 + CF_03) +
            ((1.0/24.0))*(-CS_06*max_lambda_00 + CF_06) + ((13.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_1 +
            (-(1.0/24.0)*(-CS_02*max_lambda_00 + CF_02) - (1.0/24.0)*(-CS_05*max_lambda_00 + CF_05) +
            ((7.0/24.0))*(-CS_03*max_lambda_00 + CF_03) + ((7.0/24.0))*(-CS_04*max_lambda_00 + CF_04))*omega_2 +
            Recon_0;

       beta_0 = ((547.0/960.0))*((CS_16*max_lambda_11 + CF_16)*(CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_16*max_lambda_11 + CF_16) + ((7043.0/480.0))*(CS_15*max_lambda_11 +
            CF_15))*(CS_15*max_lambda_11 + CF_15) + ((1.0/2.0))*(CS_13*max_lambda_11 +
            CF_13)*(-(1567.0/80.0)*(CS_14*max_lambda_11 + CF_14) - (309.0/80.0)*(CS_16*max_lambda_11 + CF_16) +
            ((2107.0/480.0))*(CS_13*max_lambda_11 + CF_13) + ((3521.0/240.0))*(CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14)*(-(8623.0/240.0)*(CS_15*max_lambda_11 + CF_15) +
            ((2321.0/240.0))*(CS_16*max_lambda_11 + CF_16) + ((11003.0/480.0))*(CS_14*max_lambda_11 + CF_14));

       beta_1 = ((89.0/320.0))*((CS_15*max_lambda_11 + CF_15)*(CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_15*max_lambda_11 + CF_15) + ((2843.0/480.0))*(CS_14*max_lambda_11 +
            CF_14))*(CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(1261.0/240.0)*(CS_13*max_lambda_11 + CF_13) - (247.0/240.0)*(CS_15*max_lambda_11 + CF_15) +
            ((547.0/480.0))*(CS_12*max_lambda_11 + CF_12) + ((961.0/240.0))*(CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(CS_13*max_lambda_11 + CF_13)*(-(2983.0/240.0)*(CS_14*max_lambda_11 + CF_14) +
            ((267.0/80.0))*(CS_15*max_lambda_11 + CF_15) + ((3443.0/480.0))*(CS_13*max_lambda_11 + CF_13));

       beta_2 = ((547.0/960.0))*((CS_14*max_lambda_11 + CF_14)*(CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_14*max_lambda_11 + CF_14) + ((3443.0/480.0))*(CS_13*max_lambda_11 +
            CF_13))*(CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(CS_11*max_lambda_11 +
            CF_11)*(-(247.0/240.0)*(CS_14*max_lambda_11 + CF_14) + ((89.0/160.0))*(CS_11*max_lambda_11 + CF_11) +
            ((267.0/80.0))*(CS_13*max_lambda_11 + CF_13)) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(2983.0/240.0)*(CS_13*max_lambda_11 + CF_13) - (821.0/240.0)*(CS_11*max_lambda_11 + CF_11) +
            ((961.0/240.0))*(CS_14*max_lambda_11 + CF_14) + ((2843.0/480.0))*(CS_12*max_lambda_11 + CF_12));

       beta_3 = ((2107.0/960.0))*((CS_13*max_lambda_11 + CF_13)*(CS_13*max_lambda_11 + CF_13)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_13*max_lambda_11 + CF_13) + ((547.0/480.0))*(CS_10*max_lambda_11 +
            CF_10))*(CS_10*max_lambda_11 + CF_10) + ((1.0/2.0))*(CS_11*max_lambda_11 +
            CF_11)*(-(647.0/80.0)*(CS_10*max_lambda_11 + CF_10) + ((3521.0/240.0))*(CS_13*max_lambda_11 + CF_13) +
            ((7043.0/480.0))*(CS_11*max_lambda_11 + CF_11)) + ((1.0/2.0))*(CS_12*max_lambda_11 +
            CF_12)*(-(8623.0/240.0)*(CS_11*max_lambda_11 + CF_11) - (1567.0/80.0)*(CS_13*max_lambda_11 + CF_13) +
            ((2321.0/240.0))*(CS_10*max_lambda_11 + CF_10) + ((11003.0/480.0))*(CS_12*max_lambda_11 + CF_12));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj1 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_1 = (-(23.0/24.0)*(CS_12*max_lambda_11 + CF_12) - (1.0/8.0)*(CS_10*max_lambda_11 + CF_10) +
            ((13.0/24.0))*(CS_11*max_lambda_11 + CF_11) + ((25.0/24.0))*(CS_13*max_lambda_11 + CF_13))*omega_3 +
            (-(5.0/24.0)*(CS_12*max_lambda_11 + CF_12) + ((1.0/8.0))*(CS_14*max_lambda_11 + CF_14) +
            ((1.0/24.0))*(CS_11*max_lambda_11 + CF_11) + ((13.0/24.0))*(CS_13*max_lambda_11 + CF_13))*omega_2 +
            (-(5.0/24.0)*(CS_15*max_lambda_11 + CF_15) + ((1.0/8.0))*(CS_13*max_lambda_11 + CF_13) +
            ((1.0/24.0))*(CS_16*max_lambda_11 + CF_16) + ((13.0/24.0))*(CS_14*max_lambda_11 + CF_14))*omega_0 +
            (-(1.0/24.0)*(CS_12*max_lambda_11 + CF_12) - (1.0/24.0)*(CS_15*max_lambda_11 + CF_15) +
            ((7.0/24.0))*(CS_13*max_lambda_11 + CF_13) + ((7.0/24.0))*(CS_14*max_lambda_11 + CF_14))*omega_1 + Recon_1;

       beta_0 = ((547.0/960.0))*((-CS_17*max_lambda_11 + CF_17)*(-CS_17*max_lambda_11 + CF_17)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_17*max_lambda_11 + CF_17) + ((7043.0/480.0))*(-CS_16*max_lambda_11 +
            CF_16))*(-CS_16*max_lambda_11 + CF_16) + ((1.0/2.0))*(-CS_14*max_lambda_11 +
            CF_14)*(-(1567.0/80.0)*(-CS_15*max_lambda_11 + CF_15) - (309.0/80.0)*(-CS_17*max_lambda_11 + CF_17) +
            ((2107.0/480.0))*(-CS_14*max_lambda_11 + CF_14) + ((3521.0/240.0))*(-CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15)*(-(8623.0/240.0)*(-CS_16*max_lambda_11 + CF_16) +
            ((2321.0/240.0))*(-CS_17*max_lambda_11 + CF_17) + ((11003.0/480.0))*(-CS_15*max_lambda_11 + CF_15));

       beta_1 = ((89.0/320.0))*((-CS_16*max_lambda_11 + CF_16)*(-CS_16*max_lambda_11 + CF_16)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_16*max_lambda_11 + CF_16) + ((2843.0/480.0))*(-CS_15*max_lambda_11 +
            CF_15))*(-CS_15*max_lambda_11 + CF_15) + ((1.0/2.0))*(-CS_13*max_lambda_11 +
            CF_13)*(-(1261.0/240.0)*(-CS_14*max_lambda_11 + CF_14) - (247.0/240.0)*(-CS_16*max_lambda_11 + CF_16) +
            ((547.0/480.0))*(-CS_13*max_lambda_11 + CF_13) + ((961.0/240.0))*(-CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-CS_14*max_lambda_11 + CF_14)*(-(2983.0/240.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((267.0/80.0))*(-CS_16*max_lambda_11 + CF_16) + ((3443.0/480.0))*(-CS_14*max_lambda_11 + CF_14));

       beta_2 = ((547.0/960.0))*((-CS_15*max_lambda_11 + CF_15)*(-CS_15*max_lambda_11 + CF_15)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_15*max_lambda_11 + CF_15) + ((3443.0/480.0))*(-CS_14*max_lambda_11 +
            CF_14))*(-CS_14*max_lambda_11 + CF_14) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
            CF_12)*(-(821.0/240.0)*(-CS_13*max_lambda_11 + CF_13) - (247.0/240.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((89.0/160.0))*(-CS_12*max_lambda_11 + CF_12) + ((267.0/80.0))*(-CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13)*(-(2983.0/240.0)*(-CS_14*max_lambda_11 + CF_14) +
            ((961.0/240.0))*(-CS_15*max_lambda_11 + CF_15) + ((2843.0/480.0))*(-CS_13*max_lambda_11 + CF_13));

       beta_3 = ((2107.0/960.0))*((-CS_14*max_lambda_11 + CF_14)*(-CS_14*max_lambda_11 + CF_14)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_14*max_lambda_11 + CF_14) + ((11003.0/480.0))*(-CS_13*max_lambda_11 +
            CF_13))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_11*max_lambda_11 +
            CF_11)*(-(309.0/80.0)*(-CS_14*max_lambda_11 + CF_14) + ((547.0/480.0))*(-CS_11*max_lambda_11 + CF_11) +
            ((2321.0/240.0))*(-CS_13*max_lambda_11 + CF_13)) + ((1.0/2.0))*(-CS_12*max_lambda_11 +
            CF_12)*(-(8623.0/240.0)*(-CS_13*max_lambda_11 + CF_13) - (647.0/80.0)*(-CS_11*max_lambda_11 + CF_11) +
            ((3521.0/240.0))*(-CS_14*max_lambda_11 + CF_14) + ((7043.0/480.0))*(-CS_12*max_lambda_11 + CF_12));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj1 = fmax(rj1, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_1 = (-(23.0/24.0)*(-CS_15*max_lambda_11 + CF_15) - (1.0/8.0)*(-CS_17*max_lambda_11 + CF_17) +
            ((13.0/24.0))*(-CS_16*max_lambda_11 + CF_16) + ((25.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_0 +
            (-(5.0/24.0)*(-CS_12*max_lambda_11 + CF_12) + ((1.0/8.0))*(-CS_14*max_lambda_11 + CF_14) +
            ((1.0/24.0))*(-CS_11*max_lambda_11 + CF_11) + ((13.0/24.0))*(-CS_13*max_lambda_11 + CF_13))*omega_3 +
            (-(5.0/24.0)*(-CS_15*max_lambda_11 + CF_15) + ((1.0/8.0))*(-CS_13*max_lambda_11 + CF_13) +
            ((1.0/24.0))*(-CS_16*max_lambda_11 + CF_16) + ((13.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_1 +
            (-(1.0/24.0)*(-CS_12*max_lambda_11 + CF_12) - (1.0/24.0)*(-CS_15*max_lambda_11 + CF_15) +
            ((7.0/24.0))*(-CS_13*max_lambda_11 + CF_13) + ((7.0/24.0))*(-CS_14*max_lambda_11 + CF_14))*omega_2 +
            Recon_1;

       beta_0 = ((547.0/960.0))*((CS_26*max_lambda_22 + CF_26)*(CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_26*max_lambda_22 + CF_26) + ((7043.0/480.0))*(CS_25*max_lambda_22 +
            CF_25))*(CS_25*max_lambda_22 + CF_25) + ((1.0/2.0))*(CS_23*max_lambda_22 +
            CF_23)*(-(1567.0/80.0)*(CS_24*max_lambda_22 + CF_24) - (309.0/80.0)*(CS_26*max_lambda_22 + CF_26) +
            ((2107.0/480.0))*(CS_23*max_lambda_22 + CF_23) + ((3521.0/240.0))*(CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24)*(-(8623.0/240.0)*(CS_25*max_lambda_22 + CF_25) +
            ((2321.0/240.0))*(CS_26*max_lambda_22 + CF_26) + ((11003.0/480.0))*(CS_24*max_lambda_22 + CF_24));

       beta_1 = ((89.0/320.0))*((CS_25*max_lambda_22 + CF_25)*(CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_25*max_lambda_22 + CF_25) + ((2843.0/480.0))*(CS_24*max_lambda_22 +
            CF_24))*(CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(1261.0/240.0)*(CS_23*max_lambda_22 + CF_23) - (247.0/240.0)*(CS_25*max_lambda_22 + CF_25) +
            ((547.0/480.0))*(CS_22*max_lambda_22 + CF_22) + ((961.0/240.0))*(CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(CS_23*max_lambda_22 + CF_23)*(-(2983.0/240.0)*(CS_24*max_lambda_22 + CF_24) +
            ((267.0/80.0))*(CS_25*max_lambda_22 + CF_25) + ((3443.0/480.0))*(CS_23*max_lambda_22 + CF_23));

       beta_2 = ((547.0/960.0))*((CS_24*max_lambda_22 + CF_24)*(CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_24*max_lambda_22 + CF_24) + ((3443.0/480.0))*(CS_23*max_lambda_22 +
            CF_23))*(CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(CS_21*max_lambda_22 +
            CF_21)*(-(247.0/240.0)*(CS_24*max_lambda_22 + CF_24) + ((89.0/160.0))*(CS_21*max_lambda_22 + CF_21) +
            ((267.0/80.0))*(CS_23*max_lambda_22 + CF_23)) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(2983.0/240.0)*(CS_23*max_lambda_22 + CF_23) - (821.0/240.0)*(CS_21*max_lambda_22 + CF_21) +
            ((961.0/240.0))*(CS_24*max_lambda_22 + CF_24) + ((2843.0/480.0))*(CS_22*max_lambda_22 + CF_22));

       beta_3 = ((2107.0/960.0))*((CS_23*max_lambda_22 + CF_23)*(CS_23*max_lambda_22 + CF_23)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_23*max_lambda_22 + CF_23) + ((547.0/480.0))*(CS_20*max_lambda_22 +
            CF_20))*(CS_20*max_lambda_22 + CF_20) + ((1.0/2.0))*(CS_21*max_lambda_22 +
            CF_21)*(-(647.0/80.0)*(CS_20*max_lambda_22 + CF_20) + ((3521.0/240.0))*(CS_23*max_lambda_22 + CF_23) +
            ((7043.0/480.0))*(CS_21*max_lambda_22 + CF_21)) + ((1.0/2.0))*(CS_22*max_lambda_22 +
            CF_22)*(-(8623.0/240.0)*(CS_21*max_lambda_22 + CF_21) - (1567.0/80.0)*(CS_23*max_lambda_22 + CF_23) +
            ((2321.0/240.0))*(CS_20*max_lambda_22 + CF_20) + ((11003.0/480.0))*(CS_22*max_lambda_22 + CF_22));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj2 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_2 = (-(23.0/24.0)*(CS_22*max_lambda_22 + CF_22) - (1.0/8.0)*(CS_20*max_lambda_22 + CF_20) +
            ((13.0/24.0))*(CS_21*max_lambda_22 + CF_21) + ((25.0/24.0))*(CS_23*max_lambda_22 + CF_23))*omega_3 +
            (-(5.0/24.0)*(CS_22*max_lambda_22 + CF_22) + ((1.0/8.0))*(CS_24*max_lambda_22 + CF_24) +
            ((1.0/24.0))*(CS_21*max_lambda_22 + CF_21) + ((13.0/24.0))*(CS_23*max_lambda_22 + CF_23))*omega_2 +
            (-(5.0/24.0)*(CS_25*max_lambda_22 + CF_25) + ((1.0/8.0))*(CS_23*max_lambda_22 + CF_23) +
            ((1.0/24.0))*(CS_26*max_lambda_22 + CF_26) + ((13.0/24.0))*(CS_24*max_lambda_22 + CF_24))*omega_0 +
            (-(1.0/24.0)*(CS_22*max_lambda_22 + CF_22) - (1.0/24.0)*(CS_25*max_lambda_22 + CF_25) +
            ((7.0/24.0))*(CS_23*max_lambda_22 + CF_23) + ((7.0/24.0))*(CS_24*max_lambda_22 + CF_24))*omega_1 + Recon_2;

       beta_0 = ((547.0/960.0))*((-CS_27*max_lambda_22 + CF_27)*(-CS_27*max_lambda_22 + CF_27)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_27*max_lambda_22 + CF_27) + ((7043.0/480.0))*(-CS_26*max_lambda_22 +
            CF_26))*(-CS_26*max_lambda_22 + CF_26) + ((1.0/2.0))*(-CS_24*max_lambda_22 +
            CF_24)*(-(1567.0/80.0)*(-CS_25*max_lambda_22 + CF_25) - (309.0/80.0)*(-CS_27*max_lambda_22 + CF_27) +
            ((2107.0/480.0))*(-CS_24*max_lambda_22 + CF_24) + ((3521.0/240.0))*(-CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25)*(-(8623.0/240.0)*(-CS_26*max_lambda_22 + CF_26) +
            ((2321.0/240.0))*(-CS_27*max_lambda_22 + CF_27) + ((11003.0/480.0))*(-CS_25*max_lambda_22 + CF_25));

       beta_1 = ((89.0/320.0))*((-CS_26*max_lambda_22 + CF_26)*(-CS_26*max_lambda_22 + CF_26)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_26*max_lambda_22 + CF_26) + ((2843.0/480.0))*(-CS_25*max_lambda_22 +
            CF_25))*(-CS_25*max_lambda_22 + CF_25) + ((1.0/2.0))*(-CS_23*max_lambda_22 +
            CF_23)*(-(1261.0/240.0)*(-CS_24*max_lambda_22 + CF_24) - (247.0/240.0)*(-CS_26*max_lambda_22 + CF_26) +
            ((547.0/480.0))*(-CS_23*max_lambda_22 + CF_23) + ((961.0/240.0))*(-CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-CS_24*max_lambda_22 + CF_24)*(-(2983.0/240.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((267.0/80.0))*(-CS_26*max_lambda_22 + CF_26) + ((3443.0/480.0))*(-CS_24*max_lambda_22 + CF_24));

       beta_2 = ((547.0/960.0))*((-CS_25*max_lambda_22 + CF_25)*(-CS_25*max_lambda_22 + CF_25)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_25*max_lambda_22 + CF_25) + ((3443.0/480.0))*(-CS_24*max_lambda_22 +
            CF_24))*(-CS_24*max_lambda_22 + CF_24) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
            CF_22)*(-(821.0/240.0)*(-CS_23*max_lambda_22 + CF_23) - (247.0/240.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((89.0/160.0))*(-CS_22*max_lambda_22 + CF_22) + ((267.0/80.0))*(-CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23)*(-(2983.0/240.0)*(-CS_24*max_lambda_22 + CF_24) +
            ((961.0/240.0))*(-CS_25*max_lambda_22 + CF_25) + ((2843.0/480.0))*(-CS_23*max_lambda_22 + CF_23));

       beta_3 = ((2107.0/960.0))*((-CS_24*max_lambda_22 + CF_24)*(-CS_24*max_lambda_22 + CF_24)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_24*max_lambda_22 + CF_24) + ((11003.0/480.0))*(-CS_23*max_lambda_22 +
            CF_23))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_21*max_lambda_22 +
            CF_21)*(-(309.0/80.0)*(-CS_24*max_lambda_22 + CF_24) + ((547.0/480.0))*(-CS_21*max_lambda_22 + CF_21) +
            ((2321.0/240.0))*(-CS_23*max_lambda_22 + CF_23)) + ((1.0/2.0))*(-CS_22*max_lambda_22 +
            CF_22)*(-(8623.0/240.0)*(-CS_23*max_lambda_22 + CF_23) - (647.0/80.0)*(-CS_21*max_lambda_22 + CF_21) +
            ((3521.0/240.0))*(-CS_24*max_lambda_22 + CF_24) + ((7043.0/480.0))*(-CS_22*max_lambda_22 + CF_22));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj2 = fmax(rj2, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_2 = (-(23.0/24.0)*(-CS_25*max_lambda_22 + CF_25) - (1.0/8.0)*(-CS_27*max_lambda_22 + CF_27) +
            ((13.0/24.0))*(-CS_26*max_lambda_22 + CF_26) + ((25.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_0 +
            (-(5.0/24.0)*(-CS_22*max_lambda_22 + CF_22) + ((1.0/8.0))*(-CS_24*max_lambda_22 + CF_24) +
            ((1.0/24.0))*(-CS_21*max_lambda_22 + CF_21) + ((13.0/24.0))*(-CS_23*max_lambda_22 + CF_23))*omega_3 +
            (-(5.0/24.0)*(-CS_25*max_lambda_22 + CF_25) + ((1.0/8.0))*(-CS_23*max_lambda_22 + CF_23) +
            ((1.0/24.0))*(-CS_26*max_lambda_22 + CF_26) + ((13.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_1 +
            (-(1.0/24.0)*(-CS_22*max_lambda_22 + CF_22) - (1.0/24.0)*(-CS_25*max_lambda_22 + CF_25) +
            ((7.0/24.0))*(-CS_23*max_lambda_22 + CF_23) + ((7.0/24.0))*(-CS_24*max_lambda_22 + CF_24))*omega_2 +
            Recon_2;

       beta_0 = ((547.0/960.0))*((CS_36*max_lambda_33 + CF_36)*(CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_36*max_lambda_33 + CF_36) + ((7043.0/480.0))*(CS_35*max_lambda_33 +
            CF_35))*(CS_35*max_lambda_33 + CF_35) + ((1.0/2.0))*(CS_33*max_lambda_33 +
            CF_33)*(-(1567.0/80.0)*(CS_34*max_lambda_33 + CF_34) - (309.0/80.0)*(CS_36*max_lambda_33 + CF_36) +
            ((2107.0/480.0))*(CS_33*max_lambda_33 + CF_33) + ((3521.0/240.0))*(CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34)*(-(8623.0/240.0)*(CS_35*max_lambda_33 + CF_35) +
            ((2321.0/240.0))*(CS_36*max_lambda_33 + CF_36) + ((11003.0/480.0))*(CS_34*max_lambda_33 + CF_34));

       beta_1 = ((89.0/320.0))*((CS_35*max_lambda_33 + CF_35)*(CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_35*max_lambda_33 + CF_35) + ((2843.0/480.0))*(CS_34*max_lambda_33 +
            CF_34))*(CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(1261.0/240.0)*(CS_33*max_lambda_33 + CF_33) - (247.0/240.0)*(CS_35*max_lambda_33 + CF_35) +
            ((547.0/480.0))*(CS_32*max_lambda_33 + CF_32) + ((961.0/240.0))*(CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(CS_33*max_lambda_33 + CF_33)*(-(2983.0/240.0)*(CS_34*max_lambda_33 + CF_34) +
            ((267.0/80.0))*(CS_35*max_lambda_33 + CF_35) + ((3443.0/480.0))*(CS_33*max_lambda_33 + CF_33));

       beta_2 = ((547.0/960.0))*((CS_34*max_lambda_33 + CF_34)*(CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_34*max_lambda_33 + CF_34) + ((3443.0/480.0))*(CS_33*max_lambda_33 +
            CF_33))*(CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(CS_31*max_lambda_33 +
            CF_31)*(-(247.0/240.0)*(CS_34*max_lambda_33 + CF_34) + ((89.0/160.0))*(CS_31*max_lambda_33 + CF_31) +
            ((267.0/80.0))*(CS_33*max_lambda_33 + CF_33)) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(2983.0/240.0)*(CS_33*max_lambda_33 + CF_33) - (821.0/240.0)*(CS_31*max_lambda_33 + CF_31) +
            ((961.0/240.0))*(CS_34*max_lambda_33 + CF_34) + ((2843.0/480.0))*(CS_32*max_lambda_33 + CF_32));

       beta_3 = ((2107.0/960.0))*((CS_33*max_lambda_33 + CF_33)*(CS_33*max_lambda_33 + CF_33)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_33*max_lambda_33 + CF_33) + ((547.0/480.0))*(CS_30*max_lambda_33 +
            CF_30))*(CS_30*max_lambda_33 + CF_30) + ((1.0/2.0))*(CS_31*max_lambda_33 +
            CF_31)*(-(647.0/80.0)*(CS_30*max_lambda_33 + CF_30) + ((3521.0/240.0))*(CS_33*max_lambda_33 + CF_33) +
            ((7043.0/480.0))*(CS_31*max_lambda_33 + CF_31)) + ((1.0/2.0))*(CS_32*max_lambda_33 +
            CF_32)*(-(8623.0/240.0)*(CS_31*max_lambda_33 + CF_31) - (1567.0/80.0)*(CS_33*max_lambda_33 + CF_33) +
            ((2321.0/240.0))*(CS_30*max_lambda_33 + CF_30) + ((11003.0/480.0))*(CS_32*max_lambda_33 + CF_32));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj3 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_3 = (-(23.0/24.0)*(CS_32*max_lambda_33 + CF_32) - (1.0/8.0)*(CS_30*max_lambda_33 + CF_30) +
            ((13.0/24.0))*(CS_31*max_lambda_33 + CF_31) + ((25.0/24.0))*(CS_33*max_lambda_33 + CF_33))*omega_3 +
            (-(5.0/24.0)*(CS_32*max_lambda_33 + CF_32) + ((1.0/8.0))*(CS_34*max_lambda_33 + CF_34) +
            ((1.0/24.0))*(CS_31*max_lambda_33 + CF_31) + ((13.0/24.0))*(CS_33*max_lambda_33 + CF_33))*omega_2 +
            (-(5.0/24.0)*(CS_35*max_lambda_33 + CF_35) + ((1.0/8.0))*(CS_33*max_lambda_33 + CF_33) +
            ((1.0/24.0))*(CS_36*max_lambda_33 + CF_36) + ((13.0/24.0))*(CS_34*max_lambda_33 + CF_34))*omega_0 +
            (-(1.0/24.0)*(CS_32*max_lambda_33 + CF_32) - (1.0/24.0)*(CS_35*max_lambda_33 + CF_35) +
            ((7.0/24.0))*(CS_33*max_lambda_33 + CF_33) + ((7.0/24.0))*(CS_34*max_lambda_33 + CF_34))*omega_1 + Recon_3;

       beta_0 = ((547.0/960.0))*((-CS_37*max_lambda_33 + CF_37)*(-CS_37*max_lambda_33 + CF_37)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_37*max_lambda_33 + CF_37) + ((7043.0/480.0))*(-CS_36*max_lambda_33 +
            CF_36))*(-CS_36*max_lambda_33 + CF_36) + ((1.0/2.0))*(-CS_34*max_lambda_33 +
            CF_34)*(-(1567.0/80.0)*(-CS_35*max_lambda_33 + CF_35) - (309.0/80.0)*(-CS_37*max_lambda_33 + CF_37) +
            ((2107.0/480.0))*(-CS_34*max_lambda_33 + CF_34) + ((3521.0/240.0))*(-CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35)*(-(8623.0/240.0)*(-CS_36*max_lambda_33 + CF_36) +
            ((2321.0/240.0))*(-CS_37*max_lambda_33 + CF_37) + ((11003.0/480.0))*(-CS_35*max_lambda_33 + CF_35));

       beta_1 = ((89.0/320.0))*((-CS_36*max_lambda_33 + CF_36)*(-CS_36*max_lambda_33 + CF_36)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_36*max_lambda_33 + CF_36) + ((2843.0/480.0))*(-CS_35*max_lambda_33 +
            CF_35))*(-CS_35*max_lambda_33 + CF_35) + ((1.0/2.0))*(-CS_33*max_lambda_33 +
            CF_33)*(-(1261.0/240.0)*(-CS_34*max_lambda_33 + CF_34) - (247.0/240.0)*(-CS_36*max_lambda_33 + CF_36) +
            ((547.0/480.0))*(-CS_33*max_lambda_33 + CF_33) + ((961.0/240.0))*(-CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-CS_34*max_lambda_33 + CF_34)*(-(2983.0/240.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((267.0/80.0))*(-CS_36*max_lambda_33 + CF_36) + ((3443.0/480.0))*(-CS_34*max_lambda_33 + CF_34));

       beta_2 = ((547.0/960.0))*((-CS_35*max_lambda_33 + CF_35)*(-CS_35*max_lambda_33 + CF_35)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_35*max_lambda_33 + CF_35) + ((3443.0/480.0))*(-CS_34*max_lambda_33 +
            CF_34))*(-CS_34*max_lambda_33 + CF_34) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
            CF_32)*(-(821.0/240.0)*(-CS_33*max_lambda_33 + CF_33) - (247.0/240.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((89.0/160.0))*(-CS_32*max_lambda_33 + CF_32) + ((267.0/80.0))*(-CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33)*(-(2983.0/240.0)*(-CS_34*max_lambda_33 + CF_34) +
            ((961.0/240.0))*(-CS_35*max_lambda_33 + CF_35) + ((2843.0/480.0))*(-CS_33*max_lambda_33 + CF_33));

       beta_3 = ((2107.0/960.0))*((-CS_34*max_lambda_33 + CF_34)*(-CS_34*max_lambda_33 + CF_34)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_34*max_lambda_33 + CF_34) + ((11003.0/480.0))*(-CS_33*max_lambda_33 +
            CF_33))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_31*max_lambda_33 +
            CF_31)*(-(309.0/80.0)*(-CS_34*max_lambda_33 + CF_34) + ((547.0/480.0))*(-CS_31*max_lambda_33 + CF_31) +
            ((2321.0/240.0))*(-CS_33*max_lambda_33 + CF_33)) + ((1.0/2.0))*(-CS_32*max_lambda_33 +
            CF_32)*(-(8623.0/240.0)*(-CS_33*max_lambda_33 + CF_33) - (647.0/80.0)*(-CS_31*max_lambda_33 + CF_31) +
            ((3521.0/240.0))*(-CS_34*max_lambda_33 + CF_34) + ((7043.0/480.0))*(-CS_32*max_lambda_33 + CF_32));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj3 = fmax(rj3, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_3 = (-(23.0/24.0)*(-CS_35*max_lambda_33 + CF_35) - (1.0/8.0)*(-CS_37*max_lambda_33 + CF_37) +
            ((13.0/24.0))*(-CS_36*max_lambda_33 + CF_36) + ((25.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_0 +
            (-(5.0/24.0)*(-CS_32*max_lambda_33 + CF_32) + ((1.0/8.0))*(-CS_34*max_lambda_33 + CF_34) +
            ((1.0/24.0))*(-CS_31*max_lambda_33 + CF_31) + ((13.0/24.0))*(-CS_33*max_lambda_33 + CF_33))*omega_3 +
            (-(5.0/24.0)*(-CS_35*max_lambda_33 + CF_35) + ((1.0/8.0))*(-CS_33*max_lambda_33 + CF_33) +
            ((1.0/24.0))*(-CS_36*max_lambda_33 + CF_36) + ((13.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_1 +
            (-(1.0/24.0)*(-CS_32*max_lambda_33 + CF_32) - (1.0/24.0)*(-CS_35*max_lambda_33 + CF_35) +
            ((7.0/24.0))*(-CS_33*max_lambda_33 + CF_33) + ((7.0/24.0))*(-CS_34*max_lambda_33 + CF_34))*omega_2 +
            Recon_3;

       beta_0 = ((547.0/960.0))*((CS_46*max_lambda_44 + CF_46)*(CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-(647.0/80.0)*(CS_46*max_lambda_44 + CF_46) + ((7043.0/480.0))*(CS_45*max_lambda_44 +
            CF_45))*(CS_45*max_lambda_44 + CF_45) + ((1.0/2.0))*(CS_43*max_lambda_44 +
            CF_43)*(-(1567.0/80.0)*(CS_44*max_lambda_44 + CF_44) - (309.0/80.0)*(CS_46*max_lambda_44 + CF_46) +
            ((2107.0/480.0))*(CS_43*max_lambda_44 + CF_43) + ((3521.0/240.0))*(CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44)*(-(8623.0/240.0)*(CS_45*max_lambda_44 + CF_45) +
            ((2321.0/240.0))*(CS_46*max_lambda_44 + CF_46) + ((11003.0/480.0))*(CS_44*max_lambda_44 + CF_44));

       beta_1 = ((89.0/320.0))*((CS_45*max_lambda_44 + CF_45)*(CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-(821.0/240.0)*(CS_45*max_lambda_44 + CF_45) + ((2843.0/480.0))*(CS_44*max_lambda_44 +
            CF_44))*(CS_44*max_lambda_44 + CF_44) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(1261.0/240.0)*(CS_43*max_lambda_44 + CF_43) - (247.0/240.0)*(CS_45*max_lambda_44 + CF_45) +
            ((547.0/480.0))*(CS_42*max_lambda_44 + CF_42) + ((961.0/240.0))*(CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(CS_43*max_lambda_44 + CF_43)*(-(2983.0/240.0)*(CS_44*max_lambda_44 + CF_44) +
            ((267.0/80.0))*(CS_45*max_lambda_44 + CF_45) + ((3443.0/480.0))*(CS_43*max_lambda_44 + CF_43));

       beta_2 = ((547.0/960.0))*((CS_44*max_lambda_44 + CF_44)*(CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(CS_44*max_lambda_44 + CF_44) + ((3443.0/480.0))*(CS_43*max_lambda_44 +
            CF_43))*(CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(CS_41*max_lambda_44 +
            CF_41)*(-(247.0/240.0)*(CS_44*max_lambda_44 + CF_44) + ((89.0/160.0))*(CS_41*max_lambda_44 + CF_41) +
            ((267.0/80.0))*(CS_43*max_lambda_44 + CF_43)) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(2983.0/240.0)*(CS_43*max_lambda_44 + CF_43) - (821.0/240.0)*(CS_41*max_lambda_44 + CF_41) +
            ((961.0/240.0))*(CS_44*max_lambda_44 + CF_44) + ((2843.0/480.0))*(CS_42*max_lambda_44 + CF_42));

       beta_3 = ((2107.0/960.0))*((CS_43*max_lambda_44 + CF_43)*(CS_43*max_lambda_44 + CF_43)) +
            ((1.0/2.0))*(-(309.0/80.0)*(CS_43*max_lambda_44 + CF_43) + ((547.0/480.0))*(CS_40*max_lambda_44 +
            CF_40))*(CS_40*max_lambda_44 + CF_40) + ((1.0/2.0))*(CS_41*max_lambda_44 +
            CF_41)*(-(647.0/80.0)*(CS_40*max_lambda_44 + CF_40) + ((3521.0/240.0))*(CS_43*max_lambda_44 + CF_43) +
            ((7043.0/480.0))*(CS_41*max_lambda_44 + CF_41)) + ((1.0/2.0))*(CS_42*max_lambda_44 +
            CF_42)*(-(8623.0/240.0)*(CS_41*max_lambda_44 + CF_41) - (1567.0/80.0)*(CS_43*max_lambda_44 + CF_43) +
            ((2321.0/240.0))*(CS_40*max_lambda_44 + CF_40) + ((11003.0/480.0))*(CS_42*max_lambda_44 + CF_42));

       alpha_0 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj4 = 0.027027027027027*fabs(-1.0 + 35*omega_3) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_0) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_2) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_1);

       Recon_4 = (-(23.0/24.0)*(CS_42*max_lambda_44 + CF_42) - (1.0/8.0)*(CS_40*max_lambda_44 + CF_40) +
            ((13.0/24.0))*(CS_41*max_lambda_44 + CF_41) + ((25.0/24.0))*(CS_43*max_lambda_44 + CF_43))*omega_3 +
            (-(5.0/24.0)*(CS_42*max_lambda_44 + CF_42) + ((1.0/8.0))*(CS_44*max_lambda_44 + CF_44) +
            ((1.0/24.0))*(CS_41*max_lambda_44 + CF_41) + ((13.0/24.0))*(CS_43*max_lambda_44 + CF_43))*omega_2 +
            (-(5.0/24.0)*(CS_45*max_lambda_44 + CF_45) + ((1.0/8.0))*(CS_43*max_lambda_44 + CF_43) +
            ((1.0/24.0))*(CS_46*max_lambda_44 + CF_46) + ((13.0/24.0))*(CS_44*max_lambda_44 + CF_44))*omega_0 +
            (-(1.0/24.0)*(CS_42*max_lambda_44 + CF_42) - (1.0/24.0)*(CS_45*max_lambda_44 + CF_45) +
            ((7.0/24.0))*(CS_43*max_lambda_44 + CF_43) + ((7.0/24.0))*(CS_44*max_lambda_44 + CF_44))*omega_1 + Recon_4;

       beta_0 = ((547.0/960.0))*((-CS_47*max_lambda_44 + CF_47)*(-CS_47*max_lambda_44 + CF_47)) +
            ((1.0/2.0))*(-(647.0/80.0)*(-CS_47*max_lambda_44 + CF_47) + ((7043.0/480.0))*(-CS_46*max_lambda_44 +
            CF_46))*(-CS_46*max_lambda_44 + CF_46) + ((1.0/2.0))*(-CS_44*max_lambda_44 +
            CF_44)*(-(1567.0/80.0)*(-CS_45*max_lambda_44 + CF_45) - (309.0/80.0)*(-CS_47*max_lambda_44 + CF_47) +
            ((2107.0/480.0))*(-CS_44*max_lambda_44 + CF_44) + ((3521.0/240.0))*(-CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45)*(-(8623.0/240.0)*(-CS_46*max_lambda_44 + CF_46) +
            ((2321.0/240.0))*(-CS_47*max_lambda_44 + CF_47) + ((11003.0/480.0))*(-CS_45*max_lambda_44 + CF_45));

       beta_1 = ((89.0/320.0))*((-CS_46*max_lambda_44 + CF_46)*(-CS_46*max_lambda_44 + CF_46)) +
            ((1.0/2.0))*(-(821.0/240.0)*(-CS_46*max_lambda_44 + CF_46) + ((2843.0/480.0))*(-CS_45*max_lambda_44 +
            CF_45))*(-CS_45*max_lambda_44 + CF_45) + ((1.0/2.0))*(-CS_43*max_lambda_44 +
            CF_43)*(-(1261.0/240.0)*(-CS_44*max_lambda_44 + CF_44) - (247.0/240.0)*(-CS_46*max_lambda_44 + CF_46) +
            ((547.0/480.0))*(-CS_43*max_lambda_44 + CF_43) + ((961.0/240.0))*(-CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-CS_44*max_lambda_44 + CF_44)*(-(2983.0/240.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((267.0/80.0))*(-CS_46*max_lambda_44 + CF_46) + ((3443.0/480.0))*(-CS_44*max_lambda_44 + CF_44));

       beta_2 = ((547.0/960.0))*((-CS_45*max_lambda_44 + CF_45)*(-CS_45*max_lambda_44 + CF_45)) +
            ((1.0/2.0))*(-(1261.0/240.0)*(-CS_45*max_lambda_44 + CF_45) + ((3443.0/480.0))*(-CS_44*max_lambda_44 +
            CF_44))*(-CS_44*max_lambda_44 + CF_44) + ((1.0/2.0))*(-CS_42*max_lambda_44 +
            CF_42)*(-(821.0/240.0)*(-CS_43*max_lambda_44 + CF_43) - (247.0/240.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((89.0/160.0))*(-CS_42*max_lambda_44 + CF_42) + ((267.0/80.0))*(-CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43)*(-(2983.0/240.0)*(-CS_44*max_lambda_44 + CF_44) +
            ((961.0/240.0))*(-CS_45*max_lambda_44 + CF_45) + ((2843.0/480.0))*(-CS_43*max_lambda_44 + CF_43));

       beta_3 = ((2107.0/960.0))*((-CS_44*max_lambda_44 + CF_44)*(-CS_44*max_lambda_44 + CF_44)) +
            ((1.0/2.0))*(-(1567.0/80.0)*(-CS_44*max_lambda_44 + CF_44) + ((11003.0/480.0))*(-CS_43*max_lambda_44 +
            CF_43))*(-CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(-CS_41*max_lambda_44 +
            CF_41)*(-(309.0/80.0)*(-CS_44*max_lambda_44 + CF_44) + ((547.0/480.0))*(-CS_41*max_lambda_44 + CF_41) +
            ((2321.0/240.0))*(-CS_43*max_lambda_44 + CF_43)) + ((1.0/2.0))*(-CS_42*max_lambda_44 +
            CF_42)*(-(8623.0/240.0)*(-CS_43*max_lambda_44 + CF_43) - (647.0/80.0)*(-CS_41*max_lambda_44 + CF_41) +
            ((3521.0/240.0))*(-CS_44*max_lambda_44 + CF_44) + ((7043.0/480.0))*(-CS_42*max_lambda_44 + CF_42));

       alpha_0 = 0.0285714285714286 + ((1.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.342857142857143 + ((12.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_1)*(1.0e-40 + beta_1));

       alpha_2 = 0.514285714285714 + ((18.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_2)*(1.0e-40 + beta_2));

       alpha_3 = 0.114285714285714 + ((4.0/35.0))*(fabs(-beta_3 - 3*beta_2 + 3*beta_1 + beta_0)*fabs(-beta_3 - 3*beta_2
            + 3*beta_1 + beta_0))/((1.0e-40 + beta_3)*(1.0e-40 + beta_3));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      omega_2 = alpha_2*inv_alpha_sum;

      omega_3 = alpha_3*inv_alpha_sum;

       rj4 = fmax(rj4, 0.027027027027027*fabs(-1.0 + 35*omega_0) + 0.027027027027027*fabs(-1.0 + ((35.0/4.0))*omega_3) +
            0.027027027027027*fabs(-1.0 + ((35.0/12.0))*omega_1) + 0.027027027027027*fabs(-1.0 +
            ((35.0/18.0))*omega_2));

       Recon_4 = (-(23.0/24.0)*(-CS_45*max_lambda_44 + CF_45) - (1.0/8.0)*(-CS_47*max_lambda_44 + CF_47) +
            ((13.0/24.0))*(-CS_46*max_lambda_44 + CF_46) + ((25.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_0 +
            (-(5.0/24.0)*(-CS_42*max_lambda_44 + CF_42) + ((1.0/8.0))*(-CS_44*max_lambda_44 + CF_44) +
            ((1.0/24.0))*(-CS_41*max_lambda_44 + CF_41) + ((13.0/24.0))*(-CS_43*max_lambda_44 + CF_43))*omega_3 +
            (-(5.0/24.0)*(-CS_45*max_lambda_44 + CF_45) + ((1.0/8.0))*(-CS_43*max_lambda_44 + CF_43) +
            ((1.0/24.0))*(-CS_46*max_lambda_44 + CF_46) + ((13.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_1 +
            (-(1.0/24.0)*(-CS_42*max_lambda_44 + CF_42) - (1.0/24.0)*(-CS_45*max_lambda_44 + CF_45) +
            ((7.0/24.0))*(-CS_43*max_lambda_44 + CF_43) + ((7.0/24.0))*(-CS_44*max_lambda_44 + CF_44))*omega_2 +
            Recon_4;

       Recon_0 = (((1.0/840.0))*(-533*CF_03 - 533*CF_04 - 29*CF_01 - 29*CF_06 + 3*CF_00 + 3*CF_07 + 139*CF_02 +
            139*CF_05) + Recon_0)*rj0;

       Recon_1 = (((1.0/840.0))*(-533*CF_13 - 533*CF_14 - 29*CF_11 - 29*CF_16 + 3*CF_10 + 3*CF_17 + 139*CF_12 +
            139*CF_15) + Recon_1)*rj1;

       Recon_2 = (((1.0/840.0))*(-533*CF_23 - 533*CF_24 - 29*CF_21 - 29*CF_26 + 3*CF_20 + 3*CF_27 + 139*CF_22 +
            139*CF_25) + Recon_2)*rj2;

       Recon_3 = (((1.0/840.0))*(-533*CF_33 - 533*CF_34 - 29*CF_31 - 29*CF_36 + 3*CF_30 + 3*CF_37 + 139*CF_32 +
            139*CF_35) + Recon_3)*rj3;

       Recon_4 = (((1.0/840.0))*(-533*CF_43 - 533*CF_44 - 29*CF_41 - 29*CF_46 + 3*CF_40 + 3*CF_47 + 139*CF_42 +
            139*CF_45) + Recon_4)*rj4;

       rho_RKold_B0(0,0,0) = 0.707106781186547*AVG_2_rho*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_2_rho*Recon_4*inv_AVG_a + Recon_2;

       rhou0_RKold_B0(0,0,0) = AVG_2_u0*Recon_2 - AVG_2_rho*Recon_1 +
            0.707106781186547*AVG_2_rho*AVG_2_u0*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_2_rho*AVG_2_u0*Recon_4*inv_AVG_a;

       rhou1_RKold_B0(0,0,0) = AVG_2_rho*Recon_0 + AVG_2_u1*Recon_2 +
            0.707106781186547*AVG_2_rho*AVG_2_u1*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_2_rho*AVG_2_u1*Recon_4*inv_AVG_a;

       rhou2_RKold_B0(0,0,0) = AVG_2_u2*Recon_2 + 0.707106781186547*(-AVG_2_a + AVG_2_u2)*AVG_2_rho*Recon_4*inv_AVG_a +
            0.707106781186547*(AVG_2_a + AVG_2_u2)*AVG_2_rho*Recon_3*inv_AVG_a;

       rhoE_RKold_B0(0,0,0) = (((1.0/2.0))*(AVG_2_u0*AVG_2_u0) + ((1.0/2.0))*(AVG_2_u1*AVG_2_u1) +
            ((1.0/2.0))*(AVG_2_u2*AVG_2_u2))*Recon_2 + AVG_2_rho*AVG_2_u1*Recon_0 - AVG_2_rho*AVG_2_u0*Recon_1 +
            0.707106781186547*(((AVG_2_a*AVG_2_a) + ((1.0/2.0))*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) +
            (AVG_2_u2*AVG_2_u2))*gamma_m1)*invgamma_m1 + AVG_2_a*AVG_2_u2)*AVG_2_rho*Recon_3*inv_AVG_a +
            0.707106781186547*(((AVG_2_a*AVG_2_a) + ((1.0/2.0))*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) +
            (AVG_2_u2*AVG_2_u2))*gamma_m1)*invgamma_m1 - AVG_2_a*AVG_2_u2)*AVG_2_rho*Recon_4*inv_AVG_a;

   }

   else{

      rho_RKold_B0(0,0,0) = 0.0;

      rhou0_RKold_B0(0,0,0) = 0.0;

      rhou1_RKold_B0(0,0,0) = 0.0;

      rhou2_RKold_B0(0,0,0) = 0.0;

      rhoE_RKold_B0(0,0,0) = 0.0;

   }

}

 void opensbliblock00Kernel061(const ACC<double> &D11_B0, const ACC<double> &Residual0_B0, const ACC<double>
&Residual1_B0, const ACC<double> &Residual2_B0, const ACC<double> &Residual3_B0, const ACC<double> &Residual4_B0, const
ACC<double> &kappa_B0, const ACC<double> &rhoE_RKold_B0, const ACC<double> &rho_RKold_B0, const ACC<double>
&rhou0_RKold_B0, const ACC<double> &rhou1_RKold_B0, const ACC<double> &rhou2_RKold_B0, const ACC<double> &wk0_B0, const
ACC<double> &wk1_B0, const ACC<double> &wk2_B0, const ACC<double> &wk3_B0, const ACC<double> &wk4_B0, ACC<double>
&WENO_filter_B0, ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0, ACC<double>
&rhou2_B0, const int *idx)
{
   double Grid_0 = 0.0;
   double Grid_1 = 0.0;
   double Grid_2 = 0.0;
   double Wall = 0.0;
   Grid_0 = idx[0];

   Grid_1 = idx[1];

   Grid_2 = idx[2];

   Wall = ((Grid_1 >= -6 + block0np1 || Grid_1 <= 5) ? (
   0
)
: (
   1
));

    WENO_filter_B0(0,0,0) = ((fmax(kappa_B0(-1,0,0), fmax(kappa_B0(0,-1,0), fmax(kappa_B0(0,0,2), fmax(kappa_B0(0,0,-1),
      fmax(kappa_B0(2,0,0), fmax(kappa_B0(0,2,0), fmax(kappa_B0(0,0,1), fmax(kappa_B0(1,0,0), fmax(kappa_B0(0,1,0),
      kappa_B0(0,0,0)))))))))) >= Ducros_select) ? (
   1
)
: (
   0.0
));

    rho_B0(0,0,0) = (-(-wk0_B0(-1,0,0) + wk0_B0(0,0,0))*inv_rfact0_block0 - (-rho_RKold_B0(0,0,-1) +
      rho_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual0_B0(0,-1,0) +
      Residual0_B0(0,0,0))*inv_rfact1_block0*D11_B0(0,0,0)*Wall)*dt*WENO_filter_B0(0,0,0) + rho_B0(0,0,0);

    rhou0_B0(0,0,0) = (-(-wk1_B0(-1,0,0) + wk1_B0(0,0,0))*inv_rfact0_block0 - (-rhou0_RKold_B0(0,0,-1) +
      rhou0_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual1_B0(0,-1,0) +
      Residual1_B0(0,0,0))*inv_rfact1_block0*D11_B0(0,0,0)*Wall)*dt*WENO_filter_B0(0,0,0) + rhou0_B0(0,0,0);

    rhou1_B0(0,0,0) = (-(-wk2_B0(-1,0,0) + wk2_B0(0,0,0))*inv_rfact0_block0 - (-rhou1_RKold_B0(0,0,-1) +
      rhou1_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual2_B0(0,-1,0) +
      Residual2_B0(0,0,0))*inv_rfact1_block0*D11_B0(0,0,0)*Wall)*dt*WENO_filter_B0(0,0,0) + rhou1_B0(0,0,0);

    rhou2_B0(0,0,0) = (-(-wk3_B0(-1,0,0) + wk3_B0(0,0,0))*inv_rfact0_block0 - (-rhou2_RKold_B0(0,0,-1) +
      rhou2_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual3_B0(0,-1,0) +
      Residual3_B0(0,0,0))*inv_rfact1_block0*D11_B0(0,0,0)*Wall)*dt*WENO_filter_B0(0,0,0) + rhou2_B0(0,0,0);

    rhoE_B0(0,0,0) = (-(-wk4_B0(-1,0,0) + wk4_B0(0,0,0))*inv_rfact0_block0 - (-rhoE_RKold_B0(0,0,-1) +
      rhoE_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual4_B0(0,-1,0) +
      Residual4_B0(0,0,0))*inv_rfact1_block0*D11_B0(0,0,0)*Wall)*dt*WENO_filter_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel070(ACC<double> &E_mean_B0, ACC<double> &M_mean_B0, ACC<double> &TT_mean_B0, ACC<double>
&T_mean_B0, ACC<double> &a_mean_B0, ACC<double> &mu_mean_B0, ACC<double> &p_mean_B0, ACC<double> &pp_mean_B0,
ACC<double> &rhomean_B0, ACC<double> &rhou0mean_B0, ACC<double> &rhou0u0mean_B0, ACC<double> &rhou1mean_B0, ACC<double>
&rhou1u0mean_B0, ACC<double> &rhou1u1mean_B0, ACC<double> &rhou2mean_B0, ACC<double> &rhou2u0mean_B0, ACC<double>
&rhou2u1mean_B0, ACC<double> &rhou2u2mean_B0, ACC<double> &u0mean_B0, ACC<double> &u0u0mean_B0, ACC<double> &u1mean_B0,
ACC<double> &u1u0mean_B0, ACC<double> &u1u1mean_B0, ACC<double> &u2mean_B0, ACC<double> &u2u0mean_B0, ACC<double>
&u2u1mean_B0, ACC<double> &u2u2mean_B0)
{
   rhomean_B0(0,0,0) = invniter*rhomean_B0(0,0,0);

   rhou0mean_B0(0,0,0) = invniter*rhou0mean_B0(0,0,0);

   rhou1mean_B0(0,0,0) = invniter*rhou1mean_B0(0,0,0);

   rhou2mean_B0(0,0,0) = invniter*rhou2mean_B0(0,0,0);

   E_mean_B0(0,0,0) = invniter*E_mean_B0(0,0,0);

   u0mean_B0(0,0,0) = invniter*u0mean_B0(0,0,0);

   u1mean_B0(0,0,0) = invniter*u1mean_B0(0,0,0);

   u2mean_B0(0,0,0) = invniter*u2mean_B0(0,0,0);

   u0u0mean_B0(0,0,0) = invniter*u0u0mean_B0(0,0,0);

   u1u0mean_B0(0,0,0) = invniter*u1u0mean_B0(0,0,0);

   u1u1mean_B0(0,0,0) = invniter*u1u1mean_B0(0,0,0);

   u2u0mean_B0(0,0,0) = invniter*u2u0mean_B0(0,0,0);

   u2u1mean_B0(0,0,0) = invniter*u2u1mean_B0(0,0,0);

   u2u2mean_B0(0,0,0) = invniter*u2u2mean_B0(0,0,0);

   p_mean_B0(0,0,0) = invniter*p_mean_B0(0,0,0);

   pp_mean_B0(0,0,0) = invniter*pp_mean_B0(0,0,0);

   a_mean_B0(0,0,0) = invniter*a_mean_B0(0,0,0);

   T_mean_B0(0,0,0) = invniter*T_mean_B0(0,0,0);

   TT_mean_B0(0,0,0) = invniter*TT_mean_B0(0,0,0);

   mu_mean_B0(0,0,0) = invniter*mu_mean_B0(0,0,0);

   M_mean_B0(0,0,0) = invniter*M_mean_B0(0,0,0);

   rhou0u0mean_B0(0,0,0) = invniter*rhou0u0mean_B0(0,0,0);

   rhou1u0mean_B0(0,0,0) = invniter*rhou1u0mean_B0(0,0,0);

   rhou1u1mean_B0(0,0,0) = invniter*rhou1u1mean_B0(0,0,0);

   rhou2u0mean_B0(0,0,0) = invniter*rhou2u0mean_B0(0,0,0);

   rhou2u1mean_B0(0,0,0) = invniter*rhou2u1mean_B0(0,0,0);

   rhou2u2mean_B0(0,0,0) = invniter*rhou2u2mean_B0(0,0,0);

}

#endif
