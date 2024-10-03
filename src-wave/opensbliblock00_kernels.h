#ifndef OPENSBLIBLOCK00_KERNEL_H
#define OPENSBLIBLOCK00_KERNEL_H
void opensbliblock00Kernel000(ACC<double> &phi_B0, ACC<double> &x0_B0, const int *idx) {
  x0_B0(0) = Delta0block0 * idx[0];

  phi_B0(0) = sin(2.0 * M_PI * x0_B0(0));
}

void opensbliblock00Kernel007(const ACC<double> &phi_B0, ACC<double> &phi_RKold_B0) { phi_RKold_B0(0) = phi_B0(0); }

void opensbliblock00Kernel003(const ACC<double> &phi_B0, ACC<double> &wk0_B0) {
  wk0_B0(0) =
      (-(2.0 / 3.0) * phi_B0(-1) - (1.0 / 12.0) * phi_B0(2) + ((1.0 / 12.0)) * phi_B0(-2) + ((2.0 / 3.0)) * phi_B0(1)) *
      invDelta0block0;
}

void opensbliblock00Kernel004(const ACC<double> &wk0_B0, ACC<double> &Residual0_B0) {
  Residual0_B0(0) = -c0 * wk0_B0(0);
}

void opensbliblock00Kernel009(const ACC<double> &Residual0_B0, const ACC<double> &phi_RKold_B0, ACC<double> &phi_B0,
                              const double *rknew) {
  phi_B0(0) = rknew[0] * dt * Residual0_B0(0) + phi_RKold_B0(0);
}

void opensbliblock00Kernel008(const ACC<double> &Residual0_B0, ACC<double> &phi_RKold_B0, const double *rkold) {
  phi_RKold_B0(0) = rkold[0] * dt * Residual0_B0(0) + phi_RKold_B0(0);
}

#endif
