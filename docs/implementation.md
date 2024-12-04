# Implementation Notes

The `stride` parameter is important. Without it, the output is a mess. The number of ghost zones does not seem to
matter.

```c
int strided_size[] = {block0np0 / stride[0], block0np1 / stride[1]};
int strided_base[] = {0, 0};
int strided_d_p[] = {5, 5};
int strided_d_m[] = {-5, -5};
double *dummy = NULL;
rho_B0_strided = ops_decl_dat(block, 1, strided_size, strided_base, strided_d_m, strided_d_p, stride, dummy, "double",
                              "rho_B0_strided");
```

The restrict stencil is what does the striding, it HAS to have the word "RESTRICT" in it and it is case sensitive. This
is because the OPS translator sees this name and generates different code to "regular" stencils.

```c
int stencil_point[] = {0, 0};
stencil2d_00 = ops_decl_stencil(2, 1, stencil_point, "stencil_00");
stencil2d_restrict_00 = ops_decl_restrict_stencil(2, 1, stencil_point, stride, "stencil_RESTRICT_00");
```
