#ifndef DIFF
#define DIFF
#include <string>
void differ_backward(double h, int o, int p, double c[], double x[]);
void differ_central(double h, int o, int p, double c[], double x[]);
void differ_forward(double h, int o, int p, double c[], double x[]);
double *differ_inverse(int n, double stencil[]);
double *differ_matrix(int n, double stencil[]);
double *differ_solve(int n, double stencil[], int order);
void differ_stencil(double x0, int o, int p, double x[], double c[]);
int i4_max(int i1, int i2);
int i4_min(int i1, int i2);
std::string i4_to_string(int i4);
double inverse_error(int n, double a[], double b[]);
double r8_factorial(int n);
double *r8mat_fs_new(int n, double a[], double b[]);
double *r8mat_mm_new(int n1, int n2, int n3, double a[], double b[]);
double *r8mat_mv_new(int m, int n, double a[], double x[]);
double r8mat_norm_fro(int m, int n, double a[]);
void r8mat_print(int m, int n, double a[], std::string title);
void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi,
                      int jhi, std::string title);
double *r8mat_sub_new(int m, int n, double a[], double b[]);
void r8vec_print(int n, double a[], std::string title);
double *r8vec_uniform_01_new(int n, int &seed);
void r8vec2_print(int n, double a1[], double a2[], std::string title);
void r8vm_sl(int n, double a[], double b[], int job, double x[], int &info);
double *r8vm_sl_new(int n, double a[], double b[], int job, int &info);
void timestamp();
#endif /* ifndef DIFF */
