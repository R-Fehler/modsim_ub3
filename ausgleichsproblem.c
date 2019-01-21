#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <stdlib.h>

/*********************************************************************/
/*                                                                   */
/* Es liegen Messdaten "input.dat" in Form von Wertepaaren vor.      */
/* Durch diese Punkte muss eine Ausgleichskurve gelegt werden.       */
/*                                                                   */
/*********************************************************************/
/*             Zu Beginn einige Hilfsfunktionen                      */
// --------------------------------------------------------------------
// Funktion um die Wertepaare abzuzählen

/*

md5sum input.dat
87f09479f34c660bc9034819d894a6f7

 ./ausgleichsproblem
lambda1: 1.028357e-03
lambda2: 1.351542e+00
 */

int getNumberofPoints(char *name) {
  FILE *fp;
  char *line = NULL;
  int cnt = 0;
  size_t len = 0;

  if ((fp = fopen(name, "r")) == NULL) {
    exit(EXIT_FAILURE);
  }
  while (getline(&line, &len, fp) != -1) {
    cnt++;
  }
  free(line);

  return cnt;
}
// In dieser Funktion werden die Wertepaare eingelesen und
// in Form von Arrays x[N] und y[N] übergeben.
void readFile(char *name, double x[], double y[]) {
  FILE *fp;
  char *line = NULL;
  size_t len = 0;

  if ((fp = fopen(name, "r")) == NULL) {
    exit(EXIT_FAILURE);
  }

  int cnt = 0;
  while (getline(&line, &len, fp) != -1) {
    sscanf(line, "%lf %lf", &x[cnt], &y[cnt]);
    cnt++;
  }

  free(line);
}
// --------------------------------------------------------------------
// --------------------------------------------------------------------

// Definition der Ansatzfunktionen f1 und f2
double f1(double x) { return (x / (1.0 + x)); }
double f2(double x) { return (1.0 / (2.0 + x)); }

void gsl_printmatrix(gsl_matrix *mat) {
  int r, c;
  r = mat->size1;
  c = mat->size2;
  printf("size r: %d c: %d\n", r, c);
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      printf("[%d,%d] %3.2f\t", i, j, gsl_matrix_get(mat, i, j));
    }
    printf("\n");
  }
}
void gsl_printvector(gsl_vector *vec) {
  int r;
  r = vec->size;
  printf("size r: %d\n", r);
  for (int i = 0; i < r; i++) {
    printf("[%d] %3.2f\t", i, gsl_vector_get(vec, i));
    printf("\n");
  }
}

void gsl_gauss(gsl_matrix *A, gsl_vector *b, gsl_vector *result) {
  // Declare pointer variables for two gsl vectors, a matrix,
  // and a permutation list
  gsl_permutation *p;
  int s;

  // Solve the system A x = b
  // Start by creating a vector to receive the result
  // and a permutation list to pass to the LU decomposition
  p = gsl_permutation_alloc(b->size);

  // Do the LU decomposition on A and use it to solve the system
  gsl_linalg_LU_decomp(A, p, &s);
  gsl_linalg_LU_solve(A, p, b, result);

  // Print the result
  printf("The solution:\n");
  gsl_vector_fprintf(stdout, result, "%e");
  // Clean up
  gsl_permutation_free(p);
}

int main(int argc, char *argv[]) {
  // Abzählen der Wertepaare
  int N = getNumberofPoints("input.dat");

  double x[N];  // Vektor für den Abstand der Messung
  double y[N];  // Vektor für den gemessenen Wert
  //   Einlesen der Daten
  readFile("input.dat", x, y);

  //_____________________________________________________________________
  // benötigte Variablen einlegen und initialisieren.
  gsl_matrix *A;
  gsl_matrix *A_tr;
  gsl_matrix *A_tr_A;
  gsl_vector *y_vec;
  gsl_vector *A_tr_y;
  gsl_vector *lambda;

  A = gsl_matrix_alloc(N, 2);
  A_tr = gsl_matrix_alloc(2, N);
  A_tr_A = gsl_matrix_alloc(2, 2);
  y_vec = gsl_vector_alloc(N);
  lambda = gsl_vector_alloc(2);
  A_tr_y = gsl_vector_alloc(2);

  for (int i = 0; i < N; i++) {
    gsl_matrix_set(A, i, 0, f2(x[i]));
    gsl_matrix_set(A, i, 1, f1(x[i]));
  }

  for (int i = 0; i < N; i++) {
    gsl_vector_set(y_vec, i, y[i]);
  }
  //   Berechnung von lambda1
  gsl_printmatrix(A);
  gsl_matrix_transpose_memcpy(A_tr, A);
  gsl_printmatrix(A_tr);
  // matrix Multipl.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_tr, A, 0.0, A_tr_A);
  gsl_printmatrix(A_tr_A);
  // vector Multiplikation
  gsl_blas_dgemv(CblasNoTrans, 1.0, A_tr, y_vec, 0.0, A_tr_y);
  gsl_printvector(A_tr_y);
  // gauss LGS Verfahren
  gsl_gauss(A_tr_A, A_tr_y, lambda);
  gsl_printvector(lambda);
  // Plotten wenn plotflag!=0
  long plotflag = 1;

  if (plotflag) {
    char plotbefehl[10000];
    sprintf(plotbefehl,
            " echo "
            " ' reset; set key right top Left box; \n"
            "set term jpeg;\n"
            " set xrange [0.0:100]; set xlabel \"x\"; \n"
            " set yrange [0.0:0.75]; set ylabel \"y\"; \n"
            " f(x) = %le*x*(1+x)**(-1)+%le*(2+x)**(-1);\n "
            " plot f(x) lt -1 lw 2, \"input.dat\" u 1:2 pt 7;\n"
            " ' | gnuplot -persist>output.jpg",
            gsl_vector_get(lambda, 1), gsl_vector_get(lambda, 0));

    system(plotbefehl);
  }

  return 0;
}
