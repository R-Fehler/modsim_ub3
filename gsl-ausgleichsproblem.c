#include <gsl/gsl_linalg.h>
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
void printmatrix(gsl_matrix *matB, int N, int M) {
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < M; j++)
      printf("m(%d,%d) = %g\n", i, j, gsl_matrix_get(matB, i, j));
}
// --------------------------------------------------------------------
// --------------------------------------------------------------------

// Definition der Ansatzfunktionen f1 und f2

double f1(double x) { return (x / (1 + x)); }
double f2(double x) { return (1 / 2 + x); }
int main(int argc, char *argv[]) {
  // Abzählen der Wertepaare
  int N = getNumberofPoints("input.dat");

  double x[N];  // Vektor für den Abstand der Messung
  double y[N];  // Vektor für den gemessenen Wert
  double matData[2 * N];
  //   Einlesen der Daten
  readFile("input.dat", x, y);

  double lambda1 = 0.0;  // Koeffizient für die Funktion f1
  double lambda2 = 0.0;  // Koeffizient für die Funktion f2
  double f1_values[N], f2_values[N];
  gsl_vector_view f1_vector = gsl_vector_view_array(x, N);
  gsl_vector_view f2_vector = gsl_vector_view_array(y, N);
  gsl_vector_fprintf(stdout, &f1_vector.vector, "%g");

  // for (size_t i = 0; i < N; i++) {
  //   f1_values[N] = f1(N);
  //   f2_values[N] = f2(N);
  // }
  for (size_t i = 0; i < N; i++) {
    matData[i] = x[i];
    matData[i + N] = y[i];
  }
  /////////////////////////////////////
  // gsl_matrix *m = gsl_matrix_alloc(10, 10);

  // for (i = 0; i < 10; i++)
  //   for (j = 0; j < 10; j++) gsl_matrix_set(m, i, j, sin(i) + cos(j));

  // for (j = 0; j < 10; j++) {
  //   gsl_vector_view column = gsl_matrix_column(m, j);
  //   double d;

  //   d = gsl_blas_dnrm2(&column.vector);

  //   printf("matrix column %zu, norm = %g\n", j, d);
  // }

  // gsl_matrix_free(m);
  /////////////////////
  gsl_matrix *matB = gsl_matrix_alloc(N, 2);
  gsl_matrix *matB_trans = gsl_matrix_alloc(2, N);
  gsl_matrix *matC = gsl_matrix_alloc(2, 2);
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < 2; j++) {
      gsl_matrix_set(matB, i, 0, matData[i]);
      gsl_matrix_set(matB, i, 1, matData[N + i]);
    }
  }
  // print the matrix
  printmatrix(matB, N, 2);
  gsl_matrix_transpose_memcpy(matB_trans, matB);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &matB_trans, &matB, 0.0,
                 &matC);

  gsl_matrix_view matA = gsl_matrix_view_array(matData, N, 2);
  // gsl_matrix_fprintf(stdout, &matA.matrix, "%g");
  //_____________________________________________________________________
  // benötigte Variablen einlegen und initialisieren.

  //   Berechnung von lambda1

  // Plotten wenn plotflag!=0
  long plotflag = 0;

  if (plotflag) {
    char plotbefehl[10000];
    sprintf(plotbefehl,
            " echo "
            " ' reset; set key right top Left box; \n"
            " set xrange [0.0:100]; set xlabel \"x\"; \n"
            " set yrange [0.0:0.75]; set ylabel \"y\"; \n"
            " f(x) = %le*x*(1+x)**(-1)+%le*(2+x)**(-1);\n "
            " plot f(x) lt -1 lw 2, \"input.dat\" u 1:2 pt 7;\n"
            " set term eps; ' "
            " | gnuplot -persist",
            lambda1, lambda2);

    system(plotbefehl);
  }

  return 0;
}
