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

void printmatrix(int r, int c, double mat[r][c]) {
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      printf("[%d,%d] %3f\t", i, j, mat[i][j]);
    }
    printf("\n");
  }
}

void m_m_multiply(int m1_rows, int m1_cols, double m1[m1_rows][m1_cols],
                  int m2_rows, int m2_cols, double m2[m2_rows][m2_cols],
                  int r_r, int r_c, double result[r_r][r_c]) {
  int m, n, p, q, c, d, k;
  double sum = 0.0;
  m = m1_rows;
  n = m1_cols;
  p = m2_rows;
  q = m2_cols;

  if (n != p)
    printf("The matrices can't be multiplied with each other.\n");
  else {
    for (c = 0; c < m; c++) {
      for (d = 0; d < q; d++) {
        for (k = 0; k < p; k++) {
          sum = sum + m1[c][k] * m2[k][d];
        }
        // printf("[%d,%d]summe:%f\n", c, d, sum);
        result[c][d] = sum;
        sum = 0;
      }
    }
  }
}

void gsl_gauss(int m1_rows, int m1_cols, double m1[m1_rows][m1_cols], int r_r,
               int r_c, double vector[r_r][r_c], double result[2]) {
  // Declare pointer variables for two gsl vectors, a matrix,
  // and a permutation list
  gsl_vector *b, *x;
  gsl_matrix *A;
  gsl_permutation *p;
  int s;
  // Create the vector and the matrix
  b = gsl_vector_alloc(r_c);
  A = gsl_matrix_alloc(m1_rows, m1_cols);

  // load them with values
  for (int i = 0; i < r_c; i++) {
    gsl_vector_set(b, i, vector[i][r_c]);
  }

  for (int i = 0; i < m1_rows; i++) {
    for (int j = 0; j < m1_cols; j++) {
      gsl_matrix_set(A, i, j, m1[i][j]);
    }
  }

  // Solve the system A x = b
  // Start by creating a vector to receive the result
  // and a permutation list to pass to the LU decomposition
  x = gsl_vector_alloc(r_c);
  p = gsl_permutation_alloc(r_c);

  // Do the LU decomposition on A and use it to solve the system
  gsl_linalg_LU_decomp(A, p, &s);
  gsl_linalg_LU_solve(A, p, b, x);

  // Print the result
  printf("The solution:\n");

  gsl_vector_fprintf(stdout, x, "%f");
  for (int i = 0; i < r_c; i++) {
    result[i] = gsl_vector_get(x, i);
  }
  // Clean up
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_permutation_free(p);
  gsl_matrix_free(A);
}

void transpose(int m1_rows, int m1_cols, double m1[m1_rows][m1_cols], int r_r,
               int r_c, double result[r_r][r_c]) {
  {
    int m, n, c, d;

    m = m1_rows;
    n = m1_cols;
    c = r_r;
    d = r_c;

    for (c = 0; c < m; c++)
      for (d = 0; d < n; d++) result[d][c] = m1[c][d];

    // printf("Transpose of the matrix:\n");

    // for (c = 0; c < n; c++) {
    //   for (d = 0; d < m; d++) printf("%d\t", m1[c][d]);
    //   printf("\n");
  }
}
int main(int argc, char *argv[]) {
  // Abzählen der Wertepaare
  int N = getNumberofPoints("input.dat");

  double x[N];  // Vektor für den Abstand der Messung
  double y[N];  // Vektor für den gemessenen Wert
  //   Einlesen der Daten
  readFile("input.dat", x, y);

  double lambda[2] = {0.0, 0.0};  // Koeffizient für die Funktion f1 , f2
  double mat_A[N][2];

  double mat_A_tr[2][N];
  double mat_A_t_A[2][2];
  double vector[1][2];
  double yy[N][1];

  for (int i = 0; i < N; i++) {
    mat_A[i][0] = f2(x[i]);
    mat_A[i][1] = f1(x[i]);
  }
  for (int i = 0; i < N; i++) {
    yy[i][1] = x[i];
  }
  printf("MATRIX A\n");
  printmatrix(N, 2, mat_A);
  transpose(N, 2, mat_A, 2, N, mat_A_tr);
  printf("matrix Atr\n");
  printmatrix(2, N, mat_A_tr);
  printf("matrix Atr vollständig\n");
  m_m_multiply(2, N, mat_A_tr, N, 2, mat_A, 2, 2, mat_A_t_A);
  printf("MATRIX Atr x A \n");
  printmatrix(2, 2, mat_A_t_A);
  printf("matrix a_tr x y\n");
  m_m_multiply(2, N, mat_A_tr, N, 1, yy, 1, 2, vector);
  printmatrix(1, 2, vector);
  gsl_gauss(2, 2, mat_A_t_A, 1, 2, vector, lambda);

  for (int i = 0; i < 2; i++) {
    printf("lambda %d : %f \n", i, lambda[i]);
  }

  //_____________________________________________________________________
  // benötigte Variablen einlegen und initialisieren.

  //   Berechnung von lambda1

  // Plotten wenn plotflag!=0
  long plotflag = 1;

  if (plotflag) {
    char plotbefehl[10000];
    sprintf(plotbefehl,
            " echo "
            " ' reset; set key right top Left box; \n"
            "set term jpeg;\n"
            "set output plot.jpg;\n"
            " set xrange [0.0:100]; set xlabel \"x\"; \n"
            " set yrange [0.0:0.75]; set ylabel \"y\"; \n"
            " f(x) = %le*x*(1+x)**(-1)+%le*(2+x)**(-1);\n "
            " plot f(x) lt -1 lw 2, \"input.dat\" u 1:2 pt 7;\n"
            " ' | gnuplot -persist>>test.jpg",
            lambda[0], lambda[1]);

    system(plotbefehl);
  }

  return 0;
}
