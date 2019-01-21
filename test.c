#include <stdio.h>
#include <stdlib.h>

void m_m_multiply(const *double m1, int m1_rows, int m1_cols, const *double m2,
                  int m2_rows, int m2_cols, *double result) {
  int m, n, p, q, c, d, k, sum = 0;
  m = m1_rows;
  n = m1_cols;
  p = m2_rows;
  q = m2_cols;

  int first[10][10], second[10][10], multiply[10][10];

  if (n != p)
    printf("The matrices can't be multiplied with each other.\n");
  else {
    printf("Enter elements of second matrix\n");

    for (c = 0; c < m; c++) {
      for (d = 0; d < q; d++) {
        for (k = 0; k < p; k++) {
          sum = sum + m1[c][k] * m2[k][d];
        }

        result[c][d] = sum;
        sum = 0;
      }
    }
  }
}
int main(int argc, char const *argv[]) {
  m_m_multiply();
  return 0;
}
