#include <iostream>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

using matrix = vector<vector<int>>;
using vi = vector<int>;

long long cont = 0;

FILE *out;

matrix add_matrix(matrix &a, matrix &b) {
  int n = a.size(), m = a[0].size();

  matrix c(n, vi(m));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      c[i][j] = a[i][j] + b[i][j];
    }
  }

  return c;
}

matrix substract_matrix(matrix &a, matrix &b) {
  int n = a.size(), m = a[0].size();

  matrix c(n, vi(m));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      c[i][j] = a[i][j] - b[i][j];
    }
  }

  return c;
}

matrix strassen(matrix &a, matrix &b) {
  if (a.size() == 0) return {{}};
  
  if (a.size() == 1)
  {
    return {{a[0][0] * b[0][0]}};
  }

  int md = a.size() / 2;

  matrix a00 (md, vi(md)), a01 (md, vi(md)),
         a10 (md, vi(md)), a11 (md, vi(md)),
         b00 (md, vi(md)), b01 (md, vi(md)),
         b10 (md, vi(md)), b11 (md, vi(md));

  for (int i = 0; i < md; i++) {
    for (int j = 0; j < md; j++) {
      a00[i][j] = a[i][j];
      a01[i][j] = a[i][j + md];
      a10[i][j] = a[i + md][j];
      a11[i][j] = a[i + md][j + md];
      b00[i][j] = b[i][j];
      b01[i][j] = b[i][j + md];
      b10[i][j] = b[i + md][j];
      b11[i][j] = b[i + md][j + md];
    }
  }

  matrix m1, m2, m3, m4, m5, m6, m7;
  matrix a_a03 = add_matrix(a00, a11), a_b03 = add_matrix(b00, b11),
         a_a23 = add_matrix(a10, a11), s_b13 = substract_matrix(b01, b11),
         s_b10 = substract_matrix(b10, b00), a_a01 = add_matrix(a00, a01),
         s_a10 = substract_matrix(a10, a00), a_b01 = add_matrix(b00, b01),
         s_a13 = substract_matrix(a01, a11), a_b23 = add_matrix(b10, b11);

  m1 = strassen(a_a03, a_b03);
  m2 = strassen(a_a23, b00);
  m3 = strassen(a00, s_b13);
  m4 = strassen(a11, s_b10);
  m5 = strassen(a_a01, b11);
  m6 = strassen(s_a10, a_b01);
  m7 = strassen(s_a13, a_b23);

  matrix result(2 * md, vi(2 * md));

  for (int i = 0; i < md; i++) {
    for (int j = 0; j < md; j++) {
      result[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
      result[i][j + md] = m3[i][j] + m5[i][j];
      result[i + md][j] = m2[i][j] + m4[i][j];
      result[i + md][j + md] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
    }
  }

 return result;
}

matrix create_matrix(int n) {
  matrix m(n, vi(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      m[i][j] = rand() % 15;
    }
  }
  return m;
}

void solve(int n) {
  matrix a, b;
  a = create_matrix(n);
  b = create_matrix(n);
  matrix c1 (n, vi(n, 0));

  time_t ti_s, te_s, ti, te;

  ti_s = clock();
  matrix c = strassen(a, b);
  te_s = clock() - ti_s;

  ti = clock();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        c1[i][j] += a[i][k] * b[k][j];
      }
      if (c1[i][j] != c[i][j]) {
        printf("error\nn was %d\n", n);
        fprintf(out, "%d, error, error\n", n);
        return;
      }
    }
  }
  te = clock() - ti;

  fprintf(out, "%d, %f, %f\n", n, (double)te_s, (double)te);
}

int main() {
  out = fopen("out.csv", "wt");
  for (int i = 1; i < 10; i++) {
    solve(1 << i);
  }
}
