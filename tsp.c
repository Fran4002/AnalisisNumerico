#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef unsigned long long ull;

const uint UINT_MAX = -1;

FILE * f;

uint** allocate_matrix (int n) {
  uint** mat = malloc(n * sizeof(uint*));

  for (int i = 0; i < n; i++) {
    mat[i] = malloc(n * sizeof(uint));
  }

  return mat;
}

uint** generate_matrix(int n) {
  uint **mat = allocate_matrix(n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      do {
        mat[i][j] = rand() % 1000;
      } while (mat[i][j] == 0);
    }

    mat[i][i] = 0;
  }

  return mat;
}

void cp_matrix(int n, int mat[][5], uint ***ptr) {
  uint **p = malloc(sizeof(uint*) * n);
  for (int i = 0; i < n; i++) {
    p[i] = malloc(n * sizeof(uint));
    for (int j = 0; j < n; j++) {
      p[i][j] = mat[i][j];
    }
  }
  *ptr = p;
}

void print_matrix(uint **ptr, uint n) {
  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < n; j++) {
      printf("%u, ", ptr[i][j]);
    }
    printf("\n");
  }
}

void print_set(ull bs, uint size) {
  ull mask = (ull)1 << (size - 1);
  for (uint i = size; i > 0; i--) {
    printf("%d", (mask & bs) ? 1 : 0);
    mask >>= 1;
  }
  printf("\n");
}

short next_permutation(ull *bs, uint size) {
  ull mask = 2; //...010
  short last = (1 & *bs) ? 1 : 0;
  uint ones = 0, ind = size;

  if (last)
    ones++;

  for (uint i = 1; i < size; i++, mask <<= 1) {
    short actual = (mask & *bs) ? 1 : 0;
    if (actual) {
      ones++;
    } else if (actual < last) {
      //There is a 0 after a 1
      //we can swap those to create
      //new permutation
      ind = i;
      ones--;
      *bs |= mask;
      mask >>= 1;
      break;
    }
    last = ones;
  }
  //ind didn't actualize
  //can't permute
  if (ind == size) {
    return 0;
  }

  //reset the bits from ind to the quantity of ones
  for (uint i = ind; i > ones; i--) {
    *bs &= (~mask);
    mask >>= 1;
  }

  //set the bits of the number of ones
  for (uint i = ones; i > 0; i--) {
    *bs |= mask;
    mask >>= 1;
  }

  return 1;
}

uint tsp(int n, uint **W, uint ***PATH) {
  uint **D;
  uint **P;

  //Allocate memory
  D = malloc(n * sizeof(uint *));
  P = malloc(n * sizeof(uint *));
  P[0] = malloc(((ull)1 << (n - 1)) * sizeof(uint));

  //1 << (n - 1) = 2 ^ (n - 1) = # of sets
  for (int i = 1; i < n; i++) {
    D[i] = malloc(((ull)1 << (n - 1)) * sizeof(uint));
    P[i] = malloc(((ull)1 << (n - 1)) * sizeof(uint));
    if (D[i] == NULL || P[i] == NULL) {
      printf("%d not allocated\n", i);
      exit(1);
    }
  }

  //printf("Mem alloc: OK\n");

  //Initialize tsp
  for (int i = 1; i < n; i++) {
    //Path from i to 0
    D[i][0] = W[i][0];
    P[i][0] = 0;
  }
  
  //printf("Initialization: OK\n");

  for (uint k = 1; k <= n - 2; k++) {
    //Set the first k bits on 1
    ull bs = (1 << k) - 1;
    do {
      ull imask, jmask;

      //print_set(bs, n - 1);
      
      for (uint i = 1; i < n; i++) {
        imask = (ull)1 << (i - 1);
        if (imask & bs) continue;
        
        uint dis = UINT_MAX;
        uint p = 0;
        
        for (uint j = 1; j < n; j++) {
          jmask = (ull)1 << (j - 1);

          if (jmask & bs) {
            if (dis > W[i][j] + D[j][bs - jmask]) {
              dis = W[i][j] + D[j][bs - jmask];
              p = j;
            }
          }
        }
        D[i][bs] = dis;
        P[i][bs] = p;
      }
      
    } while (next_permutation(&bs, n - 1));

    //printf("Permutations with %u: OK\n", k);
  }
  
  ull set = ((ull)1 << (n - 1)) - 1;

  uint minlength = UINT_MAX;
  for (uint j = 1; j < n; j++) {
    ull mask = (ull)1 << (j - 1);

    uint len = W[0][j] + D[j][set - mask];

    if (minlength > len) {
      minlength = len;
      P[0][set] = j;
    }
  }

  *PATH = P;

  return minlength;
}

void print_path(uint n, uint **P) {
  //Mark all the nodes (except 0) as visited
  ull set = ((ull)1 << (n - 1)) - 1;
  //Move to the second node of the path
  uint j = P[0][set];

  //Print 0 as the begin
  printf("0, ");
  
  while (set) {
    ull j_mask = (ull)1 << (j - 1);
    //Unmark the actual node
    set ^= j_mask;

    printf("%u, ", j);

    //Move to the next node
    j = P[j][set];
  }

  //Print 0 as the end
  printf("0\n");
}

void solve(int n) {
  uint **W = generate_matrix(n), **P;
  time_t start, end, elapsed;
  start = clock();
  tsp(n, W, &P);
  for (int i = 0; i < n; i++) {
    free(W[i]);
    free(P[i]);
  }
  free(W);
  free(P);
  end = clock();
  elapsed = ((double)end - (double)start) / CLOCKS_PER_SEC;
  fprintf(f, "%d, %f\n", n, (double)elapsed);
}

int main(void) {
  f = fopen("bench.csv", "a+");
  for (int i = 5; i < 20; i += 5) {
    solve(i);
  }
  return 0;
}