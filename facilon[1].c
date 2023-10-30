#include<stdio.h>
#include <stdlib.h>

typedef long long int ll;
typedef unsigned int u;

typedef struct {
    ll num, den;
} Fraction;

typedef struct {
    Fraction *ptr;
    u rows, columns;
} Matrix;

ll gcd (ll, ll);
void reduce_fraction(ll *a, ll *b);
Fraction sum_fractions (Fraction a, Fraction b);
Fraction multiply_fractions (Fraction a, Fraction b);

void sum_rows(Fraction *a, Fraction *b, u size);
void multiply_rows(Fraction *a, Fraction scalar, u size);

void initMatrix (Matrix *mat, u r, u c);
void readMatrix (Matrix *mat);
Fraction matElement(Matrix *mat, u i, u j);
void print_matrix(Matrix *mat);

void swapRow (Fraction *a, Fraction *b, u size);
int make_pivot(Matrix *mat, int index, int print_steps);
void make_zeros (Matrix *mat, int index, int print_steps);
void reduce_to_row_echelon_form (Matrix *mat, int print_steps);

int main() {
    Matrix mat;
    //Fraction *solution;

    u rows, columns;
    char option;
    int print_steps;

    printf("Rows: ");
    scanf("%u", &rows);

    printf("Columns: ");
    scanf("%u", &columns);

    initMatrix(&mat, rows, columns);

    printf("Enter the matrix row by row:\n");

    readMatrix(&mat);

    printf("Print steps? y/N: ");
    scanf("%c\n", &option);

    print_steps = 1;

    printf("\nInitial matrix:\n");
    print_matrix(&mat);

    reduce_to_row_echelon_form(&mat, print_steps);

    free(mat.ptr);

    return 0;
}

//Functions
ll gcd (ll a, ll b) {
    if (a == 0) return b;
    else return gcd(b % a, a);
}

void reduce_fraction(ll* a, ll* b) {
    int flag = 0;

    if (a < 0) *a *= -1, flag ^= 1;
    if (b < 0) *b *= - 1, flag ^= 1;

    ll gdc = gcd(*a, *b);

    *a /= gdc;
    *b /= gdc;

    if (flag) *b *= -1;
}

Fraction sum_fractions (Fraction a, Fraction b) {
    Fraction sum;

    sum.num = a.num * b.den + a.den * b.num;
    sum.den = a.den * b.den;

    reduce_fraction(&sum.num, &sum.den);

    return sum;
}

Fraction multiply_fractions (Fraction a, Fraction b) {
    Fraction product;

    product.num = a.num * b.num;
    product.den = a.den * b.den;

    reduce_fraction(&product.num, &product.den);

    return product;
}

void sum_rows(Fraction *a, Fraction *b, u size) {
    for (int i = 0; i < size; i++) {
        *(a + i) = sum_fractions(*(a + i), *(b + i));
    }
}

void multiply_rows(Fraction *a, Fraction scalar, u size) {
    for (int i = 0; i < size; i++)
        *(a + i) = multiply_fractions(*(a + i), scalar);
}

void initMatrix (Matrix *mat, u r, u c) {
    mat->ptr = (Fraction*)calloc(r * (c + 1), sizeof(Fraction));
    mat->rows = r;
    mat->columns = c + 1;
}

void readMatrix (Matrix *mat) {
    Fraction temp;

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->columns; j++) {
            scanf("%lld/%lld", &temp.num, &temp.den);

            reduce_fraction(&temp.num, &temp.den);

            *(mat->ptr + i * mat->columns + j) =  temp;
        }
    }
}

Fraction matElement(Matrix *mat, u i, u j) {
    return *(mat->ptr + i * mat->columns + j);
}

void print_matrix(Matrix *mat) {
    for (u i = 0; i < mat->rows; i++) {
        for (u j = 0; j < mat->columns - 1; j++) {
            printf("%lld", matElement(mat, i, j).num);

            if (matElement(mat, i, j).den != 1) {
                printf("/%lld", matElement(mat, i, j).den);
            }

            printf("\t");
        }
        printf(" | %lld", matElement(mat, i, mat->columns - 1).num);
        if (matElement(mat, i, mat->columns - 1).den != 1) {
            printf("/%lld", matElement(mat, i, mat->columns - 1).den);
        }
        printf("\n");
    }
}

void swapRow (Fraction *a, Fraction *b, u size) {
    Fraction tmp;

    for (u i = 0; i < size; i++) {
        tmp = *(a + i);
        *(a + i) = *(b + i);
        *(b + i) = tmp;
    }
}

int make_pivot(Matrix *mat, int index, int print_steps) {
    int has_solution = 1;
    Fraction scalar;
    ll rows = mat->rows;

    if (matElement(mat, index, index).num != 1 || matElement(mat, index, index).den != 1) {
        for (u i = index + 1; i < rows; i++) {
            if (matElement(mat, index, index).num == 1 && matElement(mat, index, index).den == 1) {

                swapRow(mat->ptr + i * mat->columns, mat->ptr + index * mat->columns, mat->columns);

                if (print_steps) {
                    printf("\nSwaps the rows %u and %u\n", index + 1, i + 1);
                }

                return 1;
            }
        }
    }

    if (matElement(mat, index, index).num == 0) {
        has_solution = 0;
        for (u i = index + 1; i < rows; i++) {
            if (matElement(mat, index, index).num != 0) {
                sum_rows(mat->ptr + index * mat->columns, mat->ptr + i * mat->columns, mat->columns);
                has_solution = 1;

                if (print_steps) {
                    printf("\nSum the rows %u and %u\n", index + 1, i + 1);
                    print_matrix(mat);
                }

                break;
            }
        }
    }

    if (has_solution) {
        scalar.num = matElement(mat, index, index).den;
        scalar.den = matElement(mat, index, index).num;
        multiply_rows(mat->ptr + index * mat->columns, scalar, mat->columns);

        if (print_steps) {
            printf("\nMultiply the row %u by %lld/%lld\n", index + 1, scalar.num, scalar.den);
            print_matrix(mat);
        }
    }
    else {
        printf("\nNo solution found\n");
        print_matrix(mat);
    }

    return has_solution;
}

void make_zeros (Matrix *mat, int index, int print_steps) {
    Fraction scalar;
    Fraction *pivot_row = (Fraction*)calloc(mat->columns + 1, sizeof(Fraction));
    u rows = mat->rows;
    u columns = mat->columns;

    for (u i = 0; i < rows; i++) {
        if (i != index && matElement(mat, i, index).num != 0) {
            scalar.num = -matElement(mat, i, index).num;
            scalar.den = matElement(mat, i, index).den;
            reduce_fraction(&scalar.num, &scalar.den);

            for (u k = 0; k < mat->columns; k++) {
                *(pivot_row + k) = matElement(mat, index, k);
            }

            multiply_rows(pivot_row, scalar, columns);
            sum_rows(mat->ptr + i * columns, pivot_row, columns);

            if (print_steps) {
                printf("\nMultiply the pivot_row (%d) by %lld", index + 1, scalar.num);

                if (scalar.den != 1) {
                    printf("/%lld", scalar.den);
                }

                printf("\n");

                for (u j = 0; j < columns; j++) {
                    Fraction temp = *(pivot_row + j);
                    printf("%lld", temp.num);

                    if (temp.den != 1) {
                        printf("/%lld", temp.den);
                    }

                    printf("\t");
                }

                printf("\n\nSum the pivot row(%d) plus the row %d\n", index + 1, i + 1);
                print_matrix(mat);
            }
        }
    }

    free(pivot_row);
}

void reduce_to_row_echelon_form (Matrix *mat, int print_steps) {
    u rows = mat->rows;
    u columns = mat->columns;

    for (int i = 0; i < rows && i < columns; i++) {
        if(!make_pivot(mat, i, print_steps))
            break;
        make_zeros(mat, i, print_steps);
    }

    printf("\nFinal matrix:\n");
    print_matrix(mat);
}
