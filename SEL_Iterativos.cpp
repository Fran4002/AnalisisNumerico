#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
using namespace std;

using ll = long long;
using dl = double long;

class Fraction {
  private:
  ll num, den;
  ll fraction_gcd (ll a, ll b) {
    if (a == 0) return b;
    else return fraction_gcd(b % a, a);
  }

  public:
  Fraction(ll num_d) {
    this->num = num_d;
    this->den = 1;
  }

  Fraction() {
    this->num = 0;
    this->den = 1;
  }

  double long decimal_value() {
    return (double long)this->num / (double long)this->den;
  }

  void reduce_fraction() {
    bool change_sign = false;

    if (num < 0) num *= -1, change_sign ^= true;
    if (den < 0) den *= -1, change_sign ^= true;

    ll gdc = fraction_gcd(num, den);

    num /= gdc;
    den /= gdc;

    if (change_sign) num *= -1;
  }

  void read_fraction() {
    char c;
    cin >> this->num >> c >> this->den;
  }

  void print_fraction() {
    cout << this->num;
    if (den != 1) cout << "/" << this->den;
  }

  Fraction operator+ (Fraction const& fraction_f) const {
    Fraction sum;
    sum.num = this->num  * fraction_f.den + fraction_f.num * this->den;
    sum.den = this->den * fraction_f.den;

    sum.reduce_fraction();

    if (sum.den == 0) {
      exit(-2);
    }

    return sum;
  }

  Fraction operator- (Fraction const &fraction_f) const {
    Fraction sustraction;
    sustraction.num = this->num * fraction_f.den - fraction_f.num * this->den;
    sustraction.den = this->den * fraction_f.den;

    sustraction.reduce_fraction();

    if (sustraction.den == 0) {
      exit(-2);
    }

    return sustraction;
  }

  Fraction operator* (Fraction const& fraction_f) const {
    Fraction product;
    product.num = fraction_f.num * this->num;
    product.den = fraction_f.den * this->den;

    product.reduce_fraction();

    if (product.den == 0) {
      exit(-2);
    }

    return product;
  }

  Fraction operator/ (Fraction const &fraction_f) const {
    Fraction qoutient, inverse;
    
    inverse.den = fraction_f.num;
    inverse.num = fraction_f.den;

    qoutient = *this * inverse;

    qoutient.reduce_fraction();

    if (qoutient.den == 0) {
      exit(-2);
    }

    return qoutient;
  }

  bool operator< (Fraction const &fraction_f) const {
    ll f1 = this->num * fraction_f.den;
    ll f2 = this->den * fraction_f.num;

    return f1 < f2;
  }

  bool operator> (Fraction const &fraction_f) const {
    ll f1 = this->num * fraction_f.den;
    ll f2 = this->den * fraction_f.num;

    return f1 > f2;
  }
};

using vf = vector<Fraction>;
using matrix = vector<vf>;

void read_matrix (matrix &mat) {
  for (ll i = 0; mat.size() > i; ++i) {
    for(ll j = 0; j < mat[i].size(); j++) {
      mat[i][j].read_fraction();
    }
  }
}

void print_matrix (matrix mat) {
  for (ll i = 0; mat.size() > i; ++i) {
    for (ll j = 0; mat[i].size() - 1 > j; ++j) {
      mat[i][j].print_fraction();
      cout << "  ";
    }

    cout << "\t| ";
    mat[i].back().print_fraction();
    cout << '\n';
  }
}

void print_row (vf row) {
  for (auto f : row) {
    //f.print_fraction();
    cout << f.decimal_value();
    cout << "  ";
  }
  cout << '\n';
}

Fraction Jacobi_calculation (matrix mat, vf solutions, ll i) {
  Fraction x, sum;
  for (ll j = 0; j < mat[i].size() - 1; j++) {
    if (j == i) continue;
    sum = sum + (mat[i][j] * solutions[j]);
    sum.reduce_fraction();
  }

  if (mat[i][i].decimal_value() == 0) {
    cout << "Couldn't calculate the " << i << "th row of Jacobi method because a_ii was 0.\n";
    exit(-1);
  }
  
  x = (mat[i].back() - sum) / mat[i][i];
  x.reduce_fraction();

  return x;
}

vf Jacobi_iterations (const matrix &mat, ll iterations, const double long t) {
  vf result(mat.size()), answer(mat.size()), best_answer(mat.size());
  Fraction aux, other;
  ll cont = 0;
  double long e, et = 1e9;
  
  while (iterations--) {
    e = 0;
    cout << "Iteration #" << cont++ << '\n';
    for (ll i = 0; i < answer.size(); i++) {
      cout << "\tX_" << i << ":\t";
      //answer[i].print_fraction();
      cout << "\t" << answer[i].decimal_value() << '\n';
    }
      
    for (ll i = 0; i < mat.size(); i++) {
      result[i] = Jacobi_calculation(mat, answer, i);

      aux = answer[i] - result[i];
      
      e = max(e, abs(aux.decimal_value()));
    }
    
    for (ll i = 0; i < mat.size(); i++) {
      answer[i] = result[i];

      if (e < et) {
        best_answer[i] = answer[i];
      }
    }

    et = min(e, et);

    if (e < t) {
      cout << "The Jacobi method converged in " << cont << " iterations.\n";
      return answer;
    }

    if (e > et * 100) break;
  }

  cout << "The Jacobi method didn't converged after " << cont << " iterations.\n";
  return best_answer;
}

vf GS_iterations (const matrix &mat, ll iterations, const double long t) {
  vf answer(mat.size()), best_answer(mat.size());
  Fraction aux;
  ll cont = 0;
  double long e, et = 1e9;

  while (iterations--) {
    e = 0;
    
    cout << "Iteration #" << cont++ << '\n';
    for (ll i = 0; i < answer.size(); i++) {
      cout << "\tX_" << i << ":\t";
      //answer[i].print_fraction();
      cout << "\t" << answer[i].decimal_value() << '\n';
    }

    for (ll i = 0; i < mat.size(); i++) {
      aux = answer[i];
      
      answer[i] = Jacobi_calculation(mat, answer, i);

      aux = answer[i] - aux;

      e = max(e, abs(aux.decimal_value()));
    }

    if (e < et) {
      for (ll i = 0; i < mat.size(); i++) {
        best_answer[i] = answer[i];
      }
      et = e;
    }

    if (e < t) {
      cout << "The Gauss-Seidel method converged in " << cont << " iterations.\n";
      return answer;
    }

    if (e > et * 100) break;
  }

  cout << "The Gauss-Seidel method didn't converged after " << cont << " iterations.\n";
  return best_answer;
}

int main() {
  ll rows, columns;
  cin >> rows >> columns;

  matrix mat(rows, vf(columns + 1));

  read_matrix(mat);

  print_matrix(mat);

  vf ans = Jacobi_iterations(mat, 500, 1e-2);
  print_row(ans);
  
  vf o_ans = GS_iterations(mat, 500, 1e-2);
  print_row(o_ans);

  //print_matrix(mat);
}
