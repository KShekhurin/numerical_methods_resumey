#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

double f(double x) {
    return pow(x, 5) - 3.2*pow(x, 3) + 9.5*pow(x, 2) - 7*x - 7.5;
}

double integrate_segment_38(double a, double h) {
    return (f(a) + 3*(f(a + h) + f(a + 2*h)) + f(a + 3*h));
}

vector<double> gen_linspace(double a, double b, int num_of_segments) {
    int n = 1 + num_of_segments * 3;
    vector<double> x(n);

    for(int i = 0; i < n; i++) {
        x[i] = a + (b - a) / (n - 1) * i;
    }

    return x;
}

double integrate_range(double a, double b, int segments_cnt) {
    auto g = gen_linspace(a, b, segments_cnt);
    int n = g.size();
    double h = g[1] - g[0];
    double res = 0;

    for(int i = 0; i < n - 2; i += 3) {
        res += integrate_segment_38(g[i], h);
    }

    return 3*h/8*res;
}

double calc_rest(double N1_val, double N2_val, int N1, int N2) {
    return fabs(N2_val - N1_val) / (pow(2, 4) - 1);
}

double integrate_runge(double a, double b, double eps, int& N) {
    N = 0;
    int N1 = 1;
    int N2 = 1;

    double N1_val = 0;
    double N2_val = integrate_range(a, b, N2);

    do {
        N1 = N2;
        N2 = N1 * 2;
        N1_val = N2_val;
        N2_val = integrate_range(a, b, N2);
        N++;
    } while (calc_rest(N1_val, N2_val, N1, N2) > eps);

    return N2_val;
}

int main() {
    const double I = -12.712372166666665;

    const double a = -0.7;
    const double b = 1.6;

    ofstream outfile("err_eps.csv");
    ofstream noutfile("n_eps.csv");
    outfile.precision(16);

    std::cout << std::fixed << std::setprecision(16);

    for (int i = 0; i < 13; i++) {
        double eps = pow(10, -1 * i);
        int N = 0;
        auto res = integrate_runge(a, b, eps, N);

        outfile << eps << ";" << fabs(res - I) << endl;
        noutfile << eps << ";" << N << endl;
    }

    outfile.close();
    noutfile.close();

    ofstream houtfile("Err_h.csv");

    double prev = 99999999999;
    for(int i = 1; i < 1024; i++) {
        auto g = gen_linspace(a, b, i);
        double h = g[1] - g[0];
        double res = integrate_range(a, b, i);

        if(fabs(log10(fabs(res - I)) - log10(prev)) >= 0.1) {
            houtfile << h << ";" << fabs(res - I) << endl;
            prev = fabs(res - I);
        }
    }

    cout << fabs(integrate_range(a, b, 1024) - I) << endl;

    houtfile.close();

    return 0;
}
