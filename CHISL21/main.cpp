#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <format>
#include <math.h>
#include <bits/ranges_algo.h>

using namespace std;
using grid = pair<vector<double>, vector<double>>;

grid read_grid(string filename) {
    vector<double> x;
    vector<double> y;
    int n;

    ifstream file(filename, ios::binary);

    file.read(reinterpret_cast<char *>(&n), sizeof(int));
    x.resize(n);
    y.resize(n);

    file.read(reinterpret_cast<char *>(x.data()), n * sizeof(double));
    file.read(reinterpret_cast<char *>(y.data()), n * sizeof(double));

    file.close();

    return {x, y};
}

double lagrange(double x, const grid& A) {
    auto& [X, Y] = A;
    int n = X.size();
    double res = 0;

    for(int i = 0; i < n; i++) {
        double addition = Y[i];
        for(int j = 0; j < n; j++) {
            if(i != j)
                addition *= (x - X[j]) / (X[i] - X[j]);
        }
        res += addition;
    }

    return res;
}

void find_max() {
    ofstream out("rpoly_RES.csv");

    for(int i = 2; i <= 100; i++) {
        string in = format("./data/rpoly{}_SRC.bin", i);
        grid A = read_grid(in);
        double max_err = 0;
        double x_index = 0;

        for(int j = 0; j < A.first.size() - 1; j++) {
            double x = (A.first[j] + A.first[j + 1]) / 2;
            if(max_err < abs(lagrange(x, A) - tan(x) - cos(x) - 0.1)) {
                max_err = abs(lagrange(x, A) - tan(x) - cos(x) - 0.1);
                x_index = x;
            }
        }

        out << x_index << ";" << max_err << endl;
    }

    out.close();
}

void interpolate_function(string in, string out, int n, pair<double, double> range) {
    double pivot = range.first;
    double step = (range.second - range.first) / n;
    grid A = read_grid(in);

    ofstream fout(out);
    fout.precision(16);

    vector<pair<double, double>> res;

    for(int i = 0; i < n; i++) {
        double y = lagrange(pivot, A);
        res.push_back({pivot, y});

        pivot += step;
    }

    for(auto& x : A.first) {
        double y = lagrange(x, A);
        res.push_back({x, y});
    }

    sort(res.begin(), res.end());

    for(auto& p : res) {
        fout << p.first << ";" << p.second << endl;
    }

    fout.close();
}

int main() {
    interpolate_function("poly7_SRC.bin",
        "poly7_RES.csv", 1000, {-0.5, 0.5});

    find_max();

    return 0;
}
