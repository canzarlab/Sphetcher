#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <iomanip>      // std::setprecision
#include <omp.h>
#include "util_functions.h"
#include <chrono>
using namespace std;

template <class T>
vector<T> HDscore(const vector<vector<T>> & distance_mat, const vector<int> & setcover, int n_sets, int n_dim, int k)
{
    vector<T> scores(n_sets);
    int i,j;
    int cover_size = (int) setcover.size();
    #pragma omp parallel for private(j)
    for (i=0; i < n_sets; ++i){
        T dmin = correlationDistance(distance_mat.at(i), distance_mat.at(setcover[0]), n_dim);
        for (j=1; j< cover_size; ++j){
            T dd = correlationDistance(distance_mat.at(i), distance_mat.at(setcover[j]), n_dim);
            if (dd < dmin){
                dmin = dd;
            }
        }
        scores[i] = dmin;
    }

    sort(scores.begin(), scores.end());
    vector<T> newscore(scores.end()-k, scores.end());
    return newscore;

}

int main(int argc, char* argv[])
{
    std::string dataname = argv[1]; // zeisel grun
    std::string dfilename = "../../../spectral_clustering_matlab/data/"+dataname+"_pca.csv";
    auto start = chrono::steady_clock::now();
    vector<vector<double>> distance_mat = parse2DCsvFile2Double(dfilename);
    auto end = chrono::steady_clock::now();
    cout << "read file in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    int n_sets = (int) distance_mat.size();
    int n_dim =  (int) distance_mat[0].size();
    cout << "# of items: " << n_sets << ", n_dim: " << n_dim << endl;

    bool indicator_solution = true;
    vector<double> hdscore_vec;
    if(indicator_solution){
        // string solfilename = "output/" + dataname+"_setcover_indicator_solutions.csv";
        string solfilename = "output/" + dataname+"_box_solutions.csv";
        vector<vector<int>> solutions = parse2DCsvFile(solfilename);
        int q = max(1, (int) (1e-4* (float) n_sets));
        for (auto indicator_sol: solutions){
            vector<int> sol;
            for (int i=0; i < (int) indicator_sol.size(); ++i){
                if (indicator_sol.at(i)){
                    sol.push_back(i);
                }
            }
            vector<double> hdscore = HDscore(distance_mat, sol, n_sets, n_dim, q);
            hdscore_vec.push_back(hdscore.front());
            cout << "subsample size: " << sol.size() << ", score: ";
            for (size_t i=0; i < hdscore.size(); ++i){
                cout << std::fixed << std::setprecision(2) << hdscore.at(i) << "  ";
            }
            cout << "\n";
        }
    }
    else{
        string solfilename = "output/" + dataname+"_setcover_solutions.csv";
        vector<vector<int>> solutions = parse2DCsvFile(solfilename);
        int q = max(1, (int) (1e-4* (float) n_sets));
        for (auto sol: solutions){
            vector<double> hdscore = HDscore(distance_mat, sol, n_sets, n_dim, q);
            hdscore_vec.push_back(hdscore.front());
            for (size_t i=0; i < hdscore.size(); ++i){
                cout << std::fixed << std::setprecision(2) << hdscore.at(i) << "  ";
            }
            cout << "\n";
        }
    }
    vector<vector<double>> hd2D_vec;
    int n_runs = 10;
    int n_pct = 6;
    for (int i=0; i < n_pct; ++i){
        vector<double> hd2D;
        for (int j=0; j < n_runs; ++j){
            hd2D.push_back(hdscore_vec.at(i + j*n_pct));
        }
        hd2D_vec.push_back(hd2D);
    }
    writeVec2File("output/" + dataname + "_box_hdscore.csv", hd2D_vec);

    return 0;
}

