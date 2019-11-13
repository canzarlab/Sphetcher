#ifndef __GRID_COVER__
#define __GRID_COVER__

#include <vector>
#include <algorithm>
#include <iostream>

template <class T>
long int map_box(vector<T> & x, T delta, int dim_x, T N)
{
    T value = 0.0;
    T Nd = 1.0;
    for (int i=0; i < dim_x; ++i){
        value += floor(x[i]/delta)*Nd;
        Nd *= N;
    }
    return ((long int) value);
}

vector<pair<long int, long int>> unique_pair(vector<pair<long int, long int>> & solution){
    vector<pair<long int, long int>> result;
    size_t ealier_index = 0;
    for (size_t i=0; i < solution.size()-1; ++i){
        if (solution[i].second < solution[i+1].second){
            size_t randspace = 0;
            if (i > ealier_index){
                randspace = rand() % (i-ealier_index);
            }
            result.push_back(solution[i-randspace]);
            ealier_index = i;
        }
    }

    return result;
}

template<class T>
vector<int> run_box(vector<vector<T>> & distance_mat, T delta, int dim_x, T N)
{
    vector<pair<long int, long int>> solution;
    int n_sets = (int) distance_mat.size();
    for (long int i=0; i < n_sets; ++i){
        solution.push_back(make_pair(i, map_box(distance_mat[i], delta, dim_x, N)));
    }
    // sort solution according to the mapping value 
    std::sort(solution.begin(), solution.end(), [](auto &left, auto &right) {
                return left.second < right.second;
                });
    auto unique_sol = unique_pair(solution);
    vector<int> result;
    for (size_t i=0; i < unique_sol.size(); ++i){
        result.push_back(unique_sol[i].first);
    }

    return result;
}

template<class T>
vector<int> run_box_binarysearch(vector<vector<T>> & distance_mat, int dim_x, T N, int L_min, int L_max)
{
    vector<T> row1;
    for (size_t i=0; i < distance_mat.size(); ++i){
        row1.push_back(distance_mat[i][0]);
    }
    auto minmax = minmax_element(row1.begin(), row1.end());
    T delta_max = (*minmax.second - *minmax.first)*2;
    T delta_min = delta_max/1000.0;
    vector<int> result;
    int L = -1; int iter = 0;
    while ((L < L_min || L > L_max) && iter < 20)
    {
        T delta = (delta_min + delta_max)/2.0;
        result = run_box(distance_mat, delta, dim_x, N);
        L = (int) result.size();
        // cout << "iter: " << iter <<  ", L_min: " << L_min << ", L_max: " << L_max << ", box size: " << L << "\n";
        if (L < L_min){
            delta_max = delta;
        }else{
            delta_min = delta;
        }
        iter++;
    }

    return result;
}


#endif
