#ifndef _TESTTEMP_H_
#define _TESTTEMP_H_
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream
#include <string>
#include <vector>
#include <stack>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <set>

using namespace std;

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::minus<T>());
    return result;
}

template <typename T>
T distance_squared(const std::vector<T>& point_a, const std::vector<T>& point_b) {
	T d_squared = T();
	for (typename std::vector<T>::size_type i = 0; i < std::min(point_a.size(), point_b.size()); ++i) {
		auto delta = point_a.at(i) - point_b.at(i);
		d_squared += delta * delta;
	}
	return d_squared;
}

template <typename T>
T distance(const std::vector<T>& point_a, const std::vector<T>& point_b) {
	return std::sqrt(distance_squared(point_a, point_b));
}

template <typename T>
void print_vec(const std::vector<T>& vec)
{
    for (auto & i: vec)
        cout << i << " ";
    cout << endl;
}

template <typename T>
void print_vec(vector<int> & vec_index, const std::vector<T>& vec)
{
    for (auto & i: vec_index)
        cout << vec.at(i) << " ";
    cout << endl;
}

template <typename T>
void writeVec2File(std::string fileName, vector<vector<T>> & vec, std::string sep = ",")
{
    ofstream myfile;
    myfile.open (fileName);
    for (size_t i=0; i < vec.size(); ++i){
        vector<T> u = vec.at(i);
        for (size_t j=0; j < u.size()-1; ++j){
            myfile << u.at(j) << sep;
        }
        myfile << u.back() << "\n";
    }
    myfile.close();
}

void swap_vec(vector<int> & v1, vector<int> & v2);

double getJaccard (vector<int> v1, vector<int> v2); // compute the jaccard of two sets of integers v1 and v2

vector<vector<int>> parse2DCsvFile(string inputFileName);

vector<vector<string>> parse2DCsvFile2String(string inputFileName);

vector<vector<double>> parse2DCsvFile2Double(string inputFileName);

vector<int> parse1DcsvFile2Int(string inputFileName, int from_line);

void switchTwoElement(vector<int> & vec, int x, int y);

void uniqueClusterPresentation(vector<int> & vec);

int cluster_index(vector<int> vec1, vector<int> vec2);

vector<int> bestMatch(vector<int> & single_cluster, vector<vector<int>> & cluster, vector<bool> & open_option);

vector<vector<int>> findOptimalOrdering(vector<vector<int>> & ref_cluster, vector<vector<int>> & cluster);

template <typename T>
int lenIntersection(vector<T> &v1, vector<T> &v2)
{
    vector<T> v3;
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
    return (int) v3.size();
}
// compute row sum and column_sum of a matrix.
vector<int> row_sum(vector<vector<int>> & matrix);
vector<int> column_sum(vector<vector<int>> & matrix);
double _comb2(int n);
// return adjusted Random index.
double adjustedRandomScore(vector<vector<int>> & labels_true, vector<vector<int>> & labels_pred);

bool inVec(vector<int> & sol, vector<vector<int>> & ref_vec);

double kMeanObj(vector<vector<double>> & data, std::vector<std::vector<double>>& means, vector<uint32_t> labels);

double MaxMCC(vector<vector<int>> & cluster1, vector<vector<int>> & cluster2);

double MCC(vector<int> & labels_true, vector<int> & labels_pred);

std::vector<std::vector<int> > transpose(const std::vector<std::vector<int> > & data);

bool is_file_exist(const std::string& name);

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

float correlationDistance(const vector<double> & X,const vector<double> &Y, int n);

float euclideanDistance(const vector<double> & X,const vector<double> &Y, int n);

#endif
