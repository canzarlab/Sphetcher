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

template <typename T>
float distanceCorrelation(const vector<T> & X, const vector<T> &Y, int n) 
{ 

    float sum_X = 0.0, sum_Y = 0.0, sum_XY = 0.0; 
    float squareSum_X = 0.0, squareSum_Y = 0.0; 
    for (int i = 0; i < n; i++) 
    { 
        // sum of elements of array X. 
        sum_X = sum_X + X[i]; 

        // sum of elements of array Y. 
        sum_Y = sum_Y + Y[i]; 

        // sum of X[i] * Y[i]. 
        sum_XY = sum_XY + X[i] * Y[i]; 

        // sum of square of array elements. 
        squareSum_X = squareSum_X + X[i] * X[i]; 
        squareSum_Y = squareSum_Y + Y[i] * Y[i]; 
    } 
    // use formula for calculating correlation coefficient. 
    float corr = (n * sum_XY - sum_X * sum_Y) 
                / sqrt((n * squareSum_X - sum_X * sum_X) 
                    * (n * squareSum_Y - sum_Y * sum_Y)); 

    return 1.0 - corr; 
} 

// function returns the rank vector of the set of observations 
template<typename T>
vector<T> rankify(const vector<T> X){
    int N = X.size(); 
    vector<T> Rank_X(N); 
    for(int i = 0; i < N; i++)  
    { 
        int r = 1, s = 1; 
        // Count no of smaller elements 
        // in 0 to i-1 
        for(int j = 0; j < i; j++) { 
            if (X[j] < X[i] ) r++; 
            if (X[j] == X[i] ) s++; 
        } 
        // Count no of smaller elements 
        // in i+1 to N-1 
        for (int j = i+1; j < N; j++) { 
            if (X[j] < X[i] ) r++; 
            if (X[j] == X[i] ) s++; 
        } 
        // Use Fractional Rank formula 
        // fractional_rank = r + (n-1)/2 
        Rank_X[i] = r + (s-1) * 0.5;         
    } 
    // Return Rank Vector 
    return Rank_X; 
}

template<typename T> 
float distanceMetrics(const vector<T> & X, const vector<T> &Y, int n, std::string metrics)
{
    if (metrics == "manhattan"){
        float sum_X = 0.0;
        for (int i = 0; i < n; i++) 
        { 
            sum_X += abs(X[i]-Y[i]);
        } 
        return sum_X;
    }
    else if (metrics == "euclidean"){
        float sum_X = 0.0;
        for (int i = 0; i < n; i++) 
        { 
            sum_X += (X[i]-Y[i])*(X[i]-Y[i]);
        } 
        return sqrt(sum_X);
    }
    else if (metrics == "correlation"){
        float sum_X = 0.0, sum_Y = 0.0, sum_XY = 0.0; 
        float squareSum_X = 0.0, squareSum_Y = 0.0; 
        for (int i = 0; i < n; i++) 
        { 
            // sum of elements of array X. 
            sum_X = sum_X + X[i]; 

            // sum of elements of array Y. 
            sum_Y = sum_Y + Y[i]; 

            // sum of X[i] * Y[i]. 
            sum_XY = sum_XY + X[i] * Y[i]; 

            // sum of square of array elements. 
            squareSum_X = squareSum_X + X[i] * X[i]; 
            squareSum_Y = squareSum_Y + Y[i] * Y[i]; 
        } 
        // use formula for calculating correlation coefficient. 
        float corr = (n * sum_XY - sum_X * sum_Y) 
                    / sqrt((n * squareSum_X - sum_X * sum_X) 
                        * (n * squareSum_Y - sum_Y * sum_Y)); 

        return 1.0 - corr; 
    }
    else if (metrics == "spearman"){
        vector<T> vec_x = rankify(X);
        vector<T> vec_y = rankify(Y);
        return distanceCorrelation(vec_x, vec_y, n);
    }
    else if(metrics == "cosine"){
        float mul = 0.0;
        float d_a = 0.0;
        float d_b = 0.0 ;

        // Prevent Division by zero
        if (X.size() < 1)
        {
            throw std::logic_error("Vector X and Vector Y are empty");
        }

        typename std::vector<T>::const_iterator Y_iter = Y.begin();
        typename std::vector<T>::const_iterator X_iter = X.begin();
        for( ; X_iter != X.end(); X_iter++ , Y_iter++ )
        {
            mul += *X_iter * *Y_iter;
            d_a += *X_iter * *X_iter;
            d_b += *Y_iter * *Y_iter;
        }
        if (d_a == 0.0f || d_b == 0.0f)
        {
            throw std::logic_error(
                    "cosine similarity is not defined whenever one or both "
                    "input vectors are zero-vectors.");
        }

        return 1.0 - mul/(sqrt(d_a) * sqrt(d_b));

        }
        else{
            cout << "metrics is undefined\n";
            return -1000.0;
        }

}

template<typename T>
vector<vector<T>> transpose(vector<vector<T> > &b)
{
    vector<vector<T> > trans_vec(b[0].size(), vector<T>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    return trans_vec;    // <--- reassign here
}


#endif
