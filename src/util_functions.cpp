#include "util_functions.h"

void swap_vec(vector<int> & v1, vector<int> & v2)
{
    vector<int> temp = v2;
    v2 = v1;
    v1 = temp;
}

double getJaccard (vector<int> v1, vector<int> v2) // compute the jaccard of two sets of integers v1 and v2
{
    vector<int> v3; // v3 is the intersection of v1 and v2
    vector<int> v4; // v4 is the union of v1 and v2
    
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
    set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v4));
    return (double) v3.size()/(double) v4.size();
}

vector<vector<int>> parse2DCsvFile(string inputFileName) {
    vector<vector<int> > _data;
    ifstream inputFile(inputFileName);
    int l = 0; // 0
    int count =0;
 
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<int> record;
            int removefist = 0;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    if (removefist >= 0)    // matrix
                    record.push_back(stof(line));
                    else
                        removefist++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
            _data.push_back(record);
//             if (count > 0)
//                 data.push_back(record);
//             else    
//                 count ++;
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return _data;
}

vector<vector<string>> parse2DCsvFile2String(string inputFileName) {
    vector<vector<string> > data;
    ifstream inputFile(inputFileName);
    int l = 0; // 0
    int count =0;
 
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<string> record;
            int removefist = 0;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    if (removefist >= 0)    
                    record.push_back(line);
                    else
                        removefist++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
            data.push_back(record);
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}

vector<int> parse1DcsvFile2Int(string inputFileName, int from_line){
    vector<int> _data;
    ifstream inputFile(inputFileName);
    int l = 0; // 0
    int count =0;
 
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            int removefist = 0;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    if (removefist >= 0 && l >= from_line)    // matrix
                        _data.push_back(stoi(line));
                    else
                        removefist++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
        }
    }
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return _data;
}

vector<vector<double>> parse2DCsvFile2Double(string inputFileName) {
    vector<vector<double> > data;
    ifstream inputFile(inputFileName);
    int l = 0; // 0
    int count =0;
 
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> record;
            int removefist = 0;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    if (removefist >= 0)    
                    record.push_back(std::stod(line));
                    else
                        removefist++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
            data.push_back(record);
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}

void switchTwoElement(vector<int> & vec, int x, int y){
    for (auto & u: vec)
    {
        if (u == x)
        {
            u = y;
        } else if (u == y){
            u = x;
        } else 
        {}
    }
}

void uniqueClusterPresentation(vector<int> & vec)
{
    int min_vec = *std::min_element(vec.begin(),vec.end());
    if (min_vec != 0)
    {
        for (auto &u: vec)
        {
            u = u - min_vec;
        }
    }
    
    int n = (int) vec.size();
    switchTwoElement(vec, vec[0], 0);
    int current = 1;
    for (int i=1; i < n; ++ i)
    {
        int max_value = vec[0];
        for (int j=0; j < i; j++)
        {
            max_value = std::max(max_value, vec[j]);
        }
        
        if (vec[i] > max_value + 1)
        {
            switchTwoElement(vec, current, vec[i]);
            current++;
        }
    }
}

int cluster_index(vector<int> vec1, vector<int> vec2)
{
    assert(vec1.size() == vec2.size());
    int result = 0;
    for (int i=0; i < (int) vec1.size(); ++i)
    {
        result += std::abs(vec1.at(i)-vec2.at(i));
    }
    return result;
}

vector<int> bestMatch(vector<int>  & single_cluster, vector<vector<int>> & cluster, vector<bool> & open_option)
{
    int n_clusters = (int) cluster.size();
    int available_slots = 0;
    for (auto u: open_option)
    {
        if (u)  available_slots++;
    }
    assert(available_slots > 0);
    
    int max_value = -1;
    int max_index;
    for (int i=0; i < n_clusters; ++i){
        if (open_option.at(i))
        {
            std::vector<int> v_intersection;
            std::set_intersection(single_cluster.begin(), single_cluster.end(),
                                cluster[i].begin(), cluster[i].end(),
                                std::back_inserter(v_intersection));
            if ((int) v_intersection.size() > max_value){
                max_value = (int) v_intersection.size();
                max_index = i;
            }
        }
    }
    open_option.at(max_index) = false;
    cout << "MAX INDEX: " << max_index << endl;
    return cluster.at(max_index);
}

vector<vector<int>> findOptimalOrdering(vector<vector<int>> & ref_cluster, vector<vector<int>> & cluster)
{    
    int n_clusters = (int) cluster.size();
    for (int i=0; i < n_clusters; ++i)
    {
        sort(ref_cluster[i].begin(), ref_cluster[i].end());
        sort(cluster[i].begin(), cluster[i].end());
    }
    
    vector<bool> open_option(n_clusters, true);
    vector<vector<int>> results;
    for (int i=0; i < n_clusters; ++i)
    {
        results.push_back(bestMatch(ref_cluster.at(i), cluster, open_option));
    }
    
    return results;
}

vector<int> row_sum(vector<vector<int>> & matrix)
{
    vector<int> result;
    for (auto & u: matrix){
        int cumsum = 0;
        for (auto & v: u){
            cumsum += v;
        }
        result.push_back(cumsum);
    }
    return result;
}

vector<int> column_sum(vector<vector<int>> & matrix){
    int s = (int) matrix[0].size();
    vector<int> results(s, 0);
    for (auto & u: matrix){
        results = results + u;
    }
    return results;
}

double _comb2(int n)
{
    return (double) (n*(n-1))/2.0;
}

double adjustedRandomScore(vector<vector<int>> & labels_true, vector<vector<int>> & labels_pred)
{
    int n_classes  = (int) labels_true.size();
    int n_clusters = (int) labels_pred.size();
    int n_samples = 0;
    for (auto & u: labels_true)
        n_samples += (int) u.size();
    
    int n_samples2 = 0;
    for (auto & u: labels_pred)
        n_samples2 += (int) u.size();
    assert(n_samples == n_samples2);

//     # Special limit cases: no clustering since the data is not split;
//     # or trivial clustering where each document is assigned a unique cluster.
//     # These are perfect matches hence return 1.0.
    if (n_classes == n_clusters && n_clusters == 1 || n_classes == n_clusters && n_clusters == 0 || 
        n_classes == n_clusters && n_clusters == n_samples)
            return 1.0;

//     # Compute the ARI using the contingency data
    vector<vector<int>> contingency;
    for (int i=0; i < n_classes; ++i){
        vector<int> row;
        for (int j=0; j < n_clusters; ++j){
            row.push_back(lenIntersection(labels_true.at(i), labels_pred.at(j)));
        }
        contingency.push_back(row);
    }
    
    vector<int> a = row_sum(contingency);
    vector<int> b = column_sum(contingency);
    double sum_comb_a = 0.0;
    for (auto & v: a){
        sum_comb_a += _comb2(v);
    }
    double sum_comb_b = 0.0;
    for (auto & v: b){
        sum_comb_b += _comb2(v);
    }
    double sum_comb = 0.0;
    for (auto & u: contingency){
        for (auto & v: u){
            sum_comb += _comb2(v);
        }
    }
    
    double prod_comb = (sum_comb_a * sum_comb_b) / _comb2(n_samples);
    double mean_comb = (sum_comb_a + sum_comb_b) / 2.0;
    return (sum_comb - prod_comb) / (mean_comb - prod_comb);
}

bool inVec(vector<int> & sol, vector<vector<int>> & ref_vec){
    for (auto & u: ref_vec)
        if (sol == u)
            return true;
    return false;
}

double kMeanObj(vector<vector<double>> & data, std::vector<std::vector<double>>& means, vector<uint32_t> labels)
{
    double obj = 0.0;
    for (size_t i=0; i < labels.size(); ++i)
    {
        obj += distance(data.at(i), means[labels.at(i)]);
    }
    return obj;
}

double MaxMCC(vector<vector<int>> & cluster1, vector<vector<int>> & cluster2)
{
    int n = (int) cluster1[0].size() + (int) cluster1[1].size();
    
    vector<int> labels_true(n, 0);
    vector<int> labels_true2(n, 1);
    for (auto & u: cluster1.at(1)){
        labels_true.at(u) = 1;
        labels_true2.at(u) = 0;
    }
    
    vector<int> labels_pred(n, 0);
    for (auto & u: cluster2.at(1)){
        labels_pred.at(u) = 1;
    }
    return std::max(MCC(labels_true, labels_pred), MCC(labels_true2, labels_pred));
}

double MCC(vector<int> & labels_true, vector<int> & labels_pred)
{
    for (std::vector<int>::size_type i=0; i != labels_true.size(); ++i){
        assert(labels_true[i] >= 0 && labels_true[i] <= 1 && labels_pred[i] >= 0 && labels_pred[i] <= 1);
    }
    assert(labels_pred.size() == labels_true.size());
    double TP=0.0, TN =0.0, FP = 0.0, FN = 0.0;
    for (std::vector<int>::size_type i=0; i != labels_true.size(); ++i){
        if (labels_true[i]==1 && labels_pred[i]==1){
            TP += 1.0;
        }
        else if (labels_true[i]==1 && labels_pred[i]==0){
            FN += 1.0;
        }
        else if (labels_true[i]==0 && labels_pred[i]==1){
            FP += 1.0;
        }
        else {
            TN += 1.0;
        }
    }
    double temp = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    if (temp == 0.0)
        return 0.0;
    else
        return (TP*TN - FP*FN)/temp;
}

std::vector<std::vector<int> > transpose(const std::vector<std::vector<int>> & data) {
    // this assumes that all inner vectors have the same size and
    // allocates space for the complete result in advance
    std::vector<std::vector<int> > result(data[0].size(),
                                          std::vector<int>(data.size()));
    for (std::vector<int>::size_type i = 0; i < data[0].size(); i++) 
        for (std::vector<int>::size_type j = 0; j < data.size(); j++) {
            result[i][j] = data[j][i];
        }
    return result;
}

bool is_file_exist(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

float correlationDistance(const vector<double> & X,const vector<double> &Y, int n) 
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

float euclideanDistance(const vector<double> & X,const vector<double> &Y, int n) 
{ 

    float sum_X = 0.0;
    for (int i = 0; i < n; i++) 
    { 
        sum_X += (X[i]-Y[i])*(X[i]-Y[i]);
    } 
    return sqrt(sum_X);
} 
