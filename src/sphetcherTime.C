#define OPENMP 1
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "parallel.h"
#include "IO.h"
#include "graphIO.h"
#include "parseCommandLine.h"
#include "setCover.h"
#include <omp.h>
#include "util_functions.h"
#include "grid_cover.h"
#include <chrono>
using namespace std;
using namespace benchIO;

template <class T>
void create_threshold_data(vector<vector<T>> & Dist, T LL, int set_size, vector<int> & instance){
    int n;
    n=(int) set_size;
    int total_elem = 0;
    int first_elem_set_idx = 0;
    instance.clear(); // clear the vector 
    instance.push_back(n);
    instance.push_back(total_elem); // will change later 
    for (int i=0; i < n; ++i){
        instance.push_back(i);
    }

    for(int i = 0; i<n; i++){
        instance[i+2] = first_elem_set_idx;
        for(int j = 0; j < n; j++){
            if(Dist[i][j] <= LL)
            {
                instance.push_back(j);
                first_elem_set_idx++;
            }
        }
    }
    instance[1] = first_elem_set_idx;
}


template <class T>
void create_threshold_data(vector<vector<T>> & Dist, T LL, int set_size, vector<int> & instance, vector<int> labels, int n_splits, int n_classes = 300){
    // hit every labels, n_splits: split each class into n_splits 
    cout << " ----------------------------------- RUN HITTING CONDITION --------------------------------\n";
    cout << "size of labels: " << labels.size() << ", l_min: " << n_splits << endl;
    int min_label = *std::min_element(labels.begin(), labels.end());
    vector<int> increase_demand(n_classes*2, 0);
    for (size_t i=0; i < labels.size(); ++i){
        labels[i] -= min_label;
    }
    vector<int> origin_labels = labels;
    for (size_t i=0; i < labels.size(); ++i){
        labels[i] *= n_splits;
    }
    for (size_t i=0; i < labels.size(); ++i){
        labels[i] += set_size + increase_demand[origin_labels[i]];
        increase_demand[origin_labels[i]]++;
        increase_demand[origin_labels[i]] = increase_demand[origin_labels[i]] % n_splits;
    }
    // print_vec(labels);

    int n;
    n=(int) set_size;
    int total_elem = 0;
    int first_elem_set_idx = 0;
    instance.clear(); // clear the vector 
    instance.push_back(n);
    instance.push_back(total_elem); // will change later 
    for (int i=0; i < n; ++i){
        instance.push_back(i);
    }

    for(int i = 0; i<n; i++){
        instance[i+2] = first_elem_set_idx;
        for(int j = 0; j < n; j++){
            if(Dist[i][j] <= LL)
            {
                instance.push_back(j);
                first_elem_set_idx++;
            }
        }
        // add hitting set
        instance.push_back(labels[i]);
        first_elem_set_idx++;
    }
    instance[1] = first_elem_set_idx;
    // instance[0] = *std::max_element(labels.begin(), labels.end());
    // cout << instance[0] << "xxx";
}

void timesetCover(graph<intT> G, int rounds, char* outFile) {
  _seq<intT> Sets = _seq<intT>(NULL,0);
  for (int i=0; i < rounds; i++) {
    if (Sets.A != NULL) free(Sets.A);
    graph<intT> GG = G.copy();
    startTime();
    Sets = setCover(GG);
    nextTimeN();
    GG.del();
  }
  cout << endl;

  if (outFile != NULL)
    writeIntArrayToFile(Sets.A, Sets.n, outFile);

  Sets.del();
  G.del();
}

template <class T>
vector<int> run_setcover_binarysearch(vector<vector<T>> & Dist, T d_min, T d_max, int n_points, int n_iters, string dataname, int n_sets, vector<int> & v, int n_splits=-1)
{
    // create the instance 
    vector<int> results;
    for (int iter=0; iter < n_iters; ++iter)
    {
        T L = (d_min + d_max)/2.0;
        vector<int> instance;
        if (n_splits == -1){
            create_threshold_data(Dist, L, n_sets, instance);
        }else{
            create_threshold_data(Dist, L, n_sets, instance, v, n_splits);
        }
        // run setcover
        int rounds = 1;
        graph<intT> G = readGraphFromVector<intT>(instance);
        string outputFile = dataname+"_setcover_sol__.csv";
        char * oFile = &outputFile[0];
        timesetCover(G, rounds, oFile);
        vector<int> set_cover =  parse1DcsvFile2Int(outputFile, 1);
        std::remove(oFile);
        int cover_size = (int) set_cover.size();
        cout << "iter: " << iter << ", sketch size: " << n_points << ", set_cover size: " << cover_size << "\n\n";
        if (cover_size < n_points){
            d_max = L;
        }
        else{
            d_min = L;
        }

        if(iter == n_iters - 1){ results =  set_cover; }
    }
    return results;
}

vector<int> convert2original_idx(const vector<int> & shuffle_idx, const vector<int> & result){
    vector<int> origin;
    for (auto v: result){
        origin.push_back(shuffle_idx[v]);
    }
    return origin;
}

vector<int> fill_in(vector<int> & solution, int n_sets, int n_points){
    int sol_size = (int) solution.size();
    if (n_points > sol_size){
       for (int i=0; i < n_points-sol_size-1; ++i){
          int key = rand() % (n_sets-1);
          if (std::find(solution.begin(), solution.end(), key) != solution.end())
              solution.push_back(key);
       }
      return solution; 
    }
    else
    {
        vector<int> newsolution(solution.begin(), solution.begin()+n_points);
        return newsolution;
    }
}

int parallel_main(int argc, char* argv[]) {
  // commandLine P(argc, argv, "[-o <outFile>] [-r <rounds>] <inFile>");
  // char* iFile = P.getArgument(0);
  // char* oFile = P.getOptionValue("-o");
  // int rounds = P.getOptionIntValue("-r",1);

    std::string dataname = "dataset";
    std::string dfilename, lfilename;
    vector<int> origin_labels;
    int n_splits = -1;
    // std::string dfilename = "/data/hoan/spectral_clustering_matlab/data/"+dataname+"_pca.csv";
    // std::string lfilename = "/data/hoan/spectral_clustering_matlab/data/"+dataname+"_pca_labels.csv";
    if (argc < 3){
        cout << "missing arguments \n";
        cout << "./sphetcher expression_matrix.csv sketch_size \n";
        cout << "./sphetcher expression_matrix.csv sketch_size class_labels.csv l_min\n";
        assert(argc >= 3);
    }
    dfilename = argv[1];
    int sketch_size = stoi(argv[2]);
    if (argc == 4){
        lfilename = argv[3];
        n_splits = 1;
        cout << "sample at least one element from each class \n";
    }
    if (argc == 5){
        lfilename = argv[3];
        n_splits = stoi(argv[4]);
        cout << "sample at least " << n_splits << " elements from each class\n";
    }
    if (n_splits != -1)
    {
        vector<vector<int>> labels_vec =  parse2DCsvFile(lfilename);
        if (labels_vec.size() == 1){
            origin_labels = labels_vec.at(0);
        }else{
            for (auto u: labels_vec){
                origin_labels.push_back(u[0]);
            }
        }
    }
    // TODO: need to remove the labels 
    // std::string dfilename = "../../../data/"+dataname+"-prepare-log_count_pca.csv";
    // if (!is_file_exist(dfilename)){
    //     dfilename = "../../../data/"+dataname+"_pca.csv";
    // }
    // if (!is_file_exist(dfilename)){
    //     dfilename = "../../../data/"+dataname+"-prepare-log_count_pca2000.csv";
    // }
    auto start = chrono::steady_clock::now();
    vector<vector<double>> distance_mat_origin = parse2DCsvFile2Double(dfilename);
    int n_sets = (int) distance_mat_origin.size();
    // shuffle data 
    vector<int> shuffle_idx;
    for (int i=0; i < n_sets; ++i){
        shuffle_idx.push_back(i);
    }
    // srand(time(0));
    vector<vector<int>> sc_solutions;
    vector<vector<int>> box_vec; 
    for(int rep=0; rep < 1; ++rep)
    {
        std::random_shuffle(shuffle_idx.begin(), shuffle_idx.end());
        vector<vector<double>> distance_mat;
        for (size_t i=0; i < distance_mat_origin.size(); ++i){
            distance_mat.push_back(distance_mat_origin[shuffle_idx[i]]);
        }
        vector<int> labels;
        if (n_splits != -1){
            for (size_t i =0; i < origin_labels.size(); ++i){
                labels.push_back(origin_labels[shuffle_idx[i]]);
            }
        }
        auto end = chrono::steady_clock::now();
        cout << "read file in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
        int n_dim =  (int) distance_mat[0].size();
        cout << "# of samples: " << n_sets << ", n_dim: " << n_dim << endl;

        int n_boxes = std::min(std::max(80000,(int) round(0.2*n_sets)), 100000);
        vector<int> n_points_vec;
        n_points_vec.push_back(sketch_size);
        if ( n_sets < n_boxes)
        {
            start = chrono::steady_clock::now();
            vector<vector<float>> Dist(n_sets, vector<float>(n_sets,0.0));
            int i,j; 
            float d_max = 0.0;
            #pragma omp parallel for private(i,j)
            for (i=0; i < n_sets; ++i){
                for (j=0; j<n_sets; ++j){
                    Dist[i][j] = euclideanDistance(distance_mat.at(i), distance_mat.at(j), n_dim);
                    if (Dist[i][j] > d_max){
                        d_max = Dist[i][j];
                    }
                }
            }
            end = chrono::steady_clock::now();
            cout << "compute distance matrix in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;

            start = chrono::steady_clock::now();
            d_max = d_max/2.0; float d_min = d_max/80.0; int n_iters = 20;
            // int n_points = min(10000, n_sets);
            // vector<int> setcover = run_setcover_binarysearch(Dist, d_min, d_max, n_points, n_iters, dataname, n_sets);
            // vector<int> fill_setcover = fill_in(setcover, n_sets, n_points);
            // sc_solutions.push_back(fill_setcover);

            // run for multiple n_points
            for (auto n_points: n_points_vec){
                cout << "--------------------------------------------------------- n_points: " << n_points << " -------------------------------------------\n";
                cout << "--------------------------------------------------------- l_min: " << n_splits << " -------------------------------------------\n";
                vector<int> setcover = run_setcover_binarysearch(Dist, d_min, d_max, n_points, n_iters, dataname, n_sets, labels, n_splits);
                vector<int> fill_setcover = fill_in(setcover, n_sets, n_points);
                // vector<int> fill_setcover = setcover;
                // convert to original 
                vector<int> origin_fill_setcover = convert2original_idx(shuffle_idx, fill_setcover);
                vector<int> sc_indicator(n_sets, 0);
                for(auto fi: origin_fill_setcover){
                    sc_indicator[fi] = 1;
                }
                sc_solutions.push_back(sc_indicator);
            }
            end = chrono::steady_clock::now();
            cout << "setcover in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
        }
        else 
        {
            // int dim_x = 10; double N = 17.0; 
            // for (auto L_min : n_points_vec){
            //     auto start2 = chrono::steady_clock::now();
            //     vector<int> box_idx = run_box_binarysearch(distance_mat, dim_x, N, L_min, L_min + 20);
            //     auto end2 = chrono::steady_clock::now();
            //     cout << "run box algorithm in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;
            //     vector<int> box_indicator(n_sets, 0);
            //     for (auto idx: box_idx){
            //         box_indicator.at(idx) = 1;
            //     }
            //     box_vec.push_back(box_indicator);
            // }
            // run grid box algorithm 
            auto start2 = chrono::steady_clock::now();
            int dim_x = std::min(10, n_dim); double N = 17.0; int L_box = n_boxes;
            int L_min = L_box-20; int L_max = min(n_sets, L_box+20);
            vector<int> box_idx = run_box_binarysearch(distance_mat, dim_x, N, L_min, L_max);
            auto end2 = chrono::steady_clock::now();
            cout << "run box algorithm in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;
            vector<int> box_indicator(n_sets, 0);
            for (auto idx: box_idx){
                box_indicator.at(idx) = 1;
            }
            box_vec.push_back(box_indicator);

            // compute distance matrix
            start = chrono::steady_clock::now();
            n_boxes = (int) box_idx.size();
            vector<vector<float>> Dist(n_boxes, vector<float>(n_boxes,0.0));
            int i,j; float d_max = 0.0;
            #pragma omp parallel for private(i,j)
            for (i=0; i < n_boxes; ++i){
                for (j=0; j< n_boxes; ++j){
                    Dist[i][j] = correlationDistance(distance_mat.at(box_idx[i]), distance_mat.at(box_idx[j]), n_dim);
                    if (Dist[i][j] > d_max){
                        d_max = Dist[i][j];
                    }
                }
            }
            end = chrono::steady_clock::now();
            cout << "compute distance matrix in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;

            // compute set cover 
            d_max = d_max/4.0; float d_min = d_max/100.0; int n_iters = 14;
            for (auto n_points: n_points_vec)
            {
                start = chrono::steady_clock::now();
                assert(n_points <= n_boxes);
                cout << "--------------------- sketch size: " << n_points << "\n";
                vector<int> setcover = run_setcover_binarysearch(Dist, d_min, d_max, n_points, n_iters, dataname, n_boxes, labels, n_splits);
                vector<int> convert_boxid2realid;
                for(size_t k=0; k < setcover.size(); ++k){
                    convert_boxid2realid.push_back(box_idx[setcover[k]]);
                }
                vector<int> fill_setcover = fill_in(convert_boxid2realid, n_sets, n_points);
                // vector<int> fill_setcover = convert_boxid2realid;
                // convert to original 
                vector<int> origin_fill_setcover = convert2original_idx(shuffle_idx, fill_setcover);
                vector<int> sc_indicator(n_sets, 0);
                for(auto fi: origin_fill_setcover){
                    sc_indicator[fi] = 1;
                }
                sc_solutions.push_back(sc_indicator);
                end = chrono::steady_clock::now();
                cout << "setcover in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
            }

        }
    }
    //write solution to file
    // writeVec2File(dataname + "_box_solutions.csv", box_vec);
    string fileName = "indicator_solutions.csv";
    writeVec2File(fileName, sc_solutions);
}
