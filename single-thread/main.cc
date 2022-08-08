/*  
    Copyright(c) 2020 Luo Yongping. THIS SOFTWARE COMES WITH NO WARRANTIES, 
    USE AT YOUR OWN RISK!
*/
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <random>
#include <algorithm>  
#include <omp.h>
#include <atomic>
#include <cstdio>
#include <omp.h>
 
#include "common.h"
#include "btree.h" 
#include "robtree.h"
#include "fptree.h"
#include "fastfair.h"
#include "rutree.h"
// #include "rutree.h"



using std::cout;
using std::endl;
using std::ifstream;
using std::string;

PMAllocator * galc;
uint32_t flush_cnt;
_key_t * keys;
typedef int64_t mykey_t; 

template <typename BTreeType>
void preload(BTreeType *tree, uint64_t load_size, ifstream & fin, int thread_cnt = 1) {
    #ifdef DEBUG
        fin.read((char *)keys, sizeof(mykey_t) * MILLION);
        for(int i = 0; i < LOADSCALE * KILO; i++) {
            mykey_t key = keys[i];
            //cout << key << endl;
            tree->insert((mykey_t)key, key);
        }
        //tree.printAll();
    #else 
        for(uint64_t t = 0; t < load_size; t++) {
            fin.read((char *)keys, sizeof(mykey_t) * MILLION);
 
            for(int i = 0; i < MILLION; i++) {
                mykey_t key = keys[i];
                tree->insert((mykey_t)key, key);
                
            }
        }
    #endif

    return ;
}

template <typename BTreeType>
void preRead(BTreeType *tree, ifstream & fin) {
    int op_id, notfound = 0;
    _key_t key; 
    _value_t val;
    for(int i = 0; fin >> op_id >> key; i++) {
        if(!tree->find(key, val) || !val) notfound++; // optimizer killer
    }
    return ;
}

template<typename BTreeType>
double run_test(BTreeType *tree, std::vector<QueryType> querys, int thread_cnt) {
    int small_noise = getRandom() % 99; // each time we run, we will insert different keys   

    // set the timer
    auto start = seconds();
    int notfount = 0;
    // #pragma omp barrier

    // start the section of parallel 
    // #pragma omp parallel num_threads(thread_cnt) 
    {   
        // if (omp_get_thread_num() == 0) {
            for(size_t i = 0; i < querys.size(); i++) {
                // get a unique query from querys
                OperationType op = querys[i].op;
                _key_t key = querys[i].key;
                _value_t val;

                switch (op) { 
                    case OperationType::INSERT: {
                        tree->insert(key + small_noise, key + small_noise);
                        break;
                    }
                    case OperationType::READ: {
                        auto r = tree->find(key, val);
                        // if (key != val) cout << obtain_pos << "key=" << key << "value=" << val << endl;
                        if (key != val)     notfount++;     //kill optimize
                        break;
                    } 
                    case OperationType::UPDATE: {
                        tree->update(key, key);
                        break;  
                    }
                    case OperationType::DELETE: {
                        tree->remove(key);
                        break;
                    }
                    default: 
                        cout << "wrong operation id" << endl;
                        break;
                }
            }
        // }
    }
    // #pragma omp barrier
    auto end = seconds();
    cout << "notfound:" << notfount << endl;

    return end - start; 
}

int main(int argc, char ** argv) {
    int opt_testid = 1;
    string opt_fname = "workload.txt";
    string opt_preRead = "";
    string dataset = "dataset.dat";
    int size = LOADSCALE; 
    int thread_cnt = 1;

    static const char* optstr = "w:c:t:p:s:hH";
    opterr = 0;
    char opt;
    while ((opt = getopt(argc, argv, optstr)) != -1) {
        switch (opt) {
            case 'w': 
                opt_fname = optarg;
                break;
            case 'c':
                opt_testid = atoi(optarg);
                break;
            case 's':
                size = atoi(optarg);
                break;
            case 't':
                thread_cnt = atoi(optarg);
                break;
            case 'p':
                opt_preRead = optarg;
                break;
            case 'h':
            case 'H':
            default:
                cout << "w: workload filename" << endl
                    << "t: thread cnt" << endl
                    << "c: test tree--  1:fastfair  2:basetree -- 3:hatree  4:rutree" << endl
                    << "s: dataset size" << endl
                    << "p: preRead filename" << endl;
                exit(-1);
                break;
        }
    }



    if (!file_exist(dataset.c_str())) {
        cout << "dataset doesnt exist" << endl;
        exit(0);
    }


    keys = new mykey_t[sizeof(_key_t) * MILLION];
    ifstream pre(dataset.c_str(), std::ios::binary);
    ifstream fin(opt_fname.c_str());
    double time = 0;

    std::vector<QueryType> querys;
    int op;
    _key_t key;
    if(!fin) {
        cout << "workload file not openned" << endl;
        exit(-1);
    }
    while(fin >> op >> key) {
        querys.push_back({(OperationType)op, key});
    }

    switch (opt_testid) { 
        case 1: {
            fast::btree* tree = new fast::btree("/mnt/pmem1/lgc/fastfair.pool", false);
            auto start = seconds();
            preload(tree, size, pre, 1); 
            auto end = seconds();
            // cout << "preload time:" << double(end - start) << endl;
            time = run_test(tree, querys, thread_cnt); 
            cout << "fast" ;
            delete tree; 
            break; 
        }
        case 2: {
            //cout << "btree" << endl;
            btree::btree* tree = new btree::btree("/mnt/pmem/lgc/basetree.pool", false);
            auto start = seconds();
            preload(tree, size, pre, 1); 
            auto end = seconds();
            // cout << "preload time:" << double(end - start) << endl;
            tree->statistics();
            time = run_test(tree, querys, thread_cnt);  
            delete tree; 
            cout << "basetree" ;
            break; 
        }   
        case 3: {
            robtree::robtree* tree = new robtree::robtree("/mnt/pmem1/lgc/robtree.pool", false);
            auto start = seconds();
            preload(tree, size, pre, 1); 
            auto end = seconds();
            //cout << "preload time:" << double(end - start) << endl;
            time = run_test(tree, querys, thread_cnt); 
            delete tree; 
            cout << "hatree";
            break; 
        }    
        case 4: {
            rutree::rutree* tree = new rutree::rutree("/mnt/pmem1/lgc/rutree.pool", false);
            auto start = seconds();
            preload(tree, size, pre, 1); 
            auto end = seconds();
            time = run_test(tree, querys, thread_cnt); 
            // #pragma omp barrier
            // #pragma omp parallel num_threads(2)
            // {
            //     if (omp_get_thread_num() == 0) {
            //         time = run_test(tree, querys, thread_cnt); 
            //     } else {
            //         // tree->flush();
            //     }
            // }
            // #pragma omp barrier
            delete tree; 
            cout << "rutree";
            break; 
        }
        default:
            cout << "Invalid tree type" << endl;
            break;
    } 
    cout << "'s Throughput : " << querys.size() / time / 1000000 << endl;
    delete keys;
    pre.close();
    fin.close();
    return time;
}