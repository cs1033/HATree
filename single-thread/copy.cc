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
 
#include "btree.h" 
#include "robtree.h"
#include "wbtree.h"
#include "fptree.h"
#include "fastfair.h"
#include "hatree1.h"
#include "hatree2.h"
#include "hatree3.h"
#include "hatree4.h"
#include "rutree.h"
//#include "../../cd "

using std::cout;
using std::endl;
using std::ifstream;
using std::string;

PMAllocator * galc;
uint32_t flush_cnt;
_key_t * keys;
typedef int64_t mykey_t; 

template <typename BTreeType>
void preload(BTreeType *tree, uint64_t load_size, ifstream & fin) {
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
double run_test(BTreeType *tree, ifstream & fin) {
    int small_noise = getRandom() % 99; // each time we run, we will insert different keys
    
    uint64_t sum = 0;
    double ret = 0;
    double time = 0; 
    
    //BTreeType tree(path, true);
     
    int op_id, notfound = 0;
    _key_t key; 
    _value_t val;
    auto start = seconds();
    for(int i = 0; fin >> op_id >> key; i++) {
         
        switch (op_id) { 
            case OperationType::INSERT: {
                //cout << key + seed << endl;
                tree->insert(key + small_noise, key + small_noise);
                break;
            }
            case OperationType::READ: {
                //cout << key << endl;
                auto r = tree->find(key,  val);
                assert(r);
                break;
            }
            case OperationType::UPDATE: 
                tree->update(key, key * 2);
                break;  
            case OperationType::DELETE:
                tree->remove(key);
                break;
            default: 
                cout << "wrong operation id" << endl;
                break;
        }
        
        //time[op_id] += double(end - start);
        if (key != val) notfound++;
    }
    auto end = seconds(); 
    if(notfound > 0) cout << notfound << " keys not found" << sum << endl; // optimizer killer
    cout << double(end - start) << endl;
    return notfound; 
}

int main(int argc, char ** argv) {
    int opt_testid = 1;
    string opt_fname = "workload.txt";
    string opt_preRead = "";
    string dataset = "dataset.dat";
    int size = LOADSCALE;

    if(argc > 1 && atoi(argv[1]) > 0) {
        opt_testid = atoi(argv[1]); 
    }
    if (argc > 2) {
        size = atoi(argv[2]);
    }
    if(argc > 3) { 
        if(file_exist(argv[3]))
            opt_fname = argv[3];
        else
            cout << "workload file "<< argv[3] << " not exist" << endl;
    }
    if (argc > 4) {
        if (file_exist(argv[4]))
            opt_preRead = argv[4];
    }


    if (!file_exist(dataset.c_str())) {
        cout << "dataset doesnt exist" << endl;
        exit(0);
    }



    keys = new mykey_t[sizeof(_key_t) * MILLION];
    ifstream pre(dataset.c_str(), std::ios::binary);
    ifstream fin(opt_fname.c_str());
    double time = 0;
    switch (opt_testid) { 
    case 2: {
        //cout << "btree" << endl;
        btree::btree* tree = new btree::btree("/mnt/pmem/lgc/btree.pool", false);
        auto start = seconds();
        preload(tree, size, pre); 
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        //tree->statistics();
        time = run_test(tree, fin); 
        delete tree; 
        break; 
    }    
    case 3: {
        //cout << "robtree" << endl;
        robtree::robtree* tree = new robtree::robtree("/mnt/pmem/lgc/robtree.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        if (opt_preRead != "") {
            ifstream read(opt_preRead, std::ios::binary);
            preRead(tree, read);
        }
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break; 
    }
    case 4: {
        //cout << "rutree" << endl;
        rutree::rutree* tree = new rutree::rutree("/mnt/pmem/lgc/rutree.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        if (opt_preRead != "") {
            ifstream read(opt_preRead, std::ios::binary);
            preRead(tree, read);
        }
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break; 
    }
    case 11: {
        //cout << "ratree" << endl;
        hatree1::hatree1* tree = new hatree1::hatree1("/mnt/pmem/lgc/hatree1.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break;
    }
    case 12: {
        //cout << "ratree" << endl;
        hatree2::hatree2* tree = new hatree2::hatree2("/mnt/pmem/lgc/hatree2.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break;
    }
    case 13: {
        //cout << "ratree" << endl;
        hatree3::hatree3* tree = new hatree3::hatree3("/mnt/pmem/lgc/hatree3.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break;
    }
    case 14: {
        //cout << "ratree" << endl;
        hatree4::hatree4* tree = new hatree4::hatree4("/mnt/pmem/lgc/hatree4.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break;
    }
    case 1: {
        //cout << "fptree" << endl;
        fptree::fptree* tree = new fptree::fptree("/mnt/pmem/lgc/fptree.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        auto end = seconds();
        //cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break;
    }
    case 5: {
        cout << "fastfair" << endl;
        fast::btree * tree = new fast::btree("/mnt/pmem/fasttree.pool", false);
        auto start = seconds();
        preload(tree, LOADSCALE, pre);
        auto end = seconds();
        cout << "preload time:" << double(end - start) << endl;
        time = run_test(tree, fin);
        delete tree;
        break;
    }
    default:
        cout << "Invalid tree type" << endl;
        break;
    } 
    
    //cout << "test time:" << time << endl;
    delete keys;
    pre.close();
    fin.close();
    return time;
}