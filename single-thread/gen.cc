/*  
    Copyright(c) 2020 Luo Yongping. THIS SOFTWARE COMES WITH NO WARRANTIES, 
    USE AT YOUR OWN RISK!
*/
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cstdlib>
#include <unistd.h>
#include <algorithm> 
 
#include "common.h"
#include "zipfian_int_distribution.h" 

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::ifstream;
 
enum DistributionType {RAND = 0, ZIPFIAN};

struct WorkloadType {
    int operations = KILO;
    float read = 1.0;
    float insert = 0;
    float update = 0;
    float remove = 0;
    DistributionType dist = RAND;
    float skewness = 0.8;
    bool valid() {
        return read + insert + update + remove == 1.0 && skewness > 0 && skewness < 1.0;
    }
    void print() {
        cout << "=========WORKLOAD TYPE=========" << endl; 
        cout << "Operations  : " << operations << endl;
        cout << "Read Ratio  : " << read << endl;
        cout << "Insert Ratio: " << insert << endl;
        cout << "Update Ratio: " << update << endl;
        cout << "Remove Ratio: " << remove << endl;
        cout << "Distribution: " << (dist == RAND ? "random" : "Zipfian") << endl;
        if (dist == ZIPFIAN) {
            cout << "Skewness " << skewness << endl;
        }
        cout << "===============================" << endl;
    }
};


class OperationGenerator {
public :
    OperationType mappings_[100];
    std::default_random_engine gen_;
    std::uniform_int_distribution<uint32_t> dist_;

    OperationGenerator(WorkloadType &w) : gen_(getRandom()) {
        int read_end = 100 * w.read;
        int insert_end = read_end + 100 * w.insert;
        int update_end = insert_end + 100 * w.update;
        int remove_end = update_end + 100 * w.remove;

        for(int i = 0; i < read_end; i++)
            mappings_[i] = OperationType::READ;

        for(int i = read_end; i < insert_end; i++)
            mappings_[i] = OperationType::INSERT;

        for(int i = insert_end; i < update_end; i++)
            mappings_[i] = OperationType::UPDATE;

        for(int i = update_end; i < remove_end; i++)
            mappings_[i] = OperationType::DELETE;
    }

    OperationType next() {
        return mappings_[dist_(gen_) % 100];
    }
};

// try not to generate duplicate keys, because some tree indices may not get used to it

void gen_dataset(int64_t *arr, int64_t scale, bool random) {
    std::mt19937 gen(10007);
    if(random) {
        uint64_t step = INT64_MAX / scale;
        std::uniform_int_distribution<int64_t> dist(0, INT64_MAX);
        for(int64_t i = 0; i < scale; i++) {
            #ifdef DEBUG
                arr[i] = i + 1; 
            #else  
                arr[i] = i * step + 1; // non duplicated keys, and no zero key 
            #endif
        }
        std::shuffle(arr, arr + scale - 1, gen);
    } else {
        std::normal_distribution<double> dist(INT64_MAX / 2, INT64_MAX / 8);
        int64_t i = 0;
        while (i < scale) {
            double val = (int64_t)dist(gen);
            if(val < 0 || val > (double) INT64_MAX) {
                continue;
            } else {
                arr[i++] = (int64_t)std::round(val);  
            }
        }
    }    
    return ;    
}   

void gen_workload(int64_t *arr, int64_t scale, QueryType * querys, WorkloadType w) {
    std::mt19937 gen(getRandom());
    std::uniform_int_distribution<int64_t> idx1_dist(0, scale);
    zipfian_int_distribution<int64_t> idx2_dist(0, scale, w.skewness);
    OperationGenerator op_gen(w);

    for(int i = 0; i < w.operations; i++) {
        OperationType op = op_gen.next();
        int idx = (w.dist == RAND ? idx1_dist(gen) : idx2_dist(gen));
        // for insert operations, we should make sure the key does not exist in the dataset
        int64_t key = (op == OperationType::INSERT ? arr[idx] + getRandom(): arr[idx]); 

        querys[i] = {op, key};
    }
}

void gen_workload_hotspot(int64_t *arr, int64_t scale, QueryType * querys, WorkloadType w, int ocur) {
    cout << "hotspot" << endl;
    std::mt19937 gen(getRandom());
    int64_t time = scale / 100;
    std::uniform_int_distribution<int64_t> idx1_dist(0, time*99);
    zipfian_int_distribution<int64_t> idx2_dist(0, time, w.skewness);
    zipfian_int_distribution<int64_t> idx3_dist(0, scale, w.skewness);
    OperationGenerator op_gen(w);
 
    for(int i = 0; i < w.operations; i++)  {
        int base = 0; 
        if (i % (time * 5) == 0)  
            base = idx1_dist(gen);  
        
        OperationType op = op_gen.next();
        int choose = idx1_dist(gen) % 10 ;
        if (choose < ocur)  {
            choose = ZIPFIAN;
        } else {
            choose = RAND;
        }
        int idx = (choose == RAND ? idx3_dist(gen) : base + idx2_dist(gen));
        // for insert operations, we should make sure the key does not exist in the dataset
        int64_t key = (op == OperationType::INSERT ? arr[idx] + getRandom(): arr[idx]); 

        querys[i] = {op, key};
        
    }
} 
 
// void gen_workload(int64_t *arr, int64_t scale, QueryType * querys, WorkloadType w) {
//     std::mt19937 gen(getRandom());
//     std::uniform_int_distribution<int64_t> idx1_dist(0, scale*9/10);
//     zipfian_int_distribution<int64_t> idx2_dist(0, scale/10, w.skewness);
//     OperationGenerator op_gen(w);

//     int cycle[1001];

//     for (int i = 0; i < 1001; ++i) {
//         cycle[i] = idx1_dist(gen);
//     }

//     for(int i = 0; i * MILLION < w.operations; i++) {
//         int base = idx1_dist(gen);
//         for (int j = 0; j < MILLION; ++j) {
//              OperationType op = op_gen.next();
//             int idx = cycle[idx1_dist(gen) % 1001];
//             //int idx = (w.dist == RAND ? idx1_dist(gen) : base + idx2_dist(gen));
//             // for insert operations, we should make sure the key does not exist in the dataset
//             int64_t key = (op == OperationType::INSERT ? arr[idx] + getRandom(): arr[idx]); 
//             querys[i*MILLION + j] = {op, key};
//         }
//     }
// }

int main(int argc, char ** argv) { 
    static const bool DATASET_RANDOM =  true;
    bool opt_zipfian  = false;
    WorkloadType w;
    std::string workloadName = "workload.txt";
    int ocur = 0;
    int datasize = LOADSCALE;

    static const char * optstr = "L:l:p:P:n:r:i:u:d:o:s:hz"; 
    opterr = 0;
    char opt;
    while((opt = getopt(argc, argv, optstr)) != -1) {
        switch(opt) {
        case 'l':
        case 'L':
            datasize = atoi(optarg);
            break;
        case 'p':
        case 'P':
            ocur = atoi(optarg);
            break;
        case 'n':
            workloadName = optarg;
            break;
        case 'o':
            w.operations = atoi(optarg) * 1000 * 1000;
            break;
        case 's':
            w.skewness = atof(optarg);
            break;
        case 'r':
            w.read = atof(optarg);
            break;
        case 'i':
            w.insert = atof(optarg);
            break;
        case 'd':
            w.remove = atof(optarg);
            break;
        case 'u':
            w.update = atof(optarg);
            break;
        case 'z':
            opt_zipfian = true;
            break;
        case '?':
        case 'h': 
        default:
            cout << "USAGE: "<< argv[0] << "[option]" << endl;
            cout << "\t -h: " << "Print the USAGE" << endl;
            cout << "\t -z: " << "Use zipfian distribution (Not specified: random distribution)" << endl;
            cout << "\t -o: " << "The number of operations" << endl;
            cout << "\t -s: " << "The skewness of query workload(0 - 1)" << endl;
            cout << "\t -r: " << "Read ratio" << endl;
            cout << "\t -i: " << "Insert ratio" << endl;
            cout << "\t -u: " << "update ratio" << endl;
            cout << "\t -d: " << "Delete ratio" << endl;
            exit(-1);
        }
    }

    if(!w.valid()) {
        cout << "Invalid workload configuration" << endl;
        exit(-1);
    }
    if(opt_zipfian == true) {
        w.dist = ZIPFIAN;
    }

    w.print();

    #ifdef DEBUG
        uint64_t scale = LOADSCALE * KILO;
    #else
        uint64_t scale = datasize * MILLION;
    #endif

    int64_t * arr = new int64_t[scale];
    if(!file_exist("dataset.dat")) {
        std::mt19937 gen(10007);
        //ifstream fin("sort.dat", std::ios::binary);
        //fin.read((char *)arr, sizeof(int64_t) * scale); 
        gen_dataset(arr, scale, DATASET_RANDOM); 
        //std::shuffle(arr, arr + scale - 1, gen);
        ofstream fout("dataset.dat", std::ios::binary);
        fout.write((char *) arr, sizeof(int64_t) * scale);
        //fin.close();
        fout.close(); 
        cout << "generate a dataset file" << endl;
    } else {
        ifstream fin("dataset.dat", std::ios::binary);
        fin.read((char *)arr, sizeof(int64_t) * scale);
        fin.close();
    }
 
    QueryType * querys = new QueryType[w.operations];
    if (ocur != 0) {
        gen_workload_hotspot(arr, scale, querys, w, ocur);
    } else {
        gen_workload(arr, scale, querys, w);
    }
    
    ofstream fout(workloadName.c_str());
    for(int i = 0; i < w.operations; i++) {
        fout << querys[i].op << " " << querys[i].key << endl;
    }
    fout.close();
    cout << "generate a query workload file" << endl;

    delete [] arr;
    delete [] querys;
    return 0;
}