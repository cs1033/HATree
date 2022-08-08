/*
    Copyright (c) Luo Yongping. THIS SOFTWARE COMES WITH NO WARRANTIES, 
    USE AT YOUR OWN RISK!
*/
#ifndef __COMMON_H__
#define __COMMON_H__

#include <sys/stat.h>
#include <cstdint>
#include <unistd.h>
#include <sys/time.h>
 
// in-persistent-memory log area
#define LOG_DATA_SIZE        48
#define LOG_AREA_SIZE        4194304
#define LE_DATA              0
#define LE_COMMIT            1

#define POOLSIZE 10 
#define LOADSCALE 100 

#define KILO 1024
#define MILLION (KILO * KILO)
#define CACHE_LINE_SIZE 64

#define DOFLUSH

using _key_t = int64_t;
using _value_t = int64_t;

struct Record {
    _key_t key;
    char * val; 
    Record(_key_t k=INT64_MAX, char * v=NULL) : key(k), val(v) {}
    bool operator<(double b) {
        return (double) key < b;
    }
};

template <typename T>
static inline bool get(T state, size_t index) {
    return state & ((T)1 << index);
}

template <typename T>
static inline void set(T& state, size_t index) {
    state |= (T)1 << index;
}

template <typename T>
static inline void unset(T& state, size_t index) {
    state &= ~((T)1 << index);
}


enum OperationType {READ = 0, INSERT, UPDATE, DELETE};

struct QueryType {
    OperationType op;
    int64_t key;
};

struct res_t { // a result type use to pass info when split and search
    bool flag; 
    Record rec;
    res_t(bool f, Record e):flag(f), rec(e) {}
};

static double seconds()
{
  timeval now;
  gettimeofday(&now, NULL);
  return now.tv_sec + now.tv_usec/1000000.0;
}

static int getRandom() {
    timeval now;
    gettimeofday(&now, NULL);
    return now.tv_usec;
}

bool file_exist(const char *pool_path) {
    struct stat buffer;
    return (stat(pool_path, &buffer) == 0);
}

static char finger_print(_key_t k) { 
// using fmix function from murmur hash
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return (k & 0xff);
}


static char finger_print2(_key_t k) {
// using for the buffer fp
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return (k & 0xff);
}


#endif //__COMMON_H__