#ifndef __PRADIX_H__
#define __PRADIX_H__

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include <string>

#include "common.h"
#include "flush.h"

#define GAPPED

namespace lptree {

constexpr int INNER_RADIX        = 14;                 // max radix value 
constexpr int MAX_KID_NUM        = (1 << INNER_RADIX); // the fanout of inner node
constexpr int PART_THRESHOLD     = 4096;               // partition threashold

using std::string;
using std::vector;
using std::cout;
using std::endl;


class LinearModelBuilder {
 public:
  double a_;
  double b_;
  LinearModelBuilder() { }

  inline void add(_key_t x, int y) {
    count_++;
    x_sum_ += static_cast<long double>(x);
    y_sum_ += static_cast<long double>(y);
    xx_sum_ += static_cast<long double>(x) * x;
    xy_sum_ += static_cast<long double>(x) * y;
    x_min_ = std::min(x, x_min_);
    x_max_ = std::max(x, x_max_);
    y_min_ = std::min<double>(y, y_min_);
    y_max_ = std::max<double>(y, y_max_);
  }

  void build() {
    if (count_ <= 1) {
      a_ = 0;
      b_ = static_cast<double>(y_sum_);
      return;
    }

    if (static_cast<long double>(count_) * xx_sum_ - x_sum_ * x_sum_ == 0) {
      // all values in a bucket have the same key.
      a_ = 0;
      b_ = static_cast<double>(y_sum_) / count_;
      return;
    }

    auto slope = static_cast<double>(
        (static_cast<long double>(count_) * xy_sum_ - x_sum_ * y_sum_) /
        (static_cast<long double>(count_) * xx_sum_ - x_sum_ * x_sum_));
    auto intercept = static_cast<double>(
        (y_sum_ - static_cast<long double>(slope) * x_sum_) / count_);
    a_ = slope;
    b_ = intercept;

    // If floating point precision errors, fit spline
    if (a_ <= 0) {
      a_ = (y_max_ - y_min_) / (x_max_ - x_min_);
      b_ = -static_cast<double>(x_min_) * a_;
    }
  }

 private:
  int count_ = 0;
  long double x_sum_ = 0;
  long double y_sum_ = 0;
  long double xx_sum_ = 0;
  long double xy_sum_ = 0;
  _key_t x_min_ = 0;
  _key_t x_max_ = INT64_MAX;
  double y_min_ = std::numeric_limits<double>::max();
  double y_max_ = std::numeric_limits<double>::lowest();
};

/*
    WARNING: do NOT use virtual function in Persistent Memory (PM)
    Virtual pointers are only valid in the life time of a process,
    you can NOT store the virtual pointer into PM
*/

struct PackedNode {
    // model parameter
    double a_;
    double b_; // loss some predict accuracy
    // seperator key and value
    _key_t key; 
    struct {
        uint64_t is_leaf  : 1;
        uint64_t capacity : 15;
        uint64_t val      : 48;
    } compound;
    
    
    PackedNode():key(INT64_MAX) {
        compound.capacity = 0;
        compound.val = 0;
    }

    bool operator!=(const PackedNode & b) {
        return key != b.key;
    }

    inline int predict(_key_t k) {
        return static_cast<int>((double)a_ * static_cast<double>(k) + b_);
    }

    inline bool isLeaf() {
        return compound.is_leaf == 1;
    }

    inline uint64_t get_capacaity() {
        return compound.capacity;
    }

    inline char * get_val() {
        uint64_t v = compound.val;
        return (char *) v;
    }

    inline void set_val(char * v) {
        compound.val = (uint64_t)v;
    }
};

struct DNode {
/* 
*  Data node of radix-like index.
*  The Data node is intended for better lookup performance, moderate insertions is supported.
*  More over, data node split/expand can also be exploited for better insertion performance.
*/  
    Record slots_[0];
    
    void clear() {
        ;
    }
};

struct INode {
/* 
* Inner node of radix-like index.
* It's immutable and built when bulkloading from sorted array. A inner node has MAX_KID_NUM 
* children, thay can be a inner node or a data node. Adjacent child pointers may point to the 
* same data node, that means they all share a same data node.
*/  
    PackedNode kids_[MAX_KID_NUM];

    void clear() {
        PackedNode * last_child = nullptr;
        for(int i = 0; i < MAX_KID_NUM; i++) {
            PackedNode *child = kids_ + i;
            
            if(child != last_child) {
                if(child->isLeaf()) {
                    ((DNode *)child->get_val())->clear();
                } else {
                    ((INode *)child->get_val())->clear();
                } 
            }
            last_child = child;
        }

        delete this; // free this node, because it is allocated by new operator
    }

    int adjust(int predict_pos, _key_t k) {
        int pos = std::min(std::max(0, predict_pos), MAX_KID_NUM - 1); // adjuct predict position to [0, MAX_KID_NUM - 1]
        /*
        * As we split the records evenly in the key space, we won't predict to the leftside of the accurate position.
        * But we may predict to its rightside, beacause in some setting, the search procedure needs locate to the lower 
        * bound record. We may predict ourself to a partition in which no records appears at our leftside, so go check
        * partition left-wards
        */
        if(pos > 0 && kids_[pos].key > k) { // current key is not leq than k, go check left
            pos -= 1;
            while(pos > 0) {
                if(kids_[pos].key <= k)
                    break;
                pos--;
            }
            return pos;
        }

        return pos; // we predict correctly
    }
};

class PRadix {
private:
    struct partition_t { // partition type for radix partition
        int level;
        int left_boundary;
        int right_boundary;
        partition_t(): level(0) {} // default constructor
        partition_t(int l, int lb, int rb): level(l), left_boundary(lb), right_boundary(rb) {}
    };

    static constexpr int LINEAR_SPAN = 32; // linear search threshold, after which we start binary search
    /*
        When the total record number is less equal than INITIAL_CAP, 
        rather than build a PRadix tree, store records in a small array
    */
    static constexpr int INITIAL_CAP = 32;

public:
    PRadix(vector<Record> & recs) { // Construct a new PRadix on Persistent Memory
        int num_keys = recs.size();

        if(num_keys <= INITIAL_CAP) {
            leaf_array_ = (Record *) malloc(sizeof(Record) * INITIAL_CAP);
            root_.set_val(NULL);
            for(int i = 0; i < num_keys; i++) {
                leaf_array_[i] = recs[i];
            }
        } else {
            leaf_array_ = NULL;
            root_ = bulk_load_node(recs.data(), num_keys);
            root_.key = INT64_MAX;
        }
        initial_size_ = num_keys;
    }
    
    char * find_lower(_key_t k) {
        if(leaf_array_ != NULL) {
            _key_t max_leqkey = leaf_array_[0].key;
            int max_leqi = 0;
            for(int i = 1; i < initial_size_; i++) {
                if(leaf_array_[i].key <= k && leaf_array_[i].key > max_leqkey) {
                    max_leqi = i;
                    max_leqkey = leaf_array_[i].key;
                }
            }

            return leaf_array_[max_leqi].val;
        }

        PackedNode cur = root_;
        while(!cur.isLeaf()) {
            int bucketID = cur.predict(k);
            INode * inode = (INode *) cur.get_val();
            bucketID = inode->adjust(bucketID, k);
            cur = inode->kids_[bucketID];
        }

        int bucketID = cur.predict(k);
        DNode * dnode = (DNode *) cur.get_val();
        Record * e = exponential_search(dnode, cur.get_capacaity(), bucketID, k);

        // TODO: if the search key is not in the valid range
        if(e == nullptr) {
            printf("%ld not in the valid range\n", k);
            return nullptr;
        }

        return e->val;
    }

    bool insert(_key_t k, char * v) {
        if(leaf_array_ != NULL) {
            leaf_array_[initial_size_++] = {k, v};
            
            if(initial_size_ == INITIAL_CAP) { // the initial buffer is full
                std::sort(leaf_array_, leaf_array_ + INITIAL_CAP, [](const Record &a, const Record & b) {return a.key < b.key;});

                root_ = bulk_load_node(leaf_array_, INITIAL_CAP);
                root_.key = INT64_MAX;

                delete leaf_array_;
                leaf_array_ = NULL;
            }
            return true;
        }

        PackedNode cur = root_;
        while(!cur.isLeaf()) {
            int bucketID = cur.predict(k);
            INode * inode = (INode *) cur.get_val();
            bucketID = inode->adjust(bucketID, k);
            cur = inode->kids_[bucketID];
        }

        int bucketID = cur.predict(k);
        DNode * dnode = (DNode *) cur.get_val();
        Record * e = exponential_search(dnode, cur.get_capacaity(), bucketID, k);

        if(e == nullptr) return false;

        // check the sibling slots to find vacuum
        Record * left = e - 1;
        Record * right= e + 1;
        if(left >= dnode->slots_ && left->key == e->key) {
            *e = {k, v};
            return true;
        } else if(right < dnode->slots_ + cur.get_capacaity() && right->key == e->key) {
            *right = {k, v};
            return true;
        } else { // no more space
            return false;
        }
    }

    void clear() {
        if (leaf_array_ != NULL)
            delete leaf_array_;
        else {
            if(root_.isLeaf()) {
                ((DNode *)root_.get_val())->clear();
            } else {
                ((INode *)root_.get_val())->clear();
            }
        }
    }

private:
    void adaptive_radix_partition(Record *records, int num_keys, vector<partition_t> & parts) 
    {
    /* 
    *  Adaptive Radix Partition separates the key space to k (power of 2) parts recursively till the
    *  record number in one partition converge to a lower bound. 
    *  Adaptive means partition size in different nodes are not equivalent, it's dynamiclly determined
    *  by the partition procedure of a node.
    */
        parts[0] = {0, 0, num_keys - 1};
        
        // caculate some immediate values beforehand
        _key_t key_space = records[num_keys - 1].key - records[0].key;
        double key_piece = (double) key_space / (MAX_KID_NUM - 1);

        for(int i = 0; i < INNER_RADIX; i++) {
            // caculate some immediate values beforehand
            int step = 1 << (INNER_RADIX - i);
            double tmp_imm = key_piece * step / 2;
            
            // do partition to each partition from last level 
            for(int j = 0; j < (1 << INNER_RADIX); j += step) {
                const partition_t partition = parts[j];
                // check the size of this partition
                if(partition.level < i || (partition.right_boundary - partition.left_boundary + 1) < \
                                PART_THRESHOLD) {
                    // this partition is small enough, stop partition it
                    continue;
                }

                double split_val = records[0].key + key_piece * j + tmp_imm;
                // find the median position of current partition
                int mid_pos = std::lower_bound(records + partition.left_boundary, 
                                                records + partition.right_boundary + 1, 
                                                split_val) - records;

                // do binary partition according to the parition position
                parts[j] = {i + 1, partition.left_boundary, std::max(partition.left_boundary, mid_pos - 1)};
                parts[j + step / 2] = {i + 1, mid_pos, partition.right_boundary};
                // now partition j is divide into two equal size sub partition
            }
        }

        return ;
    }

    PackedNode bulk_load_node(Record * records, int num_keys) {   
        // build a radix tree through partition from top down
        PackedNode node;
        if(num_keys <= PART_THRESHOLD) { // if records number in this partition is small, make this node a data node
            #ifdef GAPPED
                int node_cap = (num_keys * 2);
            #else
                int node_cap = num_keys;
            #endif
            DNode * dnode = (DNode *) malloc(sizeof(DNode) + node_cap * sizeof(Record));
            node.compound.is_leaf = 1;
            node.compound.capacity = node_cap;
            node.set_val((char *) dnode);
            
            LinearModelBuilder builder;
            #ifdef GAPPED
                // empty space is left sandwiched in the node
                for(int i = 0; i < num_keys; i++) {
                    dnode->slots_[2 * i] = records[i];
                    dnode->slots_[2 * i + 1] = records[i];
                    builder.add(records[i].key, 2 * i);
                }
            #else
                // empty space is left at the right side
                for(int i = 0; i < num_keys; i++) {
                    dnode->slots_[i] = records[i];
                    builder.add(records[i].key, i);
                } 
                for(int i = num_keys; i < dnode->capacity_; i++)
                    dnode->slots_[i] = {INT64_MAX, NULL};
            #endif
            builder.build();

            node.a_ = builder.a_;
            node.b_ = builder.b_;
        } else { // otherwise, recursively partition records

            // Do radix partition to region records[0 ~ numkey - 1]
            vector<partition_t> partitions(1 << INNER_RADIX);
            adaptive_radix_partition(records, num_keys, partitions);
            // allocate inode memory
            INode * inode = (INode *)malloc(sizeof(INode));
            node.compound.is_leaf = 0;
            node.set_val((char *) inode);
            // update the linear model
            node.a_= (double) (MAX_KID_NUM - 1) / (records[num_keys - 1].key - records[0].key);
            node.b_ = (-1.0 * node.a_ * records[0].key);

            _key_t key_space = records[num_keys - 1].key - records[0].key;
            double key_piece = (double) key_space / (MAX_KID_NUM - 1);
            // recursively build the children of inode and link it to inode
            for(int i = 0; i < (1 << INNER_RADIX); ) {
                const partition_t & part = partitions[i];
                
                // build a child node to region records[left_boundary, right_boundary]
                PackedNode child = bulk_load_node(records + part.left_boundary, 
                                            (part.right_boundary - part.left_boundary + 1));
                
                // link all children to the new (merged) child
                int step = 1 << (INNER_RADIX - part.level);
                double tmp_imm = records[0].key + key_piece * i;
                for(int j = i; j < i + step; j++) {
                    // record the split_key of child
                    #ifdef RECORDS_ARE_DATA
                        // Records are pure data, searching procedure just check for existence
                        child.key = tmp_imm; 
                    #else
                        // Records are not pure data, but seperators of data, store the accurate seperators
                        child.key = records[part.left_boundary].key;
                    #endif

                    inode->kids_[j] = child;
                }

                i += step; // jump to next partition
            }
        }
        return node;
    }

    Record * exponential_search(DNode * dnode, int capacity, int predict_pos, _key_t k) {   
        int node_cap = capacity;
        Record * cur = std::min(dnode->slots_ + std::max(0, predict_pos), dnode->slots_ + node_cap - 1);
        if(cur->key < k) { // go right // linearly search for several Record
            for(int i = 0; i < LINEAR_SPAN && cur < dnode->slots_ + node_cap; i++) {
                if(cur->key >= k) { // the first > k
                    return cur->key == k ? cur : cur - 1;
                }
                cur++;
            }
            // cur->key still less than k
            if(cur == dnode->slots_ + node_cap)
                return cur - 1;

            // search exponentially
            int step = LINEAR_SPAN * 2;
            while(cur < dnode->slots_ + node_cap - 1) {
                Record * right_bound = std::min(dnode->slots_ + node_cap - 1, cur + step);
                if(right_bound->key > k) { // right_bound->key > k and cur->key < k
                    return binary_search(cur, right_bound, k);
                } else if(right_bound->key == k) {
                    return right_bound;
                } else { 
                    step *= 2;
                    cur = right_bound;
                }
            }
            // right_bound->key < k and right_bound reach the right end
            if(cur == dnode->slots_ + node_cap - 1)
                return cur;
            
            return nullptr;
        } else if (cur->key == k){
            return cur;
        } else { // cur->key greater than k, go left
            for(int i = 0; i < LINEAR_SPAN && cur > dnode->slots_; i++) {
                if(cur->key <= k) { // the first <= k
                    return cur;
                }
                cur--;
            }
            if(cur == dnode->slots_)
                return cur;

            // search exponentially
            int step = LINEAR_SPAN * 2;
            while(cur > dnode->slots_) {
                Record * left_bound = std::max(dnode->slots_ + 0, cur - step);
                if(left_bound->key < k) {
                    return binary_search(left_bound, cur, k);
                } else if(left_bound->key == k) {
                    return left_bound;
                } else{
                    step *= 2;
                    cur = left_bound;
                }
            }

            return nullptr;
        }
    }

    Record * binary_search(Record * left, Record *right, _key_t k) {    
        Record * low = left;
        Record * high = right;
        Record * mid;
        while(low <= high) {
            mid = low + (high - low) / 2;
            if(mid->key < k) {
                low = mid + 1;
            } else if (mid->key == k){
                return mid; // if one key eqaul to k, return the exact pos
            } else {
                high = mid - 1;
            }
        }
        return low - 1; 
    }

public:
    Record * leaf_array_;
    uint64_t initial_size_; // used in bulkload
    PackedNode root_;
};

} // namespace lptree

#endif