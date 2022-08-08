/*
     Copyright (c) 2018, UNIST. All rights reserved.  The license is a free
     non-exclusive, non-transferable license to reproduce, use, modify and display
     the source code version of the Software, with or without modifications solely
     for non-commercial research, educational or evaluation purposes. The license
     does not entitle Licensee to technical support, telephone assistance,
     enhancements or updates to the Software. All rights, title to and ownership
     interest in the Software, including all intellectual property rights therein
     shall remain in UNIST. 
*/

#ifndef __FAST_H__
#define __FAST_H__

#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <cassert>
#include <climits>
#include <future>
#include <mutex>

#include "common.h"
#include "flush.h"
#include "pmallocator.h"

using namespace std;

namespace fast {

const int PAGESIZE = 256;

#define IS_FORWARD(c) (c % 2 == 0)
pthread_mutex_t print_mtx;

class page;
class btree;

typedef btree uptree_t;

struct entrance_t {
    page* root;
};

class btree {
    public:
        uint64_t GetHitCnt() {return 0;}
        size_t GetCacheUse() {return 0;}
        double GetVerTime() {return 0;}
    private:
        void btree_insert_internal(char *, _key_t, char *, uint32_t);
        void btree_delete_internal(_key_t, char *, uint32_t, _key_t *, bool *, page **);
        inline page * get_root() {
            return (page *)galc->absolute(entrance->root);
        }

    public:
        entrance_t * entrance;
        page * root_;

        btree(entrance_t *ent): entrance(ent){}

        btree(string path, bool recover = false, string id = "fastfair");
        
        btree(vector<Record> & recs);

        void setNewRoot(void *);
        void getNumberOfNodes();

        void insert(_key_t, _value_t);
        bool remove(_key_t);
        bool find(_key_t, _value_t &);
        bool update(_key_t, _value_t);
        int scan(_key_t, int, _value_t *); 
        
        void printAll();
		void print_height();

        friend class page;

    public:
        bool try_remove(_key_t, bool & need_rebuild);

        char ** find_lower(_key_t) const;

        char ** find_first();

        bool isEmpty() {
            return entrance->root == NULL;
        }

        void clear(string id = "fastfair"); // add a clear function to reclainment the persistent memory space
        void clear(void **, void *); 
};

class header{
    private:
        page* leftmost_ptr;         // 8 bytes
        page* sibling_ptr;          // 8 bytes
        uint32_t level;             // 4 bytes
        uint8_t switch_counter;     // 1 bytes
        uint8_t is_deleted;         // 1 bytes
        int16_t last_index;         // 2 bytes
        uint64_t unused;
        //Spinlock mtx;      // 8 bytes

        friend class page;
        friend class btree;

    public:
        header() {
            leftmost_ptr = NULL;  
            sibling_ptr = NULL;
            switch_counter = 0;
            last_index = -1;
            is_deleted = false;
        }

        ~header() {

        }
};

const int cardinality = (PAGESIZE-sizeof(header))/sizeof(Record);
const int count_in_line = CACHE_LINE_SIZE / sizeof(Record);

class page{
    private:
        header hdr;  // header in persistent memory, 16 bytes
        Record records[cardinality]; // slots in persistent memory, 16 bytes * n

    public:
        friend class btree;

        page(uint32_t level = 0) {
            hdr.level = level;
            records[0].val = NULL;
        }

        // this is called when tree grows
        page(page* left, _key_t key, page* right, uint32_t level = 0):hdr() {
            hdr.leftmost_ptr = (page *) galc->relative(left);  
            hdr.level = level;
            records[0].key = key;
            records[0].val = (char*) galc->relative(right);
            records[1].val = NULL;

            hdr.last_index = 0;

            clflush((char*)this, sizeof(page));
        }

        void *operator new(size_t size) {
            void * ret = galc->malloc(size);
            return ret;
        }

        inline int count() {
            uint8_t previous_switch_counter;
            int count = 0;
            do {
                previous_switch_counter = hdr.switch_counter;
                count = hdr.last_index + 1;

                while(count >= 0 && records[count].val != NULL) {
                    if(IS_FORWARD(previous_switch_counter))
                        ++count;
                    else
                        --count;
                } 

                if(count < 0) {
                    count = 0;
                    while(records[count].val != NULL) {
                        ++count;
                    }
                }

            } while(previous_switch_counter != hdr.switch_counter);

            return count;
        }

        inline bool remove_key(_key_t key) {
            // Set the switch_counter
            if(IS_FORWARD(hdr.switch_counter)) 
                ++hdr.switch_counter;

            bool shift = false;
            int i;
            for(i = 0; records[i].val != NULL; ++i) {
                if(!shift && records[i].key == key) {
                    records[i].val = (i == 0) ? 
                        (char *)hdr.leftmost_ptr : records[i - 1].val; 
                    shift = true;
                }

                if(shift) {
                    records[i].key = records[i + 1].key;
                    records[i].val = records[i + 1].val;

                    // flush
                    uint64_t records_ptr = (uint64_t)(&records[i]);
                    int remainder = records_ptr % CACHE_LINE_SIZE;
                    bool do_flush = (remainder == 0) || 
                        ((((int)(remainder + sizeof(Record)) / CACHE_LINE_SIZE) == 1) && 
                         ((remainder + sizeof(Record)) % CACHE_LINE_SIZE) != 0);
                    if(do_flush) {
                        clflush((char *)records_ptr, CACHE_LINE_SIZE);
                    }
                }
            }

            if(shift) {
                --hdr.last_index;
            }
            return shift;
        }

        bool remove(btree* bt, _key_t key, bool only_rebalance = false, bool with_lock = true) {
            //hdr.mtx.lock();

            bool ret = remove_key(key);

            //hdr.mtx.unlock();

            return ret;
        }

        /*
         * Although we implemented the rebalancing of B+-Tree, it is currently blocked for the performance.
         * Please refer to the follow.
         * Chi, P., Lee, W. C., & Xie, Y. (2014, August). 
         * Making B+-tree efficient in PCM-based main memory. In Proceedings of the 2014
         * international symposium on Low power electronics and design (pp. 69-74). ACM.
         */
        bool remove_rebalancing(btree* bt, _key_t key, bool only_rebalance = false, bool with_lock = true) {
            if(with_lock) {
                //hdr.mtx.lock();
            }

            if(hdr.is_deleted) {
                if(with_lock) {
                    //hdr.mtx.unlock();
                }
                return false;
            }

            if(!only_rebalance) { // the function not just rebalance the tree, but also remove the key from the node
                int num_entries_before = count();

                // This node is root
                page * root = (page *) galc->absolute(bt->entrance->root);
                if(this == root) {
                    if(hdr.level > 0) {
                        if(num_entries_before == 1 && !hdr.sibling_ptr) {
                            bt->entrance->root = hdr.leftmost_ptr;
                            clflush(&(bt->entrance->root), 8);

                            galc->free(root);
                        }
                    }

                    // Remove the key from this node
                    bool ret = remove_key(key);

                    if(with_lock) {
                        //hdr.mtx.unlock();
                    }
                    return true;
                }

                bool should_rebalance = true;
                // check the node utilization
                if(num_entries_before - 1 >= (int)((cardinality - 1) * 0.5)) { 
                    should_rebalance = false;
                }

                // Remove the key from this node
                bool ret = remove_key(key);

                if(!should_rebalance) {
                    if(with_lock) {
                        //hdr.mtx.unlock();
                    }
                    return (hdr.leftmost_ptr == NULL) ? ret : true;
                }
            } 

            //Remove a key from the parent node
            _key_t deleted_key_from_parent = 0;
            bool is_leftmost_node = false;
            page *left_sibling;
            bt->btree_delete_internal(key, (char *)this, hdr.level + 1,
                    &deleted_key_from_parent, &is_leftmost_node, &left_sibling);

            if(is_leftmost_node) {
                if(with_lock) {
                    //hdr.mtx.unlock();
                }
                
                page * sibling = (page *) galc->absolute(hdr.sibling_ptr);
                sibling->remove_rebalancing(bt, sibling->records[0].key, true, with_lock);//this call does not remove the key
                //hdr.sibling_ptr->remove(bt, hdr.sibling_ptr->records[0].key, true, with_lock);//just rebalance, not to remove the key
                return true;
            }

            if(with_lock) {
                //left_sibling->hdr.mtx.lock();
            }

            while(left_sibling->hdr.sibling_ptr != (page *)galc->relative(this)) {
                if(with_lock) {
                    page *t = (page *)galc->absolute(left_sibling->hdr.sibling_ptr);
                    //left_sibling->hdr.mtx.unlock();
                    left_sibling = t;
                    //left_sibling->hdr.mtx.lock();
                }
                else
                    left_sibling = (page *)galc->absolute(left_sibling->hdr.sibling_ptr);
            }

            int num_entries = count();
            int left_num_entries = left_sibling->count();

            // Merge or Redistribution
            int total_num_entries = num_entries + left_num_entries;
            if(hdr.leftmost_ptr)
                ++total_num_entries;

            _key_t parent_key;

            if(total_num_entries > cardinality - 1) { // Redistribution
                int m = (int) ceil(total_num_entries / 2);

                if(num_entries < left_num_entries) { // left -> right
                    if(hdr.leftmost_ptr == NULL){
                        for(int i=left_num_entries - 1; i>=m; i--){
                            insert_key
                                (left_sibling->records[i].key, left_sibling->records[i].val, &num_entries); 
                        } 

                        left_sibling->records[m].val = NULL;
                        clflush((char *)&(left_sibling->records[m].val), sizeof(char *));

                        left_sibling->hdr.last_index = m - 1;
                        clflush((char *)&(left_sibling->hdr.last_index), sizeof(int16_t));

                        parent_key = records[0].key; 
                    }
                    else{
                        insert_key(deleted_key_from_parent, (char*)hdr.leftmost_ptr,
                                &num_entries); 

                        for(int i=left_num_entries - 1; i>m; i--){
                            insert_key
                                (left_sibling->records[i].key, left_sibling->records[i].val, &num_entries); 
                        }

                        parent_key = left_sibling->records[m].key; 

                        hdr.leftmost_ptr = (page*)left_sibling->records[m].val; 
                        clflush((char *)&(hdr.leftmost_ptr), sizeof(page *));

                        left_sibling->records[m].val = NULL;
                        clflush((char *)&(left_sibling->records[m].val), sizeof(char *));

                        left_sibling->hdr.last_index = m - 1;
                        clflush((char *)&(left_sibling->hdr.last_index), sizeof(int16_t));
                    }

                    if(left_sibling == ((page *)bt->get_root())) {
                        page* new_root = new page(left_sibling, parent_key, this, hdr.level + 1);
                        bt->setNewRoot((void *)new_root);
                    }
                    else {
                        bt->btree_insert_internal
                            ((char *)left_sibling, parent_key, (char *)this, hdr.level + 1);
                    }
                }
                else { // from leftmost case
                    hdr.is_deleted = 1;
                    clflush((char *)&(hdr.is_deleted), sizeof(uint8_t));

                    page* new_sibling = new page(hdr.level); 
                    //new_sibling->hdr.mtx.lock();
                    new_sibling->hdr.sibling_ptr = hdr.sibling_ptr;

                    int num_dist_entries = num_entries - m;
                    int new_sibling_cnt = 0;

                    if(hdr.leftmost_ptr == NULL){
                        for(int i=0; i<num_dist_entries; i++){
                            left_sibling->insert_key(records[i].key, records[i].val,
                                    &left_num_entries); 
                        } 

                        for(int i=num_dist_entries; records[i].val != NULL; i++){
                            new_sibling->insert_key(records[i].key, records[i].val,
                                    &new_sibling_cnt, false); 
                        } 

                        clflush((char *)(new_sibling), sizeof(page));

                        left_sibling->hdr.sibling_ptr = galc->relative(new_sibling);
                        clflush(&(left_sibling->hdr), sizeof(header));

                        parent_key = new_sibling->records[0].key; 
                    }
                    else{
                        left_sibling->insert_key(deleted_key_from_parent,
                                (char*)hdr.leftmost_ptr, &left_num_entries);

                        for(int i=0; i<num_dist_entries - 1; i++){
                            left_sibling->insert_key(records[i].key, records[i].val,
                                    &left_num_entries); 
                        } 

                        parent_key = records[num_dist_entries - 1].key;

                        new_sibling->hdr.leftmost_ptr = (page*)records[num_dist_entries - 1].val;
                        for(int i=num_dist_entries; records[i].val != NULL; i++){
                            new_sibling->insert_key(records[i].key, records[i].val,
                                    &new_sibling_cnt, false); 
                        } 
                        clflush((char *)(new_sibling), sizeof(page));

                        left_sibling->hdr.sibling_ptr = galc->relative(new_sibling);
                        clflush(&(left_sibling->hdr), sizeof(header));
                    }

                    if(left_sibling == ((page *)bt->get_root())) {
                        page* new_root = new page(left_sibling, parent_key, new_sibling, hdr.level + 1);
                        bt->setNewRoot((void *)new_root);
                    }
                    else {
                        bt->btree_insert_internal
                            ((char *)left_sibling, parent_key, (char *)new_sibling, hdr.level + 1);
                    }

                    //new_sibling->hdr.mtx.unlock();
                }
            } else {
                hdr.is_deleted = 1;
                clflush((char *)&(hdr.is_deleted), sizeof(uint8_t));

                if(hdr.leftmost_ptr)
                    left_sibling->insert_key(deleted_key_from_parent, 
                            (char *)hdr.leftmost_ptr, &left_num_entries);

                for(int i = 0; records[i].val != NULL; ++i) { 
                    left_sibling->insert_key(records[i].key, records[i].val, &left_num_entries);
                }

                left_sibling->hdr.sibling_ptr = hdr.sibling_ptr;
                clflush((char *)&(left_sibling->hdr.sibling_ptr), sizeof(page *));
            }

            if(with_lock) {
                ;//left_sibling->hdr.mtx.unlock();
                //hdr.mtx.unlock();
            }

            return true;
        }

        inline void 
            insert_key(_key_t key, char * ptr, int *num_entries, bool flush = true,
                    bool update_last_index = true) {
                // update switch_counter
                if(!IS_FORWARD(hdr.switch_counter))
                    ++hdr.switch_counter;

                // FAST
                if(*num_entries == 0) {  // this page is empty
                    Record* new_entry = (Record*) &records[0];
                    Record* array_end = (Record*) &records[1];
                    new_entry->key = (_key_t) key;
                    new_entry->val = ptr;

                    array_end->val = (char*)NULL;

                    if(flush) {
                        clflush((char*) this, CACHE_LINE_SIZE);
                    }
                }
                else {
                    int i = *num_entries - 1, inserted = 0, to_flush_cnt = 0;
                    records[*num_entries+1].val = records[*num_entries].val; 
                    if(flush) {
                        if((uint64_t)&(records[*num_entries+1].val) % CACHE_LINE_SIZE == 0) 
                            clflush((char*)&(records[*num_entries+1].val), sizeof(char*));
                    }

                    // FAST
                    for(i = *num_entries - 1; i >= 0; i--) {
                        if(key < records[i].key ) {
                            records[i+1].val = records[i].val;
                            records[i+1].key = records[i].key;

                            if(flush) {
                                uint64_t records_ptr = (uint64_t)(&records[i+1]);

                                int remainder = records_ptr % CACHE_LINE_SIZE;
                                bool do_flush = (remainder == 0) || 
                                    ((((int)(remainder + sizeof(Record)) / CACHE_LINE_SIZE) == 1) 
                                     && ((remainder+sizeof(Record))%CACHE_LINE_SIZE)!=0);
                                if(do_flush) {
                                    clflush((char*)records_ptr,CACHE_LINE_SIZE);
                                    to_flush_cnt = 0;
                                }
                                else
                                    ++to_flush_cnt;
                            }
                        }
                        else{
                            records[i+1].val = records[i].val;
                            records[i+1].key = key;
                            records[i+1].val = const_cast<char *>(ptr);

                            if(flush)
                                clflush((char*)&records[i+1],sizeof(Record));
                            inserted = 1;
                            break;
                        }
                    }
                    if(inserted==0){
                        records[0].val =(char*) hdr.leftmost_ptr;
                        records[0].key = key;
                        records[0].val = const_cast<char *>(ptr);
                        if(flush)
                            clflush((char*) &records[0], sizeof(Record)); 
                    }
                }

                if(update_last_index) {
                    hdr.last_index = *num_entries;
                }
                ++(*num_entries);
            }

        // Insert a new key - FAST and FAIR
        page *store (btree* bt, char * left, _key_t key, char * right,
             bool flush, bool with_lock, page *invalid_sibling = NULL) {
                if(with_lock) {
                    //hdr.mtx.lock(); // Lock the write lock
                }
                if(hdr.is_deleted) {
                    if(with_lock) {
                        //hdr.mtx.unlock();
                    }

                    return NULL;
                }

                // If this node has a sibling node,
                if(hdr.sibling_ptr && (hdr.sibling_ptr != invalid_sibling)) {
                    // Compare this key with the first key of the sibling
                    page * sibling = (page *) galc->absolute(hdr.sibling_ptr);
                    if(key > sibling->records[0].key) {
                        if(with_lock) { 
                            //hdr.mtx.unlock(); // Unlock the write lock
                        }
                        return sibling->store(bt, NULL, key, right, 
                                true, with_lock, invalid_sibling);
                    }
                }

                int num_entries = count();

                // FAST
                if(num_entries < cardinality - 1) {
                    insert_key(key, hdr.leftmost_ptr == NULL ? (char *)right: (char *)galc->relative(right), &num_entries, flush);

                    if(with_lock) {
                        //hdr.mtx.unlock(); // Unlock the write lock
                    }

                    return this;
                }
                else {// FAIR
                    // overflow
                    // create a new node
                    page* sibling = new page(hdr.level); 
                    int m = (int) ceil(num_entries/2);
                    _key_t split_key = records[m].key;

                    // migrate half of keys into the sibling
                    int sibling_cnt = 0;
                    if(hdr.leftmost_ptr == NULL){ // leaf node
                        for(int i=m; i<num_entries; ++i){ 
                            sibling->insert_key(records[i].key, records[i].val, &sibling_cnt, false);
                        }
                    }
                    else{ // internal node
                        for(int i=m+1;i<num_entries;++i){ 
                            sibling->insert_key(records[i].key, records[i].val, &sibling_cnt, false);
                        }
                        sibling->hdr.leftmost_ptr = (page*) records[m].val;
                    }

                    sibling->hdr.sibling_ptr = hdr.sibling_ptr;
                    clflush((char *)sibling, sizeof(page));

                    hdr.sibling_ptr = (page *)galc->relative(sibling);
                    clflush((char*) &hdr, sizeof(hdr));

                    // set to NULL
                    if(IS_FORWARD(hdr.switch_counter))
                        hdr.switch_counter += 2;
                    else
                        ++hdr.switch_counter;
                    records[m].val = NULL;
                    clflush((char*) &records[m], sizeof(Record));

                    hdr.last_index = m - 1;
                    clflush((char *)&(hdr.last_index), sizeof(int16_t));

                    num_entries = hdr.last_index + 1;

                    page *ret;

                    // insert the key
                    if(key < split_key) {
                        insert_key(key, hdr.leftmost_ptr == NULL ? (char *)right: (char *)galc->relative(right), &num_entries);
                        ret = this;
                    }
                    else {
                        sibling->insert_key(key, hdr.leftmost_ptr == NULL ? (char *)right: (char *)galc->relative(right), &sibling_cnt);
                        ret = sibling;
                    }

                    // Set a new root or insert the split key to the parent
                    if(bt->get_root() == this) { // only one node can update the root ptr
                        page* new_root = new page((page*)this, split_key, sibling, 
                                hdr.level + 1);
                        bt->setNewRoot((void *)new_root);

                        if(with_lock) {
                            //hdr.mtx.unlock(); // Unlock the write lock
                        }
                    }
                    else {
                        if(with_lock) {
                            //hdr.mtx.unlock(); // Unlock the write lock
                        }
                        bt->btree_insert_internal(NULL, split_key, (char *)sibling, 
                                hdr.level + 1);
                    }

                    return ret;
                }

            }

        bool update_value(_key_t key, _value_t value) {
                assert(hdr.leftmost_ptr == NULL);

                //hdr.mtx.lock();
                int pos = 0;
                uint8_t previous_switch_counter;
                do {
                        previous_switch_counter = hdr.switch_counter;

                        if(IS_FORWARD(previous_switch_counter)) {
                                if(records[0].key == key && records[0].val != NULL) {
                                        pos = 0;
                                        continue;
                                }

                                for(int i = 1; records[i].val != NULL; ++i) { 
                                        if(key == records[i].key && records[i-1].val != records[i].val) { 
                                                pos = i;
                                                break;
                                        } // found a record whose key is great equal to key
                                }
                        } else { // search from right to left
                                for(int i = count() - 1; i >= 0; --i) {
                                        if(i == 0) {
                                                if(key == records[0].key && records[0].val != NULL) {
                                                        pos = 0;
                                                        break;
                                                }
                                        } else {
                                                if(key == records[i].key && records[i-1].val != records[i].val) {
                                                        pos = i;
                                                        break;
                                                }
                                        }
                                }
                        }
                } while(hdr.switch_counter != previous_switch_counter);

                if(key == records[pos].key) {
                        records[pos].val = (char *)value;
                        clflush((char *)&(records[pos]), sizeof(Record));

                        //hdr.mtx.unlock();
                        return true;
                }
                
                //hdr.mtx.unlock();
                page *t = (page *) galc->absolute(hdr.sibling_ptr);
                if(t != NULL && key >= t->records[0].key) { //the node has splited, go check its siblings
                    return ((page *)t)->update_value(key, value);
                } 

                return false;
        }

        // Search keys with linear search
        int linear_search_range(_key_t min, int scan_sz, _value_t *buf) {
                int i, off = 0;
                uint8_t previous_switch_counter;
                page *current = this;

                while(current) {
                    int old_off = off;
                    do {
                        previous_switch_counter = current->hdr.switch_counter;
                        off = old_off;

                        _key_t tmp_key;
                        char *tmp_ptr;

                        if(IS_FORWARD(previous_switch_counter)) {
                            if((tmp_key = current->records[0].key) > min) {
                                if(off < scan_sz) {
                                    if((tmp_ptr = current->records[0].val) != NULL) {
                                        if(tmp_key == current->records[0].key) {
                                            if(tmp_ptr) {
                                                buf[off++] = (_value_t)tmp_ptr;
                                            }
                                        }
                                    }
                                }
                                else
                                    return off;
                            }

                            for(i=1; current->records[i].val != NULL; ++i) { 
                                if((tmp_key = current->records[i].key) > min) {
                                    if(off < scan_sz) {
                                        if((tmp_ptr = current->records[i].val) != current->records[i - 1].val) {
                                            if(tmp_key == current->records[i].key) {
                                                if(tmp_ptr)
                                                    buf[off++] = (_value_t)tmp_ptr;
                                            }
                                        }
                                    }
                                    else
                                        return off;
                                }
                            }
                        }
                        else {
                            for(i=count() - 1; i > 0; --i) { 
                                if((tmp_key = current->records[i].key) > min) {
                                    if(off < scan_sz) {
                                        if((tmp_ptr = current->records[i].val) != current->records[i - 1].val) {
                                            if(tmp_key == current->records[i].key) {
                                                if(tmp_ptr)
                                                    buf[off++] = (_value_t)tmp_ptr;
                                            }
                                        }
                                    }
                                    else
                                        return off;
                                }
                            }

                            if((tmp_key = current->records[0].key) > min) {
                                if(off < scan_sz) {
                                    if((tmp_ptr = current->records[0].val) != NULL) {
                                        if(tmp_key == current->records[0].key) {
                                            if(tmp_ptr) {
                                                buf[off++] = (_value_t)tmp_ptr;
                                            }
                                        }
                                    }
                                }
                                else
                                    return off;
                            }
                        }
                    } while(previous_switch_counter != current->hdr.switch_counter);

                    current = (page *) galc->absolute(current->hdr.sibling_ptr);
                }

                return 0;
            }

        char *linear_search(_key_t key) {
            int i = 1;
            uint8_t previous_switch_counter;
            char *ret = NULL;
            char *t = NULL; 
            _key_t k;

            if(hdr.leftmost_ptr == NULL) { // Search a leaf node
                do {
                    previous_switch_counter = hdr.switch_counter;
                    ret = NULL;

                    // search from left ro right
                    if(IS_FORWARD(previous_switch_counter)) { 
                        if((k = records[0].key) == key) {
                            if((t = records[0].val) != NULL) { // NULL ptr 
                                if(k == records[0].key) {
                                    ret = t;
                                    continue;
                                }
                            }
                        }

                        for(i=1; records[i].val != NULL; ++i) { 
                            if((k = records[i].key) == key) {
                                if(records[i-1].val != (t = records[i].val)) {
                                    if(k == records[i].key) {
                                        ret = t;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    else { // search from right to left
                        for(i = count() - 1; i > 0; --i) {
                            if((k = records[i].key) == key) {
                                if(records[i - 1].val != (t = records[i].val) && t) {
                                    if(k == records[i].key) {
                                        ret = t;
                                        break;
                                    }
                                }
                            }
                        }

                        if(!ret) {
                            if((k = records[0].key) == key) {
                                if(NULL != (t = records[0].val) && t) {
                                    if(k == records[0].key) {
                                        ret = t;
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                } while(hdr.switch_counter != previous_switch_counter);

                if(ret) {
                    return ret;
                }

                page *t = (page *) galc->absolute(hdr.sibling_ptr);
                if(t != NULL && key >= t->records[0].key)
                    return (char *)t;

                return NULL;
            }
            else { // internal node
                do {
                    previous_switch_counter = hdr.switch_counter;
                    ret = NULL;

                    if(IS_FORWARD(previous_switch_counter)) {
                        if(key < (k = records[0].key)) {
                            if((t = (char *)hdr.leftmost_ptr) != records[0].val) { 
                                ret = t;
                                continue;
                            }
                        }

                        for(i = 1; records[i].val != NULL; ++i) { 
                            if(key < (k = records[i].key)) { 
                                if((t = records[i-1].val) != records[i].val) {
                                    ret = t;
                                    break;
                                }
                            }
                        }

                        if(!ret) {
                            ret = records[i - 1].val; 
                            continue;
                        }
                    }
                    else { // search from right to left
                        for(i = count() - 1; i >= 0; --i) {
                            if(key >= (k = records[i].key)) {
                                if(i == 0) {
                                    if((char *)hdr.leftmost_ptr != (t = records[i].val)) {
                                        ret = t;
                                        break;
                                    }
                                }
                                else {
                                    if(records[i - 1].val != (t = records[i].val)) {
                                        ret = t;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                } while(hdr.switch_counter != previous_switch_counter);

                if((t = (char *)hdr.sibling_ptr) != NULL) {
                    if(key >= ((page *)galc->absolute(t))->records[0].key)
                        return t;
                }

                if(ret) {
                    return ret;
                }
                else
                    return (char *)hdr.leftmost_ptr;
            }

            return NULL;
        }
        
        char ** linear_search_lower(_key_t key, _key_t & lower_key) {
            int i = 1;
            uint8_t previous_switch_counter;
            char **ret = NULL, *t;
            _key_t k;

            do {
                previous_switch_counter = hdr.switch_counter;
                ret = NULL;
                if(IS_FORWARD(previous_switch_counter)) {
                    for(i = 1; records[i].val != NULL; ++i) { 
                        if(key < records[i].key) { 
                            if((t = records[i-1].val) != records[i].val) {
                                k = records[i-1].key;
                                ret = &(records[i-1].val);
                                break;
                            }
                        } // found a record whose key is great equal to key
                    }
                    if(!ret) {
                        k = records[i - 1].key;
                        ret = &(records[i - 1].val); 
                        continue;
                    }
                } else { // search from right to left
                    for(i = count() - 1; i >= 0; --i) {
                        if(key >= (k = records[i].key)) {
                            if(i > 0 && records[i - 1].val != (t = records[i].val)) {
                                ret = &(records[i].val);
                                break;
                            } else if(i == 0) {
                                ret = &(records[0].val);
                            }
                        }
                    }
                }
            } while(hdr.switch_counter != previous_switch_counter);

            if(ret) {
                lower_key = k;
                return ret;
            }

            page *tt = (page *) galc->absolute(hdr.sibling_ptr);
            if(tt != NULL && key >= tt->records[0].key)
                return (char **)tt;

            return NULL; // should never be here
        }
        // print a node 
        void print() {
            if(hdr.leftmost_ptr == NULL) 
                printf("[%d] leaf %lx \n", this->hdr.level, (long unsigned int)this);
            else 
                printf("[%d] internal %lx \n", this->hdr.level, (long unsigned int)this);
            printf("last_index: %d\n", hdr.last_index);
            printf("switch_counter: %d\n", hdr.switch_counter);
            printf("search direction: ");
            if(IS_FORWARD(hdr.switch_counter))
                printf("->\n");
            else
                printf("<-\n");

            if(hdr.leftmost_ptr!=NULL) 
                printf("%lx ",(long unsigned int)hdr.leftmost_ptr);

            for(int i=0;records[i].val != NULL;++i)
                printf("%ld,%lx ",records[i].key,(long unsigned int)records[i].val);

            printf("%lx ", (long unsigned int)hdr.sibling_ptr);

            printf("\n");
        }

        void printAll() {
            if(hdr.leftmost_ptr==NULL) {
                printf("printing leaf node: ");
                print();
            }
            else {
                printf("printing internal node: ");
                print();
                ((page*) hdr.leftmost_ptr)->printAll();
                for(int i=0;records[i].val != NULL;++i){
                    ((page*) records[i].val)->printAll();
                }
            }
        }
        
};

/*
 * class btree
 */
btree::btree(string path, bool recover, string id){
    if(recover == false) {
        galc = new PMAllocator(path.c_str(), false, id.c_str());

        // add a entrance to the persistent tree
        entrance = (entrance_t *) galc->get_root(sizeof(entrance_t));
        entrance->root = NULL;
        clflush(entrance, sizeof(entrance_t));

        // allocate the root node
        page * root = (page *)new page(); // persistent new
        
        entrance->root = galc->relative(root);
        clflush(&(entrance->root), 8);
    } else {
        galc = new PMAllocator(path.c_str(), true, id.c_str());
        
        // get the entrance from the pmem file
        entrance = (entrance_t *) galc->get_root(sizeof(entrance_t));

        if(entrance->root == NULL) {
            printf("the tree is empty\n");
            exit(-1);
        }
    }
}

void btree::setNewRoot(void *new_root) {
    entrance->root = (page *)galc->relative(new_root);
    clflush(&(entrance->root), 8);
}

bool btree::find(_key_t key, _value_t & val){
    page* p = (page *) galc->absolute(entrance->root);

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)galc->absolute(p->linear_search(key));
    }

    page *t;
    while((t = (page *)p->linear_search(key)) == p->hdr.sibling_ptr) {
        p = t;
        if(!p) {
            break;
        }
    }

    if(!t) {
        //printf("NOT FOUND %lu, t = %x\n", key, t);
        return false;
    }

    val = (_value_t)t;
    return true;
}

// insert the key in the leaf node
void btree::insert(_key_t key, _value_t right){ //need to be string
    page* p = (page *) galc->absolute(entrance->root);

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)galc->absolute(p->linear_search(key));
    }
    if(!p->store(this, NULL, key, (char *)right, true, true)) { // store 
        insert(key, right);
    }
}

// store the key into the node at the given level 
void btree::btree_insert_internal(char *left, _key_t key, char *right, uint32_t level) {
    page* p = (page *) galc->absolute(entrance->root);

    if(level > p->hdr.level)
        return;

    while(p->hdr.level > level) 
        p = (page *)galc->absolute(p->linear_search(key));

    if(!p->store(this, NULL, key, right, true, true)) {
        btree_insert_internal(left, key, right, level);
    }
}

bool btree::remove(_key_t key) {
    page* p = (page *) galc->absolute(entrance->root);

    while(p->hdr.leftmost_ptr != NULL){
        p = (page *)galc->absolute(p->linear_search(key));
    }

    page *t;
    while((t = (page *)p->linear_search(key)) == p->hdr.sibling_ptr) {
        p = t;
        if(!p)
            break;
    }

    if(p) {
        //Modified by ypluo, eliminate the dead loop on delete key
        return p->remove_rebalancing(this, key, false, true);

        if(p->hdr.is_deleted == 1) {
            galc->free(p);
        }
    
    } else {
        //printf("not found the key to delete %lu\n", key);
        return false;
    }
}

void btree::btree_delete_internal(_key_t key, char *ptr, uint32_t level, _key_t *deleted_key, 
 bool *is_leftmost_node, page **left_sibling) {
    page* p = (page *) galc->absolute(entrance->root);
    
    if(level > p->hdr.level)
        return;


    while(p->hdr.level > level) {
        p = (page *)galc->absolute(p->linear_search(key));
    }

    //p->hdr.mtx.lock();

    if((char *)p->hdr.leftmost_ptr == (char *)galc->relative(ptr)) {
        *is_leftmost_node = true;
        //p->hdr.mtx.unlock();
        return;
    }

    *is_leftmost_node = false;

    for(int i=0; p->records[i].val != NULL; ++i) {
        if(p->records[i].val == (char *)galc->relative(ptr)) {
            if(i == 0) {
                if((char *)p->hdr.leftmost_ptr != p->records[i].val) {
                    *deleted_key = p->records[i].key;
                    *left_sibling = (page *)galc->absolute(p->hdr.leftmost_ptr);
                    
                    p->remove_rebalancing(this, *deleted_key, false, false); // Modified by ypluo
                    //p->remove(this, *deleted_key, false, false);
                    break;
                }
            }
            else {
                if(p->records[i - 1].val != p->records[i].val) {
                    *deleted_key = p->records[i].key;
                    *left_sibling = (page *)galc->absolute(p->records[i - 1].val);
                    p->remove_rebalancing(this, *deleted_key, false, false); // Modified by ypluo
                    //p->remove(this, *deleted_key, false, false);
                    break;
                }
            }
        }
    }

    //p->hdr.mtx.unlock();
}

// insert the key in the leaf node
bool btree::update(_key_t key, _value_t value){ //need to be string
        page* p = (page *) galc->absolute(entrance->root);

        while(p->hdr.leftmost_ptr != NULL) {
            p = (page *)galc->absolute(p->linear_search(key));
        }

        return p->update_value(key, value);
}

// Function to search keys from "min" to "max"
int btree::scan(_key_t min, int scan_sz, _value_t * buf) {
    page* p = (page *) galc->absolute(entrance->root);

    while(p) {
        if(p->hdr.leftmost_ptr != NULL) {
            // The current page is internal
            p = (page *)galc->absolute(p->linear_search(min));
            //p = (page *)p->linear_search(min);
        }
        else {
            // Found a leaf
            return p->linear_search_range(min, scan_sz, buf);
        }
    }

    return 0;
}

void btree::printAll(){
    int total_keys = 0;
    page *leftmost = (page *) galc->absolute(entrance->root);
    
	printf("root: %lx\n", (long unsigned int)leftmost);
    do {
        page *sibling = leftmost;
        while(sibling) {
            if(sibling->hdr.level == 0) {
                total_keys += sibling->hdr.last_index + 1;
            }
            sibling->print();
            sibling = (page *) galc->absolute(sibling->hdr.sibling_ptr);
        }
        printf("-----------------------------------------\n");
        leftmost = (page *) galc->absolute(leftmost->hdr.leftmost_ptr);
    } while(leftmost);

    printf("total number of keys: %d\n", total_keys);
}

// additional function to iterater the tree
char ** btree::find_lower(_key_t key) const{
    page* p = (page *) galc->absolute(entrance->root);

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)galc->absolute(p->linear_search(key));
    }

    page *t = NULL;
    _key_t discard_key;
    while((t = (page *)p->linear_search_lower(key, discard_key)) == p->hdr.sibling_ptr) { 
        p = (page *)galc->absolute(t);
        if(!p) {
            break;
        }
    }
    return (char **)t;
}

bool btree::try_remove(_key_t key, bool & need_rebuild) {
    page* p = (page *) galc->absolute(entrance->root);
    while(p->hdr.leftmost_ptr != NULL){
        p = (page *)galc->absolute(p->linear_search(key));
    }
    page *t;
    _key_t new_key = -1;
    while((t = (page *)p->linear_search_lower(key, new_key)) == p->hdr.sibling_ptr) {
        p = (page *)galc->absolute(t);
        if(!p)
            break;
    }

    p->remove_rebalancing(this, new_key);
    if(p->hdr.is_deleted == 1) {
        galc->free(p);
    }

    need_rebuild = false;
    return true;
    return new_key;
}

void btree::print_height() {
	page* p = (page *) galc->absolute(entrance->root);			
	printf("Height: %d\n", p->hdr.level + 1);
}

btree::btree(vector<Record> & rec) {
    entrance = (entrance_t *) galc->malloc(sizeof(entrance_t)); // this constructor is called by tlbtree only, so use malloc but not get_root
    entrance->root = NULL;
    clflush(entrance, sizeof(entrance_t));

    // allocate the root node
    page * root = (page *)new page();
    
    root_ = root;
    entrance->root = galc->relative(root);

    for(auto r : rec) {
        insert(r.key, (_value_t)r.val);
    }
}

char ** btree::find_first() {
    return NULL;
}

entrance_t * build_uptree(vector<Record> & recs, uptree_t * & new_tree) {
    new_tree = new btree(recs);

    return new_tree->entrance;
}

void free_uptree(entrance_t * ent) {
    // fastfair will never need to rebuild
    return ;
}

}; // fast

#endif // __FAST_H__
