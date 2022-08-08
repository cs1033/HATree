/*
    Read Optimized B+-tree (hatree3)
    Copyright(c) Luo Yongping All rights reserved!
*/

#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <cstdio>
#include <queue>
#include <algorithm>
#include <iostream>

#include "pmallocator.h"
#include "flush.h"
 
using std::string;
using std::cout;
using std::endl;

namespace hatree3 {

const uint16_t INNER_NODE_SIZE = 30;
const uint16_t SEARCH_K = 3;                     // Submit hotspot every 5 queries
const uint16_t NODE_SIZE = 14;                  // INode size 256B
const int INNER_UNDERFLOW_CARD = INNER_NODE_SIZE / 2;
const int LEAF_UNDERFLOW_CARD = NODE_SIZE / 3;  // Underflow limit
const uint16_t INNER_MAX = INNER_NODE_SIZE ; // one quarter of the slots are reserved as read buffer
const uint16_t BUFFER_SIZE = 8;
const uint16_t BUFFER2_SIZE = 24;

struct LNode;
struct INode;



struct LNode { // leaf nodes of btree, allocated on Optane
    // union state_t { // an 8 bytes states type
    //     //1
    //     uint64_t pack;
    //     //2
    //     struct unpack_t {
    //         uint64_t bitmap         : 56;
    //         uint64_t sib_version    : 1;
    //         uint64_t node_version   : 7;
    //     } unpack;

    //     inline uint8_t count() {
    //         return (uint8_t)_mm_popcnt_u64(unpack.bitmap);
    //     }

    //     inline bool read(int8_t idx) {
    //         return (unpack.bitmap & ((uint64_t)0x8000000000000000 >> (idx + 8))) > 0;
    //     }

    //     inline int8_t alloc() {
    //         uint64_t tmp = ((uint64_t)0xff00000000000000 | unpack.bitmap);
    //         return __builtin_ia32_lzcnt_u64(~tmp) - 8;
    //     }

    //     inline uint64_t add(int8_t idx) {
    //         return unpack.bitmap + ((uint64_t)0x8000000000000000 >> (idx + 8));
    //     }

    //     inline uint64_t free(int8_t idx) {
    //         return unpack.bitmap - ((uint64_t)0x8000000000000000 >> (idx + 8));
    //     }

    //     inline uint8_t get_sibver() {
    //         return (uint8_t) unpack.sib_version;
    //     }
    // };
    // 64 B header
    union state_t { // an 2 bytes states type
        //1
        uint16_t pack;
        //2
        struct unpack_t {
            uint16_t bitmap         : 14;
            uint16_t sib_version    : 1;
            uint16_t node_version   : 1;
        } unpack;

        inline uint8_t count() {
            return (uint8_t)_mm_popcnt_u32(unpack.bitmap);
        }

        inline bool read(int8_t idx) {
            return (unpack.bitmap & ((uint16_t)0x8000 >> (idx + 2))) > 0;
        }

        inline int8_t alloc() {
            uint32_t tmp = ((uint64_t)0xFFFFC000 | unpack.bitmap);
            return __builtin_ia32_lzcnt_u32(~tmp) - 18;
        }

        inline uint16_t add(int8_t idx) {
            return unpack.bitmap + ((uint16_t)0x8000 >> (idx + 2));
        }

        inline uint16_t free(int8_t idx) {
            return unpack.bitmap - ((uint16_t)0x8000 >> (idx + 2));
        }

        inline uint8_t get_sibver() {
            return (uint8_t) unpack.sib_version;
        }
    };
    state_t state_;
    char finger_prints_[NODE_SIZE]; 
    
    Record recs_[NODE_SIZE]; // all empty slots
    char* sibs_[2];         // two shadow siblings for the sake of failure-atomicity

    LNode() {
        sibs_[0] = sibs_[1] = NULL;
        state_.pack = 0;
    }

    void * operator new(size_t size) {
        return galc->malloc(size);
    }

    _key_t get_median() {
        std::priority_queue<_key_t, std::vector<_key_t>> q; // max heap

        for(int i = 0; i <= NODE_SIZE / 2; i++) {
            q.push(recs_[i].key);
        }

        for(int i = NODE_SIZE / 2 + 1; i < NODE_SIZE; i++) {
            if(q.top() > recs_[i].key) {
                q.pop();
                q.push(recs_[i].key);
            }
        }

        // now in the priority queue, store the first half of the keys
        return q.top();
    }

    char * get_child(_key_t k) {
        char fp = finger_print(k);
        for(int i = 0; i < NODE_SIZE; i++) {
            if (state_.read(i) && finger_prints_[i] == fp && recs_[i].key == k) {
                return recs_[i].val;
            }
        }

        return NULL;
    }

    bool store(_key_t k, _value_t v, _key_t &split_k, INode * &split_node_out) {
        if(state_.count() == NODE_SIZE) {
            LNode * split_node = new LNode;
            int8_t m = 0;
            state_t new_state = state_;
            
            split_k = get_median(); // get the median of current node
            split_node_out = (INode *) split_node;
            // TODO: avoid extreme cases that many duplicate keys in a node
            for(int i = 0; i < NODE_SIZE; i++) {
                if (state_.read(i) && recs_[i].key >= split_k) {
                    split_node->append(recs_[i].key, recs_[i].val);
                    new_state.unpack.bitmap = new_state.free(i);
                    m++;
                }
            }
            // do the split: modifying free space
            split_node->sibs_[0] = sibs_[state_.get_sibver()];
            uint64_t free_sibver = (state_.get_sibver() + 1) % 2;
            sibs_[free_sibver] = (char *)galc->relative(split_node);
            // flush data back into PM
            clwb(split_node, sizeof(LNode));
            clwb(&sibs_[free_sibver], sizeof(char *));
            mfence();

            new_state.unpack.sib_version = free_sibver;
            new_state.unpack.node_version += 1;
            state_.pack = new_state.pack;
            clwb(this, 8);

            if(split_k > k) {
                insert(k, (char *)v);
            } else {
                split_node->insert(k, (char *)v);
            }
            return true;
        } else {
            insert(k, (char *)v);
            return false;
        }
    }

    bool remove(_key_t k) {
        int8_t slotid = -1;
        char fp = finger_print(k);
        for(int i = 0; i < NODE_SIZE; i++) {
            if (state_.read(i) && finger_prints_[i] == fp && recs_[i].key == k) {
                slotid = i;
                break;
            }
        }

        if(slotid == -1) return false;

        state_t new_state = state_;
        new_state.unpack.bitmap = state_.free(slotid);
        new_state.unpack.node_version += 1;
        state_.pack = new_state.pack;
        clwb(this, 64);
        return true;
    }

    void update(_key_t k, _value_t v) {
        char fp = finger_print(k);
        for(int i = 0; i < NODE_SIZE; i++) {
            if (state_.read(i) && finger_prints_[i] == fp && recs_[i].key == k) {
                recs_[i].val = (char *)v;
                clwb(&recs_[i], sizeof(Record));
                return ;
            }
        }
    }

    void insert(_key_t k, char * v) {
        int8_t slotid = state_.alloc();

        finger_prints_[slotid] = finger_print(k);
        recs_[slotid] = {k, (char *) v};
        clwb(&recs_[slotid], sizeof(Record));
        mfence();

        state_t new_state = state_;
        new_state.unpack.bitmap = state_.add(slotid);
        new_state.unpack.node_version += 1;
        state_.pack = new_state.pack;
        clwb(this, 64);
    }

    // void print(string prefix) {
    //     printf("%s (%lx, %d)[", prefix.c_str(), state_.unpack.bitmap, count());

    //     for(int i = 0; i < NODE_SIZE; i++) {
    //         if(state_.read(i))
    //             printf("(%ld %ld) ", recs_[i].key, (uint64_t)recs_[i].val);
    //     }
    //     printf("]\n");
    // }

    int8_t count() {
        return state_.count();
    }

    void append(_key_t k, char * v) {
        int8_t slotid = state_.alloc();
        // append the key val into this node
        finger_prints_[slotid] = finger_print(k);
        recs_[slotid] = {k, (char *) v};

        state_.unpack.bitmap = state_.add(slotid);
    }

    static void merge(LNode * left, LNode * right, _key_t merge_key) {
        state_t new_state = left->state_;
        for(int i = 0; i < NODE_SIZE; i++) {
            if(right->state_.read(i)) {
                int8_t slotid = new_state.alloc();
                left->recs_[slotid] = right->recs_[i];
                new_state.unpack.bitmap = left->state_.add(slotid);
            }
        }
        clwb(left->recs_, sizeof(Record) * NODE_SIZE);
        mfence();
        
        // install sibling of right node to the left node
        uint64_t free_sibver = (left->state_.get_sibver() + 1) % 2;
        left->sibs_[free_sibver] = right->sibs_[right->state_.get_sibver()];
        // update the state of the left node
        new_state.unpack.node_version += 1;
        new_state.unpack.sib_version = free_sibver;
        left->state_.pack = new_state.pack;
        clwb(left, 64);

        galc->free(right);
    }
};


struct INode { // inner node, allocated on DRAM 
    char * leftmost_ptr_; // represents the leftmost child of current node
    char * sibling_ptr_;  // the sibling nodes of current 
    uint8_t count_;       // total record number in current node
    bool is_parent_of_leaf_;
    uint8_t turn_  = 0;
    uint64_t bitmap = 0;
    //char pend[3];      //pend to 16B
    
    
    //uint16_t local_version = 0;
    char finger_prints_[BUFFER_SIZE + BUFFER2_SIZE] = {0};
    //char pend_[13];      //pend to 16B
    
    
    
    Record recs_[INNER_NODE_SIZE + BUFFER2_SIZE];

    void insert(_key_t k, _value_t v) {
        uint64_t i;
        for(i = 0; i < count_; i++) {
            if(recs_[i].key > k) {
                break;
            }
        }
        // recs_[i].key <= key
        memmove(&recs_[i + 1], &recs_[i], sizeof(Record) * (count_ - i));
        recs_[i] = {k, (char *) v};
        count_ += 1;
        if (count_ > INNER_NODE_SIZE - BUFFER_SIZE) {
            bitmap  &= ~(3ULL<<(count_ - 1 + BUFFER_SIZE - INNER_NODE_SIZE)*2);
        }
    }

    INode (bool parent_of_leaf = false): leftmost_ptr_(NULL), sibling_ptr_(NULL), 
                                        count_(0), is_parent_of_leaf_(parent_of_leaf) {}
    
    ~INode() {
        if(is_parent_of_leaf_ == false) { // do not free leaf nodes
            if(leftmost_ptr_ != NULL) {
                delete (INode *)leftmost_ptr_;
            }
            for(int i = 0; i < count_; i++) {
                delete (INode *)recs_[i].val;
            }
        }
    }

    void * operator new (size_t size) { // make the allocation 64 B aligned
        void * ret;
        if(posix_memalign(&ret, 64, size) != 0) 
            exit(-1);
        return ret;
    }

    void operator delete(void * ptr) {    
        free(ptr);
    }

    void insert_hotspot(_key_t key, _value_t val) {
        char fp = finger_print2(key);
        bool flag = count_ > (INNER_NODE_SIZE - BUFFER_SIZE);
        uint8_t start = flag ? count_ - (INNER_NODE_SIZE - BUFFER_SIZE) : 0;
        uint64_t newBitmap = bitmap;
        if (_mm_popcnt_u64(bitmap) == BUFFER_SIZE - start + BUFFER2_SIZE) { // cache is full
            uint8_t i = BUFFER_SIZE - start + BUFFER2_SIZE; //Real hotspot cache size
            uint8_t index = std::max(turn_, start);         //begin index
            while (i--) {
                uint8_t state = (newBitmap >> (2*index)) & 0x3;
                if (state == 1) {
                    finger_prints_[index] = fp;
                    recs_[INNER_MAX  - BUFFER_SIZE + index] = {key, (char *) val};
                    turn_ = (index + 1) % (BUFFER_SIZE + BUFFER2_SIZE);
                    break;
                } else if (state == 2){
                    newBitmap |= 1ULL << (2*index);
                    newBitmap &= ~(1ULL << (2*index + 1));
                }
                
                index = std::max(start, (uint8_t)((index + 1) % (BUFFER2_SIZE + BUFFER_SIZE)));
            } 
         } else { // cache is not full
            for (int i = start; i < BUFFER_SIZE + BUFFER2_SIZE; ++i) {
                uint8_t state = (bitmap >> (2*i)) & 0x3;
                if (state == 0) {
                    finger_prints_[i] = fp;
                    recs_[INNER_MAX  - BUFFER_SIZE + i] = {key, (char *) val};
                    newBitmap |= 1ULL << (2*i);
                    break;
                }
            }
        }

        bitmap = newBitmap;

    }

    char * get_child(_key_t key) { // find the record whose key is the last one that is less equal to key
        uint64_t i;
        for(i = 0; i < count_; i++) {
            if(recs_[i].key > key) {
                break;
            }
        }

        if(i == 0)
            return leftmost_ptr_;
        else // recs_[i - 1].key <= key
            return recs_[i - 1].val;
    }

    char * pln_update_child(_key_t key, _value_t val) { // find the record whose key is the last one that is less equal to key
        char fp = finger_print2(key);
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE, 0, 1);
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE*2, 0, 1);
        for (int i = 0; i < BUFFER_SIZE + BUFFER2_SIZE; ++i) {
            // if ( ((bitmap >> (2*i)) & 0x3) && finger_prints_[i] == fp  ) {
            //     bitmap &= ~(3ULL<<(2*i));
            // }
            if ( ((bitmap >> (2*i)) & 0x3ULL) && finger_prints_[i] == fp && key == recs_[i + INNER_NODE_SIZE - BUFFER_SIZE].key ) {
                recs_[i + INNER_NODE_SIZE - BUFFER_SIZE].val = (char *)val;
                break ;
            } 
        }
        //bitmap &= 0ULL;

        uint64_t i;
        for(i = 0; i < count_; i++) {
            if(recs_[i].key > key) {
                break; 
            }
        }

        if(i == 0)
            return leftmost_ptr_;
        else // recs_[i - 1].key <= key
            return recs_[i - 1].val;
    }

    char * pln_delete_child(_key_t key) { // find the record whose key is the last one that is less equal to key
        char fp = finger_print2(key);
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE, 0, 1);
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE*2, 0, 1);
        for (int i = 0; i < BUFFER_SIZE + BUFFER2_SIZE; ++i) {
            // if ( ((bitmap >> (2*i)) & 0x3) && finger_prints_[i] == fp  ) {
            //     bitmap &= ~(3ULL<<(2*i));
            // }
            if ( ((bitmap >> (2*i)) & 0x3ULL) && finger_prints_[i] == fp && key == recs_[i + INNER_NODE_SIZE - BUFFER_SIZE].key ) {
                bitmap &= ~(3ULL<<(2*i));
                
                //recs_[i + INNER_NODE_SIZE - BUFFER_SIZE].val = (char *)val;
                break ;
            } 
        }
        //bitmap &= 0ULL;

        uint64_t i;
        for(i = 0; i < count_; i++) {
            if(recs_[i].key > key) {
                break; 
            }
        }

        if(i == 0)
            return leftmost_ptr_;
        else // recs_[i - 1].key <= key
            return recs_[i - 1].val;
    }

    char * pln_get_child(_key_t key, bool &hit) { // find the record from pln
        char fp = finger_print2(key);
        bool flag = count_ > (INNER_NODE_SIZE - BUFFER_SIZE);
        uint8_t start = flag ? count_ - (INNER_NODE_SIZE - BUFFER_SIZE) : 0;
        
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE, 0, 1);
        uint64_t i = start;
        for (; i < BUFFER_SIZE + BUFFER2_SIZE; ++i) {
            uint8_t state = ((bitmap >> (2*i)) & 0x3ULL);
            if ( state  && finger_prints_[i] == fp && recs_[i + INNER_MAX - BUFFER_SIZE].key == key) {
                hit = true; 
                if (state == 1) {
                    bitmap |= 1ULL << (2*i + 1);
                    bitmap &= ~(1ULL << (2*i));
                }
                return recs_[i + INNER_MAX - BUFFER_SIZE].val;
            }
        }
        
        hit = false;
        for(i = 0; i < count_; i++) {
            if(recs_[i].key > key) {
                break;
            }
        }

        if(i == 0)
            return leftmost_ptr_;
        else // recs_[i - 1].key <= key
            return recs_[i - 1].val;
    }

    bool store(_key_t k, _value_t v, _key_t & split_k, INode * & split_node) {
        if(count_ == INNER_MAX) {
            split_node = new INode(is_parent_of_leaf_);

            uint64_t m = count_ / 2;
            split_k = recs_[m].key;
            // move half records into the new node
            split_node->leftmost_ptr_ = recs_[m].val;
            split_node->count_ = count_ - m - 1;
            memcpy(&(split_node->recs_[0]), &(recs_[m + 1]), sizeof(Record) * (split_node->count_));
            count_ = m;

            // set bitmap to 0
            // bitmap &= 0U;

            // update sibling pointer
            split_node->sibling_ptr_ = sibling_ptr_;
            sibling_ptr_ = (char *) split_node;

            if(split_k > k) {
                insert(k, v);
            } else {
                split_node->insert(k, v);
            }
            return true;
        } else {
            insert(k, v);
            return false;
        }
    }

    bool remove(_key_t k) { // remove k from current node
        int pos = -1;
        for(int i = 0; i < count_; i++) {
            if(recs_[i].key == k){
                pos = i;
                break;
            }
        }
        if(pos >= 0) { // we found k in this node
            memmove(&recs_[pos], &recs_[pos + 1], sizeof(Record) * (count_ - pos - 1));
            count_ -= 1;
            return true;
        } 
        return false;
    }

    int get_lrchild(_key_t k, INode * & left, INode * & right) {
        int16_t i = 0;
        for( ; i < count_; i++) {
            if(recs_[i].key > k)
                break;
        }

        if(i == 0) {
            left = NULL;
        } else if(i == 1) {
            left = (INode *)leftmost_ptr_;
        } else {
            left = (INode *)recs_[i - 2].val;
        }

        if(i == count_) {
            right = NULL;
        } else {
            right = (INode *)recs_[i].val;
        }
        return i;
    }

    // void print(string prefix) {
    //     printf("%s[(%u) ", prefix.c_str(), count());
    //     for(int i = 0; i < count_; i++) {
    //         printf("(%ld, %ld) ", recs_[i].key, (int64_t)recs_[i].val);
    //     }
    //     printf("]\n");

    //     if(is_parent_of_leaf_ == false) {
    //         INode * child = (INode *)leftmost_ptr_;
    //         child->print(prefix + "    ");

    //         for(int i = 0; i < count_; i++) {
    //             INode * child = (INode *)recs_[i].val;
    //             child->print(prefix + "    ");
    //         }
    //     } else {
    //         LNode * child = (LNode *)leftmost_ptr_;
    //         child->print(prefix + "    ");

    //         for(int i = 0; i < count_; i++) {
    //             LNode * child = (LNode *)recs_[i].val;
    //             child->print(prefix + "    ");
    //         }
    //     }
    // }

    int8_t count() {
        return count_;
    }

    static void merge(INode * left, INode * right, _key_t merge_key) {
        left->recs_[left->count_++] = {merge_key, right->leftmost_ptr_}; 
        for(int i = 0; i < right->count_; i++) {
            left->recs_[left->count_++] = right->recs_[i];
        }
        //left->bitmap &= 0;
        free((void *)right);
    }
} __attribute__((aligned(CACHE_LINE_SIZE))) ;


class hatree3 {
    public:
        hatree3(string path, bool recover) {
            if(recover == false) {
                galc = new PMAllocator(path.c_str(), false, "hatree3");
                root_ = new INode(true);
                //entrance_ = (hatree3_entrance_t *) galc->get_root(sizeof(hatree3_entrance_t));
                //entrance_->root = galc->relative(root_);
                LNode * leaf = (LNode *)galc->get_root(sizeof(LNode));
                root_->leftmost_ptr_ = (char *)leaf;
            } else {
                galc = new PMAllocator(path.c_str(), true, "hatree3");
                LNode * leaf = (LNode *)galc->get_root(sizeof(LNode));
                root_->leftmost_ptr_ = (char *)leaf;
            }
        }

        ~hatree3() {
            std::cout << "hit cnt:" << hit_cnt_ << std::endl;
            delete root_;
            delete galc;
        }

        bool find(_key_t key, _value_t &val) {
            INode * cur = root_;
            bool hit = false;
            while(!cur->is_parent_of_leaf_) { // no prefetch here
                char * child_ptr = cur->get_child(key);
                cur = (INode *)child_ptr;
            }
            //cout << "test1" << endl;
            LNode* leaf = (LNode * )cur->pln_get_child(key, hit);

            if (hit) {
                hit_cnt_++;
                val = (_value_t) leaf;
            } else {
                val = (_value_t) leaf->get_child(key);
            }
            
            //cout << "test2" << endl;
            if (search_cnt_ == SEARCH_K) {
                if (!hit) {
                    cur->insert_hotspot(key, val);
                } 
                search_cnt_ = 0;
            } else {
                search_cnt_++;
            }
            // search_cnt_++;
            // if (search_cnt_ % SEARCH_K == 0 && !hit) {
            //     
            // }
            
            /*
            LNode * leaf = (LNode * )cur->get_child(key);

            val = (_value_t) leaf->get_child(key);
            */

            return true;
        }

        void insert(_key_t key, _value_t val) {
            _key_t split_k;
            INode * split_node;
            bool splitIf = insert_recursive(root_, key, val, split_k, split_node, false);

            if(splitIf) {
                INode *new_root = new INode(false);
                new_root->leftmost_ptr_ = (char *)root_;
                new_root->recs_[0].val = (char *)split_node;
                new_root->recs_[0].key = split_k;
                new_root->count_ = 1;
                root_ = new_root;
            }
        }
        
        bool update(_key_t key, _value_t value) { // TODO: not implemented
            INode * cur = root_;
            while(!cur->is_parent_of_leaf_) { // no prefetch here
                char * child_ptr = cur->get_child(key);
                cur = (INode *)child_ptr;
            }

            LNode * leaf = nullptr;
            if ( 1) {
                leaf = (LNode * )cur->pln_update_child(key, value);
            } else {
                leaf = (LNode * )cur->get_child(key);
            }
            

            leaf->update(key, value);

            return true;
        }

        bool remove(_key_t key) { 
            INode * child = (INode *) root_->get_child(key);

            bool shouldMrg = remove_recursive(child, key, root_->is_parent_of_leaf_);

            if(shouldMrg) {
                INode *leftsib = NULL, *rightsib = NULL;
                int pos = root_->get_lrchild(key, leftsib, rightsib);

                int8_t lcnt, cnt, rcnt;
                if(root_->is_parent_of_leaf_) {
                    lcnt = leftsib == NULL ? 0 : ((LNode *) leftsib)->count();
                    cnt = ((LNode *) child)->count();
                    rcnt = rightsib == NULL ? 0 : ((LNode *) rightsib)->count();
                } else {
                    lcnt = leftsib == NULL ? 0 : leftsib->count();
                    cnt = child->count();
                    rcnt = rightsib == NULL ? 0 : rightsib->count();
                }

                if(leftsib != NULL && lcnt + cnt < INNER_MAX) {
                    // merge with left node
                    _key_t merge_key = root_->recs_[pos - 1].key;
                    root_->remove(merge_key);
                    if(root_->is_parent_of_leaf_)
                        LNode::merge((LNode *)leftsib, (LNode *)child, merge_key);
                    else
                        INode::merge(leftsib, child, merge_key);
                } 
                else if (rightsib != NULL && cnt + rcnt < INNER_MAX) {
                    // merge with right node
                    _key_t merge_key = root_->recs_[pos].key;
                    root_->remove(merge_key);
                    if(root_->is_parent_of_leaf_)
                        LNode::merge((LNode *)child, (LNode *)rightsib, merge_key);
                    else 
                        INode::merge(child, rightsib, merge_key);

                }

                if(root_->count() == 0 && root_->is_parent_of_leaf_ == false) { // the root_ is empty
                    INode * old_root = root_;
                    root_ = (INode *)root_->leftmost_ptr_;
                    
                    old_root->leftmost_ptr_ = NULL;

                    delete old_root;
                }
            }

            return false;
        }

        // void printAll() {
        //     root_->print(string(""));
        // }

    private:
        bool insert_recursive(INode * n, _key_t k, _value_t v, _key_t &split_k, INode * &split_node, bool insert_into_leaf) {
            if(insert_into_leaf) {
                return ((LNode *)n)->store(k, v, split_k, split_node);
            } else {
                INode * child = (INode *)n->get_child(k);
                
                _key_t split_k_child;
                INode * split_node_child;
                bool splitIf = insert_recursive(child, k, v, split_k_child, split_node_child, n->is_parent_of_leaf_);

                if(splitIf) { 
                    return n->store(split_k_child, (_value_t)split_node_child, split_k, split_node);
                } 
                return false;
            }
        }

        bool remove_recursive(INode * n, _key_t k, bool delete_from_leaf) {
            if(delete_from_leaf) {
                ((LNode *)n)->remove(k);
                return ((LNode *)n)->count() <= LEAF_UNDERFLOW_CARD;
            }
            else {
                INode * child = nullptr; 
                if ((INode *)n->is_parent_of_leaf_ ) {
                    child = (INode *) n->pln_delete_child(k);
                } else {
                    child = (INode *) n->get_child(k);
                }

                bool shouldMrg = remove_recursive(child, k, n->is_parent_of_leaf_);

                if(shouldMrg) {
                    INode *leftsib = NULL, *rightsib = NULL;
                    int pos = n->get_lrchild(k, leftsib, rightsib); // pos > 0 or left_sib == NULL
                    int8_t lcnt, cnt, rcnt;
                    if(n->is_parent_of_leaf_) {
                        lcnt = leftsib == NULL ? 0 : ((LNode *) leftsib)->count();
                        cnt = ((LNode *) child)->count();
                        rcnt = rightsib == NULL ? 0 : ((LNode *) rightsib)->count();
                    } else {
                        lcnt = leftsib == NULL ? 0 : leftsib->count();
                        cnt = child->count();
                        rcnt = rightsib == NULL ? 0 : rightsib->count();
                    }
                    if(leftsib != NULL && lcnt + cnt < INNER_MAX) {
                        // merge with left node
                        _key_t merge_key = n->recs_[pos - 1].key;
                        n->remove(merge_key);
                        if(n->is_parent_of_leaf_)
                            LNode::merge((LNode *)leftsib, (LNode *)child, merge_key);
                        else
                            INode::merge(leftsib, child, merge_key);
                        
                        return n->count() <= INNER_UNDERFLOW_CARD;
                    } else if (rightsib != NULL && cnt + rcnt < INNER_MAX) {
                        // merge with right node
                        _key_t merge_key = n->recs_[pos].key;
                        n->remove(merge_key);
                        if(n->is_parent_of_leaf_)
                            LNode::merge((LNode *)child, (LNode *)rightsib, merge_key);
                        else 
                            INode::merge(child, rightsib, merge_key);
                        
                        return n->count() <= INNER_UNDERFLOW_CARD;
                    }
                }
                return false;
            }
        }

    private:
        INode * root_;
        //uint16_t  global_version_ = 0;
        uint8_t  search_cnt_ = 0;
        uint32_t   hit_cnt_ = 0;
};

} // namespace hatree3

