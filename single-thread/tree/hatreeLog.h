#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <cstdio>
#include <queue>
#include <algorithm>

#include "pmallocator.h"
#include "flush.h"

#define eADR
 
using std::string;
using std::cout;
using std::endl;

namespace hatreeLog {

const uint16_t INNER_NODE_SIZE = 30;
const uint16_t SEARCH_K = 3;                     // Submit hotspot every 5 queries
const uint16_t NODE_SIZE = 14;                  // INode size 256B
const int INNER_UNDERFLOW_CARD = INNER_NODE_SIZE / 2;
const int LEAF_UNDERFLOW_CARD = NODE_SIZE / 3;  // Underflow limit
const uint16_t INNER_MAX = INNER_NODE_SIZE ; // one quarter of the slots are reserved as read buffer
const uint16_t BUFFER_SIZE = 8;
const uint64_t LOG_ENTRY_SIZE = LOADSCALE * 1024 * 128 - 1;
const uint16_t XPLINE_SIZE = 256;

struct LNode;
struct INode;

struct log_area {
    uint64_t index;
    uint64_t next;
    Record recs[LOG_ENTRY_SIZE];

    log_area() {
        index = 0;
    }

    void add_entry(Record e) {
        recs[index] = e;
    #ifndef eADR
        clwb(&(recs[index]), sizeof(Record));
        mfence();
    #endif
        index = (index + 1) % LOG_ENTRY_SIZE;
    #ifndef eADR
        clwb(&index, 8);
    #endif
    }
}__attribute__((aligned(XPLINE_SIZE)));


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
        #ifndef eADR
            clwb(split_node, sizeof(LNode));
            clwb(&sibs_[free_sibver], sizeof(char *));
            mfence();
        #endif

            new_state.unpack.sib_version = free_sibver;
            new_state.unpack.node_version += 1;
            state_.pack = new_state.pack;
        #ifndef eADR
            clwb(this, 8);
        #endif

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
    #ifndef eADR
        clwb(this, 64);
    #endif
        return true;
    }

    void update(_key_t k, _value_t v) {
        char fp = finger_print(k);
        for(int i = 0; i < NODE_SIZE; i++) {
            if (state_.read(i) && finger_prints_[i] == fp && recs_[i].key == k) {
                recs_[i].val = (char *)v;
            #ifndef eADR
                clwb(&recs_[i], sizeof(Record));
            #endif
                return ;
            }
        }
    }

    void insert(_key_t k, char * v) {
        int8_t slotid = state_.alloc();

        finger_prints_[slotid] = finger_print(k);
        recs_[slotid] = {k, (char *) v};
    #ifndef eADR
        clwb(&recs_[slotid], sizeof(Record));
        mfence();
    #endif
        state_t new_state = state_;
        new_state.unpack.bitmap = state_.add(slotid);
        new_state.unpack.node_version += 1;
        state_.pack = new_state.pack;
    #ifndef eADR
        clwb(this, 64);
    #endif
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
    #ifndef eADR
        clwb(left->recs_, sizeof(Record) * NODE_SIZE);
        mfence();
    #endif
        // install sibling of right node to the left node
        uint64_t free_sibver = (left->state_.get_sibver() + 1) % 2;
        left->sibs_[free_sibver] = right->sibs_[right->state_.get_sibver()];
        // update the state of the left node
        new_state.unpack.node_version += 1;
        new_state.unpack.sib_version = free_sibver;
        left->state_.pack = new_state.pack;
    #ifndef eADR
        clwb(left, 64);
    #endif
        galc->free(right);
    }
}__attribute__((aligned(XPLINE_SIZE)));


struct INode { // inner node, allocated on DRAM 
    char * leftmost_ptr_; // represents the leftmost child of current node
    char * sibling_ptr_;  // the sibling nodes of current 
    uint8_t count_;       // total record number in current node
    bool is_parent_of_leaf_;
    uint8_t turn_  = 0;
    //uint8_t search_cnt_ = 0;
    //uint8_t ocur_;
    uint8_t visit_ = 0;
    uint8_t bitmap_ = 0;
    uint8_t dirty_ = 0;
    char pend[2];      //pend to 16B
    
    
    //uint16_t local_version = 0;
    char finger_prints_[BUFFER_SIZE] = {0};
    //char pend[16];      //pend to 16B
    
    
    
    Record recs_[INNER_NODE_SIZE];

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
            size_t index = count_ - 1 + BUFFER_SIZE - INNER_NODE_SIZE; 
            if (get(bitmap_, index) && get(dirty_, index)) {
                _key_t key = recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].key;
                _value_t value = (_value_t)recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].val;
                LNode* leaf = (LNode*)get_child(key);
                leaf->update(key, value);
            }
            bitmap_  &= ~(1ULL<<(count_ - 1 + BUFFER_SIZE - INNER_NODE_SIZE));
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
        if (count_ == INNER_NODE_SIZE)  return ;
        char fp = finger_print2(key);
        bool flag = count_ > (INNER_NODE_SIZE - BUFFER_SIZE);
        uint8_t start = flag ? count_ - (INNER_NODE_SIZE - BUFFER_SIZE) : 0;
        
        if (_mm_popcnt_u32(bitmap_) == BUFFER_SIZE - start) {
            uint8_t i = BUFFER_SIZE - start;
            uint8_t index = std::max(turn_, start);
            while (i--) {
                if (get(visit_, index) == false) {
                    if (get(dirty_, index)) {
                        _key_t k = recs_[INNER_NODE_SIZE - BUFFER_SIZE + index].key;
                        _value_t v = (_value_t)recs_[INNER_NODE_SIZE - BUFFER_SIZE + index].val;
                        LNode* leaf = (LNode*) get_child(k);
                        leaf->update(k, v);
                    }
                    finger_prints_[index] = fp;
                    recs_[INNER_MAX  - BUFFER_SIZE + index] = {key, (char *) val};
                    turn_ = (index + 1) % BUFFER_SIZE;
                    break;
                } else {
                    unset(visit_, index);
                }
                index = std::max(start, (uint8_t)((index + 1) % BUFFER_SIZE));
            }
        } else {
            for (int i = start; i < BUFFER_SIZE; ++i) {
                if (get(bitmap_, i) == false) {
                    finger_prints_[i] = fp;
                    recs_[INNER_MAX  - BUFFER_SIZE + i] = {key, (char *) val};
                    set(bitmap_, i);
                    break;
                }
            }
        }
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


    char * pln_delete_child(_key_t key) { // find the record whose key is the last one that is less equal to key
        char fp = finger_print2(key);
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE, 0, 1);
        __builtin_prefetch((&count_)+CACHE_LINE_SIZE*2, 0, 1);
        for (int i = 0; i < BUFFER_SIZE; ++i) {
            if ( get(bitmap_, i) && finger_prints_[i] == fp && key == recs_[i + INNER_NODE_SIZE - BUFFER_SIZE].key ) {
                unset(bitmap_, i);
                break ;
            } 
        }

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

    void search_hotspot(_key_t key, bool &hit, size_t &index) {
        char fp = finger_print2(key);
        bool flag = count_ > (INNER_NODE_SIZE - BUFFER_SIZE);
        uint8_t start = flag ? count_ - (INNER_NODE_SIZE - BUFFER_SIZE) : 0;
        
        size_t i = start;
        for (; i < BUFFER_SIZE; ++i) {
            auto state = get(bitmap_, i);
            if ( state  && finger_prints_[i] == fp && recs_[i + INNER_MAX - BUFFER_SIZE].key == key) {
                hit = true; 
                index = i;
                if (get(visit_, i) == false) {
                    set(visit_, i);
                }
                return ;
            }
        }
        
        hit = false;
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
        // flush dirty entries
        for (size_t i = 0; i < BUFFER_SIZE; ++i) {
            if (get(left->bitmap_, i) && get(left->dirty_, i)) {
                _key_t key = left->recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].key;
                _value_t value = (_value_t)left->recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].val;
                LNode* leaf = (LNode*)left->get_child(key);
                leaf->update(key, value);
            }
        }
        for (size_t i = 0; i < BUFFER_SIZE; ++i) {
            if (get(right->bitmap_, i) && get(right->dirty_, i)) {
                _key_t key = right->recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].key;
                _value_t value = (_value_t)right->recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].val;
                LNode* leaf = (LNode*)right->get_child(key);
                leaf->update(key, value);
            }
        }

        left->recs_[left->count_++] = {merge_key, right->leftmost_ptr_}; 
        for(int i = 0; i < right->count_; i++) {
            left->recs_[left->count_++] = right->recs_[i];
        }
        //left->bitmap &= 0;
        free((void *)right);
    }
} __attribute__((aligned(CACHE_LINE_SIZE))) ;


class hatreeLog {
    public:
        hatreeLog(string path, bool recover) {
            if(recover == false) {
                galc = new PMAllocator(path.c_str(), false, "hatreeLog");
                root_ = new INode(true);
                //entrance_ = (hatreeLog_entrance_t *) galc->get_root(sizeof(hatreeLog_entrance_t));
                //entrance_->root = galc->relative(root_);
                LNode * leaf = (LNode *)galc->get_root(sizeof(LNode));
                root_->leftmost_ptr_ = (char *)leaf;
                
                //init log_
                log_ = (log_area *) galc->malloc(sizeof(log_area));
                log_->index = 0;
            #ifndef eADR
                clwb(log_, 8);
            #endif
                cout << "the size of log area is " << sizeof(log_area) / 1024.0 /1024 <<  "MB" << endl;
            } else {
                galc = new PMAllocator(path.c_str(), true, "hatreeLog");
                LNode * leaf = (LNode *)galc->get_root(sizeof(LNode));
                root_->leftmost_ptr_ = (char *)leaf;
            }
        }

        ~hatreeLog() {
            //std::cout << "hit cnt:" << hit_cnt_ << std::endl;
            delete root_;
            delete galc;
        }

        bool find(_key_t key, _value_t &val) {
            search_cnt_ = (search_cnt_ + 1) % SEARCH_K;

            //get the pln node
            INode * pln = root_;
            while(!pln->is_parent_of_leaf_) { // no prefetch here
                char * child_ptr = pln->get_child(key);
                pln = (INode *)child_ptr;
            }

            //prefetch the pln node
            __builtin_prefetch(&(pln->count_), 0, 1);
            __builtin_prefetch(&(pln->count_)+CACHE_LINE_SIZE, 0, 1);
            __builtin_prefetch(&(pln->count_)+CACHE_LINE_SIZE*2, 0, 1);

            //search pln cache
            size_t index = 0;
            bool hit = false;
            pln->search_hotspot(key, hit, index);


            if (hit) {
                val = (_value_t)pln->recs_[index + INNER_NODE_SIZE - BUFFER_SIZE].val;
            } else {
                LNode* leaf = (LNode*) pln->get_child(key);
                val = (_value_t) leaf->get_child(key);
                if (search_cnt_ == 0) {
                    pln->insert_hotspot(key, val);
                }
            }
            
            if (val) 
                return true;
            else 
                return false;
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
            search_cnt_ = (search_cnt_ + 1) % SEARCH_K;
            
            INode * pln = root_;
            while(!pln->is_parent_of_leaf_) { // no prefetch here
                char * child_ptr = pln->get_child(key);
                pln = (INode *)child_ptr;
            }

            //prefetch the pln node
            __builtin_prefetch(&(pln->count_), 0, 1);
            __builtin_prefetch(&(pln->count_)+CACHE_LINE_SIZE, 0, 1);
            __builtin_prefetch(&(pln->count_)+CACHE_LINE_SIZE*2, 0, 1);

            //search pln cache
            size_t index = 0;
            bool hit = false;
            pln->search_hotspot(key, hit, index);

            if (hit) {
                pln->recs_[index + INNER_NODE_SIZE - BUFFER_SIZE].val = (char*)value;
                if (update_with_log_) {
                    if (get(pln->dirty_, index) == false) {
                        set(pln->dirty_, index);
                    }
                    add_log(key, value);
                } else {
                    // update the leaf node
                    LNode * leaf = (LNode * )pln->get_child(key);
                    leaf->update(key, value);
                }
            } else {
                // update the leaf node
                LNode * leaf = (LNode * )pln->get_child(key);
                leaf->update(key, value);
                if (search_cnt_ == 0) {
                    pln->insert_hotspot(key, value);
                }
            }

            if (!update_with_log_) {
                flush_one_node();
            }
            

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
        void add_log(_key_t key, _value_t value) {
            log_->add_entry({key, (char*)value});
            if (log_->index == LOG_ENTRY_SIZE) {
                update_with_log_ = false;
                cur_node_ = root_;
                while(!cur_node_->is_parent_of_leaf_) { // no prefetch here
                    cur_node_ = (INode*)cur_node_->leftmost_ptr_;
                }
            }
        }

        void flush_one_node() {
            if (cur_node_ == nullptr) {
                update_with_log_ = true;
                log_->index = 0;
            #ifndef eADR
                clwb(&(log->index), 8);
            #endif
            } else {
                for (size_t i = 0; i < BUFFER_SIZE; ++i) {
                    if (get(cur_node_->bitmap_, i) && get(cur_node_->dirty_, i)) {
                        _key_t key = cur_node_->recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].key;
                        _value_t value = (_value_t)cur_node_->recs_[INNER_NODE_SIZE - BUFFER_SIZE + i].val;
                        LNode* leaf = (LNode*) cur_node_->get_child(key);
                        leaf->update(key, value);
                    }
                }
            }
        }

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
                        else {
                            if (!update_with_log_ && cur_node_ == child) {
                                flush_one_node();
                            }
                            INode::merge(leftsib, child, merge_key);
                        }

                        return n->count() <= INNER_UNDERFLOW_CARD;
                    } else if (rightsib != NULL && cnt + rcnt < INNER_MAX) {
                        // merge with right node
                        _key_t merge_key = n->recs_[pos].key;
                        n->remove(merge_key);
                        if(n->is_parent_of_leaf_)
                            LNode::merge((LNode *)child, (LNode *)rightsib, merge_key);
                        else  {
                            if (!update_with_log_ && cur_node_ == rightsib) {
                                flush_one_node();
                            }
                            INode::merge(child, rightsib, merge_key);
                        }
                        
                        return n->count() <= INNER_UNDERFLOW_CARD;
                    }
                }
                return false;
            }
        }

    private:
        INode * root_;
        INode * cur_node_;
        log_area *log_;
        bool update_with_log_ = true;
        uint8_t  search_cnt_ = 0;
        uint32_t   hit_cnt_ = 0;
};

} // namespace hatreeLog
