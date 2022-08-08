/*
    An implementation of FP-Tree (2016 SIGMOD FPTree: A Hybrid SCM-DRAM Persistent 
    and Concurrent B-Tree for Storage Class Memory)
    Copyright(c) Luo Yongping All rights reserved!
*/

#ifndef __FPTREE__
#define __FPTREE__

#include <cstdint>
#include <cstring>
#include <string>
#include <cstdio>
#include <queue>

#include "pmallocator.h"
#include "flush.h"

namespace fptree {
using std::string; 

const int INNER_NODE_SIZE = 30;
const int INNER_MAX = INNER_NODE_SIZE;
const int NODE_SIZE = 14    ; // Node size 256B
const int INNER_UNDERFLOW_CARD = INNER_NODE_SIZE / 2;
const int LEAF_UNDERFLOW_CARD = NODE_SIZE / 3;  // Underflow limit


static char finger_print(_key_t k) { 
// using fmix function from murmur hash
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return (k & 0xff);
}

struct log_entry{
    unsigned int size;
    unsigned char type;
    void *addr;
    char data[LOG_DATA_SIZE];
} ;
struct log_area{
    log_entry *next_offset;
    char log_data[LOG_AREA_SIZE];
};

struct fptree_entrance_t{
    void * root;
    void * log_;
};

class fptree;
static inline void AddLog(fptree * ftree, void *addr, unsigned int size, unsigned char type);

/* We abandoned virtual Node class as it takes 8 bytes extra space to store the virtual pointer */
struct LNode { // leaf node of fptree
    struct state_t { // slots occupancy info
        uint64_t bitmap;
        inline int8_t alloc() {
            for(int i = 0; i < 7; i++) {
                unsigned char c = (bitmap >> (8 * i)) & 0xff;
                if(c != 0xff) {
                    for(int j = 0; j < 8; j++) {
                        if ((c & (0x01 << j)) == 0) return 8 * i + j;
                    }
                }
            }
            return NODE_SIZE;
        }
        inline bool read(int8_t idx) {return (bitmap & (0x1LL << idx)) > 0;}
        inline uint64_t add(int8_t idx) {return bitmap + (0x1LL << idx);}
        inline uint64_t clear(int8_t idx) {return bitmap - (0x1LL << idx);}
    };

    // 1 + 1 + 4 + 4 + 14 + 8 = 32 Bytes
    uint8_t isleaf_;
    uint8_t count_;
    //char dummy[2];
    state_t state_;
    char finger_prints_[NODE_SIZE]; // 14 bytes finger print
    LNode * sibling_ptr_;

    Record recs_[NODE_SIZE];

    LNode(): isleaf_(1), count_(0), sibling_ptr_(NULL) {
        state_.bitmap = 0;
    }

    void * operator new(size_t size) {
        return galc->malloc(size);
    }

    uint8_t count() {
        return count_;
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
            if (state_.read(i) && finger_print(k) == fp && recs_[i].key == k) {
                return recs_[i].val;
            }
        }

        return NULL;
    }

    bool store(fptree * tree, _key_t k, _value_t v, _key_t &split_k, LNode * &split_node) {
        if(count_ == NODE_SIZE) {
            //std::cout << "3" << endl;
            AddLog(tree, this, sizeof(LNode), LE_DATA);

            //std::cout << "3" << endl;
            split_node = new LNode;
            split_k = get_median(); // get the median of current node

            int8_t m = 0;
            state_t new_state = state_;
            
            //std::cout << "3" << endl;
            // TODO: Avoid extreme cases that many duplicate keys in one node
            for(int i = 0; i < NODE_SIZE; i++) {
                if (state_.read(i) && recs_[i].key >= split_k) {
                    split_node->insert(recs_[i].key, recs_[i].val);
                    new_state.bitmap = new_state.clear(i);
                    m++;
                }
            }
            //std::cout << "3" << endl;
            state_.bitmap = new_state.bitmap;
            count_ -= m; 
            split_node->sibling_ptr_ = this->sibling_ptr_;
            this->sibling_ptr_ = galc->relative(split_node);

            //std::cout << "3" << endl;
            AddLog(tree, NULL, 0, LE_COMMIT);

            if(split_k > k) {
                insert(k, (char *)v);
            } else {
                split_node->insert(k, (char *)v);
            }
            //std::cout << "3" << endl;
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
            if (state_.read(i) && finger_print(k) == fp && recs_[i].key == k) {
                slotid = i;
                break;
            }
        }

        if(slotid == -1) return false;

        state_.bitmap = state_.clear(slotid);
        count_ -= 1;
        clwb(this, 64);
        return true;
    }

    void update(_key_t k, _value_t v) {
        char fp = finger_print(k);
        for(int i = 0; i < NODE_SIZE; i++) {
            if (state_.read(i) && finger_print(k) == fp && recs_[i].key == k) {
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

        state_.bitmap = state_.add(slotid);
        count_ += 1;
        clwb(this, 64);
    }

    void print(string prefix) {
        //printf("%s (%x, %d)[", prefix.c_str(), state_.bitmap, count_);

        for(int i = 0; i < NODE_SIZE; i++) {
            if(state_.read(i))
                printf("(%ld %ld) ", recs_[i].key, (uint64_t)recs_[i].val);
        }
        printf("]\n");
    }

    void append(_key_t k, char * v) {
        int8_t slotid = state_.alloc();
        // append the key val into this node
        finger_prints_[slotid] = finger_print(k);
        recs_[slotid] = {k, (char *) v};

        state_.bitmap = state_.add(slotid);
    }

    static void merge(fptree * tree, LNode * left, LNode * right, _key_t merge_key) {
        AddLog(tree, left, sizeof(LNode), LE_DATA);
        
        for(int i = 0; i < NODE_SIZE; i++) {
            if(right->state_.read(i))
                left->append(right->recs_[i].key, right->recs_[i].val);
        }

        left->count_ += right->count_;
        left->sibling_ptr_ = right->sibling_ptr_;
        clwb(left, sizeof(LNode));

        AddLog(tree, NULL, 0, LE_COMMIT);

        galc->free(right);
    }
};


struct INode { // inner node, allocated on DRAM
    char * leftmost_ptr_; // represents the leftmost child of current node
    char * sibling_ptr_;  // the sibling nodes of current 
    uint8_t count_;       // total record number in current node
    bool is_parent_of_leaf_;
    //char finger_prints_[NODE_SIZE / 2];
    char pend_[14];
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

    void print(string prefix) {
        printf("%s[(%u) ", prefix.c_str(), count());
        for(int i = 0; i < count_; i++) {
            printf("(%ld, %ld) ", recs_[i].key, (int64_t)recs_[i].val);
        }
        printf("]\n");

        if(is_parent_of_leaf_ == false) {
            INode * child = (INode *)leftmost_ptr_;
            child->print(prefix + "    ");

            for(int i = 0; i < count_; i++) {
                INode * child = (INode *)recs_[i].val;
                child->print(prefix + "    ");
            }
        } else {
            LNode * child = (LNode *)leftmost_ptr_;
            child->print(prefix + "    ");

            for(int i = 0; i < count_; i++) {
                LNode * child = (LNode *)recs_[i].val;
                child->print(prefix + "    ");
            }
        }
    }

    int8_t count() {
        return count_;
    }

    static void merge(INode * left, INode * right, _key_t merge_key) {
        left->recs_[left->count_++] = {merge_key, right->leftmost_ptr_}; 
        for(int i = 0; i < right->count_; i++) {
            left->recs_[left->count_++] = right->recs_[i];
        }
        free((void *)right);
    }
};



class fptree {
public:
    fptree(string path, bool recover, string id = "fptree") {
        if(recover == false) {
            galc = new PMAllocator(path.c_str(), false, "btree");
            root_ = new INode(true);
            //entrance_ = (btree_entrance_t *) galc->get_root(sizeof(btree_entrance_t));
            //entrance_->root = galc->relative(root_);
            LNode * leaf = (LNode *)galc->get_root(sizeof(LNode));
            root_->leftmost_ptr_ = (char *)leaf;
            
            entrance_ = (fptree_entrance_t *) galc->get_root(sizeof(fptree_entrance_t));
            entrance_->log_ = NULL;
            clflush(entrance_, sizeof(fptree_entrance_t));
        
            log_ = (log_area *) galc->malloc(sizeof(log_area));
            log_->next_offset = (log_entry *)galc->relative(log_->log_data);
            clflush(log_, sizeof(log_entry *));
            
            
            entrance_->log_ = galc->relative(log_);
            clflush(entrance_, 16);
            
        } else {
            galc = new PMAllocator(path.c_str(), true, id.c_str());

            entrance_ = (fptree_entrance_t *) galc->get_root(sizeof(fptree_entrance_t));

            if(entrance_->root == NULL) {
                printf("The tree is empty\n");
                exit(-1);
            }

            //root_ = (INode *) galc->absolute(entrance_->root);
            log_ = (log_area *) galc->absolute(entrance_->log_);
        }
    }

    ~fptree() {
        delete galc;
        galc = nullptr;
    }

    bool find(_key_t key, _value_t &val) {
        INode * cur = root_;
        while(!cur->is_parent_of_leaf_) { // no prefetch here
            char * child_ptr = cur->get_child(key);
            cur = (INode *)child_ptr;
        }
        LNode * leaf = (LNode * )cur->get_child(key);

        val = (_value_t) leaf->get_child(key);

        return true;
    }

    void insert(_key_t key, _value_t val) {
        _key_t split_k;
        INode * split_node;
        //std::cout << "1" << std::endl;
        bool splitIf = insert_recursive(root_, key, val, split_k, split_node, false);
        //std::cout << "2" << std::endl;
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
        LNode * leaf = (LNode * )cur->get_child(key);

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
                lcnt = leftsib == NULL ? 0 : ((LNode *) leftsib)->count_;
                cnt = ((LNode *) child)->count_;
                rcnt = rightsib == NULL ? 0 : ((LNode *) rightsib)->count_;
            } else {
                lcnt = leftsib == NULL ? 0 : leftsib->count_;
                cnt = child->count();
                rcnt = rightsib == NULL ? 0 : rightsib->count_;
            }

            if(leftsib != NULL && lcnt + cnt < INNER_MAX) {
                // merge with left node
                _key_t merge_key = root_->recs_[pos - 1].key;
                root_->remove(merge_key);
                if(root_->is_parent_of_leaf_)
                    LNode::merge(this, (LNode *)leftsib, (LNode *)child, merge_key);
                else
                    INode::merge(leftsib, child, merge_key);
            } 
            else if (rightsib != NULL && cnt + rcnt < INNER_MAX) {
                // merge with right node
                _key_t merge_key = root_->recs_[pos].key;
                root_->remove(merge_key);
                if(root_->is_parent_of_leaf_)
                    LNode::merge(this, (LNode *)child, (LNode *)rightsib, merge_key);
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


    bool scan(_key_t key, int count, Record * buff) {
        return false;
    }
public:
    bool insert_recursive(INode * n, _key_t k, _value_t v, _key_t &split_k, INode * &split_node, bool insert_into_leaf) {
        if(insert_into_leaf) {
            return ((LNode *)n)->store(this, k, v, split_k, (LNode * &)split_node);
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
            return ((LNode *)n)->count_ <= LEAF_UNDERFLOW_CARD;
        }
        else {
            INode * child = (INode *) n->get_child(k);

            bool shouldMrg = remove_recursive(child, k, n->is_parent_of_leaf_);

            if(shouldMrg) {
                INode *leftsib = NULL, *rightsib = NULL;
                int pos = n->get_lrchild(k, leftsib, rightsib); // pos > 0 or left_sib == NULL
                int8_t lcnt, cnt, rcnt;
                if(n->is_parent_of_leaf_) {
                    lcnt = leftsib == NULL ? 0 : ((LNode *) leftsib)->count_;
                    cnt = ((LNode *) child)->count();
                    rcnt = rightsib == NULL ? 0 : ((LNode *) rightsib)->count_;
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
                        LNode::merge(this, (LNode *)leftsib, (LNode *)child, merge_key);
                    else
                        INode::merge(leftsib, child, merge_key);
                    
                    return n->count() <= INNER_UNDERFLOW_CARD;
                } else if (rightsib != NULL && cnt + rcnt < INNER_MAX) {
                    // merge with right node
                    _key_t merge_key = n->recs_[pos].key;
                    n->remove(merge_key);
                    if(n->is_parent_of_leaf_)
                        LNode::merge(this, (LNode *)child, (LNode *)rightsib, merge_key);
                    else 
                        INode::merge(child, rightsib, merge_key);
                    
                    return n->count() <= INNER_UNDERFLOW_CARD;
                }
            }
            return false;
        }
    }

    void add_log_entry(void *addr, unsigned int size, unsigned char type)
    {
        log_entry *log;
        int i, remain_size;

        remain_size = size - ((size / LOG_DATA_SIZE) * LOG_DATA_SIZE);

        if ((char *) galc->absolute(log_->next_offset) == (log_->log_data + LOG_AREA_SIZE))
            log_->next_offset = (log_entry *)galc->relative(log_->log_data);

        if (size <= LOG_DATA_SIZE) {
            log = (log_entry *)galc->absolute(log_->next_offset);
            log->size = size;
            log->type = type;
            log->addr = galc->relative(addr);
            memcpy(log->data, addr, size);

            if (type == LE_DATA)
                clflush(log, sizeof(log_entry), false);
            else
                clflush(log, sizeof(log_entry), true);

            log_->next_offset = log_->next_offset + 1;
        } else {
            void *next_addr = addr;

            for (i = 0; i < size / LOG_DATA_SIZE; i++) {
                log = (log_entry *)galc->absolute(log_->next_offset);
                log->size = LOG_DATA_SIZE;
                log->type = type;
                log->addr = galc->relative(next_addr);
                memcpy(log->data, next_addr, LOG_DATA_SIZE);

                clflush(log, sizeof(log_entry), false);

                log_->next_offset = log_->next_offset + 1;
                if ((char *)galc->absolute(log_->next_offset) == (log_->log_data + LOG_AREA_SIZE))
                    log_->next_offset = (log_entry *)galc->relative(log_->log_data);

                next_addr = (char *)next_addr + LOG_DATA_SIZE;
            }

            if (remain_size > 0) {
                log = (log_entry *)galc->absolute(log_->next_offset);
                log->size = LOG_DATA_SIZE;
                log->type = type;
                log->addr = next_addr;
                memcpy(log->data, next_addr, remain_size);

                clflush(log, sizeof(log_entry), false);
                
                log_->next_offset = log_->next_offset + 1;
            }
        }
    }

    
private:
    fptree_entrance_t *entrance_;
    INode * root_;
    log_area *log_;
}; // class fptree

void AddLog(fptree * tree, void *addr, unsigned int size, unsigned char type) {
    tree->add_log_entry(addr, size, type);
}

}; // namespace fptree

#endif // __FPTREE__