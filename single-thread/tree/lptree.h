/*
    Learned Persistent B+tree (LPTree)
    Copyright(c) Luo Yongping All rights reserved!
*/

#ifndef __LPTREE_H__
#define __LPTREE_H__

#include <string>
#include <queue>
#include <vector>
#include <thread>
#include <atomic>

#include "pmallocator.h"
#include "pradix.h"

using std::string;
using std::vector;

namespace lptree {
    struct wbtree_entrance_t{
        void * root;
    };

    struct LNode {        // leaf nodes with fingerprints
        static constexpr int MAX_SLOTS = 28; // Node size is 512 Bytes
        struct state_t {  // slots occupancy info
            uint32_t bitmap;
            inline int8_t alloc() {
                for(int i = 0; i < 4; i++) {
                    unsigned char c = (bitmap >> (8 * i)) & 0xff;
                    if(c != 0xff) {
                        for(int j = 0; j < 8; j++) {
                            if ((c & (0x01 << j)) == 0) return 8 * i + j;
                        }
                    }
                }
                return MAX_SLOTS;
            }
            inline size_t count() {return __builtin_popcount(bitmap & 0x7fffffff);}
            inline bool read(int8_t idx) {return (bitmap & (0x1 << idx)) > 0;}
            inline uint32_t add(int8_t idx) {return bitmap + (0x1 << idx);}
            inline uint32_t clear(int8_t idx) {return bitmap - (0x1 << idx);}
            inline size_t get_version() {return (bitmap & 0x80000000) > 0 ? 1 : 0;}
            inline uint32_t flip_version() {return bitmap ^ 0x80000000;}
        };

        state_t state_;      // 4 Bytes
        char finger_prints_[MAX_SLOTS]; // 28 bytes finger print
        Record siblings_[2]; // 32 Bytes

        Record recs_[MAX_SLOTS];

        LNode() {
            state_.bitmap = 0;
            siblings_[0] = {INT64_MAX, NULL};
        }

        void * operator new(size_t size) {
            return galc->malloc(size);
        }

        _key_t get_median() { // retrieve the median key of current node
            std::priority_queue<_key_t, std::vector<_key_t>> q; // max heap

            for(int i = 0; i <= MAX_SLOTS / 2; i++) {
                q.push(recs_[i].key);
            }

            for(int i = MAX_SLOTS / 2 + 1; i < MAX_SLOTS; i++) {
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
            for(int i = 0; i < MAX_SLOTS; i++) {
                if (state_.read(i) && finger_prints_[i] == fp && recs_[i].key == k) {
                    return recs_[i].val;
                }
            }

            return NULL;
        }

        bool store(_key_t k, _value_t v, _key_t &split_k, LNode * &split_node) {
            if(state_.count() == MAX_SLOTS) {
                split_node = new LNode;
                split_k = get_median(); // get the median of current LNode

                state_t new_state = state_;
                int8_t m = 0;
                // TODO: Avoid extreme cases that many duplicate keys in one LNode
                for(int i = 0; i < MAX_SLOTS; i++) {
                    if (state_.read(i) && recs_[i].key >= split_k) {
                        split_node->insert(recs_[i].key, recs_[i].val);
                        new_state.bitmap = new_state.clear(i);
                        m++;
                    }
                }
                // record the next sibling into free sibling slot
                size_t sib_version = this->state_.get_version();
                split_node->siblings_[0] = this->siblings_[sib_version];
                this->siblings_[(sib_version + 1) % 2] = {split_k, (char *)galc->relative(split_node)};

                // update state field by a Atomic write
                new_state.bitmap = new_state.flip_version();
                this->state_.bitmap = new_state.bitmap;
                clwb(this, 64);
                mfence();

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
            for(int i = 0; i < MAX_SLOTS; i++) {
                if (state_.read(i) && finger_prints_[i] == fp && recs_[i].key == k) {
                    slotid = i;
                    break;
                }
            }

            if(slotid == -1) return false;

            state_.bitmap = state_.clear(slotid);
            clwb(this, 64);
            
            return state_.count() < MAX_SLOTS / 3; // the node is under flow
        }

        void update(_key_t k, _value_t v) {
            char fp = finger_print(k);
            for(int i = 0; i < MAX_SLOTS; i++) {
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

            state_.bitmap = state_.add(slotid);
            clwb(this, 64);
        }

        void print(string prefix) {
            printf("%s (%u)[", prefix.c_str(), state_.bitmap);

            for(int i = 0; i < MAX_SLOTS; i++) {
                if(state_.read(i))
                    printf("(%ld %ld) ", recs_[i].key, (uint64_t)recs_[i].val);
            }
            printf("]\n");
        }

        void append(_key_t k, char * v) {
            int8_t slotid = state_.alloc();
            // append the key val into this LNode
            finger_prints_[slotid] = finger_print(k);
            recs_[slotid] = {k, (char *) v};

            state_.bitmap = state_.add(slotid);
        }

        inline Record get_sibling() {
            return siblings_[state_.get_version()];
        }

        static void merge(LNode * left, LNode * right, _key_t merge_key) {
            for(int i = 0; i < MAX_SLOTS; i++) {
                if(right->state_.read(i))
                    left->append(right->recs_[i].key, right->recs_[i].val);
            }
            clwb((void *) ((uint64_t)left + 64), sizeof(LNode) - 64);
            size_t left_free_version = (left->state_.get_version() + 1) % 2;
            left->siblings_[left_free_version] = right->siblings_[right->state_.get_version()];
            
            // persist the header
            mfence();
            clwb((void *)left, 64);
            mfence();
            
            // free the right LNode
            galc->free(right);
        }
    };

    struct lptree_entrance_t {
        char * head_node;
    };

    class LPTree {
    private:
        PRadix * ist_; // inner search tree
        LNode * head_; // the first leaf node
        lptree_entrance_t * entrance_;

        // statistic values, need not to be strictly thread-safe.
        std::atomic_bool in_rebuild;
        uint64_t access_num;
        double average_length;
        
        static constexpr double LAZY_THRESHOLD  = 1.8;
        static constexpr int MINIMUM_ACCESS_NUM = 200;
        typedef LPTree SelfType;

    public:
        LPTree(string filename, bool recover, string id = "lptree") {
            if(recover == false) {
                galc = new PMAllocator(filename.c_str(), false, id.c_str());
                entrance_ = (lptree_entrance_t *) galc->get_root(sizeof(lptree_entrance_t));

                // build an empty LPTree: only one Leaf nodes
                head_ = new LNode();
                // store important infos into entrance_ (the root entry of this LPTree)
                entrance_->head_node = (char *) galc->relative(head_);
                clwb(entrance_, sizeof(lptree_entrance_t));
                mfence();

                vector<Record> tmp = {Record((_key_t)0, (char *)head_)};
                ist_ = new PRadix(tmp);
            }

            in_rebuild = false;
            access_num = 0;
            average_length = 1;
        }

        ~LPTree() {
            delete ist_;
            delete galc;
        }

        bool find(_key_t key, _value_t & value) {
            LNode * root = (LNode *) ist_->find_lower(key);

            int length = 1;
            Record sib = root->get_sibling();
            while(key >= sib.key) {
                root = (LNode *) galc->absolute(sib.val);
                sib = root->get_sibling();
                length += 1;
            }
            average_length = (((double) average_length * access_num + length) / (access_num + 1));

            value = (_value_t) root->get_child(key);

            return (char *) value != NULL;
        }

        void insert(_key_t key, _value_t value) {
            if(key == 2101)
                int a = 0;
            LNode * root = (LNode *) ist_->find_lower(key);

            int length = 1;
            Record sib = root->get_sibling();
            while(key >= sib.key) {
                root = (LNode *) galc->absolute(sib.val);
                sib = root->get_sibling();
                length += 1;
            }
            average_length = (((double) average_length * access_num + length) / (++access_num));

            _key_t split_k;
            LNode * split_n;
            bool split = root->store(key, value, split_k, split_n);
            if(split) {
                bool success = ist_->insert(split_k, (char *)split_n); // this may fail as currently the inner search tree can not grow bigger
            }

            if(average_length >= LPTree::LAZY_THRESHOLD && access_num > MINIMUM_ACCESS_NUM) {
                bool rebuilding = false;
                if(in_rebuild.compare_exchange_weak(rebuilding, true, std::memory_order_relaxed)) {
                    std::thread rebuild_thread(&SelfType::rebuild, this);
                    rebuild_thread.detach();
                }
            }

            return ;
        }

        void update(_key_t key, _value_t value) {
            LNode * root = (LNode *) ist_->find_lower(key);

            int length = 1;
            Record sib = root->get_sibling();
            while(key >= sib.key) {
                root = (LNode *) galc->absolute(sib.val);
                sib = root->get_sibling();
                length += 1;
            }
            average_length = (((double) average_length * access_num + length) / (access_num + 1));

            root->update(key, value);

            return ;
        }

        bool remove(_key_t key) {
            LNode * root = (LNode *) ist_->find_lower(key);

            int length = 1;
            Record sib = root->get_sibling();
            while(key >= sib.key) {
                root = (LNode *) galc->absolute(sib.val);
                sib = root->get_sibling();
                length += 1;
            }
            average_length = (((double) average_length * access_num + length) / (access_num + 1));

            bool underflow = root->remove(key);
            
            // TODO: merge leaf nodes and delete the router from inner search tree
            return true;
        }

        void printAll() {
            LNode * cur = head_;

            Record sib = cur->get_sibling();
            while (cur != NULL) {
                cur->print("    ");

                Record sib = cur->get_sibling();
                cur = (LNode *) galc->absolute(sib.val);
            }

            printf("\n\n");
        }

    private:
        void rebuild() { 
        /* Rebuilding is done in a background manner
         * The rebuilding procedure scans the leaf nodes and retrieves all the infos 
         * to build a inner search tree. Succesive operations are not blocked.
         */
            vector<Record> subroots({Record(0, (char *)head_)});

            LNode * cur = head_;
            Record sib = head_->get_sibling();
            while(sib.key < INT64_MAX) {
                subroots.push_back({sib.key, galc->absolute(sib.val)});
                cur = (LNode *) galc->absolute(sib.val);
                sib = cur->get_sibling();
            }

            // rebuild the inner search tree
            PRadix * new_ist = new PRadix(subroots);
            // swap the old inner search tree with a new one
            PRadix * old_ist = ist_;
            ist_ = new_ist; // this swap is atomic

            usleep(100); // TODO: wait for 100ms, try use a epoch-based method to corrdinate R/W with rebuilding ?
            
            // free the old inner search tree
            old_ist->clear();
            delete old_ist;

            access_num = 0;
            average_length = 1;
            in_rebuild = false; // only one thread in this context
            return ;
        }
    };

} // lptree

#endif // __LPTREE_H__