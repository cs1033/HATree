/*
    Copyright(c) SeKwonLee (github user)
    Modified by Luo Yongping. For Research Use Only.
    Tranform it into CPP style and Add complete remove op
*/

#ifndef __WBTREE_BITMAP__
#define __WBTREE_BITMAP__

#include <cstdint>
#include <cstring>
#include <string>
#include <cstdio>

#include "pmallocator.h"
#include "flush.h" 

using std::string;

namespace wbtreelog {

const int PAGE_SIZE = 256;
const int NODE_SIZE = (PAGE_SIZE - 17) / 17; // 25 bytes for header, one byte for slot[0], each Record takes 17 bytes
const int SLOT_SIZE = NODE_SIZE + 1;
const int MIN_LIVE_ENTRIES = NODE_SIZE / 2;
const int BITMAP_SIZE = NODE_SIZE + 1;

static uint64_t findEmptyBit(const uint64_t bitmap)
{
    for(int i = 0; i < (BITMAP_SIZE + 7) / 8; i++) {
        unsigned char c = (bitmap >> (56 - 8 * i)) & 0xff;
        if(c != 0xff) {
            for(int j = 0; j < 8; j++) {
                if ((c & (0x80 >> j)) == 0) return 8 * i + j;
            }
        }
    }
    return BITMAP_SIZE;
}

static inline uint64_t setBit(uint64_t bitmap, uint64_t loc) {
    return bitmap + (1UL << (63 - loc));
}

static inline uint64_t clearBit(uint64_t bitmap, uint64_t loc) {
    return bitmap - (1UL << (63 - loc));
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

struct wbtree_entrance_t{
    void * root;
    void * start_log;
};

class wbtree;

class Node {
    private:
        uint64_t bitmap;
        Node *leftmostPtr;
        uint8_t isleaf;
        char slot[SLOT_SIZE];
        struct Record entries[NODE_SIZE];
    
    public:
        Node():leftmostPtr(NULL) {
            memset(slot, 0, sizeof(slot));
            bitmap = 1;
            isleaf = 1;
        }

        void * operator new(size_t size) {
            return galc->malloc(size);
        }
        
        friend class wbtree;

        int insertone_noflush(uint64_t key, void *value) {
            bitmap = bitmap - 1;

            int loc, mid, j;
            loc = Append(key, value);
            mid = BinSearch(key);

            for (j = slot[0]; j >= mid; j--)
                slot[j + 1] = slot[j];

            slot[mid] = loc;
            slot[0] = slot[0] + 1;
            bitmap = setBit(bitmap, loc) + 1;

            return mid;
        }

        int insertone(uint64_t key, void *value) {
            int loc, mid, j;
            loc = Append(key, value);
            clflush(&(entries[loc]), sizeof(Record), false);

            mid = BinSearch(key);

            bitmap = bitmap - 1;
            mfence();

            for (j = slot[0]; j >= mid; j--)
                slot[j + 1] = slot[j];
            slot[mid] = loc;

            slot[0] = slot[0] + 1;
            mfence();
            
            bitmap = setBit(bitmap, loc) + 1;
            clflush(&bitmap, 64, true);

            return mid;
        }

        bool removeone_leaf(uint64_t key) {
            int mid, j;
            mid = BinSearch(key);
            int loc = slot[mid];

            if(entries[loc].key != key) { // not found
                return false;
            }

            bitmap = bitmap - 1;
            mfence();

            for(j = mid; j < slot[0]; j++)
                slot[j] = slot[j + 1];
            slot[0] = slot[0] - 1;
            mfence();
            
            bitmap = clearBit(bitmap, loc) + 1;
            clflush(&bitmap, 64, true);

            return slot[0] < MIN_LIVE_ENTRIES;
        }

        bool removeone_inner(uint64_t key) {
            int mid = BinSearch(key);
            if(mid > slot[0]) {
                mid = slot[0];
            } else if(entries[slot[mid]].key == key) {
                mid = mid;
            } else {
                mid = mid - 1;
            }
            int loc = slot[mid];

            bitmap = bitmap - 1;
            mfence();

            for(int j = mid; j < slot[0]; j++)
                slot[j] = slot[j + 1];
            slot[0] = slot[0] - 1;
            mfence();
            
            bitmap = clearBit(bitmap, loc) + 1;
            clflush(&bitmap, 64, true);

            return slot[0] < MIN_LIVE_ENTRIES;
        }

        int BinSearch(uint64_t key) {
            int low = 1, mid = 1;
            int high = slot[0];

            while (low <= high){
                mid = (low + high) / 2;
                if (entries[slot[mid]].key > key)
                    high = mid - 1;
                else if (entries[slot[mid]].key < key)
                    low = mid + 1;
                else
                    break;
            }

            if (low > mid) // make sure the key in the middle is exactly biggger than param key 
                mid = low;

            return mid;
        }

        int Append(uint64_t key, void *value) {
            int errval = -1;
            uint64_t index;

            index = findEmptyBit(bitmap);
            if (index == BITMAP_SIZE)
                return errval;

            entries[index].key = key;
            entries[index].val = (char *)value;
            return index;
        }
        
        void print(bool recursively, string prefix="") {
            printf("%s[%lx(%d) ", prefix.c_str(), bitmap, slot[0]);
            for(int i = 1; i <= slot[0]; i++) {
                printf("(%ld 0x%lx) ", entries[slot[i]].key, (uint64_t)entries[slot[i]].val);
            }
            printf("]\n");

            if(isleaf == 0 && recursively) {
                if(leftmostPtr != NULL) {
                    Node * child = (Node *)galc->absolute(leftmostPtr);
                    child->print(recursively, prefix + "    ");
                }

                for(int i = 1; i <= slot[0]; i++) {
                    Node * child = (Node *)galc->absolute(entries[slot[i]].val);
                    child->print(recursively, prefix + "    ");
                }
            }
        }

        void get_sibling(Node * parent, Node * &left, Node * &right) {
            if(parent == NULL) { // current node is root
                left = right = NULL;
            } 
            else {
                auto splitkey = entries[slot[1]].key;
                int loc = parent->BinSearch(splitkey);

                if (loc == 1) {
                    if(parent->entries[parent->slot[1]].key == splitkey)
                        left = (Node *)galc->absolute(parent->leftmostPtr);
                    else
                        left = NULL; 
                } else if(loc == 2) {
                    if(parent->entries[parent->slot[2]].key == splitkey)
                        left = (Node *)galc->absolute(parent->entries[parent->slot[1]].val);
                    else
                        left = (Node *)galc->absolute(parent->leftmostPtr);
                } else {
                    if(parent->entries[parent->slot[loc]].key == splitkey)
                        left = (Node *)galc->absolute(parent->entries[parent->slot[loc - 1]].val);
                    else
                        left = (Node *)galc->absolute(parent->entries[parent->slot[loc - 2]].val);
                }

                if(loc > parent->slot[0])  {
                    right = NULL; 
                } else if(loc == parent->slot[0]) {
                    if(parent->entries[parent->slot[loc]].key == splitkey)
                        right = NULL;
                    else
                        right = (Node *)galc->absolute(parent->entries[parent->slot[loc]].val);
                } else {
                    if(parent->entries[parent->slot[loc]].key == splitkey)
                        right = (Node *)galc->absolute(parent->entries[parent->slot[loc + 1]].val);
                    else
                        right = (Node *)galc->absolute(parent->entries[parent->slot[loc]].val);
                }
            }
        }
};

class wbtree {
    private:
        wbtree_entrance_t *entrance;
        Node * root;
        log_area *start_log;
        string tree_id;

	public:
		uint64_t GetHitCnt() {
			return	0; 
		}
		
		double GetVerTime() {
			return 0;
		}
		size_t GetCacheUse() {
			return 0;
		}
    public:
        wbtree(string path, bool recover, string id = "wbtree") {
            if(recover == false) {
                galc = new PMAllocator(path.c_str(), false, id.c_str());

                entrance = (wbtree_entrance_t *) galc->get_root(sizeof(wbtree_entrance_t));
                entrance->root = NULL;
                entrance->start_log = NULL;
                clflush(entrance, sizeof(wbtree_entrance_t), true);
                
                root = new Node;
                entrance->root = galc->relative(root);
                clflush(&(entrance->root), 8);

                start_log = (log_area *) galc->malloc(sizeof(log_area));
                start_log->next_offset = (log_entry *)galc->relative(start_log->log_data);
                clflush(start_log, sizeof(log_entry *), true);
                
                entrance->start_log = galc->relative(start_log);
            } else {
                galc = new PMAllocator(path.c_str(), true, id.c_str());

                entrance = (wbtree_entrance_t *) galc->get_root(sizeof(wbtree_entrance_t));

                if(entrance->root == NULL) {
                    printf("The tree is empty\n");
                    exit(-1);
                }

                root = (Node *) galc->absolute(entrance->root);
                start_log = (log_area *) galc->absolute(entrance->start_log);
            }

            tree_id = id;
        }

        ~wbtree() {
            delete galc;
            galc = nullptr;
        }

        bool find(uint64_t key, _value_t &val) {
            Node *leaf = find_leaf(key);
            
            int loc = leaf->BinSearch(key);

            if (leaf->entries[leaf->slot[loc]].key != key || loc > leaf->slot[0])
                return false;

            val = (_value_t)leaf->entries[leaf->slot[loc]].val;

            return true;
        }

        void insert(_key_t key, _value_t val) {
            // if tree level in the threshold, return false, else return the splited new root
            _key_t split_k;
            Node * split_node;
            bool splitIf = insert_recursive(root, key, val, split_k, split_node);

            if(splitIf) { // splitting cascades to the root node
                Node *new_root = new Node;
                new_root->isleaf = 0;
                new_root->leftmostPtr = (Node *)galc->relative(root);
                new_root->bitmap = 1 + setBit(0, 0);
                new_root->entries[0].val = (char *)galc->relative(split_node);
                new_root->entries[0].key = split_k;
                new_root->slot[1] = 0;
                new_root->slot[0] = 1;
                clflush(new_root, sizeof(Node), false);

                entrance->root = galc->relative(new_root);
                clflush(&(entrance->root), 8, false);

                root = new_root;
            }
        }

        bool update(_key_t key, _value_t val) {
            Node *leaf = find_leaf(key);

            int loc = leaf->BinSearch(key);
            if (leaf->entries[leaf->slot[loc]].key != key || loc > leaf->slot[0])
                return false;

            leaf->entries[leaf->slot[loc]].val = (char *)val;
            clflush(&leaf->entries[leaf->slot[loc]].val, 8, true);

            return true;
        }
        
        bool remove(_key_t key) {
            if(root->isleaf) {
                root->removeone_leaf(key);
                return true;
            } else {
                Node * child = NULL;
                int loc = root->BinSearch(key);
                if (loc > root->slot[0]) 
                    child = (Node *)root->entries[root->slot[loc - 1]].val;
                else if (root->entries[root->slot[loc]].key <= key) 
                    child = (Node *)root->entries[root->slot[loc]].val;
                else if (loc == 1) 
                    child = (Node *)root->leftmostPtr;
                else 
                    child = (Node *)root->entries[root->slot[loc - 1]].val;
                child = (Node *)galc->absolute(child);

                bool underflowIf = remove_recursive(child, key);

                if(underflowIf) {
                    Node *leftsib = NULL, *rightsib = NULL;
                    child->get_sibling(root, leftsib, rightsib);

                    if(leftsib != NULL && (child->slot[0] + leftsib->slot[0]) <= NODE_SIZE) {
                        // merge with left node
                        add_log_entry(&(root->slot[0]), sizeof(SLOT_SIZE), LE_DATA);
                        root->removeone_inner(child->entries[child->slot[1]].key);

                        merge(leftsib, child);
                    } 
                    else if (rightsib != NULL && (child->slot[0] + rightsib->slot[0]) <= NODE_SIZE) {
                        // merge with right node
                        add_log_entry(&(root->slot[0]), sizeof(SLOT_SIZE), LE_DATA);
                        root->removeone_inner(rightsib->entries[rightsib->slot[1]].key);

                        merge(child, rightsib);
                    }

                    if(root->slot[0] == 0) { // the root is empty
                        Node * old_root = root;

                        entrance->root = root->leftmostPtr;
                        clflush(&(entrance->root), 8);

                        root = (Node *)galc->absolute(root->leftmostPtr);

                        galc->free(old_root);
                    }
                }

                return true;
            }
        }

    void printAll() {
        root->print(true);
    }

    void print_height() {
        int height = 1;
        Node * curr = root;
        _key_t key = 0;
        while(!curr->isleaf) {
            int loc = curr->BinSearch(key);

            if (loc > curr->slot[0])
                curr = (Node *)curr->entries[curr->slot[loc - 1]].val;
            else if (curr->entries[curr->slot[loc]].key <= key)
                curr = (Node *)curr->entries[curr->slot[loc]].val;
            else if (loc == 1)
                curr = (Node *)curr->leftmostPtr;
            else
                curr = (Node *)curr->entries[curr->slot[loc - 1]].val;
            curr = (Node *)galc->absolute(curr);
            height++;
        }
        printf("height: %d\n", height);
    }

    private:
        bool insert_recursive(Node * n, _key_t k, _value_t v, _key_t &split_k, Node * &split_node) {
            if(n->isleaf) {
                return store(n, k, v, split_k, split_node);
            } else {
                Node * child = NULL;
                int loc = n->BinSearch(k);
                if (loc > n->slot[0]) 
                    child = (Node *)n->entries[n->slot[loc - 1]].val;
                else if (n->entries[n->slot[loc]].key <= k) 
                    child = (Node *)n->entries[n->slot[loc]].val;
                else if (loc == 1) 
                    child = (Node *)n->leftmostPtr;
                else 
                    child = (Node *)n->entries[n->slot[loc - 1]].val;
                child = (Node *)galc->absolute(child);
                
                _key_t split_k_child;
                Node * split_node_child;
                bool splitIf = insert_recursive(child, k, v, split_k_child, split_node_child);

                if(splitIf) { 
                    return store(n, split_k_child, (_value_t)galc->relative(split_node_child), split_k, split_node);
                } 
                return false;
            }
        }

        bool remove_recursive(Node * n, _key_t key) {
            if(n->isleaf) {
                n->removeone_leaf(key);
                return n->slot[0] < MIN_LIVE_ENTRIES;
            } else {
                Node * child = NULL;
                int loc = n->BinSearch(key);
                if (loc > n->slot[0]) 
                    child = (Node *)n->entries[n->slot[loc - 1]].val;
                else if (n->entries[n->slot[loc]].key <= key) 
                    child = (Node *)n->entries[n->slot[loc]].val;
                else if (loc == 1) 
                    child = (Node *)n->leftmostPtr;
                else 
                    child = (Node *)n->entries[n->slot[loc - 1]].val;
                child = (Node *)galc->absolute(child);

                bool underflowIf = remove_recursive(child, key);

                if(underflowIf) {
                    Node *leftsib = NULL, *rightsib = NULL;
                    child->get_sibling(n, leftsib, rightsib);

                    if(leftsib != NULL && (child->slot[0] + leftsib->slot[0]) <= NODE_SIZE) {
                        // merge with left node
                        add_log_entry(&(n->slot[0]), sizeof(SLOT_SIZE), LE_DATA);
                        n->removeone_inner(child->entries[child->slot[1]].key);

                        merge(leftsib, child);

                        return n->slot[0] < MIN_LIVE_ENTRIES;
                    } 
                    else if (rightsib != NULL && (child->slot[0] + rightsib->slot[0]) <= NODE_SIZE) {
                        // merge with right node
                        add_log_entry(&(n->slot[0]), sizeof(SLOT_SIZE), LE_DATA);
                        n->removeone_inner(rightsib->entries[rightsib->slot[1]].key);

                        merge(child, rightsib);

                        return n->slot[0] < MIN_LIVE_ENTRIES;
                    }
                    else { // can not merge
                        return false;
                    }
                }

                return false;
            }
        }

        Node * find_leaf(_key_t key) {
            Node *curr = root;
            while(!curr->isleaf) {
                int loc = curr->BinSearch(key);

                if (loc > curr->slot[0]) 
                    curr = (Node *)curr->entries[curr->slot[loc - 1]].val;
                else if (curr->entries[curr->slot[loc]].key <= key) 
                    curr = (Node *)curr->entries[curr->slot[loc]].val;
                else if (loc == 1) 
                    curr = (Node *)curr->leftmostPtr;
                else 
                    curr = (Node *)curr->entries[curr->slot[loc - 1]].val;

                curr = (Node *)galc->absolute(curr);
            }

            return curr;
        }
        
        void add_log_entry(void *addr, unsigned int size, unsigned char type)
        {
            log_entry *log;
            int i, remain_size;

            remain_size = size - ((size / LOG_DATA_SIZE) * LOG_DATA_SIZE);

            if ((char *) galc->absolute(start_log->next_offset) == (start_log->log_data + LOG_AREA_SIZE))
                start_log->next_offset = (log_entry *)galc->relative(start_log->log_data);

            if (size <= LOG_DATA_SIZE) {
                log = (log_entry *)galc->absolute(start_log->next_offset);
                log->size = size;
                log->type = type;
                log->addr = galc->relative(addr);
                memcpy(log->data, addr, size);

                if (type == LE_DATA)
                    clflush(log, sizeof(log_entry), false);
                else
                    clflush(log, sizeof(log_entry), true);

                start_log->next_offset = start_log->next_offset + 1;
            } else {
                void *next_addr = addr;

                for (i = 0; i < size / LOG_DATA_SIZE; i++) {
                    log = (log_entry *)galc->absolute(start_log->next_offset);
                    log->size = LOG_DATA_SIZE;
                    log->type = type;
                    log->addr = galc->relative(next_addr);
                    memcpy(log->data, next_addr, LOG_DATA_SIZE);

                    clflush(log, sizeof(log_entry), false);

                    start_log->next_offset = start_log->next_offset + 1;
                    if ((char *)galc->absolute(start_log->next_offset) == (start_log->log_data + LOG_AREA_SIZE))
                        start_log->next_offset = (log_entry *)galc->relative(start_log->log_data);

                    next_addr = (char *)next_addr + LOG_DATA_SIZE;
                }

                if (remain_size > 0) {
                    log = (log_entry *)galc->absolute(start_log->next_offset);
                    log->size = LOG_DATA_SIZE;
                    log->type = type;
                    log->addr = next_addr;
                    memcpy(log->data, next_addr, remain_size);

                    clflush(log, sizeof(log_entry), false);
                    
                    start_log->next_offset = start_log->next_offset + 1;
                }
            }
        }
        
        bool store(Node * n, _key_t key, _value_t val, _key_t & split_k, Node * & splitNode) {
            if(n->slot[0] == NODE_SIZE) {
                add_log_entry(n, sizeof(Node), LE_DATA);
                splitNode = new Node;
                int j, loc, cp = n->slot[0];
                
                splitNode->isleaf = n->isleaf;
                if(n->isleaf)
                    splitNode->leftmostPtr = n->leftmostPtr;

                //overflown Node
                for (j = MIN_LIVE_ENTRIES; j > 0; j--) {
                    loc = splitNode->Append(n->entries[n->slot[cp]].key, n->entries[n->slot[cp]].val);
                    splitNode->slot[j] = loc;
                    splitNode->slot[0]++;

                    splitNode->bitmap = setBit(splitNode->bitmap, loc);

                    n->bitmap = clearBit(n->bitmap, n->slot[cp]);
                    cp--;
                }

                n->slot[0] -= MIN_LIVE_ENTRIES;

                if (splitNode->entries[splitNode->slot[1]].key > key) {
                    loc = n->insertone_noflush(key, (void *)val);
                    clflush(&(n->entries[loc]), sizeof(Record), false);
                } else {
                    splitNode->insertone_noflush(key, (void *)val);
                }

                split_k = splitNode->entries[splitNode->slot[1]].key;

                if(n->isleaf) {
                    n->leftmostPtr = galc->relative(splitNode);
                    clflush(&(n->leftmostPtr), 8);
                }
                
                clflush(n->slot, (char *)n->entries - (char *)n->slot, false);

                add_log_entry(NULL, 0, LE_COMMIT);

                return true;
            }
            else{
                n->insertone(key, (void *)val);

                return false;
            }
        }
        
        void merge(Node *lnode, Node * rnode) {
            add_log_entry(lnode, sizeof(Node), LE_DATA); // keep the node in the log first

            for(int i = 1; i <= rnode->slot[0]; i++) {
                int slotid = rnode->slot[i];
                lnode->insertone_noflush(rnode->entries[slotid].key, rnode->entries[slotid].val);
            }

            if(lnode->isleaf) lnode->leftmostPtr = rnode->leftmostPtr;

            clflush(lnode, sizeof(Node), true);

            add_log_entry(NULL, 0, LE_COMMIT);

            galc->free(rnode);
        }
}; // class wbtree

}; // namespace wbtree

#endif // __WBTREE_BITMAP__