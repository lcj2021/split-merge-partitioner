#ifndef HEP_MIN_HEAP_HPP
#define HEP_MIN_HEAP_HPP

#include "common.hpp"

template<typename ValueType, typename KeyType, typename IdxType = vid_t>
class HepMinHeap {
private:
    IdxType n;
    std::vector<std::pair<ValueType, KeyType>> heap;
    std::vector<IdxType> key2idx;

public:
    HepMinHeap() : n(0), heap(), key2idx() { }

    std::vector<std::pair<ValueType, KeyType>> getHeap(){
        return heap;
    }

    IdxType getSize(){
        return n;
    }

    IdxType shift_up(IdxType cur) {
        if (cur == 0) return 0;
        IdxType p = (cur-1) / 2;

        if (heap[cur].first < heap[p].first) {
            std::swap(heap[cur], heap[p]);
            std::swap(key2idx[heap[cur].second], key2idx[heap[p].second]);
            return shift_up(p);
        }
        return cur;
    }

    void shift_down(IdxType cur) {
        IdxType l = cur*2 + 1;
        IdxType r = cur*2 + 2;

        if (l >= n)
            return;

        IdxType m = cur; // m is the minimum child
        if (heap[l].first < heap[cur].first) {
            m = l;
        }

        if (r < n && heap[r].first < heap[m].first){
            m = r;
        }

        if (m != cur) {
            std::swap(heap[cur], heap[m]);
            std::swap(key2idx[heap[cur].second], key2idx[heap[m].second]);
            shift_down(m);
        }
    }

    void insert(ValueType value, KeyType key) {
        heap[n] = std::make_pair(value, key);
        key2idx[key] = n++;
        IdxType cur = shift_up(n-1);
        shift_down(cur);

    }

    bool contains(KeyType key) {
        return key2idx[key] < n && heap[key2idx[key]].second == key;
    }

    void decrease_key(KeyType key, ValueType d = 1, ValueType old_value = 0) {

        if (!contains(key)){
            insert(old_value, key);
        }

        IdxType cur = key2idx[key];


        CHECK(cur < n && heap[cur].second == key) << "key not found " << key << std::endl;

        CHECK_GE(heap[cur].first, d) << "value cannot be negative at key " << key << std::endl;
        heap[cur].first -= d;
        shift_up(cur);
    }

    bool remove(KeyType key) {
        IdxType cur = key2idx[key];
        if (cur >= n || heap[cur].second != key)
            return false;

        n--;
        if (n > 0) {
            heap[cur] = heap[n];
            key2idx[heap[cur].second] = cur;
            cur = shift_up(cur);
            shift_down(cur);
        }
        return true;
    }

    bool get_min(ValueType& value, KeyType& key) {
        if (n > 0) {
            value = heap[0].first;
            key = heap[0].second;
            return true;
        } else
            return false;
    }



    void reserve(IdxType nelements) {
        n = 0;
        heap.resize(nelements);
        key2idx.resize(nelements);
    }

    void clear() {
        n = 0;
    }

};

#endif