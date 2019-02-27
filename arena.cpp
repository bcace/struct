#include "arena.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>


Arena::Arena(int _capacity) {
    capacity = _capacity;
    data = (char *)malloc(capacity);
    taken = 0;
    locked = false;
}

Arena::~Arena() {
    free(data);
}

void Arena::clear() {
    taken = 0;
    locked = false;
}

void Arena::unlock() {
    locked = false;
}

char *Arena::alloc_bytes(int bytes) {
    assert(taken + bytes < capacity);
    assert(!locked);
    char *mem = data + taken;
    taken += bytes;
    return mem;
}

char *Arena::lock_bytes(int bytes) {
    assert(taken + bytes < capacity);
    assert(!locked);
    locked = true;
    return data + taken;
}
