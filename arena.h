#ifndef arena_h
#define arena_h


struct Arena {
    char *data;
    int capacity;
    int taken;
    bool locked;

    Arena(int _capacity);
    ~Arena();

    void clear();
    void unlock();
    char *alloc_bytes(int bytes);
    char *lock_bytes(int bytes);

    template<typename T>
    T *alloc(int count=1) {
        return (T *)alloc_bytes(sizeof(T) * count);
    }

    template<typename T>
    T *lock(int count=1) {
        return (T *)lock_bytes(sizeof(T) * count);
    }

    template<typename T>
    T *rest() {
        return (T *)(data + taken);
    }
};

#endif
