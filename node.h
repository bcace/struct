#ifndef node_h
#define node_h

#include "vec.h"

#define CONSTR_DOF_EQ 0xffffffff


enum DofFlags {
    DOF_X   = 0x1,
    DOF_Y   = 0x2,
    DOF_Z   = 0x4,
    DOF_RX  = 0x8,
    DOF_RY  = 0x10,
    DOF_RZ  = 0x20,
    DOF_ALL = 0x3f,
};


struct Equations {
    union {
        struct {
            unsigned x, y, z, rx, ry, rz;
        };
        unsigned v[6];
    };
};


struct Node {
    dvec3 pos;
    Equations eqs;    // node dof to global matrix dof mapping, if eq set to -1 corresponding dof is constrained
};


struct Nodes {
    int cap, count;
    Node *nodes;

    Nodes(int _cap);
    ~Nodes();

    void clear();
    void add_node(dvec3 pos, unsigned dofs);
    void add_load(int node_index, dvec3 force, double *F);
    int index_node_eqs();
};

#endif
