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
    Equations eqs; /* node dof to global matrix dof mapping, if eq set to -1 corresponding dof is constrained */
};

/*
USAGE:
1. add *all* nodes to array first using node_add(), then
2. index equations using node_index_node_eqs(), and then
3. clear loads array using node_clear_loads(), then
4. add loads using node_add_force() or node_add_moment()
because node_add_load() uses equation indices to fill
force vector correctly, and node_index_node_eqs() can
index node equations only after all nodes are created.
*/

void node_add(Node *nodes, int &nodes_count, dvec3 pos, unsigned dofs);
int node_index_node_eqs(Node *nodes, int nodes_count);
void node_clear_loads(double *loads, int eq_count);
void node_add_force(Node *nodes, int node_index, double *loads, dvec3 force);
void node_add_moment(Node *nodes, int node_index, double *loads, dvec3 moment);

#endif
