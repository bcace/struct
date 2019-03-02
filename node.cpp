#include "node.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>


void node_add(Node *nodes, int &nodes_count, dvec3 pos, unsigned dofs) {
    Node &n = nodes[nodes_count++];
    n.pos = pos;
    n.eqs.x =  ((dofs & DOF_X)  == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.y =  ((dofs & DOF_Y)  == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.z =  ((dofs & DOF_Z)  == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.rx = ((dofs & DOF_RX) == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.ry = ((dofs & DOF_RY) == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.rz = ((dofs & DOF_RZ) == 0) ? CONSTR_DOF_EQ : 0;
}

int node_index_node_eqs(Node *nodes, int nodes_count) {
    int eq_count = 0;
    for (int i = 0; i < nodes_count; ++i) {
        Node &n = nodes[i];
        if (n.eqs.x != CONSTR_DOF_EQ)
            n.eqs.x = eq_count++;
        if (n.eqs.y != CONSTR_DOF_EQ)
            n.eqs.y = eq_count++;
        if (n.eqs.z != CONSTR_DOF_EQ)
            n.eqs.z = eq_count++;
        if (n.eqs.rx != CONSTR_DOF_EQ)
            n.eqs.rx = eq_count++;
        if (n.eqs.ry != CONSTR_DOF_EQ)
            n.eqs.ry = eq_count++;
        if (n.eqs.rz != CONSTR_DOF_EQ)
            n.eqs.rz = eq_count++;
    }
    return eq_count;
}

void node_clear_loads(double *loads, int eq_count) {
    memset(loads, 0, sizeof(double) * eq_count);
}

void node_add_force(Node *nodes, int node_index, double *loads, dvec3 force) {
    Node &n = nodes[node_index];
    if (n.eqs.x != CONSTR_DOF_EQ)
        loads[n.eqs.x] = force.x;
    if (n.eqs.y != CONSTR_DOF_EQ)
        loads[n.eqs.y] = force.y;
    if (n.eqs.z != CONSTR_DOF_EQ)
        loads[n.eqs.z] = force.z;
}

void node_add_moment(Node *nodes, int node_index, double *loads, dvec3 moment) {
    Node &n = nodes[node_index];
    if (n.eqs.rx != CONSTR_DOF_EQ)
        loads[n.eqs.rx] = moment.x;
    if (n.eqs.ry != CONSTR_DOF_EQ)
        loads[n.eqs.ry] = moment.y;
    if (n.eqs.rz != CONSTR_DOF_EQ)
        loads[n.eqs.rz] = moment.z;
}
