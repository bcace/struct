#include "node.h"
#include <math.h>
#include <stdlib.h>


int node_add(Node *nodes, int nodes_count, dvec3 pos, unsigned dofs) {
    Node &n = nodes[nodes_count];
    n.pos = pos;
    n.eqs.x =  ((dofs & DOF_X)  == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.y =  ((dofs & DOF_Y)  == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.z =  ((dofs & DOF_Z)  == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.rx = ((dofs & DOF_RX) == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.ry = ((dofs & DOF_RY) == 0) ? CONSTR_DOF_EQ : 0;
    n.eqs.rz = ((dofs & DOF_RZ) == 0) ? CONSTR_DOF_EQ : 0;
    return nodes_count++;
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

void node_add_load(Node *nodes, int node_index, double *forces_vector, dvec3 force) {
    Node &n = nodes[node_index];
    if (n.eqs.x != CONSTR_DOF_EQ)
        forces_vector[n.eqs.x] = force.x;
    if (n.eqs.y != CONSTR_DOF_EQ)
        forces_vector[n.eqs.y] = force.y;
    if (n.eqs.z != CONSTR_DOF_EQ)
        forces_vector[n.eqs.z] = force.z;
}


// Nodes::Nodes(int _cap) : cap(_cap), count(0), nodes((Node *)malloc(sizeof(Node) * _cap)) {}

// Nodes::~Nodes() {
//     free(nodes);
// }

// void Nodes::clear() {
//     count = 0;
// }

// void Nodes::add_node(dvec3 pos, unsigned dofs) {
//     Node &n = nodes[count];
//     n.pos = pos;
//     n.eqs.x =  ((dofs & DOF_X)  == 0) ? CONSTR_DOF_EQ : 0;
//     n.eqs.y =  ((dofs & DOF_Y)  == 0) ? CONSTR_DOF_EQ : 0;
//     n.eqs.z =  ((dofs & DOF_Z)  == 0) ? CONSTR_DOF_EQ : 0;
//     n.eqs.rx = ((dofs & DOF_RX) == 0) ? CONSTR_DOF_EQ : 0;
//     n.eqs.ry = ((dofs & DOF_RY) == 0) ? CONSTR_DOF_EQ : 0;
//     n.eqs.rz = ((dofs & DOF_RZ) == 0) ? CONSTR_DOF_EQ : 0;
//     ++count;
// }

// void Nodes::add_load(int node_index, dvec3 force, double *F) {
//     Node &n = nodes[node_index];
//     if (n.eqs.x != CONSTR_DOF_EQ)
//         F[n.eqs.x] = force.x;
//     if (n.eqs.y != CONSTR_DOF_EQ)
//         F[n.eqs.y] = force.y;
//     if (n.eqs.z != CONSTR_DOF_EQ)
//         F[n.eqs.z] = force.z;
// }

/*
  Map node dofs to equation indices in the global matrix.
  Assumes:
  - nodes are not shuffled after this function call
  - node dofs are added into the global matrix in the same order
*/
// int Nodes::index_node_eqs() {
//     int eq_count = 0;
//     for (int i = 0; i < count; ++i) {
//         Node &n = nodes[i];
//         if (n.eqs.x != CONSTR_DOF_EQ)
//             n.eqs.x = eq_count++;
//         if (n.eqs.y != CONSTR_DOF_EQ)
//             n.eqs.y = eq_count++;
//         if (n.eqs.z != CONSTR_DOF_EQ)
//             n.eqs.z = eq_count++;
//         if (n.eqs.rx != CONSTR_DOF_EQ)
//             n.eqs.rx = eq_count++;
//         if (n.eqs.ry != CONSTR_DOF_EQ)
//             n.eqs.ry = eq_count++;
//         if (n.eqs.rz != CONSTR_DOF_EQ)
//             n.eqs.rz = eq_count++;
//     }
//     return eq_count;
// }
