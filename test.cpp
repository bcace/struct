#include "test.h"
#include "arena.h"
#include <stdio.h>


void test_pinched_cylinder(Arena &arena) {
    arena.clear();

    // Node *nodes = arena.alloc<Node>(50);
    // int nodes_count = 0;
    // node_add(nodes, nodes_count, dvec3(0, 0, 0), DOF_X | DOF_Y | DOF_Z);
    // node_add(nodes, nodes_count, dvec3(5, 0, 0), DOF_RX | DOF_RY | DOF_RZ);
    // node_add(nodes, nodes_count, dvec3(5, 5, 0), DOF_X | DOF_Y | DOF_Z);
    // node_add(nodes, nodes_count, dvec3(0, 5, 0), DOF_RX | DOF_RY | DOF_RZ);

    // int eqs_count = node_index_node_eqs(nodes, nodes_count);

    // double *interface = arena.alloc<double>(eqs_count);
    // node_clear_loads(interface, eqs_count);

    // SparseMatrix K(eqs_count, arena.alloc<double>(eqs_count * eqs_count));

    // Mitc4Properties props(0.1, 0.1, 0.1);

    // mitc4_add_to_global(K, nodes[0], nodes[1], nodes[2], nodes[3], props);


    printf("test_pinched_cylinder finished.\n");
}
