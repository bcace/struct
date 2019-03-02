#include "node.h"
#include "beam.h"
#include "arena.h"
#include "sparse.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>


double _error(double a, double b) {
    return fabs(a - b) / fabs(a);
}

void _assert_double_eq(double a, double b, double max_error=0.00001) {
    double error = _error(a, b);
    if (error > max_error) {
        fprintf(stderr, "error: %g\n", error);
        assert(false);
    }
}

void beam_test(Arena &arena) {
    arena.clear();

    Node *nodes = arena.alloc<Node>(3);
    int nodes_count = 0;
    node_add(nodes, nodes_count, dvec3(0, 0, 0), 0);
    node_add(nodes, nodes_count, dvec3(1, 0, 0), DOF_X | DOF_Y);
    node_add(nodes, nodes_count, dvec3(0, 1, 0), DOF_Y);

    int eqs_count = node_index_node_eqs(nodes, nodes_count);

    double *interface = arena.alloc<double>(eqs_count);
    node_add_load(nodes, 1, interface, dvec3(0.0, -1000.0, 0.0));

    BeamProperties beam_props;
    beam_props.set_rectangular(0.02, 0.1, 0.1, 0.1);

    SparseMatrix K(eqs_count, arena.alloc<double>(eqs_count * eqs_count));

    beam_add_to_global(K, nodes[0], nodes[1], dvec3(0, 0, 1), beam_props);
    beam_add_to_global(K, nodes[0], nodes[2], dvec3(0, 0, 1), beam_props);
    beam_add_to_global(K, nodes[1], nodes[2], dvec3(0, 0, 1), beam_props);

    K.solve(interface);
    // K.print();

    _assert_double_eq(interface[0], -1.42857e-07);

    printf("beam_test done.\n");
}


int main() {
    Arena arena(500000);

    beam_test(arena);

    return 0;
}
