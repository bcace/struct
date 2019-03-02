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

void _print_interface(double *interface, int count) {
    printf("interface: ");
    for (int i = 0; i < count; ++i)
        printf("%g ", interface[i]);
    printf("\n");
}

void beam_test(Arena &arena) {
    arena.clear();

    Node *nodes = arena.alloc<Node>(50);
    int nodes_count = 0;
    node_add(nodes, nodes_count, dvec3(0, 0, 0), DOF_X | DOF_Y | DOF_Z);
    node_add(nodes, nodes_count, dvec3(5, 0, 0), DOF_RX | DOF_RY | DOF_RZ);

    int eqs_count = node_index_node_eqs(nodes, nodes_count);

    double *interface = arena.alloc<double>(eqs_count);
    node_add_force(nodes, 0, interface, dvec3(25, 0, 25));
    node_add_moment(nodes, 1, interface, dvec3(0, 100, 0));

    BeamProperties beam_props;
    beam_props.set_rectangular(30e+6, 0.3, 0.5, 0.5);

    SparseMatrix K(eqs_count, arena.alloc<double>(eqs_count * eqs_count));

    beam_add_to_global(K, nodes[0], nodes[1], dvec3(0, 0, -1), beam_props);

    K.solve(interface);
    // K.print();

    _print_interface(interface, eqs_count);

    /* ux= 1.667 × 10−5, uz= 0.01467 at node A, ϕy= 5.2 × 10−3 at node B. */

    printf("beam_test done.\n");
}


int main() {
    Arena arena(500000);

    beam_test(arena);

    return 0;
}
