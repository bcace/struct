#include "node.h"
#include "beam.h"
#include "mitc4.h"
#include "arena.h"
#include "sparse.h"
#include "test.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>


void _assert_double_eq(double a, double b) {
    double error = fabs(a - b) / fabs(a);
    if (error > 0.00001) {
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
    node_clear_loads(interface, eqs_count);
    node_add_force(nodes, 0, interface, dvec3(25, 0, 25));
    node_add_moment(nodes, 1, interface, dvec3(0, 100, 0));

    BeamProperties beam_props;
    beam_props.set_rectangular(30e+6, 0.3, 0.5, 0.5);

    SparseMatrix K(eqs_count, arena.alloc<double>(eqs_count * eqs_count));

    beam_add_to_global(K, nodes[0], nodes[1], dvec3(0, 0, -1), beam_props);

    K.solve(interface);

    _assert_double_eq(interface[0], 1.66667e-05);
    _assert_double_eq(interface[1], 0.0);
    _assert_double_eq(interface[2], 0.0146667);
    _assert_double_eq(interface[3], 0.0);
    _assert_double_eq(interface[4], 0.0052);
    _assert_double_eq(interface[5], 0.0);

    printf("beam_test done.\n");
}

void mitc4_test(Arena &arena) {
    arena.clear();

    Node *nodes = arena.alloc<Node>(50);
    int nodes_count = 0;
    node_add(nodes, nodes_count, dvec3(0, 0, 0), DOF_X | DOF_Y | DOF_Z);
    node_add(nodes, nodes_count, dvec3(5, 0, 0), DOF_RX | DOF_RY | DOF_RZ);
    node_add(nodes, nodes_count, dvec3(5, 5, 0), DOF_X | DOF_Y | DOF_Z);
    node_add(nodes, nodes_count, dvec3(0, 5, 0), DOF_RX | DOF_RY | DOF_RZ);

    int eqs_count = node_index_node_eqs(nodes, nodes_count);

    double *interface = arena.alloc<double>(eqs_count);
    node_clear_loads(interface, eqs_count);

    SparseMatrix K(eqs_count, arena.alloc<double>(eqs_count * eqs_count));

    Mitc4Properties props(0.1, 0.1, 0.1);

    mitc4_add_to_global(K, nodes[0], nodes[1], nodes[2], nodes[3], props);

    printf("mitc4_test done.\n");
}


int main() {
    Arena arena(10000000);

    beam_test(arena);
    mitc4_test(arena);
    test_pinched_cylinder(arena);

    return 0;
}
