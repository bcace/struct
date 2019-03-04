#include "test.h"
#include "node.h"
#include "mitc4.h"
#include "sparse.h"
#include "arena.h"
#include <math.h>
#include <stdio.h>


void __print_interface(double *interface, int count) {
    printf("interface: ");
    for (int i = 0; i < count; ++i)
        printf("%g ", interface[i]);
    printf("\n");
}

void test_pinched_cylinder(Arena &arena) {
    arena.clear();

    double R = 0.3;
    double L = 0.3; /* half of cylinder's length */
    double h = 0.003; /* thickness */

    double E = 3.0e+6;
    double v = 0.3;
    double F = 0.25;

    int side_elements_count = 10;
    int side_nodes_count = side_elements_count + 1;
    double da = M_PI_2 / side_elements_count;
    double dy = L / side_elements_count;

    Node *nodes = arena.alloc<Node>(side_nodes_count * side_nodes_count);
    int nodes_count = 0;

    for (int i = 0; i < side_nodes_count; ++i) {
        double y = dy * i;

        unsigned i_bcs = 0; /* constrained dofs */
        if (i == 0)
            i_bcs = DOF_X | DOF_Y | DOF_RY; // Tx = Ty = Ry = 0
        else if (i == side_elements_count)
            i_bcs = DOF_Y | DOF_RX | DOF_RZ; // Ty = Rx = Rz = 0

        for (int j = 0; j < side_nodes_count; ++j) {
            double a = da * j;
            double x = cos(a) * R;
            double z = sin(a) * R;

            unsigned bcs = i_bcs;
            if (j == 0)
                bcs |= DOF_Z | DOF_RX | DOF_RY; // Tz = Rx = Ry = 0
            else if (j == side_elements_count)
                bcs |= DOF_X | DOF_RY | DOF_RZ; // Tx = Ry = Rz = 0

            node_add(nodes, nodes_count, dvec3(x, y, z), (~bcs) & DOF_ALL);
        }
    }

    int eqs_count = node_index_node_eqs(nodes, nodes_count);

    printf("%d\n", eqs_count);

    double *interface = arena.alloc<double>(eqs_count);
    node_clear_loads(interface, eqs_count);
    node_add_force(nodes, nodes_count - 1, interface, dvec3(0, 0, -F));

    SparseMatrix K(eqs_count, arena.alloc<double>(eqs_count * eqs_count));

    Mitc4Properties props(E, v, h);

    for (int i = 0; i < side_elements_count; ++i) {
        int base = i * side_nodes_count;
        int next_base = base + side_nodes_count;
        for (int j = 0; j < side_elements_count; ++j) {
            mitc4_add_to_global(K,
                                nodes[base + j],
                                nodes[next_base + j],
                                nodes[next_base + j + 1],
                                nodes[base + j + 1],
                                props);
        }
    }

    K.solve(interface);

    __print_interface(interface, eqs_count);

    printf("test_pinched_cylinder finished.\n");
}
