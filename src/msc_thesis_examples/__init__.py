import itertools
import json
from typing import Any, cast

import galois
import numpy as np
from numpy.typing import NDArray

from ortools.sat.python import cp_model


GF2: type[galois.FieldArray] = galois.GF(2)


class Solution:
    A: NDArray[np.uint8]
    B: NDArray[np.uint8]
    x_stab: NDArray[np.uint8]
    target: galois.FieldArray
    name: str

    def __init__(self, sol: dict[str, Any]) -> None:
        self.A, self.B = np.asarray(sol["sol"]["A"], dtype=np.uint8), np.asarray(
            sol["sol"]["B"], dtype=np.uint8
        )
        self.x_stab = np.asarray(sol["x_stabilizers"], dtype=np.uint8)
        self.target = GF2(self.x_stab).null_space().row_reduce()
        self.name = sol["name"]

        # dynamically assert the dimensions of A and B.
        assert self.A.shape[0] == self.B.shape[0]

        # assert the weight of each row which is a very basic requirement
        assert np.all(np.sum(self.A, axis=1) == 2)
        assert np.all(np.sum(self.B, axis=1) == 1)

    def _verify_non_overlapping_detecting_regions(self) -> None:
        # we want to make sure there is a non-overlapping basis of the kernel of A
        # each vector in the kernel represents a (red) detecting region.
        # This can be done by first find the nullspace and do the row reduction on that as a matrix
        kernel: galois.FieldArray = GF2(self.A).null_space().row_reduce()
        assert np.sum(np.asarray(kernel, dtype=np.uint8)) == self.A.shape[1]

    def _verify_code_consistency(self) -> None:
        # we need to ensure that the zero-rows of A after reduction give the same row space as the targeted ones.
        kernel: galois.FieldArray = GF2(self.A.transpose()).null_space()
        # RREF is a good canonical representation of the basis.
        my_basis: galois.FieldArray = (kernel @ GF2(self.B)).row_reduce()
        # now check if row space are the same
        assert np.array_equal(my_basis, self.target)

    def _verify_fault_tolerance(self) -> None:
        # within each non-overlapping detecting region, check arbitrary weight-2 faults would be pushed out without increasing weight.
        # first find the parity check matrix
        M: NDArray[np.uint8] = np.concat((self.A, self.B), axis=1)
        H: galois.FieldArray = GF2(M).null_space()
        # now x in rowspan(M) <=> H x = 0

        # build a table of the "syndrome" for the B part when x pushes out to weight 0, 1, 2
        detecting_regions: galois.FieldArray = GF2(self.A).null_space().row_reduce()
        low_weight_B_syndrome: set[tuple[np.uint8, ...]] = {
            tuple(np.zeros(H.shape[1], dtype=np.uint8))
        }
        for w in (1, 2):
            for loc in itertools.combinations(
                range(self.A.shape[1], self.A.shape[1] + self.B.shape[1]), w
            ):
                low_weight_B_syndrome.add(
                    tuple(np.asarray(H[:, loc].sum(axis=1), dtype=np.uint8))
                )

        # now for every weight 2 fault in each detecting region, check its syndrome is gonna match any existing low weight B-syndrome.
        for region in detecting_regions:
            region = cast(
                list[int], [i for i in range(0, len(region)) if region[i] == 1]
            )
            for loc in itertools.combinations(region, 2):
                assert (
                    tuple(np.asarray(H[:, loc].sum(axis=1), dtype=np.uint8))
                    in low_weight_B_syndrome
                )

    def verify(self) -> None:
        self._verify_non_overlapping_detecting_regions()
        self._verify_code_consistency()
        self._verify_fault_tolerance()

    """
    Circuit extraction using SAT solver
    """

    def extract(self) -> cp_model.CpSolver:
        m, k, n = self.A.shape[1], self.A.shape[0], self.B.shape[1]
        # number of vertices + number of edges
        max_time = m + k + n + 3 * k
        mdl = cp_model.CpModel()

        # First build the edge and node variables from the incidence matrices described by self.A and self.B

        # for each edge, it has the following variables
        # - exactly one of "red -> green", "green -> red", "vertical"
        # - an integer "time" if it's vertical

        # for each spider,
        # - start: if it's the start of a qubit line
        # - end: if it's the end of a qubit line (the boundary spiders are automatically end of a qubit line)
        # - t_lb, t_ub which specifies the time of the CNOT that the spider is involved in.

        # all the other variables can be deduced from the above.
        # we use the following convention for vertices labelling:
        # - 0, ..., m -1: the deepest green
        # - m, ..., m + k - 1: the intermediate red
        # - m + k, ..., m + k + n - 1: the outermost green
        # and each edge is labeled from green to red

        # node variables
        start = [mdl.new_bool_var(f"start_{v}") for v in range(m + k + n)]
        end = [mdl.new_bool_var(f"end_{v}") for v in range(m + k + n)]
        t_lb = [mdl.new_int_var(0, max_time, f"t_lb_{v}") for v in range(m + k + n)]
        t_ub = [mdl.new_int_var(0, max_time, f"t_ub_{v}") for v in range(m + k + n)]

        ## auxiliary variables
        incoming = [[] for v in range(m + k + n)]
        outgoing = [[] for v in range(m + k + n)]

        # edge variables
        cnot = {}
        forward = {}
        backward = {}
        time = {}

        def append_edge(g, r):
            cnot[g, r] = mdl.new_bool_var(f"cnot_{g}_{r}")
            forward[g, r] = mdl.new_bool_var(f"forward_{g}_{r}")
            backward[g, r] = mdl.new_bool_var(f"backward_{g}_{r}")
            time[g, r] = mdl.new_int_var(0, max_time, f"time_{g}_{r}")

            # enforce the exactly one condition
            mdl.add(cnot[g, r] + forward[g, r] + backward[g, r] == 1)

            # record the auxiliary
            incoming[r].append(forward[g, r])
            outgoing[g].append(forward[g, r])

            incoming[g].append(backward[g, r])
            outgoing[r].append(backward[g, r])

            # now based on the type of the edge, enforce the condition on time
            ## if it's a cnot
            mdl.add(t_lb[g] <= time[g, r]).only_enforce_if(cnot[g, r])
            mdl.add(time[g, r] <= t_ub[g]).only_enforce_if(cnot[g, r])

            mdl.add(t_lb[r] <= time[g, r]).only_enforce_if(cnot[g, r])
            mdl.add(time[g, r] <= t_ub[r]).only_enforce_if(cnot[g, r])

            ## if it's a horizontal edge then enforce the ordering of timing interval
            mdl.add(t_lb[r] >= t_ub[g] + 1).only_enforce_if(forward[g, r])
            mdl.add(t_lb[g] >= t_ub[r] + 1).only_enforce_if(backward[g, r])

        # start with matrix A
        for i in range(k):
            for j in range(m):
                if self.A[i, j] == 1:
                    append_edge(j, m + i)

        # then matrix B
        for i in range(k):
            for j in range(n):
                if self.B[i, j] == 1:
                    append_edge(m + k + j, m + i)

        # boundary spiders should have end = 1
        for v in range(m + k, m + k + n):
            mdl.add(end[v] == 1)

        # we further require that the end spiders should all be green
        for v in range(m, m + k):
            mdl.add(end[v] == 0)

        # enforce the number of incoming/outgoing condition
        for v in range(m + k + n):
            mdl.add(sum(incoming[v]) + start[v] == 1)
            mdl.add(sum(outgoing[v]) + end[v] == 1)
            mdl.add(t_lb[v] <= t_ub[v])

        # minimize the number of CNOTs
        mdl.minimize(sum(cnot.values()))

        solver = cp_model.CpSolver()
        status = solver.solve(mdl)

        print(self.name, solver.status_name(status), solver.objective_value)

        return solver


if __name__ == "__main__":
    with open("./data/zoo_codes.json", "r") as f:
        for code in json.load(f):
            sol = Solution(code)
            sol.verify()
            sol.extract()
