import itertools
import json
from typing import Any, cast

import galois
import numpy as np
from numpy.typing import NDArray


GF2: type[galois.FieldArray] = galois.GF(2)


class Solution:
    A: NDArray[np.uint8]
    B: NDArray[np.uint8]
    x_stab: NDArray[np.uint8]
    target: galois.FieldArray

    def __init__(self, sol: dict[str, Any]) -> None:
        self.A, self.B = np.asarray(sol["sol"]["A"], dtype=np.uint8), np.asarray(
            sol["sol"]["B"], dtype=np.uint8
        )
        self.x_stab = np.asarray(sol["x_stabilizers"], dtype=np.uint8)
        self.target = GF2(self.x_stab).null_space().row_reduce()

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


if __name__ == "__main__":
    with open("./data/zoo_codes.json", "r") as f:
        for code in json.load(f):
            sol = Solution(code)
            sol.verify()
