from . import Solution
import json
from ortools.sat.python import cp_model


def test_all() -> None:
    with open("./data/zoo_codes.json", "r") as f:
        for code in json.load(f):
            sol = Solution(code)
            sol.verify()
            extracted = sol.extract()
            assert int(extracted.objective_value) == code["optimal_cnot"]
