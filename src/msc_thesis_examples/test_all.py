from . import Solution
import json


def test_all() -> None:
    with open("./data/zoo_codes.json", "r") as f:
        for code in json.load(f):
            sol = Solution(code)
            sol.verify()
