import argparse
import hashlib
import sys

import pytest

#parser = argparse.ArgumentParser()
#parser.add_argument("--e", required=True)
#parser.add_argument("--r", required=True)
#rgs = parser.parse_args()
#result = args.r
#exp_result = args.e


def get_files():
    with open("expected_simple_test.tsv", "r") as expected_file:
        exp_r = expected_file.read()
    with open("simple_test.tsv", "r") as res_file:
        r = res_file.read()
    return exp_r, r


def test_answer():
    assert result == exp_result
    # assert exp_hashed == result_hashed


result, exp_result = get_files()
exp_hashed = hashlib.sha256(exp_result.encode('utf-8')).hexdigest()
result_hashed = hashlib.sha256(result.encode('utf-8')).hexdigest()
print(f"expected: {hashlib.sha256(exp_result.encode('utf-8')).hexdigest()}")
print(f"result: {hashlib.sha256(result.encode('utf-8')).hexdigest()}")
mark = test_answer()
# if mark is not True:
#     sys.exit()
print("Tests passed!")

# if __name__ == "__main__":
#     # print(f"expected: {args.e}")
#     # print(f"result: {args.r}")
#

# This is a test push
