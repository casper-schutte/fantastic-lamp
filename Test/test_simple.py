import argparse
import hashlib

parser = argparse.ArgumentParser()
parser.add_argument("--e", required=True)
parser.add_argument("--r", required=True)
args = parser.parse_args()
result = args.r
exp_result = args.e

if __name__ == "__main__":
    # print(f"expected: {args.e}")
    # print(f"result: {args.r}")
    exp_hashed = hashlib.sha256(exp_result.encode('utf-8')).hexdigest()
    result_hashed = hashlib.sha256(result.encode('utf-8')).hexdigest()
    print(f"expected: {hashlib.sha256(exp_result.encode('utf-8')).hexdigest()}")
    print(f"result: {hashlib.sha256(result.encode('utf-8')).hexdigest()}")
    assert result == exp_result
    assert exp_hashed == result_hashed
    print("Tests passed!")

# This is a test push
