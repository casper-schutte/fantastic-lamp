import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--e", required=True)
parser.add_argument("--r", required=True)
args = parser.parse_args()
result = args.r
exp_result = args.e

if __name__ == "__main__":
    # print(f"expected: {args.e}")
    # print(f"result: {args.r}")
    assert result == exp_result
    print("Tests passed!")

# This is a test push
