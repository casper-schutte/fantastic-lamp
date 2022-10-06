import hashlib


def get_files():
    with open("expected_simple_test.tsv", "r") as expected_file:
        exp_r = expected_file.read()
    with open("simple_test.tsv", "r") as res_file:
        r = res_file.read()
    return exp_r, r


def test_answer():
    assert result == exp_result


result, exp_result = get_files()
exp_hashed = hashlib.sha256(exp_result.encode('utf-8')).hexdigest()
result_hashed = hashlib.sha256(result.encode('utf-8')).hexdigest()
print(f"expected: {hashlib.sha256(exp_result.encode('utf-8')).hexdigest()}")
print(f"result: {hashlib.sha256(result.encode('utf-8')).hexdigest()}")
test_answer()

print("Tests passed!")
