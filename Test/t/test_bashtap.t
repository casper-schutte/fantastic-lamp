#!/usr/bin/env bash

#. $(dirname $0)/bash-tap-bootstrap


BASH_TAP_ROOT=../bash-tap/
. ../bash-tap/bash-tap-bootstrap

#PATH=../:$PATH

exp_result=$(head ../expected_simple_test.tsv)
result=$(head ../simple_test.tsv)

plan tests 1

is "$exp_result" "$result" "Test passed (bash-tap)"
