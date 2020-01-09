import os
import pytest
from simutator import utils

modules_dir = os.path.dirname(os.path.abspath(utils.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "utils")


def test_syscall_get_stdout():
    """test syscall with no error"""
    got = utils.syscall("echo testing 123")
    assert got.stdout == "testing 123\n"
    assert got.stderr == ""


def test_syscall_get_stderr():
    """test syscall with no error"""
    got = utils.syscall("echo testing 123 >&2")
    assert got.stdout == ""
    assert got.stderr == "testing 123\n"


def test_syscall_with_error():
    """test syscall when there is an error"""
    with pytest.raises(RuntimeError):
        utils.syscall("notacommandunlessyoumadeitone")
