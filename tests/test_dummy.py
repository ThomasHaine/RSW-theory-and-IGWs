import pytest
from rsw_theory_and_igws.dummy_module import dummy_foo


def test_dummy():
    assert dummy_foo(4) == (4 + 4)
