# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


def test_dummy():
    from dummy import dummy
    assert dummy() == 1
