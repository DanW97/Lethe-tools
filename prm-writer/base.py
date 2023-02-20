#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File     :   base.py
# Time     :   14/02/2023
# Author   :   Daniel Weston
# Version  :   0.1.0
# Contact  :   dtw545@bham.ac.uk

from pathlib import Path

# from itertools import product


# control prm details for all the core stuff that is used everywhere in Lethe
class PRMBase:
    def __init__(self, prm_file: str | Path = None, **kwargs):
        pass

    # write each subsection to file
    def write_prm(self, destination: str | Path):
        pass

    # write multiple prm files based on varied parameters
    # this one writes stuff in lockstep, terminating at the shortest iterable provided
    def write_prm_lockstep(
        self, pattern: str | Path, suffix: list[str] | None, *args: dict
    ):
        pass

    # write multiple prm files based on all combinations of varied parameters
    def write_prm_cartesian_product(
        self, pattern: str | Path, suffix: list[str] | None, *args: dict
    ):
        pass

    # TODO have a better name
    # convert prm file into dictionary
    def _serialise(self, prm_file: str | Path):
        pass


# general stuff
# box refinement
# dimensionality
# linear solver control
# mesh
# non-linear solver control
# restart
# simulation control
# timer


# Each class has an init method and a validate method, as well as getters and setters
class BoxRefinement:
    def __init__(self, parameters: dict = None):
        pass

    # only called if parameters is not none
    def _validate(self, parameters: dict):
        pass
