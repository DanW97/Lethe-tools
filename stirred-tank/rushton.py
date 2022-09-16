#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File    :   rushton.py
# Time    :   16/09/2022
# Author  :   Daniel Weston
# Version :   0.1.0
# Contact :   dtw545@student.bham.ac.uk

import gc
import gmsh
import sys
import numpy as np

# Global constants - these form the "dim" section of the dimtag parlance used in gmsh
POINT = 0
CURVE = 1
SURFACE = 2
VOLUME = 3


class RushtonTurbineBuilder:
    def __init__(
        self,
        shaft_height=0.12333,
        shaft_radius=0.00375,
        connector_height=0.00956,
        connector_radius=0.00625,
        hub_height=0.00058,
        hub_radius=0.025,
        blade_height=0.01333,
        blade_width=0.00089,
        blade_depth=0.01667,
        axis_alignment=1,
        offset=0.060,
        nblades=6,
        filepath="rushton",
        view = True
    ) -> None:
        self.shaft_height = shaft_height
        self.shaft_radius = shaft_radius
        self.connector_height = connector_height
        self.connector_radius = connector_radius
        self.hub_height = hub_height
        self.hub_radius = hub_radius
        self.blade_height = blade_height
        self.blade_width = blade_width
        self.blade_depth = blade_depth
        self.axis_alignment = axis_alignment
        self.offset = offset  # height of impeller from bottom of blade to tank bottom
        self.nblades = nblades
        self.filepath = filepath
        self.view = view

    def draw(self):
        gmsh.initialize(sys.argv)
        gmsh.model.add("Rushton Turbine")
        gm = gmsh.model.occ
        # this assembly consists of 3 cylinders and n boxes, all mated
        # draw cylindrical parts
        shaft = gm.addCylinder(x = 0, y = 0, z = 0, dx = 0, dy = 0, dz = self.shaft_height, r = self.shaft_radius)
        connector = gm.addCylinder(x = 0, y = 0, z = 0, dx = 0, dy = 0, dz = self.connector_height, r = self.connector_radius)
        hub = gm.addCylinder(x = 0, y = 0, z = 0, dx = 0, dy = 0, dz = self.hub_height, r = self.hub_radius)
        # draw blades, these start with a center at 0, 0, 0, the thin edge aligned to x, height to z and the less thin edge to y
        blades = [gm.addBox(x = -self.blade_depth/2, y = -self.blade_width/2, z = 0, dx = self.blade_depth, dy = self.blade_width, dz = self.blade_height) for _ in range(self.nblades)]
        angles = np.linspace(0, 2 * np.pi, self.nblades, endpoint=False)
        for (blade, angle) in zip(blades, angles):
            dimtags = [(VOLUME, blade)]
            # rotate each blade
            gm.rotate(dimtags, x = 0, y = 0, z = 0, ax = 0, ay = 0, az = 1, angle = angle)
            # translate to hub edge
            gm.translate(dimtags, dx = self.hub_radius*np.cos(angle), dy = self.hub_radius*np.sin(angle), dz = 0)
        gm.synchronize()
        # align everything to the correct axis and 
        # raise each object by their relative offsets
        hub_offset = self.offset + self.blade_height/2
        connector_offset = hub_offset + self.hub_height
        shaft_offset = connector_offset + self.connector_height
        match self.axis_alignment:
            case 0:
                for dim, tag in gm.getEntities(VOLUME):
                    gm.rotate(
                        [(dim, tag)],
                        x=0,
                        y=0,
                        z=0,
                        ax=0,
                        ay=1,
                        az=0,
                        angle=-np.pi / 2,
                    )
                gm.translate([(VOLUME, hub)], dx = hub_offset, dy = 0, dz = 0)
                gm.translate([(VOLUME, connector)], dx = connector_offset, dy = 0, dz = 0)
                gm.translate([(VOLUME, shaft)], dx = shaft_offset, dy = 0, dz = 0)
                for blade in blades:
                    gm.translate([(VOLUME, blade)], dx = self.offset, dy = 0, dz = 0)
            case 1:
                for dim, tag in gm.getEntities(VOLUME):
                    gm.rotate(
                        [(dim, tag)],
                    x=0,
                    y=0,
                    z=0,
                    ax=1,
                    ay=0,
                    az=0,
                    angle=-np.pi / 2,
                )
                gm.translate([(VOLUME, hub)], dx = 0, dy = hub_offset, dz = 0)
                gm.translate([(VOLUME, connector)], dx = 0, dy = connector_offset, dz = 0)
                gm.translate([(VOLUME, shaft)], dx = 0, dy = shaft_offset, dz = 0)
                for blade in blades:
                    gm.translate([(VOLUME, blade)], dx = 0, dy = self.offset, dz = 0)


            case 2:  # do nothing as already aligned
                gm.translate([(VOLUME, hub)], dx = 0, dy = 0, dz = hub_offset)
                gm.translate([(VOLUME, connector)], dx = 0, dy = 0, dz = connector_offset)
                gm.translate([(VOLUME, shaft)], dx = 0, dy = 0, dz = shaft_offset)
                for blade in blades:
                    gm.translate([(VOLUME, blade)], dx = 0, dy = 0, dz = self.offset)

            case _:
                raise ValueError("Axis alignment value is *only* 0, 1 or 2!")
        gm.synchronize()
        gm.fragment(gm.getEntities(SURFACE), gm.getEntities(VOLUME))
        # set physical groups
        gmsh.model.addPhysicalGroup(
            SURFACE, [tag for _, tag in gmsh.model.getEntities(SURFACE)], name="surf"
        )
        gmsh.model.addPhysicalGroup(
            VOLUME, [tag for _, tag in gmsh.model.getEntities(VOLUME)], name="vol"
        )
        gm.synchronize()

    def export(self):
        gm = gmsh.model.occ
        msh = gmsh.model.mesh
        gm.synchronize()
        gmsh.option.setNumber("Mesh.Smoothing", 100)
        gmsh.option.setNumber("Mesh.SaveAll", 0)  # only save physical groups
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        gmsh.option.setNumber("General.NumThreads", 8)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1e-3)
        gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 10)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        raw_file = f"{self.filepath}.geo_unrolled"
        new_file = f"{self.filepath}.geo"
        gmsh.write(raw_file)
        msh.generate(VOLUME)
        gmsh.write(f"{self.filepath}.msh")
        import os

        os.replace(f"{raw_file}", f"{new_file}")
        if self.view:
            self._view()       

    def _view(self):
        if "-nopopup" not in sys.argv:
            gmsh.fltk.run()

turbine = RushtonTurbineBuilder(axis_alignment=1)
turbine.draw()
turbine.export()