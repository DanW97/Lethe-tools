#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File    :   packing_generator.py
# Date    :   11/03/2022
# Author  :   Daniel Weston
# Contact :   dtw545@student.bham.ac.uk


import gmsh
import numpy as np
import sys

# Global constants - these form the "dim" section of the dimtag parlance used in gmsh
POINT = 0
CURVE = 1
SURFACE = 2
VOLUME = 3

# using parameterised functions from: https://en.wikipedia.org/wiki/Superquadrics
# auxilliary functions
def f(w, m):
    return np.sign(np.sin(w)) * np.power(np.abs(np.sin(w)), m)


def g(w, m):
    return np.sign(np.cos(w)) * np.power(np.abs(np.cos(w)), m)


# parameterised functions that return make_points along a superquadric
# -pi/2<=v<=pi/2, -pi<=u<=pi
def xpts(u, v, scaleX, exponentX):
    return scaleX * g(v, 2 / exponentX) * g(u, 2 / exponentX)


def ypts(u, v, scaleY, exponentY):
    return scaleY * g(v, 2 / exponentY) * f(u, 2 / exponentY)


def zpts(v, scaleZ, exponentZ):
    return scaleZ * f(v, 2 / exponentZ)


class Superquadric:
    """Class handling the construction, transformation, meshing and export of a superquadric object that has points on:
    abs(x/A)**r + abs(y/B)**s + abs(z/C)**t = 1

    This class can handle values of r, s, t from 1 and above.

    Attributes
    ----------
    indices : list, optional
        Exponents used to define superquadric shape, by default [8, 8, 8]
    scale : list, optional
        Scale factors used to determine superquadric size in each axis, by default [1, 1, 1]
    filepath : str, optional
        Descriptor for saving .geo and .msh files. Extensions are not required as these are set by
        the export method, by default 'quad'
    rotation : list, optional
        Rotations about the x, y and z asix, respectively. These are in DEGREES and converted into
        radians by the rotate method, by default [0, 0, 0]
    npts : int, optional
        Number of points used to draw splines representing the superquadric, by default 20
    gr : int, optional
        Mesh size at each point, by default 1
    view : bool, optional
        Boolean to check whether the mesh is visualised, by default False
    init : bool, optional
        Ignore this, by default True
    ib : bool, optional
        Ignore this, by default True
    """

    def __init__(
        self,
        indices=[8, 8, 8],
        scale=[1, 1, 1],
        filepath="quad",
        rotation=[0, 0, 0],
        npts=20,
        gr=1,
        view=True,
        init=True,
        ib=True,
    ) -> None:
        # checking
        assert len(scale) == 3
        assert len(rotation) == 3
        assert len(indices) == 3
        if sum(np.array(rotation) != 0) > 0:
            self.rotatable = True
        else:
            self.rotatable = False
        self.scale = scale
        self.indices = indices
        self.rotation = rotation
        self.umin = np.pi / 2
        self.umax = np.pi
        self.vmin = 0
        self.vmax = np.pi / 2
        self.filepath = filepath
        self.npts = npts
        self.gr = gr
        self.view = view
        self.init = init
        self.ib = ib

    def draw(self):
        """Draw the superquadric object.

        Points are generated utilising the parameterised functions at the top of this file.

        The following conventions are used to define the orientation of the splines
        formed:
        1) N - the vertical spline parallel with x
        2) E - the vertical spline parallel with y
        3) NE - the horizontal spline parallel to the xy plane
        ... etc
        """
        if self.init:
            gmsh.initialize(sys.argv)
            gmsh.model.add("Superquadric")
        gm = gmsh.model.geo  # "with" doesn't work
        # Define coordinates needed for the 3 spline lines
        # North points
        vN = np.linspace(self.vmin, self.vmax, self.npts)
        uN = np.ones_like(vN) * self.umin
        xN = xpts(uN, vN, self.scale[0], self.indices[0])
        yN = ypts(uN, vN, self.scale[1], self.indices[1])
        zN = zpts(vN, self.scale[2], self.indices[2])
        # East points
        vE = np.linspace(self.vmin, self.vmax, self.npts)
        uE = np.ones_like(vE) * self.umax
        xE = xpts(uE, vE, self.scale[0], self.indices[0])
        yE = ypts(uE, vE, self.scale[1], self.indices[1])
        zE = zpts(vE, self.scale[2], self.indices[2])
        # South points
        vS = np.linspace(self.vmin, self.vmax, self.npts)
        uS = np.ones_like(vS) * self.umin - np.pi
        xS = xpts(uS, vS, self.scale[0], self.indices[0])
        yS = ypts(uS, vS, self.scale[1], self.indices[1])
        zS = zpts(vS, self.scale[2], self.indices[2])
        # West points
        vW = np.linspace(self.vmin, self.vmax, self.npts)
        uW = np.ones_like(vW) * self.umax - np.pi
        xW = xpts(uW, vW, self.scale[0], self.indices[0])
        yW = ypts(uW, vW, self.scale[1], self.indices[1])
        zW = zpts(vW, self.scale[2], self.indices[2])
        # North -> East points
        uNE = np.linspace(self.umin, self.umax, self.npts)
        vNE = np.zeros_like(uNE)
        xNE = xpts(uNE, vNE, self.scale[0], self.indices[0])
        yNE = ypts(uNE, vNE, self.scale[1], self.indices[1])
        zNE = zpts(vNE, self.scale[2], self.indices[2])
        # East -> South points
        uES = np.linspace(self.umin + np.pi / 2, self.umax + np.pi / 2, self.npts)
        vES = np.zeros_like(uES)
        xES = xpts(uES, vES, self.scale[0], self.indices[0])
        yES = ypts(uES, vES, self.scale[1], self.indices[1])
        zES = zpts(vES, self.scale[2], self.indices[2])
        # South -> West points
        uSW = np.linspace(self.umin + np.pi, self.umax + np.pi, self.npts)
        vSW = np.zeros_like(uSW)
        xSW = xpts(uSW, vSW, self.scale[0], self.indices[0])
        ySW = ypts(uSW, vSW, self.scale[1], self.indices[1])
        zSW = zpts(vSW, self.scale[2], self.indices[2])
        # West -> North points
        uWN = np.linspace(
            self.umin + np.pi * 3 / 2, self.umax + np.pi * 3 / 2, self.npts
        )
        vWN = np.zeros_like(uWN)
        xWN = xpts(uWN, vWN, self.scale[0], self.indices[0])
        yWN = ypts(uWN, vWN, self.scale[1], self.indices[1])
        zWN = zpts(vWN, self.scale[2], self.indices[2])
        # Draw splines
        # North spline
        Nspline_pts = []
        pts_append = Nspline_pts.append
        for i, (x, y, z) in enumerate(zip(xN, yN, zN)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        Nspline = gm.addSpline(Nspline_pts)
        # East spline
        Espline_pts = []
        pts_append = Espline_pts.append
        for i, (x, y, z) in enumerate(zip(xE, yE, zE)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        Espline_pts[-1] = Nspline_pts[-1]  # this spline ends at the start of Nspline
        Espline = gm.addSpline(Espline_pts)
        # South spline
        Sspline_pts = []
        pts_append = Sspline_pts.append
        for i, (x, y, z) in enumerate(zip(xS, yS, zS)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        Sspline_pts[-1] = Nspline_pts[-1]  # this spline ends at the start of Nspline
        Sspline = gm.addSpline(Sspline_pts)
        # West spline
        Wspline_pts = []
        pts_append = Wspline_pts.append
        for i, (x, y, z) in enumerate(zip(xW, yW, zW)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        Wspline_pts[-1] = Nspline_pts[-1]  # this spline ends at the start of Nspline
        Wspline = gm.addSpline(Wspline_pts)
        # North -> East spline
        NEspline_pts = []
        pts_append = NEspline_pts.append
        for i, (x, y, z) in enumerate(zip(xNE, yNE, zNE)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        NEspline_pts[0] = Nspline_pts[0]  # it starts at the start of the north spline
        NEspline_pts[-1] = Espline_pts[0]  # it finishes at the start of the east spline
        NEspline = gm.addSpline(NEspline_pts)
        # East -> South spline
        ESspline_pts = []
        pts_append = ESspline_pts.append
        for i, (x, y, z) in enumerate(zip(xES, yES, zES)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        ESspline_pts[0] = Espline_pts[0]  # it starts at the start of the east spline
        ESspline_pts[-1] = Sspline_pts[
            0
        ]  # it finishes at the start of the south spline
        ESspline = gm.addSpline(ESspline_pts)
        # South -> West spline
        SWspline_pts = []
        pts_append = SWspline_pts.append
        for i, (x, y, z) in enumerate(zip(xSW, ySW, zSW)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        SWspline_pts[0] = Sspline_pts[0]  # it starts at the start of the south spline
        SWspline_pts[-1] = Wspline_pts[0]  # it finishes at the start of the west spline
        SWspline = gm.addSpline(SWspline_pts)
        # West -> North spline
        WNspline_pts = []
        pts_append = WNspline_pts.append
        for i, (x, y, z) in enumerate(zip(xWN, yWN, zWN)):
            pts_append(gm.addPoint(x, y, z, self.gr))
        WNspline_pts[0] = Wspline_pts[0]  # it starts at the start of the west spline
        WNspline_pts[-1] = Nspline_pts[
            0
        ]  # it finishes at the start of the north spline
        WNspline = gm.addSpline(WNspline_pts)
        # Draw lower splines
        # North spline
        Nspline_pts_lower = []
        pts_append = Nspline_pts_lower.append
        for i, (x, y, z) in enumerate(zip(xN, yN, zN)):
            pts_append(gm.addPoint(x, y, -z, self.gr))
        Nspline_pts_lower[0] = NEspline_pts[0]
        Nspline_lower = gm.addSpline(Nspline_pts_lower)
        # East spline
        Espline_pts_lower = []
        pts_append = Espline_pts_lower.append
        for i, (x, y, z) in enumerate(zip(xE, yE, zE)):
            pts_append(gm.addPoint(x, y, -z, self.gr))
        Espline_pts_lower[0] = NEspline_pts[
            -1
        ]  # this spline starts at the end of NEspline
        Espline_pts_lower[-1] = Nspline_pts_lower[
            -1
        ]  # this spline ends at the end of lower Nspline
        Espline_lower = gm.addSpline(Espline_pts_lower)
        # South spline
        Sspline_pts_lower = []
        pts_append = Sspline_pts_lower.append
        for i, (x, y, z) in enumerate(zip(xS, yS, zS)):
            pts_append(gm.addPoint(x, y, -z, self.gr))
        Sspline_pts_lower[0] = ESspline_pts[
            -1
        ]  # this spline starts at the end of ESspline
        Sspline_pts_lower[-1] = Nspline_pts_lower[
            -1
        ]  # this spline ends at the end of lower Nspline
        Sspline_lower = gm.addSpline(Sspline_pts_lower)
        # West spline
        Wspline_pts_lower = []
        pts_append = Wspline_pts_lower.append
        for i, (x, y, z) in enumerate(zip(xW, yW, zW)):
            pts_append(gm.addPoint(x, y, -z, self.gr))
        Wspline_pts_lower[0] = SWspline_pts[
            -1
        ]  # this spline starts at the end of SWspline
        Wspline_pts_lower[-1] = Nspline_pts_lower[
            -1
        ]  # this spline ends at the end of lower Nspline
        Wspline_lower = gm.addSpline(Wspline_pts_lower)
        # Draw faces
        curve_loopNE = gm.addCurveLoop([-Nspline, NEspline, Espline])
        NEface = gm.addSurfaceFilling([curve_loopNE])
        curve_loopES = gm.addCurveLoop([-Espline, ESspline, Sspline])
        ESface = gm.addSurfaceFilling([curve_loopES])
        curve_loopSW = gm.addCurveLoop([-Sspline, SWspline, Wspline])
        SWface = gm.addSurfaceFilling([curve_loopSW])
        curve_loopWN = gm.addCurveLoop([-Wspline, WNspline, Nspline])
        WNface = gm.addSurfaceFilling([curve_loopWN])
        curve_loopNE_lower = gm.addCurveLoop(
            [-Nspline_lower, NEspline, Espline_lower], reorient=True
        )
        NEface_lower = gm.addSurfaceFilling([curve_loopNE_lower])
        curve_loopES_lower = gm.addCurveLoop(
            [-Espline_lower, ESspline, Sspline_lower], reorient=True
        )
        ESface_lower = gm.addSurfaceFilling([curve_loopES_lower])
        curve_loopSW_lower = gm.addCurveLoop(
            [-Sspline_lower, SWspline, Wspline_lower], reorient=True
        )
        SWface_lower = gm.addSurfaceFilling([curve_loopSW_lower])
        curve_loopWN_lower = gm.addCurveLoop(
            [-Wspline_lower, WNspline, Nspline_lower], reorient=True
        )
        WNface_lower = gm.addSurfaceFilling([curve_loopWN_lower])
        sl1 = gm.addSurfaceLoop(
            [
                NEface,
                NEface_lower,
                ESface,
                ESface_lower,
                SWface,
                SWface_lower,
                WNface,
                WNface_lower,
            ]
        )
        if self.ib:
            v1 = gm.addVolume([sl1])
            # Apply any rotations requested
            if self.rotatable:
                self._rotate(VOLUME, v1)
            gm.addPhysicalGroup(VOLUME, [v1])
        elif self.rotatable:
            self._rotate(SURFACE, sl1)

        gm.synchronize()

    def export(self):
        """Export .msh and .geo_unrolled (converted to .geo manually) files."""
        msh = gmsh.model.mesh
        gmsh.model.geo.synchronize()
        raw_file = f"{self.filepath}.geo_unrolled"
        new_file = f"{self.filepath}.geo"
        gmsh.write(raw_file)
        gmsh.option.setNumber("Mesh.Smoothing", 100)
        gmsh.option.setNumber("Mesh.SaveAll", 0)  # only save physical groups
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        gmsh.option.setNumber("General.NumThreads", 8)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
        # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1e-3)
        # gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 10)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.geo.synchronize()
        msh.generate(VOLUME)
        gmsh.model.geo.synchronize()
        gmsh.write(f"{self.filepath}.msh")
        import os

        os.replace(f"{raw_file}", f"{new_file}")
        if not self.view:
            gmsh.finalize()
        else:
            self._view_mesh()

    def _rotate(self, dim, tag):
        # TODO take dimtag as argument instead
        """Perform rotations around each axis as requested

        Parameters
        ----------
        volume : int
             Tag for the quadric volume. This is automatically populated by the call to rotate within self.draw().

        Notes
        -----
        Angles are in DEGREES and converted within this method."""
        # origin points
        x0 = 0
        y0 = 0
        z0 = 0
        gm = gmsh.model.geo
        for axis, angle in enumerate(self.rotation):
            # default axis values
            ax = 0
            ay = 0
            az = 0
            if axis == 0 and angle != 0:  # rotate around x-axis
                ax = 1
                gm.rotate([(dim, tag)], x0, y0, z0, ax, ay, az, angle * np.pi / 180)
            elif axis == 1 and angle != 0:  # rotate around y-axis
                ay = 1
                gm.rotate([(dim, tag)], x0, y0, z0, ax, ay, az, angle * np.pi / 180)
            elif axis == 2 and angle != 0:  # rotate around z-axis
                az = 1
                gm.rotate([(dim, tag)], x0, y0, z0, ax, ay, az, angle * np.pi / 180)

    def _view_mesh(self):
        """Visualise the mesh file"""
        if "-nopopup" not in sys.argv:
            gmsh.fltk.run()
        gmsh.finalize()


quad = Superquadric(indices=[8, 8, 8], scale=[6.5 / 1000, 6.5 / 1000, 3 / 1000])
quad.draw()
quad.export()
