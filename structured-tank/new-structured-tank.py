#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File    :   new-structured-tank.py
# Time    :   05/05/2023
# Author  :   Daniel Weston
# Version :   0.1.0
# Contact :   dtw545@student.bham.ac.uk

import gmsh
import sys
import numpy as np

# Global constants - these form the "dim" section of the dimtag parlance used in gmsh
POINT = 0
CURVE = 1
SURFACE = 2
VOLUME = 3


class StructuredStirredTank:
    def __init__(
        self,
        radius: float = 0.1,
        height: float = 0.1,
        axis_alignment: int = 2,
        baffle_width: float = 0.005,
        baffle_depth: float = 0.02,
        filepath: str = "stirred_tank",
        view: bool = True,
        height_spacing: float = 0.001,
        radial_coef: float = 1.0,
        height_coef: float = 1.0,
        radial_spacing: float = 0.005,
    ):
        self.filepath = filepath
        self.radius = radius
        self.height = height
        self.axis_alignment = axis_alignment
        self.baffle_angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        self.baffle_origins = np.zeros((3, 8))
        self.baffle_coordinates = np.zeros((3, 4, 4))
        self.baffle_width = baffle_width
        self.baffle_depth = baffle_depth
        self.view = view

    def draw(self):
        gmsh.initialize(sys.argv)
        gmsh.model.add("Stirred Tank")
        gm = gmsh.model.geo
        # define all the control points in the centre
        centre = gm.addPoint(0, 0, 0)
        # upper_centre = gm.addPoint(0, 0, self.height)
        self.baffle_origins[0, :] += self.radius * np.cos(self.baffle_angles)
        self.baffle_origins[1, :] += self.radius * np.sin(self.baffle_angles)
        # now assign 5 pts corresponding to the 4 corners of the baffle and its midpoint
        # from 0 rad:
        # 1 ----------- 0
        #   |
        #   |
        #   |
        # 2 ----------- 3
        baffle_pts = []
        outer_split_circle_pts = []
        inner_split_circle_pts = []
        centre_square_pts = []
        inside_baffle_pts = []
        outside_baffle_pts = []
        outer_square_pts = []
        for i in range(4):
            baffle = []
            unit_vec = self.baffle_origins[:, 2 * i] / np.sqrt(
                sum(self.baffle_origins[:, 2 * i] ** 2)
            )
            other_unit_vec = self.baffle_origins[:, 2 * i + 1] / np.sqrt(
                sum(self.baffle_origins[:, 2 * i + 1] ** 2)
            )
            normal_unit_vec = np.array([-unit_vec[1], unit_vec[0]])  # only need xy part
            self.baffle_coordinates[0:2, i, 0] = (
                self.baffle_origins[0:2, 2 * i]
                + 0.5 * self.baffle_width * normal_unit_vec
            )
            self.baffle_coordinates[0:2, i, 1] = (
                self.baffle_origins[0:2, 2 * i]
                + 0.5 * self.baffle_width * normal_unit_vec
                - unit_vec[0:2] * self.baffle_depth
            )
            self.baffle_coordinates[0:2, i, 2] = (
                self.baffle_origins[0:2, 2 * i]
                - 0.5 * self.baffle_width * normal_unit_vec
                - unit_vec[0:2] * self.baffle_depth
            )
            self.baffle_coordinates[0:2, i, 3] = (
                self.baffle_origins[0:2, 2 * i]
                - 0.5 * self.baffle_width * normal_unit_vec
            )
            outer_split_circle_pts.append(
                gm.addPoint(*self.baffle_origins[:, 2 * i + 1])
            )
            self.baffle_origins[0:2, 2 * i + 1] -= (
                other_unit_vec[0:2] * self.baffle_depth
            )
            inner_split_circle_pts.append(
                gm.addPoint(*self.baffle_origins[:, 2 * i + 1])
            )
            # assign points
            for j in range(4):
                baffle.append(
                    gm.addPoint(
                        self.baffle_coordinates[0, i, j],
                        self.baffle_coordinates[1, i, j],
                        self.baffle_coordinates[2, i, j],
                    )
                )
            baffle_pts.append(baffle)
            self.baffle_coordinates[0:2, i, 1] = (
                self.baffle_origins[0:2, 2 * i]
                + 0.5 * self.baffle_width * normal_unit_vec
                - unit_vec[0:2] * self.baffle_depth
            )
            centre_square_pts.append(
                gm.addPoint(
                    x=self.baffle_width
                    * 0.5
                    * np.sqrt(2)
                    * np.cos(self.baffle_angles[2 * i + 1]),
                    y=self.baffle_width
                    * 0.5
                    * np.sqrt(2)
                    * np.sin(self.baffle_angles[2 * i + 1]),
                    z=0,
                )
            )
            square_coords = self.baffle_origins[:, 2 * i + 1] / 2
            outer_square_pts.append(gm.addPoint(*(square_coords)))
            if np.tan(self.baffle_angles[2 * i + 1]) > 0:
                outside_baffle_pts.append(
                    gm.addPoint(
                        x=square_coords[0],
                        y=self.baffle_width
                        * 0.5
                        * np.sqrt(2)
                        * np.sin(self.baffle_angles[2 * i + 1]),
                        z=0,
                    )
                )
                inside_baffle_pts.append(
                    gm.addPoint(
                        x=square_coords[0],
                        y=self.baffle_width
                        * 0.5
                        * np.sqrt(2)
                        * np.sin(self.baffle_angles[(2 * i - 1) % 8]),
                        z=0,
                    )
                )
            else:
                outside_baffle_pts.append(
                    gm.addPoint(
                        x=self.baffle_width
                        * 0.5
                        * np.sqrt(2)
                        * np.cos(self.baffle_angles[2 * i + 1]),
                        y=square_coords[1],
                        z=0,
                    )
                )
                inside_baffle_pts.append(
                    gm.addPoint(
                        x=self.baffle_width
                        * 0.5
                        * np.sqrt(2)
                        * np.cos(self.baffle_angles[(2 * i - 1) % 8]),
                        y=square_coords[1],
                        z=0,
                    )
                )

        # curve drawing
        outer_curves = []
        inner_curves = []
        outer_connecting_curves = []
        inner_connecting_curves = []
        baffle_inside_edges = []
        baffle_outside_edges = []
        baffle_extended_edges = []
        centre_square_edges = []
        outer_inside_straight_edges = []
        outer_outside_straight_edges = []
        inner_inside_straight_edges = []
        inner_outside_straight_edges = []
        outer_square_vertical_edges = []
        outer_square_horizontal_edges = []
        middle_straight_edges = []
        for i in range(4):
            # draw outer curves
            outer_curves.append(
                gm.addCircleArc(baffle_pts[i][0], centre, outer_split_circle_pts[i])
            )
            outer_curves.append(
                gm.addCircleArc(
                    outer_split_circle_pts[i], centre, baffle_pts[(i + 1) % 4][3]
                )
            )
            # draw inner curves
            inner_curves.append(
                gm.addCircleArc(baffle_pts[i][1], centre, inner_split_circle_pts[i])
            )
            inner_curves.append(
                gm.addCircleArc(
                    inner_split_circle_pts[i], centre, baffle_pts[(i + 1) % 4][2]
                )
            )
            # draw curves connecting outer and inner curve midpoints
            outer_connecting_curves.append(
                gm.addLine(outer_split_circle_pts[i], inner_split_circle_pts[i])
            )
            inner_connecting_curves.append(
                gm.addLine(inner_split_circle_pts[i], outer_square_pts[i])
            )
            # draw inside edge baffle lines
            baffle_inside_edges.append(gm.addLine(baffle_pts[i][3], baffle_pts[i][2]))
            # draw outside edge baffle lines
            baffle_outside_edges.append(gm.addLine(baffle_pts[i][1], baffle_pts[i][0]))
            # draw extended baffle edge
            baffle_extended_edges.append(gm.addLine(baffle_pts[i][2], baffle_pts[i][1]))
            centre_square_edges.append(
                gm.addLine(centre_square_pts[i], centre_square_pts[(i + 1) % 4])
            )
            outer_inside_straight_edges.append(
                gm.addLine(baffle_pts[i][2], inside_baffle_pts[i])
            )
            outer_outside_straight_edges.append(
                gm.addLine(baffle_pts[i][1], outside_baffle_pts[i])
            )
            inner_inside_straight_edges.append(
                gm.addLine(centre_square_pts[(i - 1) % 4], inside_baffle_pts[i])
            )
            inner_outside_straight_edges.append(
                gm.addLine(centre_square_pts[i], outside_baffle_pts[i])
            )
            outer_square_vertical_edges.append(
                gm.addLine(outside_baffle_pts[i], outer_square_pts[i])
            )
            outer_square_horizontal_edges.append(
                gm.addLine(inside_baffle_pts[i], outer_square_pts[(i - 1) % 4])
            )
            middle_straight_edges.append(
                gm.addLine(inside_baffle_pts[i], outside_baffle_pts[i])
            )

        # surfaces
        outer_curved_surfaces = []
        inner_curved_surfaces = []
        outer_straight_surfaces = []
        inner_straight_surfaces = []
        square_surfaces = []
        inner_square_surfaces = []
        for i in range(4):
            outer_loop = gm.addCurveLoop(
                [
                    baffle_outside_edges[i],
                    inner_curves[2 * i],
                    outer_connecting_curves[i],
                    outer_curves[2 * i],
                ],
                reorient=True,
            )
            outer_curved_surfaces.append(gm.addPlaneSurface([outer_loop]))
            outer_loop = gm.addCurveLoop(
                [
                    outer_connecting_curves[i],
                    inner_curves[2 * i + 1],
                    baffle_inside_edges[(i + 1) % 4],
                    outer_curves[2 * i + 1],
                ],
                reorient=True,
            )
            outer_curved_surfaces.append(gm.addPlaneSurface([outer_loop]))

            inner_loop = gm.addCurveLoop(
                [
                    outer_outside_straight_edges[i],
                    outer_square_vertical_edges[i],
                    inner_connecting_curves[i],
                    inner_curves[2 * i],
                ],
                reorient=True,
            )
            inner_curved_surfaces.append(gm.addPlaneSurface([inner_loop]))
            inner_loop = gm.addCurveLoop(
                [
                    inner_connecting_curves[i],
                    outer_square_horizontal_edges[(i + 1) % 4],
                    outer_inside_straight_edges[(i + 1) % 4],
                    inner_curves[2 * i + 1],
                ],
                reorient=True,
            )
            inner_curved_surfaces.append(gm.addPlaneSurface([inner_loop]))

            outer_loop = gm.addCurveLoop(
                [
                    baffle_extended_edges[i],
                    outer_inside_straight_edges[i],
                    middle_straight_edges[i],
                    outer_outside_straight_edges[i],
                ],
                reorient=True,
            )
            outer_straight_surfaces.append(gm.addPlaneSurface([outer_loop]))

            inner_loop = gm.addCurveLoop(
                [
                    middle_straight_edges[i],
                    inner_inside_straight_edges[i],
                    centre_square_edges[(i - 1) % 4],
                    inner_outside_straight_edges[i],
                ],
                reorient=True,
            )
            inner_straight_surfaces.append(gm.addPlaneSurface([inner_loop]))

            square_loop = gm.addCurveLoop(
                [
                    inner_outside_straight_edges[i],
                    inner_inside_straight_edges[(i + 1) % 4],
                    outer_square_horizontal_edges[(i + 1) % 4],
                    outer_square_vertical_edges[i],
                ],
                reorient=True,
            )
            square_surfaces.append(gm.addPlaneSurface([square_loop]))

        inner_square_surfaces.append(
            gm.addPlaneSurface(
                [gm.addCurveLoop([line for line in centre_square_edges], reorient=True)]
            )
        )

        gm.synchronize()
        msh = gm.mesh
        # transfinite definitions
        for _, line in gmsh.model.getEntities(CURVE):
            msh.setTransfiniteCurve(line, nPoints=5)

        for _, surface in gmsh.model.getEntities(SURFACE):
            msh.setTransfiniteSurface(surface)
            msh.setRecombine(SURFACE, surface)

    def export(self):
        gmsh.model.geo.synchronize()

        # gmsh.model.mesh.generate(VOLUME)
        # gmsh.model.mesh.setOrder(2)
        gmsh.model.geo.synchronize()
        gmsh.write(f"{self.filepath}.msh")
        if self.view:
            self._view()

    def _view(self):
        if "-nopopup" not in sys.argv:
            gmsh.fltk.run()


tank = StructuredStirredTank()
tank.draw()
tank.export()
