#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File    :   structured-tank.py
# Time    :   03/03/2023
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
        self.baffle_angles = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        self.baffle_origins = np.zeros((3, 4))
        self.baffle_coordinates = np.zeros((3, 4, 4))
        self.baffle_width = baffle_width
        self.baffle_depth = baffle_depth
        self.view = view
        self.radial_coef = radial_coef
        self.height_coef = height_coef
        # set number of points based on spacing, all values must end up >= 1 otherwise Bad Things(tm) occur
        self.outer_curve_npts = (
            int((0.5 * np.pi * radius - baffle_width) / radial_spacing)
            if int((0.5 * np.pi * radius - baffle_width) / radial_spacing) > 0
            else 1
        )
        self.baffle_npts = (
            int(baffle_depth / radial_spacing)
            if int(baffle_depth / radial_spacing) > 0
            else 1
        )
        half_angle = np.pi * 135 / (2 * 180)
        self.baffle_to_octagon_npts = (
            int(
                (
                    self.radius
                    - self.baffle_depth
                    - np.tan(half_angle) * self.baffle_width / 2
                )
                / radial_spacing
            )
            if int(
                (
                    self.radius
                    - self.baffle_depth
                    - np.tan(half_angle) * self.baffle_width / 2
                )
                / radial_spacing
            )
            > 0
            else 1
        )
        self.baffle_extended_edge_npts = (
            int(baffle_width / radial_spacing)
            if int(baffle_width / radial_spacing) > 0
            else 1
        )
        self.vertical_npts = (
            int(height / height_spacing) if int(height / height_spacing) > 0 else 1
        )
        self.octagon_segment_npts = (
            int(np.tan(half_angle) * self.baffle_width / 2)
            if int(np.tan(half_angle) * self.baffle_width / 2) > 0
            else 1
        )

    def draw(self):
        gmsh.initialize(sys.argv)
        gmsh.model.add("Stirred Tank")
        gm = gmsh.model.geo
        # define all the control points in the centre
        centre = gm.addPoint(0, 0, 0)
        upper_centre = gm.addPoint(0, 0, self.height)
        self.baffle_origins[0, :] += self.radius * np.cos(self.baffle_angles)
        self.baffle_origins[1, :] += self.radius * np.sin(self.baffle_angles)
        # now assign 4 pts corresponding to the 4 corners of the baffle
        # from 0 rad:
        # 1 ----------- 0
        #   |
        #   |
        # 2 ----------- 3
        baffle_pts = []
        for i in range(4):
            baffle = []
            unit_vec = self.baffle_origins[:, i] / np.sqrt(
                sum(self.baffle_origins[:, i] ** 2)
            )
            normal_unit_vec = np.array([-unit_vec[1], unit_vec[0]])  # only need xy part
            self.baffle_coordinates[0:2, i, 0] = (
                self.baffle_origins[0:2, i] + 0.5 * self.baffle_width * normal_unit_vec
            )
            self.baffle_coordinates[0:2, i, 1] = (
                self.baffle_origins[0:2, i]
                + 0.5 * self.baffle_width * normal_unit_vec
                - unit_vec[0:2] * self.baffle_depth
            )
            self.baffle_coordinates[0:2, i, 2] = (
                self.baffle_origins[0:2, i]
                - 0.5 * self.baffle_width * normal_unit_vec
                - unit_vec[0:2] * self.baffle_depth
            )
            self.baffle_coordinates[0:2, i, 3] = (
                self.baffle_origins[0:2, i] - 0.5 * self.baffle_width * normal_unit_vec
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

        # draw central octagon
        # the points are 1 and 2 from the baffles, located baffle_width * tan(67.5) / 2 from
        # the centre - this ensures the octagon is regular
        half_angle = np.pi * 135 / (2 * 180)
        octagon_pts = []
        for i in range(4):
            octagon_pts.append(
                gm.addPoint(
                    x=self.baffle_coordinates[0, i, 2]
                    - (
                        self.radius
                        - self.baffle_depth
                        - np.tan(half_angle) * self.baffle_width / 2
                    )
                    * np.cos(i * 0.5 * np.pi),
                    y=self.baffle_coordinates[1, i, 2]
                    - (
                        self.radius
                        - self.baffle_depth
                        - np.tan(half_angle) * self.baffle_width / 2
                    )
                    * np.sin(i * 0.5 * np.pi),
                    z=self.baffle_coordinates[2, i, 2],
                )
            )
            octagon_pts.append(
                gm.addPoint(
                    x=self.baffle_coordinates[0, i, 1]
                    - (
                        self.radius
                        - self.baffle_depth
                        - np.tan(half_angle) * self.baffle_width / 2
                    )
                    * np.cos(i * 0.5 * np.pi),
                    y=self.baffle_coordinates[1, i, 1]
                    - (
                        self.radius
                        - self.baffle_depth
                        - np.tan(half_angle) * self.baffle_width / 2
                    )
                    * np.sin(i * 0.5 * np.pi),
                    z=self.baffle_coordinates[2, i, 1],
                )
            )

        # curve drawing
        # sides of octagon
        # octagon_edges = []
        # for i, _ in enumerate(octagon_pts):
        #     octagon_edges.append(gm.addLine(octagon_pts[i], octagon_pts[(i + 1) % 8]))
        # rest of surface
        outer_curves = []
        inner_curves = []
        baffle_extended_edges = []
        baffle_inside_edges = []
        baffle_outside_edges = []
        # baffle_to_octagon_inside = []
        # baffle_to_octagon_outside = []
        # octagon_segments = []
        kite_lines = []
        for i in range(4):
            # draw outer curves
            outer_curves.append(
                gm.addCircleArc(baffle_pts[i][0], centre, baffle_pts[(i + 1) % 4][3])
            )
            # draw inner curves
            inner_curves.append(
                gm.addCircleArc(baffle_pts[i][1], centre, baffle_pts[(i + 1) % 4][2])
            )
            # draw extended baffle lines
            baffle_extended_edges.append(gm.addLine(baffle_pts[i][2], baffle_pts[i][1]))
            # draw inside edge baffle lines
            baffle_inside_edges.append(gm.addLine(baffle_pts[i][3], baffle_pts[i][2]))
            # draw outside edge baffle lines
            baffle_outside_edges.append(gm.addLine(baffle_pts[i][1], baffle_pts[i][0]))
            # # draw inside edge connecting lines to octagon
            # baffle_to_octagon_inside.append(gm.addLine(baffle_pts[i][2], octagon_pts[2 * i]))
            # # draw outside edge connecting lines to octagon
            # baffle_to_octagon_outside.append(gm.addLine(baffle_pts[i][1], octagon_pts[2 * i + 1]))
            # # draw octagon segments
            # octagon_segments.append(gm.addLine(octagon_pts[2 * i], centre))
            kite_lines.append(gm.addLine(baffle_pts[i][2], centre))

        # surfaces
        outer_surfaces = []
        inner_surfaces = []
        straight_surfaces = []
        octagon_surfaces = []
        kite_surfaces = []
        for i in range(4):
            # curved surfaces between baffles
            outer_loop = gm.addCurveLoop(
                [
                    baffle_outside_edges[i],
                    outer_curves[i],
                    baffle_inside_edges[(i + 1) % 4],
                    inner_curves[i],
                ],
                reorient=True,
            )
            outer_surfaces.append(gm.addPlaneSurface([outer_loop]))
            # inner surface
            # inner_loop = gm.addCurveLoop([octagon_edges[2 * i + 1], baffle_to_octagon_inside[(i + 1) % 4], inner_curves[i], baffle_to_octagon_outside[i]], reorient = True)
            # inner_surfaces.append(gm.addPlaneSurface([inner_loop]))
            # straight surfaces from baffle to octagon
            # straight_loop = gm.addCurveLoop([baffle_to_octagon_inside[i], octagon_edges[2 * i], baffle_to_octagon_outside[i], baffle_extended_edges[i]], reorient= True)
            # straight_surfaces.append(gm.addPlaneSurface([straight_loop]))
            # inner octagon surfaces
            # octagon_loop = gm.addCurveLoop([octagon_edges[2 * i], octagon_segments[i], octagon_segments[(i + 1)%4], octagon_edges[2 * i + 1]], reorient= True)
            # octagon_surfaces.append(gm.addPlaneSurface([octagon_loop]))
            kite_loop = gm.addCurveLoop(
                [
                    baffle_extended_edges[i],
                    kite_lines[i],
                    kite_lines[(i + 1) % 4],
                    inner_curves[i],
                ],
                reorient=True,
            )
            kite_surfaces.append(gm.addPlaneSurface([kite_loop]))

        gm.synchronize()

        # # copy and move points to +z direction
        # upper_baffle_pts = []
        # for baffle in baffle_pts:
        #     upper_baffle = []
        #     for pt in baffle:
        #         new_pt = gm.copy([(POINT, pt)])
        #         gm.translate(new_pt, dx = 0, dy = 0, dz = self.height)
        #         upper_baffle.append(new_pt[0][1])
        #     upper_baffle_pts.append(upper_baffle)
        # upper_octagon_pts = []
        # for pt in octagon_pts:
        #     new_pt = gm.copy([(POINT, pt)])
        #     gm.translate(new_pt, dx = 0, dy = 0, dz = self.height)
        #     upper_octagon_pts.append(new_pt[0][1])
        # # draw curves
        # # sides of octagon
        # upper_octagon_edges = []
        # for i, _ in enumerate(octagon_pts):
        #     upper_octagon_edges.append(gm.addLine(upper_octagon_pts[i], upper_octagon_pts[(i + 1) % 8]))
        # # rest of surface
        # upper_outer_curves = []
        # upper_inner_curves = []
        # upper_baffle_extended_edges = []
        # upper_baffle_inside_edges = []
        # upper_baffle_outside_edges = []
        # upper_baffle_to_octagon_inside = []
        # upper_baffle_to_octagon_outside = []
        # upper_octagon_segments = []
        # for i in range(4):
        #     # draw outer curves
        #     upper_outer_curves.append(gm.addCircleArc(upper_baffle_pts[i][0], upper_centre, upper_baffle_pts[(i + 1) % 4][3]))
        #     # draw inner curves
        #     upper_inner_curves.append(gm.addCircleArc(upper_baffle_pts[i][1], upper_centre, upper_baffle_pts[(i + 1) % 4][2]))
        #     # draw extended baffle lines
        #     upper_baffle_extended_edges.append(gm.addLine(upper_baffle_pts[i][2], upper_baffle_pts[i][1]))
        #     # draw inside edge baffle lines
        #     upper_baffle_inside_edges.append(gm.addLine(upper_baffle_pts[i][3], upper_baffle_pts[i][2]))
        #     # draw outside edge baffle lines
        #     upper_baffle_outside_edges.append(gm.addLine(upper_baffle_pts[i][1], upper_baffle_pts[i][0]))
        #     # draw inside edge connecting lines to octagon
        #     upper_baffle_to_octagon_inside.append(gm.addLine(upper_baffle_pts[i][2], upper_octagon_pts[2 * i]))
        #     # draw outside edge connecting lines to octagon
        #     upper_baffle_to_octagon_outside.append(gm.addLine(upper_baffle_pts[i][1], upper_octagon_pts[2 * i + 1]))
        #     # draw octagon segments
        #     upper_octagon_segments.append(gm.addLine(upper_octagon_pts[2 * i], upper_centre))

        # # draw surfaces
        # upper_outer_surfaces = []
        # upper_inner_surfaces = []
        # upper_straight_surfaces = []
        # upper_octagon_surfaces = []
        # for i in range(4):
        #     # curved surfaces between baffles
        #     outer_loop = gm.addCurveLoop([upper_baffle_outside_edges[i], upper_outer_curves[i], upper_baffle_inside_edges[(i + 1) % 4], upper_inner_curves[i]], reorient= True)
        #     upper_outer_surfaces.append(gm.addPlaneSurface([outer_loop]))
        #     # inner surface
        #     inner_loop = gm.addCurveLoop([upper_octagon_edges[2 * i + 1], upper_baffle_to_octagon_inside[(i + 1) % 4], upper_inner_curves[i], upper_baffle_to_octagon_outside[i]], reorient = True)
        #     upper_inner_surfaces.append(gm.addPlaneSurface([inner_loop]))
        #     # straight surfaces from baffle to octagon
        #     straight_loop = gm.addCurveLoop([upper_baffle_to_octagon_inside[i], upper_octagon_edges[2 * i], upper_baffle_to_octagon_outside[i], upper_baffle_extended_edges[i]], reorient= True)
        #     upper_straight_surfaces.append(gm.addPlaneSurface([straight_loop]))
        #     # inner octagon surfaces
        #     octagon_loop = gm.addCurveLoop([upper_octagon_edges[2 * i], upper_octagon_segments[i], upper_octagon_segments[(i + 1)%4], upper_octagon_edges[2 * i + 1]], reorient= True)
        #     upper_octagon_surfaces.append(gm.addPlaneSurface([octagon_loop]))

        # draw lines to connect -z and +z faces
        # for i in range(4):

        # draw connecting surfaces

        # draw volumes

        # define transfinite characteristics
        # TODO investigate sensible coef values for "bump" mode
        msh = gmsh.model.geo.mesh
        # TODO logically separate nPoints for the different size scales present
        for i in range(4):
            # curved lines
            msh.setTransfiniteCurve(outer_curves[i], nPoints=self.outer_curve_npts)
            msh.setTransfiniteCurve(inner_curves[i], nPoints=self.outer_curve_npts)

            # msh.setTransfiniteCurve(upper_outer_curves[i], nPoints=self.outer_curve_npts)
            # msh.setTransfiniteCurve(upper_inner_curves[i], nPoints=self.outer_curve_npts)
            # straight sections
            msh.setTransfiniteCurve(
                baffle_extended_edges[i], nPoints=self.outer_curve_npts
            )
            msh.setTransfiniteCurve(baffle_inside_edges[i], nPoints=self.baffle_npts)
            msh.setTransfiniteCurve(baffle_outside_edges[i], nPoints=self.baffle_npts)
            # msh.setTransfiniteCurve(baffle_to_octagon_inside[i], nPoints=self.baffle_to_octagon_npts)
            # msh.setTransfiniteCurve(baffle_to_octagon_outside[i], nPoints=self.baffle_to_octagon_npts)
            # msh.setTransfiniteCurve(octagon_edges[2*i+1], nPoints=self.outer_curve_npts)
            # msh.setTransfiniteCurve(octagon_edges[2*i], nPoints=self.outer_curve_npts)

            # msh.setTransfiniteCurve(upper_baffle_extended_edges[i], nPoints=self.outer_curve_npts)
            # msh.setTransfiniteCurve(upper_baffle_inside_edges[i], nPoints=self.baffle_npts)
            # msh.setTransfiniteCurve(upper_baffle_outside_edges[i], nPoints=self.baffle_npts)
            # msh.setTransfiniteCurve(upper_baffle_to_octagon_inside[i], nPoints=self.baffle_to_octagon_npts)
            # msh.setTransfiniteCurve(upper_baffle_to_octagon_outside[i], nPoints=self.baffle_to_octagon_npts)
            # msh.setTransfiniteCurve(upper_octagon_edges[2*i+1], nPoints=self.outer_curve_npts)
            # msh.setTransfiniteCurve(upper_octagon_edges[2*i], nPoints=self.outer_curve_npts)
            # octagon segments
            # msh.setTransfiniteCurve(octagon_segments[i], nPoints=self.outer_curve_npts)

            # msh.setTransfiniteCurve(upper_octagon_segments[i], nPoints=self.outer_curve_npts)

            # connecting lines

        for i in range(4):
            # -z face surfaces
            msh.setTransfiniteSurface(outer_surfaces[i])
            msh.setRecombine(SURFACE, outer_surfaces[i])
            # msh.setTransfiniteSurface(inner_surfaces[i])
            # msh.setRecombine(SURFACE, inner_surfaces[i])
            # msh.setTransfiniteSurface(straight_surfaces[i])
            # msh.setRecombine(SURFACE, straight_surfaces[i])
            # msh.setTransfiniteSurface(octagon_surfaces[i])
            # msh.setRecombine(SURFACE, octagon_surfaces[i])
            # +z face surfaces
            # msh.setTransfiniteSurface(upper_outer_surfaces[i])
            # msh.setRecombine(SURFACE, upper_outer_surfaces[i])
            # msh.setTransfiniteSurface(upper_inner_surfaces[i])
            # msh.setRecombine(SURFACE, upper_inner_surfaces[i])
            # msh.setTransfiniteSurface(upper_straight_surfaces[i])
            # msh.setRecombine(SURFACE, upper_straight_surfaces[i])
            # msh.setTransfiniteSurface(upper_octagon_surfaces[i])
            # msh.setRecombine(SURFACE, upper_octagon_surfaces[i])
            # connecting surfaces

        # define physical surfaces and volume

    def export(self):
        gmsh.model.geo.synchronize()
        raw_file = f"{self.filepath}.geo_unrolled"
        new_file = f"{self.filepath}.geo"
        gmsh.write(raw_file)
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
