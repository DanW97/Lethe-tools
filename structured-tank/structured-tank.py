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
    """
    Class responsible for generating a structured stirred tank mesh from given
    dimensions.
    """

    def __init__(
        self,
        radius: float = 0.1,
        height: float = 0.1,
        axis_alignment: int = 2,
        baffle_width: float = 0.005,
        baffle_depth: float = 0.02,
        filepath: str = "stirred_tank",
        view: bool = True,
        height_spacing: float = 0.005,
        radial_spacing: float = 0.005,
    ):
        """
        Create a stirred tank with 4 solid baffles based on user defined dimensions.

        *ALL DIMENSIONS ARE ASSUMED TO BE IN METRES*

        Parameters
        ----------
        radius : float, optional
            Tank inner radius, by default 0.1
        height : float, optional
            *Fill* height (unless you are doing free surface simulations,
            in which case then this can be the tank height), by default 0.1
        axis_alignment : int, optional
            Axis that the tank is aligned to, numbered from 0-2 for x-z respectively,
            by default 2
        baffle_width : float, optional
            Baffle thickness, by default 0.005
        baffle_depth : float, optional
            Extension of baffle towards tank centre,
            by default 0.02
        filepath : str, optional
            Path to save output files to, by default "stirred_tank"
        view : bool, optional
            Boolean to check if mesh is visualised after creation,
            by default True
        height_spacing : float, optional
            Vertical cell spacing, by default 0.005
        radial_spacing : float, optional
            Horizontal cell spacing, by default 0.005
        """
        self.filepath = filepath
        self.radius = radius
        self.height = height
        self.axis_alignment = axis_alignment
        self.baffle_angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        self.baffle_origins = np.zeros((3, 8))
        self.baffle_origins[0, :] += self.radius * np.cos(self.baffle_angles)
        self.baffle_origins[1, :] += self.radius * np.sin(self.baffle_angles)
        self.baffle_coordinates = np.zeros((3, 4, 4))
        self.baffle_width = baffle_width
        self.baffle_depth = baffle_depth
        self.view = view
        self.height_npts = (
            int(height / height_spacing) + 1 if int(height / height_spacing) > 0 else 1
        )
        self.baffle_front_spacing = (
            int(baffle_width / radial_spacing) + 1
            if int(baffle_width / radial_spacing) > 0
            else 1
        )
        # (2*pi*r / 4 - baffle_width) / 2: arc length of each of the 8 segments
        segment_length = 0.5 * (0.5 * radius * np.pi - baffle_width)
        self.arc_spacing = (
            int(segment_length / radial_spacing) + 1
            if int(segment_length / radial_spacing) > 0
            else 1
        )
        self.baffle_edge_spacing = (
            int(baffle_depth / radial_spacing) + 1
            if int(baffle_depth / radial_spacing) > 0
            else 1
        )
        connector_length = (
            radius - baffle_depth - np.sum((self.baffle_origins[:, 1] / 2) ** 2)
        )
        self.straight_edge_spacing = (
            int(connector_length / radial_spacing)
            if int(connector_length / radial_spacing) > 0
            else 1
        )

    def draw(self):
        gmsh.initialize(sys.argv)
        gmsh.model.add("Stirred Tank")
        gm = gmsh.model.geo
        msh = gm.mesh
        centre = gm.addPoint(0, 0, 0)
        # now assign 4 pts corresponding to the 4 corners of the baffle
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
                    -baffle_outside_edges[i],
                    inner_curves[2 * i],
                    -outer_connecting_curves[i],
                    -outer_curves[2 * i],
                ],
            )
            outer_curved_surfaces.append(gm.addPlaneSurface([outer_loop]))
            outer_loop = gm.addCurveLoop(
                [
                    outer_connecting_curves[i],
                    inner_curves[2 * i + 1],
                    -baffle_inside_edges[(i + 1) % 4],
                    -outer_curves[2 * i + 1],
                ],
            )
            outer_curved_surfaces.append(gm.addPlaneSurface([outer_loop]))

            inner_loop = gm.addCurveLoop(
                [
                    outer_outside_straight_edges[i],
                    outer_square_vertical_edges[i],
                    -inner_connecting_curves[i],
                    -inner_curves[2 * i],
                ],
            )
            inner_curved_surfaces.append(gm.addPlaneSurface([inner_loop]))
            inner_loop = gm.addCurveLoop(
                [
                    inner_connecting_curves[i],
                    -outer_square_horizontal_edges[(i + 1) % 4],
                    -outer_inside_straight_edges[(i + 1) % 4],
                    -inner_curves[2 * i + 1],
                ],
            )
            inner_curved_surfaces.append(gm.addPlaneSurface([inner_loop]))

            outer_loop = gm.addCurveLoop(
                [
                    -baffle_extended_edges[i],
                    outer_inside_straight_edges[i],
                    middle_straight_edges[i],
                    -outer_outside_straight_edges[i],
                ],
            )
            outer_straight_surfaces.append(gm.addPlaneSurface([outer_loop]))

            inner_loop = gm.addCurveLoop(
                [
                    -middle_straight_edges[i],
                    -inner_inside_straight_edges[i],
                    centre_square_edges[(i - 1) % 4],
                    inner_outside_straight_edges[i],
                ],
            )
            inner_straight_surfaces.append(gm.addPlaneSurface([inner_loop]))

            square_loop = gm.addCurveLoop(
                [
                    -inner_outside_straight_edges[i],
                    inner_inside_straight_edges[(i + 1) % 4],
                    outer_square_horizontal_edges[(i + 1) % 4],
                    -outer_square_vertical_edges[i],
                ],
            )
            square_surfaces.append(gm.addPlaneSurface([square_loop]))

        inner_square_surfaces.append(
            gm.addPlaneSurface(
                [gm.addCurveLoop([line for line in centre_square_edges])]
            )
        )

        gm.synchronize()

        # test hypothesis that the transfinite behaviour is mirrored on extrude
        # there are up to 4 unique values for nPoints for these curves:
        # 1. across baffle front, and centre square
        # 2. across circle arcs and outer squares
        # 3. along lines connecting inner and outer circle arcs
        # 4. along lines connecting inner circle arcs and outer square, affecting
        #    the straight outer sections too
        # inner straight sections, baffle edges and central square are spoken for
        # already by the above constraints
        for i in range(4):
            # implement 1...
            #     for baffle front
            gm.mesh.setTransfiniteCurve(
                baffle_extended_edges[i], nPoints=self.baffle_front_spacing
            )
            #     for straight section divider
            gm.mesh.setTransfiniteCurve(
                middle_straight_edges[i], nPoints=self.baffle_front_spacing
            )
            #     for central square
            gm.mesh.setTransfiniteCurve(
                centre_square_edges[i], nPoints=self.baffle_front_spacing
            )

            # implement 2...
            #     for outer circle arcs
            gm.mesh.setTransfiniteCurve(outer_curves[2 * i], nPoints=self.arc_spacing)
            gm.mesh.setTransfiniteCurve(
                outer_curves[2 * i + 1], nPoints=self.arc_spacing
            )
            #     for inner circle arcs
            gm.mesh.setTransfiniteCurve(inner_curves[2 * i], nPoints=self.arc_spacing)
            gm.mesh.setTransfiniteCurve(
                inner_curves[2 * i + 1], nPoints=self.arc_spacing
            )
            #     for outer square edges
            gm.mesh.setTransfiniteCurve(
                outer_square_vertical_edges[i], nPoints=self.arc_spacing
            )
            gm.mesh.setTransfiniteCurve(
                outer_square_horizontal_edges[i], nPoints=self.arc_spacing
            )
            #     for inner straight sections
            gm.mesh.setTransfiniteCurve(
                inner_inside_straight_edges[i], nPoints=self.arc_spacing
            )
            gm.mesh.setTransfiniteCurve(
                inner_outside_straight_edges[i], nPoints=self.arc_spacing
            )

            # implement 3...
            #     for circle arc connectors
            gm.mesh.setTransfiniteCurve(
                outer_connecting_curves[i], nPoints=self.baffle_edge_spacing
            )
            #     for baffle edges
            gm.mesh.setTransfiniteCurve(
                baffle_inside_edges[i], nPoints=self.baffle_edge_spacing
            )
            gm.mesh.setTransfiniteCurve(
                baffle_outside_edges[i], nPoints=self.baffle_edge_spacing
            )
            # implement 4...
            #     for circle arc - square connectors
            gm.mesh.setTransfiniteCurve(
                inner_connecting_curves[i], nPoints=self.straight_edge_spacing
            )
            #     for outer straight sections
            gm.mesh.setTransfiniteCurve(
                outer_inside_straight_edges[i], nPoints=self.straight_edge_spacing
            )
            gm.mesh.setTransfiniteCurve(
                outer_outside_straight_edges[i], nPoints=self.straight_edge_spacing
            )

        # rotate if required
        # gm.synchronize()
        if self.axis_alignment == 0:
            gm.rotate(
                gmsh.model.getEntities(SURFACE),
                x=0,
                y=0,
                z=0,
                ax=0,
                ay=1,
                az=0,
                angle=-np.pi / 2,
            )
            dx = self.height
            dy = 0
            dz = 0
        elif self.axis_alignment == 1:
            gm.rotate(
                gmsh.model.getEntities(SURFACE),
                x=0,
                y=0,
                z=0,
                ax=1,
                ay=0,
                az=0,
                angle=-np.pi / 2,
            )
            dx = 0
            dy = self.height
            dz = 0
        elif self.axis_alignment == 2:
            dx = 0
            dy = 0
            dz = self.height
        else:
            raise ValueError("Valid axis_alignment values are 0, 1 or 2 only!")

        # now extrude them

        outer_curved_extrude = [
            gm.extrude(
                [(SURFACE, surface)],
                dx=dx,
                dy=dy,
                dz=dz,
                numElements=[self.height_npts],
                recombine=True,
            )
            for surface in outer_curved_surfaces
        ]
        inner_curved_extrude = [
            gm.extrude(
                [(SURFACE, surface)],
                dx=dx,
                dy=dy,
                dz=dz,
                numElements=[self.height_npts],
                recombine=True,
            )
            for surface in inner_curved_surfaces
        ]
        outer_straight_extrude = [
            gm.extrude(
                [(SURFACE, surface)],
                dx=dx,
                dy=dy,
                dz=dz,
                numElements=[self.height_npts],
                recombine=True,
            )
            for surface in outer_straight_surfaces
        ]
        inner_straight_extrude = [
            gm.extrude(
                [(SURFACE, surface)],
                dx=dx,
                dy=dy,
                dz=dz,
                numElements=[self.height_npts],
                recombine=True,
            )
            for surface in inner_straight_surfaces
        ]
        square_extrude = [
            gm.extrude(
                [(SURFACE, surface)],
                dx=dx,
                dy=dy,
                dz=dz,
                numElements=[self.height_npts],
                recombine=True,
            )
            for surface in square_surfaces
        ]
        inner_square_extrude = gm.extrude(
            [(SURFACE, surface) for surface in inner_square_surfaces],
            dx=dx,
            dy=dy,
            dz=dz,
            numElements=[self.height_npts],
            recombine=True,
        )

        gm.synchronize()

        # transfinite surface definitions
        for _, surface in gmsh.model.getEntities(SURFACE):
            msh.setTransfiniteSurface(surface)
            msh.setRecombine(SURFACE, surface)

        # we got 2 physical surfaces - noslip on tank surfaces, and weakly imposed
        # noslip at top of tank. To find the surfaces from the extrudes, they are
        # in the same order as the curve loops are defined, given that the curves
        # have been manually oriented.

        # id 1, noslip, this is relevant for the outer_curved and outer_straight
        # extrudes, as well as all the pre-extruded surfaces
        tank_surfaces = [inner_square_surfaces[0]]
        for i in range(4):
            # bottom of tank
            tank_surfaces.append(outer_curved_surfaces[2 * i])
            tank_surfaces.append(outer_curved_surfaces[2 * i + 1])
            tank_surfaces.append(inner_curved_surfaces[2 * i])
            tank_surfaces.append(inner_curved_surfaces[2 * i + 1])
            tank_surfaces.append(outer_straight_surfaces[i])
            tank_surfaces.append(inner_straight_surfaces[i])
            tank_surfaces.append(square_surfaces[i])
            # lateral surfaces - curves: curves are last in the loops defining
            # outer_curved_surfaces, so should be the last elements
            tank_surfaces.append(outer_curved_extrude[2 * i][-1][1])
            tank_surfaces.append(outer_curved_extrude[2 * i + 1][-1][1])
            # lateral surfaces - baffle sides: outer_curved_extrude 2*i has
            # outside edge, and 2*i + 1 the inside edge
            tank_surfaces.append(outer_curved_extrude[2 * i][2][1])
            tank_surfaces.append(outer_curved_extrude[2 * i + 1][4][1])
            # lateral surfaces - baffle front: outer_straight extrude has
            # these at the 1st surface past the volume in the returned list
            tank_surfaces.append(outer_straight_extrude[i][2][1])

        gm.addPhysicalGroup(SURFACE, tank_surfaces)

        # id 2, weakly imposed noslip, this is the first element of each extrude
        top_surfaces = [inner_square_extrude[0][1]]
        for i in range(4):
            top_surfaces.append(outer_curved_extrude[2 * i][0][1])
            top_surfaces.append(outer_curved_extrude[2 * i + 1][0][1])
            top_surfaces.append(inner_curved_extrude[2 * i][0][1])
            top_surfaces.append(inner_curved_extrude[2 * i + 1][0][1])
            top_surfaces.append(outer_straight_extrude[i][0][1])
            top_surfaces.append(inner_straight_extrude[i][0][1])
            top_surfaces.append(square_extrude[i][0][1])

        gm.addPhysicalGroup(SURFACE, top_surfaces)

        # physical volume
        gm.addPhysicalGroup(VOLUME, [tag for _, tag in gmsh.model.getEntities(VOLUME)])

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


diameter = 17.4 / 100

# assume baffle width is same as rushton blade width,
# therefore width = Tl/Wl * Wl/Dl * Dl/Dt
tank = StructuredStirredTank(
    radius=diameter / 2,
    height=diameter,
    baffle_depth=diameter / 10,
    baffle_width=diameter * 0.333 * 0.2 * 0.1,
)
tank.draw()
tank.export()
