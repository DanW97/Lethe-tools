#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File    :   stirred_tank.py
# Time    :   09/09/2022
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


class StirredTankBuilder:
    """
        Class that generates hex meshes of stirred tanks with N number of baffles.

        Parameters
        ----------
        radius : float, optional
            Tank radius, by default 0.1
        height : float, optional
            Tank height, by default 0.2
        nbaffles : int, optional
            Number of baffles, by default 4
        pts_mode : str, optional
            Mode for number of vertical points calculation, choices are "number" which is the user-specified number of nodes or "physical" which takes a physical spacing and calculates the number of nodes to achieve this. For the latter, the spacing *must* be an integer multiple of the tank height, by default "number"
        height_spacing : int, optional
            Number of nodes if pts_mode == "number", or spacing if pts_mode == "physical", by default 20
        filepath : str, optional
            Name of msh and geo files created by this class, by default "stirred_tank"
        axis_alignment : int, optional
            Axis alignment of the impeller, either 0, 1 or 2 for x, y or z respectively, by default 1
        view : bool, optional
            Display the resultant mesh in gmsh, by default True
        baffle_width : float, optional
            Baffle thickness, by default 0.00569
        baffle_depth : float, optional
            Baffle depth, or protrusion from the wall, by default 0.01996

        Raises
        ------
        ValueError
            Exception raised if pts_mode == "physical" and the height_spacing value is *not* an integer multiple of the height value.
        """
    def __init__(
        self,
        radius=0.1,
        height=0.2,
        nbaffles=4,
        pts_mode="number",
        height_spacing=20,
        filepath="stirred_tank",
        axis_alignment=1,
        view=True,
        baffle_width=0.00569,
        baffle_depth=0.01996,
    ) -> None:
        if pts_mode == "number":
            self.height_spacing = height_spacing
        elif pts_mode == "physical":
            if height % height_spacing == 0:
                self.height_spacing = int(height / height_spacing)
            else:
                raise ValueError(
                    "The height must be wholly divisable by the height spacing!"
                )
        self.filepath = filepath
        self.radius = radius
        self.height = height
        self.nbaffles = nbaffles
        self.axis_alignment = axis_alignment
        self.baffle_angles = np.linspace(0, 2 * np.pi, self.nbaffles, endpoint=False)
        self.baffle_origins = np.zeros((3, self.nbaffles))
        self.baffle_coordinates = np.zeros((3, self.nbaffles, 4))
        self.baffle_width = baffle_width
        self.baffle_depth = baffle_depth
        self.view = view

    def draw(self):
        """
        Draw the tank

        Raises
        ------
        ValueError
            Raised if incorrect axis_alignment parameter given.
        """
        gmsh.initialize(sys.argv)
        gmsh.model.add("Stirred Tank")
        gm = gmsh.model.geo
        # draw baffles
        # we start at 0 rad and draw baffles anti-clockwise
        # list of lists so we can maintain structure of baffle points to aid drawing
        baffle_pts = []
        curves = []
        # assign a point to the tank centre
        center = gm.addPoint(0.0, 0.0, 0.0)
        # set origin coordinates for baffles
        self.baffle_origins[0, :] += self.radius * np.cos(self.baffle_angles)
        self.baffle_origins[1, :] += self.radius * np.sin(self.baffle_angles)
        # now assign 4 pts corresponding to the 4 corners of the baffle
        # from 0 rad:
        # 1 ----------- 0
        #   |
        #   |
        # 2 ----------- 3
        for i in range(self.nbaffles):
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
        curves = []
        for i in range(self.nbaffles):
            # for each baffle, join 3-2 2-1 1-0 and then curve from 0 to 3 on next baffle
            lines = [
                gm.addLine(baffle_pts[i][3 - j], baffle_pts[i][2 - j]) for j in range(3)
            ]
            if i < self.nbaffles - 1:
                arc_03 = gm.addCircleArc(baffle_pts[i][0], center, baffle_pts[i + 1][3])
            else:
                arc_03 = gm.addCircleArc(baffle_pts[i][0], center, baffle_pts[0][3])
            curves.extend([*lines, arc_03])
        curve_loop = gm.addCurveLoop(curves)
        surface = gm.addPlaneSurface([curve_loop])
        gm.synchronize()
        match self.axis_alignment:
            case 0:
                gm.rotate(
                    [(SURFACE, surface)],
                    x=0,
                    y=0,
                    z=0,
                    ax=0,
                    ay=1,
                    az=0,
                    angle=-np.pi / 2,
                )
                gm.extrude(
                    [(SURFACE, surface)],
                    dx=self.height,
                    dy=0,
                    dz=0,
                    numElements=[self.height_spacing],
                    recombine=True,
                )
            case 1:
                gm.rotate(
                    [(SURFACE, surface)],
                    x=0,
                    y=0,
                    z=0,
                    ax=1,
                    ay=0,
                    az=0,
                    angle=-np.pi / 2,
                )
                gm.extrude(
                    [(SURFACE, surface)],
                    dx=0,
                    dy=self.height,
                    dz=0,
                    numElements=[self.height_spacing],
                    recombine=True,
                )
            case 2:  # only extrude
                gm.extrude(
                    [(SURFACE, surface)],
                    dx=0,
                    dy=0,
                    dz=self.height,
                    numElements=[self.height_spacing],
                    recombine=True,
                )
            case _:
                raise ValueError("Axis alignment value is *only* 0, 1 or 2!")
        gm.synchronize()
        gm.addPhysicalGroup(
            SURFACE, [tag for _, tag in gmsh.model.getEntities(SURFACE)], name="1"
        )
        gm.addPhysicalGroup(
            VOLUME, [tag for _, tag in gmsh.model.getEntities(VOLUME)], name="1"
        )
        gm.synchronize()

    def export(self):
        """Export .msh and .geo_unrolled (converted to .geo manually) files. 
           """
        gm = gmsh.model.geo
        msh = gmsh.model.mesh
        gm.synchronize()
        gmsh.option.setNumber("Mesh.Smoothing", 100)
        gmsh.option.setNumber("Mesh.SaveAll", 0)  # only save physical groups
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        gmsh.option.setNumber("General.NumThreads", 8)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
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


tank = StirredTankBuilder()
tank.draw()
tank.export()
