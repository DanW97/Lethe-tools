#!/usr/bin/env python
# -*-coding:utf-8 -*-
#File    :   packing_generator.py
#Date    :   03/03/2022
#Author  :   Daniel Weston
#Contact :   dtw545@student.bham.ac.uk

import gmsh
import sys


# Global constants - these form the "dim" section of the dimtag parlance used in gmsh
POINT = 0
CURVE = 1
SURFACE = 2
VOLUME = 3


class Cylinder:
    """Class that constructs a structured mesh of a cylinder, oriented along the z-axis

    Attributes
    ----------
    radius : int, optional
            Cylinder radius in mm, by default 49
    height : int, optional
            Cylinder height in mm, by default 2000
    radial_npts : int, optional
            Number of points on transfinite curves in x and y, by default 5
    height_npts : int, optional
            Number of points on transfinite curves in z, by default 20
    filepath : str, optional
            Path to save .geo and .msh files - extensions added automatically by this class, by default "cylinder"
    radial_coef : float , optional
            Geometric progression coefficient for transfinite curves defining circle segments in x and y, by default 1.
    height_coef : float , optional
            Geometric progression coefficient for transfinite lines in z, by default 1.
    diamond_coef : float , optional
            Geometric progression coefficient for transfinite curves defining center diamond in x and y, by default 1.
    diamond_ratio : float , optional
            Ratio of diamond half-length (half the diagonal length of a square) to cylinder radius
    """

    def __init__(self, radius=0.049, height=2, radial_npts=5, height_npts=20, filepath="cylinder", radial_coef=1.0, height_coef=1.0, diamond_coef=1.0, diamond_ratio=0.5) -> None:
        """Class constructor for Cylinder

        Parameters
        ----------
        radius : int, optional
                Cylinder radius in m, by default 0.049
        height : int, optional
                Cylinder height in m, by default 2
        radial_npts : int, optional
                Number of points on transfinite curves in x and y, by default 5
        height_npts : int, optional
                Number of points on transfinite curves in z, by default 20
        filepath : str, optional
                Path to save .geo and .msh files - extensions added automatically by this class, by default "cylinder"
        radial_coef : float , optional
                Geometric progression coefficient for transfinite curves defining circle segments in x and y, by default 1.
        height_coef : float , optional
                Geometric progression coefficient for transfinite lines in z, by default 1.
        diamond_coef : float , optional
                Geometric progression coefficient for transfinite curves defining center diamond in x and y, by default 1.
        diamond_ratio : float , optional
                Ratio of diamond half-length (half the diagonal length of a square) to cylinder radius
        """
        # Parameter check
        assert diamond_ratio < 1 and diamond_ratio > 0
        self.radius = radius/1000
        self.height = height/1000
        self.radial_npts = radial_npts
        self.height_npts = height_npts
        self.filepath = filepath
        self.radial_coef = radial_coef
        self.height_coef = height_coef
        self.diamond_coef = diamond_coef
        self.diamond_ratio = diamond_ratio

    def draw(self):
        """Draw a cylinder.

        This method draws 4 quadrants of the cylinder as separate volumes, and eventually defines the boundary IDs as follows:
        id     entity 
        1      walls
        2      -z face
        3      +z face
        """
        gmsh.initialize(sys.argv)
        gmsh.model.add("Cylinder")
        gm = gmsh.model.geo  # "with" doesn't work
        # Points defining segment boundaries
        center_point = gm.addPoint(0, 0, 0)
        arc_north_point = gm.addPoint(0, self.radius, 0)
        arc_east_point = gm.addPoint(self.radius, 0, 0)
        arc_south_point = gm.addPoint(0, -self.radius, 0)
        arc_west_point = gm.addPoint(-self.radius, 0, 0)
        # Points defining diamond edges
        diamond_north_point = gm.addPoint(0, self.radius*self.diamond_ratio, 0)
        diamond_east_point = gm.addPoint(self.radius*self.diamond_ratio, 0, 0)
        diamond_south_point = gm.addPoint(0, -self.radius*self.diamond_ratio, 0)
        diamond_west_point = gm.addPoint(-self.radius*self.diamond_ratio, 0, 0)
        # Lines from diamond edge to segment points
        trans_line_N = gm.addLine(diamond_north_point, arc_north_point)
        trans_line_E = gm.addLine(diamond_east_point, arc_east_point)
        trans_line_S = gm.addLine(diamond_south_point, arc_south_point)
        trans_line_W = gm.addLine(diamond_west_point, arc_west_point)
        # Draw the diamond
        diamond_line_NE = gm.addLine(diamond_north_point, diamond_east_point)
        diamond_line_ES = gm.addLine(diamond_east_point, diamond_south_point)
        diamond_line_SW = gm.addLine(diamond_south_point, diamond_west_point)
        diamond_line_WN = gm.addLine(diamond_west_point, diamond_north_point)
        diamond_loop    = gm.addCurveLoop([diamond_line_NE, diamond_line_ES, diamond_line_SW, diamond_line_WN])
        center_face_base= gm.addPlaneSurface([diamond_loop])
        extruded_center = gm.extrude([(SURFACE, center_face_base)], dx = 0, dy = 0, dz = self.height, numElements=[self.height_npts], recombine=True)
        center_face_top = extruded_center[0][1]
        center_volume   = extruded_center[1][1]
        # Draw segment 1: North - East
        arcNE = gm.addCircleArc(arc_north_point, center_point, arc_east_point)
        segmentNE = gm.addCurveLoop([trans_line_N, arcNE, -trans_line_E, -diamond_line_NE])
        faceNE_base = gm.addPlaneSurface([segmentNE])
        extrudedNE = gm.extrude([(SURFACE, faceNE_base)], dx=0, dy=0, dz=self.height, numElements=[self.height_npts], recombine=True)
        faceNE_top = extrudedNE[0][1]
        volumeNE = extrudedNE[1][1]
        curved_faceNE = extrudedNE[3][1]
        # Draw segment 2: East - South
        arcES_base = gm.addCircleArc(
            arc_east_point, center_point, arc_south_point)
        segment_ES_base = gm.addCurveLoop([trans_line_E, arcES_base, -trans_line_S, -diamond_line_ES])
        face_ES_base = gm.addPlaneSurface([segment_ES_base])
        extrudedES = gm.extrude([(SURFACE, face_ES_base)], dx=0, dy=0, dz=self.height, numElements=[self.height_npts], recombine=True)
        face_ES_top = extrudedES[0][1]
        volumeES = extrudedES[1][1]
        curved_faceES = extrudedES[3][1]
        # Draw segment 3: South - West
        arcSW_base = gm.addCircleArc(arc_south_point, center_point, arc_west_point)
        segment_SW_base = gm.addCurveLoop([trans_line_S, arcSW_base, -trans_line_W, -diamond_line_SW])
        face_SW_base = gm.addPlaneSurface([segment_SW_base])
        extrudedSW = gm.extrude([(SURFACE, face_SW_base)], dx=0, dy=0, dz=self.height, numElements=[self.height_npts], recombine=True)
        face_SW_top = extrudedSW[0][1]
        volumeSW = extrudedSW[1][1]
        curved_faceSW = extrudedSW[3][1]
        # Draw segment 4: West - North
        arcWN_base = gm.addCircleArc(arc_west_point, center_point, arc_north_point)
        segment_WN_base = gm.addCurveLoop([trans_line_W, arcWN_base, -trans_line_N, -diamond_line_WN])
        face_WN_base = gm.addPlaneSurface([segment_WN_base])
        extrudedWN = gm.extrude([(SURFACE, face_WN_base)], dx=0, dy=0, dz=self.height, numElements=[self.height_npts], recombine=True)
        face_WN_top = extrudedWN[0][1]
        volumeWN = extrudedWN[1][1]
        curved_faceWN = extrudedWN[3][1]
        gm.synchronize()
        # Define transfinite behaviour
        msh = gm.mesh
        # Shared lines
        msh.setTransfiniteCurve(
            trans_line_N, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteCurve(
            trans_line_E, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteCurve(
            trans_line_S, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteCurve(
            trans_line_W, self.radial_npts, coef=self.radial_coef)
        # Segment 1: North - East
        msh.setTransfiniteCurve(
            arcNE, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteSurface(faceNE_base)
        msh.setRecombine(SURFACE, faceNE_base)
        # Segment 2: East - South
        msh.setTransfiniteCurve(
            arcES_base, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteSurface(face_ES_base)
        msh.setRecombine(SURFACE, face_ES_base)
        # Segment 3: South - West
        msh.setTransfiniteCurve(
            arcSW_base, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteSurface(face_SW_base)
        msh.setRecombine(SURFACE, face_SW_base)
        # Segment 4: West - North
        msh.setTransfiniteCurve(
            arcWN_base, self.radial_npts, coef=self.radial_coef)
        msh.setTransfiniteSurface(face_WN_base)
        msh.setRecombine(SURFACE, face_WN_base)
        # Center diamond
        msh.setTransfiniteCurve(diamond_line_NE, self.radial_npts, coef=self.diamond_coef)
        msh.setTransfiniteCurve(diamond_line_ES, self.radial_npts, coef=self.diamond_coef)
        msh.setTransfiniteCurve(diamond_line_SW, self.radial_npts, coef=self.diamond_coef)
        msh.setTransfiniteCurve(diamond_line_WN, self.radial_npts, coef=self.diamond_coef)
        msh.setTransfiniteSurface(center_face_base)
        msh.setRecombine(SURFACE, center_face_base)
        # Define physical groups for lethe
        gr1 = gm.addPhysicalGroup(
            SURFACE, [curved_faceNE, curved_faceES, curved_faceSW, curved_faceWN])
        gr2 = gm.addPhysicalGroup(
            SURFACE, [faceNE_base, face_ES_base, face_SW_base, face_WN_base, center_face_base])
        gr3 = gm.addPhysicalGroup(
            SURFACE, [faceNE_top, face_ES_top, face_SW_top, face_WN_top, center_face_top])
        gr4 = gm.addPhysicalGroup(
            VOLUME, [volumeNE, volumeES, volumeSW, volumeWN, center_volume])
        gm.synchronize()

        self.export()


    def export(self):
        """unclean way to add all the hard-coded stuff in, and rename the file"""
        gmsh.model.geo.synchronize()
        raw_file = f"{self.filepath}.geo_unrolled"
        new_file = f"{self.filepath}.geo"
        gmsh.write(raw_file)
        gmsh.model.mesh.generate(VOLUME)
        gmsh.model.geo.synchronize()
        types = [i + 1 for i in range(31)]
        types.append(92)
        types.append(93)
        # Should only see types 1, 3, 5 and 15
        for type in types:
            els = gmsh.model.mesh.getElementsByType(type)
            print(type, els)
        gmsh.write(f"{self.filepath}.msh")
        import os
        os.system(f"mv {raw_file} {new_file}")

    def view(self):
        """Visualise the geo file"""
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()


cyl = Cylinder(height_npts=30, radial_npts=7)
cyl.draw()
cyl.view()
