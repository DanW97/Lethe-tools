# Lethe-tools

Collection of geometry creation python classes and occasionally .prm file writers for using programmatically generated (and hex-meshed) geometries in [lethe](https://github.com/lethe-cfd/lethe).

[Gmsh](https://gitlab.onelab.info/gmsh/gmsh/) is used exclusively to generate the geometries and meshes, however, the versions on conda-forge lack HXT support, so you are advised to install gmsh through `pip`.

```python
pip install gmsh
```

In each directory is a separate instance of geometry builders and .prm writers that form a base to modify as required. Currently in this repo:

- `extruded-cylinder`
  - code to generate an extruded cylinder, with finer control over the meshing than the deal.ii primitive cylinder.
- `packed-spheres`
  - code to pack spheres into a box to simulate flow over an obstacle with varying volume fraction and Reynolds number. Requires [coexist](https://github.com/uob-positron-imaging-centre/Coexist) and [a LIGGGHTS Python wrapper](https://github.com/uob-positron-imaging-centre/PICI-LIGGGHTS) for particle insertion.
- `stirred-tank`
  - code to generate an n-baffled tank and Rushton turbine.
- `superquadric-obstacles`
  - code to generate a superquadric object(s) to simulate flow over an obstacle.
