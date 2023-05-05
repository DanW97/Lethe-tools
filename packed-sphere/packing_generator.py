#!/usr/bin/env python
# -*-coding:utf-8 -*-
# File    :   packing_generator.py
# Date    :   03/03/2022
# Author  :   Daniel Weston
# Contact :   dtw545@student.bham.ac.uk

from __future__ import annotations

from pathlib import Path
import coexist

# from coexist import LiggghtsSimulation
import numpy as np
from scipy.stats import qmc
from sympy import nextprime
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


class PackingGenerator:
    """
    Class to generate packings of particles with given volume fraction and PSD.
    This class also builds Lethe .prm files that simulate a range of Reynolds numbers
    to calculate the effect that Reynolds number and packing fraction have on the
    drag coefficient.

        Parameters
        ----------
        m : int, optional
            Exponent for Sobol sampling, 2^m samples will be drawn, by default 4
        Re_max : int, optional
            Maximum Reynolds number, by default 100
        Vx : float, optional
            Length of insertion region in x direction in m, by default 0.005
        Vy : float, optional
            Length of insertion region in y direction in m, by default 0.005
        Vz : float, optional
            Length of insertion region in z direction in m, by default 0.005
        liggghts_filename : str, optional
            Filename pattern for each set of particle positions, by default "packingrun"
        liggghts_folder : str, optional
            Filename pattern for each set of particle positions,
            by default "packingsims"
        sim_filename : str, optional
            Filename pattern for each simulation, by default "packed_spheres"
        sim_folder : str, optional
            Destination folder for .prm files, by default "lethefiles"
    """

    def __init__(
        self,
        m=4,
        Re_max=100,
        Vx=0.005,
        Vy=0.005,
        Vz=0.005,
        liggghts_filename="packingrun",
        liggghts_folder="packingsims",
        sim_filename="packed_spheres",
        sim_folder="lethefiles",
    ) -> None:
        self.m = m
        self.phi_max = 0.3
        self.phi_min = 0.01
        self.Re_max = Re_max
        self.Re_min = 0.1
        self.packing_header_file = "packing_header.sim"
        self.sim_header_file = "packed_spheres_header.prm"
        self.particle_file = "particles"
        self.Vx = Vx
        self.Vy = Vy
        self.Vz = Vz
        self.liggghts_filename = liggghts_filename
        self.liggghts_folder = liggghts_folder
        self.sim_filename = sim_filename
        self.sim_folder = sim_folder

    def psd(
        self,
        discrete: bool = True,
        mu: float = 0.001,
        sigma: float = 0.0005,
        N: int = 30,
    ):
        """Create the required PSD for LIGGGHTS to pack particles.

        This method sets the *radius* for each particle.

        Parameters
        ----------
        discrete : bool, optional
            True if a discrete psd is required, False if a Gaussian psd is
            required, by default True
        mu : float, optional
            Mean of the psd, only required if discrete == False, by default 0.001
        sigma : float, optional
            Standard deviation of the psd, only required if discrete == False,
            by default 0.0005
        N : int, optional
            Number of sizes to draw from the psd, only required
            if discrete == False, by default 30
        """
        if discrete:
            self.radii = np.array(
                [4.58e-4, 5.2e-4, 5.91e-4, 6.72e-4, 7.63e-4, 8.67e-4, 9.85e-4]
            )
            self.size_fractions = np.array(
                [0.0002, 0.0147, 0.1934, 0.6517, 0.1322, 0.0077, 0.0001]
            )
            assert np.abs(np.sum(self.size_fractions) - 1) < 1e-10
        else:  # pull from normal distribution
            self.radii = np.random.normal(loc=mu, scale=sigma, size=(N,))
            self.size_fractions = np.array([1 / n for n in range(N)])
            assert np.abs(np.sum(self.size_fractions) - 1) < 1e-10

    def generate_run_parameters(self):
        """
        Create the list of phi, Re from sobol sampling to form the required .prm files.

        The created setpoints' first column is phi, the other is Re.
        """
        sampler = qmc.Sobol(d=2, scramble=False)
        sample = sampler.random_base2(self.m)
        self.setpoints = qmc.scale(
            sample, [self.phi_min, self.Re_min], [self.phi_max, self.Re_max]
        )

    def generate_packing_scripts(self, clean: bool = True):
        """Create LIGGGHTS packing scripts

        Parameters
        ----------
        clean : bool, optional
            Remove previously generated files if True, by default True
        """
        if clean is True:
            self._clean(self.liggghts_folder, self.liggghts_filename)
        nprimes = len(self.radii)  # generate prime seeds for liggghts insertion
        primes = [nextprime(10000)]  # liggghts requires seeds > 10000
        for i in range(
            nprimes
        ):  # generate nprimes more, so nprimes-1 goes to the particles and the last
            # for insertion
            primes.append(nextprime(primes[i]))
        self.end_prime = primes[-1]
        self.packingpaths = []
        futures = []
        with ThreadPoolExecutor() as exe:
            for i, phi in enumerate(self.setpoints[:, 0]):
                filepath = Path(
                    self.liggghts_folder, self.liggghts_filename + f"_{i}.sim"
                )
                self.packingpaths.append(filepath)
                future = exe.submit(
                    self._write_liggghts_sim, i, phi, nprimes, primes  # type: ignore
                )
                futures.append(future)
            for f in futures:
                f.result()

    def _write_liggghts_sim(self, i: int, phi: float, nprimes: int, primes: list[int]):
        """Private method called by concurrent.futures for each packing value.

        Parameters
        ----------
        i : int
            Index within setpoints array.
        phi : float
            Target packing fraction.
        nprimes : int
            Total number of primes.
        primes : int
            List of primes for LIGGGHTS seeds.
        """
        with open(self.packing_header_file, "r") as f:
            header = f.readlines()
            metadata = header[3:7]
            header = header[8:-1]
            # domain setup
            domain = [
                (
                    "region domain block"
                    f" {-self.Vx} {self.Vx} {-self.Vy} {self.Vy} {-self.Vz} {self.Vz} "
                    "units box\n"
                ),
                "create_box 1 domain\n",
                f"neighbor {np.min(self.radii)} bin\n",
            ]
            # particle info
            particle_templates = []
            app = (
                particle_templates.append
            )  # cut some lookup time if we have some crazy normal distribution scenario
            for j, (radius, seed) in enumerate(zip(self.radii, primes)):
                app(
                    f"fix pts{j} all particletemplate/sphere {seed} atom_type 1 density"
                    f" constant 1 radius constant {radius}\n"
                )
            liggghts_psd = [
                f"pts{j} {frac} " for j, frac in enumerate(self.size_fractions)
            ]
            psd_cmd = (
                "fix pdd1 all particledistribution/discrete"
                f" {self.end_prime} {nprimes} "
                + "".join(psd for psd in liggghts_psd)
                + "\n"
            )
            particle_info = [particle for particle in particle_templates]
            particle_info.append(psd_cmd)
            # particle insertion
            # Inserting by volume fraction is further from target than
            # numerical insertion
            insertion = [
                (
                    "fix ins all insert/pack seed"
                    f" {nextprime(self.end_prime)} distributiontemplate pdd1"
                    " insert_every once overlapcheck yes all_in yes region domain"
                    f" volumefraction_region {phi} ntry_mc 10000000\n"
                ),
                "run 2\n",
            ]
            with open(self.packingpaths[i], "w") as sim:
                sim.writelines(
                    "#================================================================="
                    "================================\n"
                )
                sim.writelines(
                    "# LIGGGHTS script to pack particles with volume fraction target"
                    f" {phi}\n"
                )
                sim.writelines(metadata)
                sim.writelines(
                    "#================================================================="
                    "================================\n\n"
                )
                sim.writelines(header)
                sim.writelines(domain)
                sim.writelines(particle_info)
                sim.writelines(insertion)

    def generate_packings(self):
        """Run LIGGGHTS for each packing script"""
        self.real_phi = np.empty((np.shape(self.setpoints)[0]))
        futures = []
        with ProcessPoolExecutor() as exe:
            for i, path in enumerate(self.packingpaths):
                future = exe.submit(self._liggghts_parallel, path, i)
                futures.append(future)
            for f in futures:
                f.result()

    def _liggghts_parallel(self, path: str, i: int):
        sim = coexist.LiggghtsSimulation(path)
        pos = sim.positions()
        rad = sim.radii()
        real_phi = sum(4 / 3 * np.pi * rad**3) / (self.Vx * self.Vy * self.Vz * 8)
        self.real_phi[i] = real_phi
        self._write_particle_file(pos, rad, i)
        self._write_prm(real_phi, i)

    def _write_particle_file(self, pos: np.ndarray, rad: np.ndarray, suffix: int):
        """Write Lethe particle file.

        Parameters
        ----------
        pos : np.ndarray
            3D particle position.
        rad : np.ndarray
            Particle radius.
        suffix : int
            Numbered identifier for corresponding sim.
        """
        with open(self.particle_file, "w") as f:
            f.write(
                "type shape_argument_0 p_x p_y p_z v_x v_y v_z omega_x omega_y"
                " omega_z\n"
            )
            for i, (position, radius) in enumerate(zip(pos, rad)):
                f.write(
                    f"sphere, {radius}, {position[0]}, {position[1]}, {position[2]}, 0,"
                    " 0, 0, 0, 0, 0\n"
                )

    def _write_prm(self, real_phi: float, suffix: int, clean: bool = True):
        """Write Lethe .prm files.

        Parameters
        ----------
        real_phi : float
            Actual volume fraction, this may be different to self.phi[suffix] due to
            limitations in the packing algorithm.
        suffix : int
            File suffix.
        clean : bool, optional
            If True, remove previously generated files, by default True
        """
        if suffix == 0 and clean is True:
            self._clean(self.sim_folder, self.sim_filename)
        kinematic_viscosity = 1e-6  # water
        with open(self.sim_header_file, "r") as f:
            header = f.readlines()
            header[15] = (
                f"# Volume fraction = {real_phi}, Reynolds number ="
                f" {self.setpoints[suffix,1]}\n"
            )
            header[36] = (
                f"    set output path                     = ./results_{suffix}/        "
                "    # Output directory\n"
            )
            # make box dimensions 10x the packed cube
            header[106] = (
                "    set grid arguments                  ="
                f" {-2*self.Vx},{-2*self.Vy},{-2*self.Vz} :"
                f" {2*self.Vx},{2*self.Vy},{2*self.Vz} : true\n"
            )
            # calculate inlet velocity
            cross_section = self.Vy * self.Vz
            reynolds_number = self.setpoints[suffix, 1]
            hydraulic_diameter = 4 * cross_section / (2 * (self.Vy + self.Vz))
            inlet_velocity = kinematic_viscosity * reynolds_number / hydraulic_diameter
            header[
                177
            ] = f"        set Function expression         = {inlet_velocity};0;0;0\n"
            header[191] = f"          set Function expression = {inlet_velocity}\n"
            with open(
                Path(self.sim_folder, self.sim_filename + f"_{suffix}.prm"), "w"
            ) as sim:
                sim.writelines(header)

    def _clean(self, target_dir: str, file_pattern: str):
        """Remove files in target_dir that start with file_pattern.

        Parameters
        ----------
        target_dir : str
            Folder with files to remove.
        file_pattern : str
            Pattern required for glob expression.
        """
        path = Path(Path(__file__).parent, target_dir)
        for file in path.glob(file_pattern + "*"):
            file.unlink()


a = PackingGenerator(m=5)
a.psd()
a.generate_run_parameters()
a.generate_packing_scripts()
a.generate_packings()
