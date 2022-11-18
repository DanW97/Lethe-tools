#!/usr/bin/env python
# -*-coding:utf-8 -*-
#File    :   packing_generator.py
#Date    :   03/03/2022
#Author  :   Daniel Weston
#Contact :   dtw545@student.bham.ac.uk

"""
Workflow as follows:

1) use sobol technique to generate list of packings and Re
2) use liggghts (no gravity) to generate packed particles following grinding media psd
3) write .prms for these immersed boundary particles 
"""

import coexist
#from coexist import LiggghtsSimulation
import numpy as np 
from scipy.stats import qmc
import os
from sympy import nextprime

class PackingGenerator:
    """Generate packings of particles with given volume fraction and PSD"""
    def __init__(self, dim = 2, m = 4, Re_max = 100, Vx = 0.005, Vy = 0.005, Vz = 0.005, savename = "packingsims/packingrun", sim_filename = "packed_spheres", sim_folder = "lethefiles") -> None:
        self.dim = dim 
        self.m = m #exponent for the 2^m samples drawn
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
        self.savename = savename
        self.sim_filename = sim_filename
        self.sim_folder = sim_folder
    
    def psd(self, discrete: bool = True, mu: float = None, sigma: float = None, N: float = None):
        """Create the required PSD for LIGGGHTS to pack particles"""
        if discrete:
            self.radii = np.array([4.58e-4, 5.2e-4, 5.91e-4, 6.72e-4, 7.63e-4, 8.67e-4, 9.85e-4])
            self.size_fractions = np.array([0.0002, 0.0147, 0.1934, 0.6517, 0.1322, 0.0077, 0.0001])
            assert np.abs(np.sum(self.size_fractions)-1) < 1e-10
        else: #pull from normal distribution
            self.radii = np.random.normal(loc = mu, scale = sigma, size = (N,))
            self.size_fractions = np.array([1/n for n in range(N)])
            assert np.abs(np.sum(self.size_fractions)-1) < 1e-10
    
    def generate_run_parameters(self):
        """
        Create the list of phi, Re from sobol sampling to form the required .prm files. The created setpoints' first column is phi, the other is Re.        
        """
        sampler = qmc.Sobol(d=2, scramble=False)
        sample = sampler.random_base2(self.m)
        self.setpoints = qmc.scale(sample, [self.phi_min, self.Re_min], [self.phi_max, self.Re_max])

    def generate_packing_scripts(self, clean: bool = True):
        """This method generates all packing scripts"""
        if clean:
            os.remove("packingsims/*")
        nprimes = len(self.radii) #generate prime seeds for liggghts insertion
        primes = [nextprime(10000)] #liggghts requires seeds > 10000
        for i in range(nprimes): #generate nprimes more, so nprimes-1 goes to the particles and the last for insertion
            primes.append(nextprime(primes[i]))
        self.end_prime = primes[-1]
        self.packingpaths = []
        with open(self.packing_header_file, "r") as f:
            header = f.readlines()
            metadata = header[3:7]
            header = header[8:-1]
            for i, phi in enumerate(self.setpoints[:,0]):
                #domain setup
                domain = [f"region domain block {-self.Vx} {self.Vx} {-self.Vy} {self.Vy} {-self.Vz} {self.Vz} units box\n",
                          "create_box 1 domain\n",
                          f"neighbor {np.min(self.radii)} bin\n"]
                #particle info
                particle_templates = []
                app = particle_templates.append #cut some lookup time if we have some crazy normal distribution scenario
                for j, (radius, seed) in enumerate(zip(self.radii, primes)):
                    app(f"fix pts{j} all particletemplate/sphere {seed} atom_type 1 density constant 1 radius constant {radius}\n")
                liggghts_psd = [f"pts{j} {frac} " for j, frac in enumerate(self.size_fractions) ]
                psd_cmd = f"fix pdd1 all particledistribution/discrete {self.end_prime} {nprimes} " + ''.join(psd for psd in liggghts_psd) + "\n"
                particle_info = [particle for particle in particle_templates]
                particle_info.append(psd_cmd)
                #particle insertion
                # Inserting by volume fraction is further from target than numerical insertion
                insertion = [f"fix ins all insert/pack seed {nextprime(self.end_prime)} distributiontemplate pdd1 insert_every once overlapcheck yes all_in yes region domain volumefraction_region {phi} ntry_mc 10000000\n",
                "run 2\n"]
                filepath = self.savename + f"_{i}.sim"
                self.packingpaths.append(filepath)
                with open(filepath, 'w') as sim:
                    sim.writelines("#=================================================================================================\n")
                    sim.writelines(f"# LIGGGHTS script to pack particles with volume fraction target {phi}\n")
                    sim.writelines(metadata)
                    sim.writelines("#=================================================================================================\n\n")
                    sim.writelines(header)
                    sim.writelines(domain)
                    sim.writelines(particle_info)
                    sim.writelines(insertion)

    def generate_packings(self):
        """Run LIGGGHTS for each packing script"""
        self.real_phi = np.empty((np.shape(self.setpoints)[0]))
        for i, path in enumerate(self.packingpaths):
            sim = coexist.LiggghtsSimulation(path)
            pos = sim.positions()
            rad = sim.radii()
            real_phi = sum(4/3*np.pi*rad**3)/(self.Vx*self.Vy*self.Vz*8)
            self.real_phi[i] = real_phi
            self.write_particle_file(pos, rad)
            self.write_prm(real_phi, i)

    def write_particle_file(self, pos, rad):
        """Write Lethe particles file"""
        with open(self.particle_file, "w") as f:
            f.write("type shape_argument_0 p_x p_y p_z v_x v_y v_z omega_x omega_y omega_z\n")
            for i, (position, radius) in enumerate(zip(pos, rad)):
                f.write(f"sphere, {radius}, {position[0]}, {position[1]}, {position[2]}, 0, 0, 0, 0, 0, 0\n")
    
    def write_prm(self, real_phi, suffix):
        """Write Lethe .prm file for each simulation"""
        kinematic_viscosity = 1e-6 #water
        with open(self.sim_header_file, 'r') as f:
            header = f.readlines()
            header[15] = f"# Volume fraction = {real_phi}, Reynolds number = {self.setpoints[suffix,1]}\n"
            header[36] = f"    set output path                     = ./results_{suffix}/            # Output directory\n"
            # make box dimensions 10x the packed cube
            L = max([self.Vx, self.Vy, self.Vz])
            header[106] = f"    set grid arguments                  = {-2*self.Vx},{-2*self.Vy},{-2*self.Vz} : {2*self.Vx},{2*self.Vy},{2*self.Vz} : true\n"
            #calculate inlet velocity
            cross_section = self.Vy*self.Vz 
            reynolds_number = self.setpoints[suffix,1]
            hydraulic_diameter = 4*cross_section/(2*(self.Vy + self.Vz))
            inlet_velocity = kinematic_viscosity*reynolds_number/hydraulic_diameter
            header[177] = f"        set Function expression         = {inlet_velocity};0;0;0\n"
            header[191] = f"          set Function expression = {inlet_velocity}\n"

a = PackingGenerator(m=5)
a.psd()
a.generate_run_parameters()
a.generate_packing_scripts(clean=False)
a.generate_packings()