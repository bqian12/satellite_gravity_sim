# Satellite Simulation with EGM2008

## Overview

This project simulates the orbit of a satellite using the Earth Gravitational Model 2008 (EGM2008) and performs numerical integration using the 4th order Runge-Kutta method. The simulation reads input from a configuration file, which specifies mass, initial conditions, time parameters, and the path to the EGM2008 coefficients file. The EGM2008 gravitational potential model coefficients were obtained from ICGEM.

## Environment

This simulation requires Rust to be installed on the system, along with the `ndarray` crate.

## Running the Simulation

To run the simulation, simply `cargo build` and `cargo run`. An example input file is included called `input.txt`. Note that the simulation is currently slow due to the amount of computation. Improvements should be made in the future to performance.
