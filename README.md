# Satellite Simulation with EGM2008

## Overview

This project simulates the orbit of a satellite and performs numerical integration using the 4th order Runge-Kutta method. The simulation reads input from a configuration file, which specifies satellite mass, initial conditions, time parameters, and earth parameters. Previously, the EGM2008 model was used but issues were encountered uploading the coefficients to Github.

## Environment

This simulation requires Rust to be installed on the system.

## Running the Simulation

To run the simulation, simply `cargo build` and `cargo run` from the top level `satellite_gravity_sim` directory. An example input file is included called `input.txt`. The user will be prompted for an input file, but the simulation will also by default use the provided `input.txt` file.
