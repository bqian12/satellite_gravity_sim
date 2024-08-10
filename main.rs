use ndarray::{Array1, Array};
use std::fs::File;
use std::io::{BufRead, BufReader, Result};

/// Constants
const G: f64 = 6.67430e-11;         // Gravitational constant (m^3 kg^-1 s^-2)
const EARTH_MASS: f64 = 5.972e24;   // Mass of Earth (kg)
const EARTH_RADIUS: f64 = 6371e3;   // Radius of Earth (m)
const TIME_STEP: f64 = 1.0;         // Simulation time step (s)

/// Struct for a satellite
struct Satellite {
    position: Array1<f64>, // Position vector [x, y, z] in meters
    velocity: Array1<f64>, // Velocity vector [vx, vy, vz] in meters per second
    mass: f64,             // Mass of the satellite in kilograms
}

struct EGM2008Coefficient {
    n: usize,
    m: usize,
    cnm: f64,
    snm: f64,
    sigma_c: f64,
    sigma_s: f64,
}

#[derive(Debug)]
struct EGM2008Header {
    product_type: String,
    modelname: String,
    earth_gravity_constant: f64,
    radius: f64,
    max_degree: usize,
    errors: String,
    norm: String,
    tide_system: String,
    url: String,
}

struct EGM2008 {
    header: EGM2008Header,
    coefficients: Vec<EGM2008Coefficient>,
}

impl EGM2008 {
    /// Create a new EGM2008 instance by loading from the file
    pub fn new(file_path: &str) -> Result<Self> {
        let file = File::open(file_path)?;
        let mut reader = BufReader::new(file);

        let header = Self::parse_header(&mut reader)?;
        let coefficients = Self::parse_coefficients(&mut reader)?;

        Ok(EGM2008 { header, coefficients })
    }

    /// Parse the header section from the file
    fn parse_header(reader: &mut dyn BufRead) -> Result<EGM2008Header> {
        let mut header = EGM2008Header {
            product_type: String::new(),
            modelname: String::new(),
            earth_gravity_constant: 0.0,
            radius: 0.0,
            max_degree: 0,
            errors: String::new(),
            norm: String::new(),
            tide_system: String::new(),
            url: String::new(),
        };

        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            let parts: Vec<&str> = line.trim().split_whitespace().collect();
            if parts.is_empty() {
                break;
            }

            match parts[0] {
                "product_type" => header.product_type = parts[1].to_string(),
                "modelname" => header.modelname = parts[1].to_string(),
                "earth_gravity_constant" => header.earth_gravity_constant = parts[1].parse().unwrap(),
                "radius" => header.radius = parts[1].parse().unwrap(),
                "max_degree" => header.max_degree = parts[1].parse().unwrap(),
                "errors" => header.errors = parts[1].to_string(),
                "norm" => header.norm = parts[1].to_string(),
                "tide_system" => header.tide_system = parts[1].to_string(),
                "url" => header.url = parts[1].to_string(),
                _ => {}
            }
            line.clear();
        }

        Ok(header)
    }

    /// Parse the coefficients section from the file
    fn parse_coefficients(reader: &mut dyn BufRead) -> Result<Vec<EGM2008Coefficient>> {
        let mut coefficients = Vec::new();

        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            if line.starts_with("end_of_head") {
                line.clear();
                break;
            }
            line.clear();
        }

        while reader.read_line(&mut line)? != 0 {
            let parts: Vec<&str> = line.trim().split_whitespace().collect();
            if parts.len() < 6 {
                continue;
            }
            let n: usize = parts[1].parse().unwrap();
            let m: usize = parts[2].parse().unwrap();
            let cnm: f64 = parts[3].parse().unwrap();
            let snm: f64 = parts[4].parse().unwrap();
            let sigma_c: f64 = parts[5].parse().unwrap();
            let sigma_s: f64 = parts[6].parse().unwrap();

            coefficients.push(EGM2008Coefficient { n, m, cnm, snm, sigma_c, sigma_s });
            line.clear();
        }

        Ok(coefficients)
    }

    /// Compute the gravitational potential using EGM2008 coefficients
    pub fn gravitational_potential(&self, latitude: f64, longitude: f64, altitude: f64) -> f64 {
        let radius = self.header.radius; // Mean Earth radius from header
        let mut potential = 0.0;

        for coeff in &self.coefficients {
            let (n, m) = (coeff.n, coeff.m);
            if n == 0 {
                continue;
            }

            let lat_rad = latitude.to_radians();
            let lon_rad = longitude.to_radians();
            let alt_radius = radius + altitude;
            let sin_lat = lat_rad.sin();
            let cos_lat = lat_rad.cos();
            let sin_lon = lon_rad.sin();
            let cos_lon = lon_rad.cos();

            let r_p = alt_radius.powi(-(n as i32));
            let legendre = self.legendre_polynomial(n, m, sin_lat);
            let potential_term = (coeff.cnm * cos_lon + coeff.snm * sin_lon) * legendre * r_p;

            potential += potential_term;
        }

        potential
    }

    /// Compute the associated Legendre polynomial
    fn legendre_polynomial(&self, n: usize, m: usize, x: f64) -> f64 {
        let mut pmm = 1.0;
        if m > 0 {
            let mut fact = 1.0;
            for i in 0..m {
                pmm *= -fact * (2 * i + 1) as f64 / (2 * i + 2) as f64;
                fact += 1.0;
            }
        }
        if n == m {
            return pmm;
        }

        let mut pmmp1 = x * (2.0 * m as f64 + 1.0) * pmm;
        if n == m + 1 {
            return pmmp1;
        }

        let mut pnm = 0.0;
        for i in m + 2..=n {
            pnm = (x * (2.0 * i as f64 - 1.0) * pmmp1 - (i + m - 1) as f64 * pmm) / (i - m) as f64;
            pmm = pmmp1;
            pmmp1 = pnm;
        }

        pnm
    }
}

/// Implementation of a satellite object that propagates dynamics forward using RK4 (Runge-Kutta) integration
impl Satellite {
    /// Calculate the derivatives for the Runge-Kutta method
    fn derivatives(&self, force: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
        let acceleration = force / self.mass;
        let velocity = self.velocity.clone();
        (velocity, acceleration)
    }

    /// Update the position and velocity of the satellite using the 4th-order Runge-Kutta method
    fn update(&mut self, force_fn: impl Fn(&Satellite) -> Array1<f64>, dt: f64) {
        // Calculate k1
        let force1 = force_fn(self);
        let (k1_vel, k1_acc) = self.derivatives(&force1);

        // Calculate k2
        let temp_sat = Satellite {
            position: &self.position + &(&k1_vel * (dt / 2.0)),
            velocity: &self.velocity + &(&k1_acc * (dt / 2.0)),
            ..*self
        };
        let force2 = force_fn(&temp_sat);
        let (k2_vel, k2_acc) = temp_sat.derivatives(&force2);

        // Calculate k3
        let temp_sat = Satellite {
            position: &self.position + &(&k2_vel * (dt / 2.0)),
            velocity: &self.velocity + &(&k2_acc * (dt / 2.0)),
            ..*self
        };
        let force3 = force_fn(&temp_sat);
        let (k3_vel, k3_acc) = temp_sat.derivatives(&force3);

        // Calculate k4
        let temp_sat = Satellite {
            position: &self.position + &(&k3_vel * dt),
            velocity: &self.velocity + &(&k3_acc * dt),
            ..*self
        };
        let force4 = force_fn(&temp_sat);
        let (k4_vel, k4_acc) = temp_sat.derivatives(&force4);

        // Combine to get the final update
        self.position = &self.position
            + dt / 6.0
                * (&k1_vel + &(2.0 * &k2_vel) + &(2.0 * &k3_vel) + &k4_vel);
        self.velocity = &self.velocity
            + dt / 6.0
                * (&k1_acc + &(2.0 * &k2_acc) + &(2.0 * &k3_acc) + &k4_acc);
    }
}

/// Compute the norm of an Array1, helper function for satellite to compute r
fn array1_norm(array: &Array1<f64>) -> f64 {
    array.mapv(|x| x.powi(2)).sum().sqrt()
}

/// Calculate the gravitational force using the EGM2008 model
fn gravitational_force(egm2008: &EGM2008, satellite: &Satellite) -> Array1<f64> {
    let lat = satellite.position[1].atan2((satellite.position[0].powi(2) + satellite.position[2].powi(2)).sqrt());
    let lon = satellite.position[2].atan2(satellite.position[0]);
    let r = array1_norm(&satellite.position);

    let potential = egm2008.gravitational_potential(lat, lon, r);

    // The force is the negative gradient of the potential
    let force_magnitude = -1.0 * potential / r;
    -&satellite.position * force_magnitude / r
}

// Helper structs for parsing input file
#[derive(Debug)]
struct InitialConditions {
    position: Array1<f64>,
    velocity: Array1<f64>,
}

#[derive(Debug)]
struct TimeParameters {
    time_step: f64,
    duration: f64,
}

#[derive(Debug)]
struct GravitationalModel {
    file_path: String,
}

#[derive(Debug)]
struct SimulationInput {
    mass: f64,
    initial_conditions: InitialConditions,
    time_parameters: TimeParameters,
    gravitational_model: GravitationalModel,
}

impl SimulationInput {
    pub fn new(file_path: &str) -> Result<Self> {
        let file = File::open(file_path)?;
        let mut reader = BufReader::new(file);

        let mass = Self::parse_mass(&mut reader)?;
        let initial_conditions = Self::parse_initial_conditions(&mut reader)?;
        let time_parameters = Self::parse_time_parameters(&mut reader)?;
        let gravitational_model = Self::parse_gravitational_model(&mut reader)?;

        Ok(SimulationInput {
            mass,
            initial_conditions,
            time_parameters,
            gravitational_model,
        })
    }

    fn parse_mass(reader: &mut dyn BufRead) -> Result<f64> {
        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            if line.trim().starts_with("mass:") {
                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                let mass = parts[1].replace(',', ".").parse().unwrap();
                return Ok(mass);
            }
            line.clear();
        }
        Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Mass not found"))
    }

    fn parse_initial_conditions(reader: &mut dyn BufRead) -> Result<InitialConditions> {
        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            if line.trim().starts_with("initial_conditions:") {
                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                let position: Array1<f64> = Array::from_vec(parts[1..4].iter().map(|&x| x.parse().unwrap()).collect());
                let velocity: Array1<f64> = Array::from_vec(parts[4..7].iter().map(|&x| x.parse().unwrap()).collect());
                return Ok(InitialConditions { position, velocity });
            }
            line.clear();
        }
        Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Initial conditions not found"))
    }

    fn parse_time_parameters(reader: &mut dyn BufRead) -> Result<TimeParameters> {
        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            if line.trim().starts_with("time_step:") {
                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                let time_step = parts[1].parse().unwrap();
                line.clear();
                reader.read_line(&mut line)?;
                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                let duration = parts[1].parse().unwrap();
                return Ok(TimeParameters { time_step, duration });
            }
            line.clear();
        }
        Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Time parameters not found"))
    }

    fn parse_gravitational_model(reader: &mut dyn BufRead) -> Result<GravitationalModel> {
        let mut line = String::new();
        let mut model = GravitationalModel {
            file_path: String::new(),
        };
        while reader.read_line(&mut line)? != 0 {
            if line.trim().starts_with("file path:") {
                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                model.file_path = parts[2].to_string();
                return Ok(model);
            }
            line.clear();
        }
        Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Gravitational model file path not found"))
    }
}

fn main() -> Result<()>{
    let input_file_path = "input.txt";
    println!("Reading in simulation input file.");
    let input = SimulationInput::new(input_file_path)?;
    
    println!("Reading in EGM2008 gravity model at path: {}. May take a while.",&input.gravitational_model.file_path);
    let egm2008 = EGM2008::new(&input.gravitational_model.file_path)?;

    let mut position = input.initial_conditions.position;
    let mut velocity = input.initial_conditions.velocity;
    let time_step = input.time_parameters.time_step;
    let duration = input.time_parameters.duration;
    let mass = input.mass;

    let mut satellite = Satellite {
        position,
        velocity,
        mass,
    };

    let num_steps = duration / time_step;
    
    println!(
        "Initial Conditions: Time: {:.2} s, Position: [{:.2}, {:.2}, {:.2}] m, Velocity: [{:.2}, {:.2}, {:.3}] m/s",
        0.0,
        satellite.position[0],
        satellite.position[1],
        satellite.position[2],
        satellite.velocity[0],
        satellite.velocity[1],
        satellite.velocity[2]
    );
 
    for step in 0..(num_steps as i64) {
        satellite.update(|sat| gravitational_force(&egm2008, sat), time_step);

        if step % 10 == 0 {
            println!(
                "Time: {:.2} s, Position: [{:.2}, {:.2}, {:.2}] m, Velocity: [{:.2}, {:.2}, {:.3}] m/s",
                step as f64 * time_step,
                satellite.position[0],
                satellite.position[1],
                satellite.position[2],
                satellite.velocity[0],
                satellite.velocity[1],
                satellite.velocity[2]
            );
        }
    }
    Ok(())
}
