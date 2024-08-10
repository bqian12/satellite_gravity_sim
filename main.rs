use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// Struct for 3D vector and implementation of associated operations
#[derive(Debug, Clone, Copy)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Vector3 { x, y, z }
    }

    fn add(&self, other: &Vector3) -> Self {
        Vector3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    fn sub(&self, other: &Vector3) -> Self {
        Vector3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    fn scale(&self, scalar: f64) -> Self {
        Vector3 {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }

    fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    fn normalize(&self) -> Self {
        let mag = self.magnitude();
        self.scale(1.0 / mag)
    }
}

/// Struct for satellite simulation
#[derive(Debug)]
struct Satellite {
    mass: f64,
    position: Vector3,
    velocity: Vector3,
}

impl Satellite {
    fn new(mass: f64, position: Vector3, velocity: Vector3) -> Self {
        Satellite { mass, position, velocity }
    }

    fn gravitational_force(&self, earth_mass: f64, earth_position: Vector3) -> Vector3 {
        let gravitational_constant = 6.67430e-11; // m^3 kg^-1 s^-2
        let distance_vector = earth_position.sub(&self.position);
        let distance = distance_vector.magnitude();

        // Gravitational force: F = G * (m1 * m2) / r^2
        let force_magnitude = gravitational_constant * (self.mass * earth_mass) / (distance * distance);
        let force_direction = distance_vector.normalize();
        force_direction.scale(force_magnitude)
    }

    fn rk4_update(&mut self, time_step: f64, earth_mass: f64, earth_position: Vector3) {
        // k1 values
        let k1_velocity = self.gravitational_force(earth_mass, earth_position).scale(1.0 / self.mass);
        let k1_position = self.velocity;

        // k2 values
        let temp_velocity = self.velocity.add(&k1_velocity.scale(time_step / 2.0));
        let temp_position = self.position.add(&k1_position.scale(time_step / 2.0));
        let k2_velocity = self.gravitational_force(earth_mass, earth_position.sub(&temp_position)).scale(1.0 / self.mass);
        let k2_position = temp_velocity;

        // k3 values
        let temp_velocity = self.velocity.add(&k2_velocity.scale(time_step / 2.0));
        let temp_position = self.position.add(&k2_position.scale(time_step / 2.0));
        let k3_velocity = self.gravitational_force(earth_mass, earth_position.sub(&temp_position)).scale(1.0 / self.mass);
        let k3_position = temp_velocity;

        // k4 values
        let temp_velocity = self.velocity.add(&k3_velocity.scale(time_step));
        let temp_position = self.position.add(&k3_position.scale(time_step));
        let k4_velocity = self.gravitational_force(earth_mass, earth_position.sub(&temp_position)).scale(1.0 / self.mass);
        let k4_position = temp_velocity;

        // Final position and velocity update
        self.velocity = self.velocity.add(
            &k1_velocity
                .add(&k2_velocity.scale(2.0))
                .add(&k3_velocity.scale(2.0))
                .add(&k4_velocity)
                .scale(time_step / 6.0),
        );

        self.position = self.position.add(
            &k1_position
                .add(&k2_position.scale(2.0))
                .add(&k3_position.scale(2.0))
                .add(&k4_position)
                .scale(time_step / 6.0),
        );
    }
}

/// Struct for Satellite Simulation configuration and input
#[derive(Debug)]
struct SimulationConfig {
    mass: f64,
    initial_position: Vector3,
    initial_velocity: Vector3,
    time_step: f64,
    duration: f64,
    earth_mass: f64,
}

impl SimulationConfig {

    // Read config from file
    fn from_file(file_path: &str) -> Result<Self, io::Error> {
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);

        let mut mass = 0.0;
        let mut initial_position = Vector3::new(0.0, 0.0, 0.0);
        let mut initial_velocity = Vector3::new(0.0, 0.0, 0.0);
        let mut time_step = 0.0;
        let mut duration = 0.0;
        let mut earth_mass = 0.0;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let clean_line = line.split('#').next().unwrap().trim(); // Remove comments

            if clean_line.contains("satellite_mass:") {
                mass = clean_line.split(':').nth(1).unwrap().trim().parse().unwrap();
            } else if clean_line.contains("initial_conditions:") {
                let values: Vec<f64> = clean_line.split(':').nth(1).unwrap().trim()
                    .split_whitespace()
                    .map(|s| s.parse().unwrap())
                    .collect();
                println!("{}, {}, {}, {}, {}, {}", values[0], values[1], values[2], values[3], values[4], values[5]);
                initial_position = Vector3::new(values[0], values[1], values[2]);
                initial_velocity = Vector3::new(values[3], values[4], values[5]);
            } else if clean_line.contains("time_step:") {
                time_step = clean_line.split(':').nth(1).unwrap().trim().parse().unwrap();
            } else if clean_line.contains("duration:") {
                duration = clean_line.split(':').nth(1).unwrap().trim().parse().unwrap();
            } else if clean_line.contains("earth_mass:") {
                earth_mass = clean_line.split(':').nth(1).unwrap().trim().parse().unwrap();
            }
        }

        Ok(SimulationConfig {
            mass,
            initial_position,
            initial_velocity,
            time_step,
            duration,
            earth_mass,
        })
    }
}

fn main() {

    // Get the command-line arguments
    let args: Vec<String> = env::args().collect();

    // Determine the input file: either from command line or prompt the user
     let input_file = if args.len() > 1 {
        args[1].clone()
    } else {
        // Prompt user for input file
        print!("Please enter the input file name (default is 'input.txt'): ");
        io::stdout().flush().unwrap(); 

        let mut input = String::new();
        io::stdin().read_line(&mut input).unwrap();

        let trimmed_input = input.trim();
        if trimmed_input.is_empty() {
            "input.txt".to_string()
        } else {
            trimmed_input.to_string()
        }
    };

    let config = SimulationConfig::from_file(&input_file).expect("Failed to read input file");

    // Earth properties
    let earth_mass = config.earth_mass;
    let earth_position = Vector3::new(0.0, 0.0, 0.0); // Assuming Earth is at the origin

    // Satellite initial conditions
    let satellite_mass = config.mass;
    let initial_position = config.initial_position; 
    let initial_velocity = config.initial_velocity; 

    let mut satellite = Satellite::new(satellite_mass, initial_position, initial_velocity);

    // Simulation parameters
    let time_step = config.time_step;
    
    let total_time = config.duration;

    let total_steps = (total_time / time_step) as usize;
    for step in 0..(total_steps) {
        satellite.rk4_update(time_step, earth_mass, earth_position);
        if step % 10 == 0 { // Print every 10 timesteps
            println!(
                "Time: {}s | Position: ({:.2}, {:.2}, {:.2}) m | Velocity: ({:.2}, {:.2}, {:.2}) m/s",
                (step as f64*time_step),
                satellite.position.x,
                satellite.position.y,
                satellite.position.z,
                satellite.velocity.x,
                satellite.velocity.y,
                satellite.velocity.z,
            );
        }
    }

}
