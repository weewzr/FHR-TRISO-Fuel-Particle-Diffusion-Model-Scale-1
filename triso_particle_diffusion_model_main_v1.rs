use uom::si::f64::*;
use uom::si::diffusion_coefficient::square_meter_per_second;
use uom::si::molar_concentration::mole_per_cubic_meter;
use uom::si::length::meter;
use uom::si::time::second;
use csv::Writer;
use plotters::prelude::*;

#[derive(Clone, Copy, PartialEq)]
#[allow(dead_code)]
enum SurfaceBoundary {
    Insulated, //zero-flux surface, the paper's source-free baseline
    Dirichlet, //fixed surface concentration, used when coupling to the coating layers
}

const SURFACE_BC: SurfaceBoundary = SurfaceBoundary::Insulated; //choose the surface boundary condition here

fn main() {
    central_difference_method();
}

pub fn central_difference_method() {
    let diffusivity = DiffusionCoefficient::new::<square_meter_per_second>(1.0e-9); //placeholder value
    let delta_radius = Length::new::<meter>(0.001); //placeholder value
    let timestep = Time::new::<second>(50.0); //placeholder value

    let number_of_nodes: usize = 5; //layers in the TRISO fuel particle, node indices 0 to 4
    let number_of_steps: usize = 200; //number of timesteps, raise for the full long-time run
    let last = number_of_nodes - 1; //index of the surface node (4)

    let mut concentration = vec![MolarConcentration::new::<mole_per_cubic_meter>(0.0); number_of_nodes]; //BC symmetry = 0, interior starts empty
    concentration[last] = MolarConcentration::new::<mole_per_cubic_meter>(1.0); //BC surface = 1

    let mut current_time = Time::new::<second>(0.0); //start the clock at zero

    let mut wtr = Writer::from_path("ficks_second_law.csv").unwrap(); //open the csv file for writing
    wtr.write_record(&["time_s", "node_0", "node_1", "node_2", "node_3", "node_4"]).unwrap(); //write the header row

    let mut history: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(number_of_steps); number_of_nodes]; //per-node time history for the concentration vs time plot
    let mut snapshots: Vec<(f64, Vec<f64>)> = Vec::new(); //radial-profile snapshots for the profiles vs radius plot
    let snapshot_every = (number_of_steps / 8).max(1); //store a profile roughly eight times across the run
    let read_profile = |c: &Vec<MolarConcentration>| -> Vec<f64> { c.iter().map(|x| x.get::<mole_per_cubic_meter>()).collect() }; //pull the raw numbers out of one concentration vector
    snapshots.push((0.0, read_profile(&concentration))); //store the initial profile at t = 0

    for step in 0..number_of_steps {
        let mut concentration_next_timestep = concentration.clone(); //concentration at the next timestep uses duplicates of the current timestep - Explicit Euler Scheme

        let lambda = (diffusivity * timestep / (delta_radius * delta_radius)).value; // lambda = D * dt / (dr)^2

        concentration_next_timestep[0] = concentration[0] + 6.0 * lambda * (concentration[1] - concentration[0]); //concentration at node 0, the centre of the sphere, at the next timestep (factor 6 from the L'Hopital limit at r = 0)

        for i in 1..last { //concentration at the interior nodes 1 to 3, between the centre and the surface, at the next timestep
            let i_f = i as f64; //node index as a float, node i sits at r = i*dr

            let term_1 = 1.0 - 2.0 * lambda; //term 1 in the equation, the coefficient of C_i
            let term_2 = lambda * (1.0 + 1.0 / i_f); //term 2 in the equation, the coefficient of C_(i+1), curvature weight 1/i
            let term_3 = lambda * (1.0 - 1.0 / i_f); //term 3 in the equation, the coefficient of C_(i-1), curvature weight 1/i

            concentration_next_timestep[i] = concentration[i] * term_1 + concentration[i + 1] * term_2 + concentration[i - 1] * term_3; //concentration at node i, between the centre and the surface, at the next timestep
        }

        match SURFACE_BC {
            SurfaceBoundary::Insulated => {
                concentration_next_timestep[last] = concentration[last] + 2.0 * lambda * (concentration[last - 1] - concentration[last]); //concentration at node 4, the surface of the sphere, at the next timestep (insulated, ghost node C_5 = C_3, factor 2 since r = R is not singular)
            }
            SurfaceBoundary::Dirichlet => {
                //fixed surface concentration, leave node 4 untouched as the clone already holds the value
            }
        }

        concentration = concentration_next_timestep; //update the concentration vector to the next timestep
        current_time = current_time + timestep; //update the current time to the next timestep

        let t = current_time.get::<second>(); //current time as a number
        let mut row = vec![t.to_string()]; //create a new row for the csv file with the current time
        for i in 0..number_of_nodes { //loop through the concentration vector and add each node to the row
            let c = concentration[i].get::<mole_per_cubic_meter>(); //concentration at node i as a number
            row.push(c.to_string()); //add the concentration at node i to the row
            history[i].push((t, c)); //store the value for the concentration vs time plot
        }
        wtr.write_record(&row).unwrap(); //write the row to the csv file

        if (step + 1) % snapshot_every == 0 || step + 1 == number_of_steps { //every so often, and on the final step
            snapshots.push((t, read_profile(&concentration))); //store a radial-profile snapshot
        }
    }
    wtr.flush().unwrap(); //make sure everything is written to disk

    plot_history(&history, number_of_nodes).expect("failed to write time plot"); //draw concentration vs time
    plot_profiles(&snapshots, delta_radius.get::<meter>()).expect("failed to write profile plot"); //draw radial profiles
    println!("Done. Wrote ficks_second_law.csv, ficks_second_law.png and ficks_radial_profiles.png"); //tell the user it finished
}

fn plot_history(history: &[Vec<(f64, f64)>], number_of_nodes: usize) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("ficks_second_law.png", (960, 640)).into_drawing_area(); //the image to draw on
    root.fill(&WHITE)?; //white background

    let t_max = history[0].last().map(|p| p.0).unwrap_or(1.0); //last time value, sets the x range
    let y_max = history.iter().flatten().map(|p| p.1).fold(0.0_f64, f64::max).max(1.0); //largest concentration, sets the y range

    let mut chart = ChartBuilder::on(&root)
        .caption("Nodal concentration vs time", ("sans-serif", 26)) //title
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(0.0..t_max, 0.0..(y_max * 1.05))?; //axes ranges

    chart.configure_mesh().x_desc("Time (s)").y_desc("Concentration (mol/m^3)").draw()?; //axis labels and grid

    for i in 0..number_of_nodes { //one line per node
        let colour = Palette99::pick(i).mix(1.0); //pick a colour for this node
        chart.draw_series(LineSeries::new(history[i].iter().cloned(), &colour))? //draw the line
            .label(format!("Node {}", i)) //legend entry
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], colour)); //legend swatch
    }

    chart.configure_series_labels().background_style(WHITE.mix(0.85)).border_style(BLACK).draw()?; //draw the legend box
    root.present()?; //save the file
    Ok(())
}

fn plot_profiles(snapshots: &[(f64, Vec<f64>)], delta_radius_m: f64) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("ficks_radial_profiles.png", (960, 640)).into_drawing_area(); //the image to draw on
    root.fill(&WHITE)?; //white background

    let number_of_nodes = snapshots[0].1.len(); //how many nodes per profile
    let r_max = (number_of_nodes - 1) as f64 * delta_radius_m; //outer radius, sets the x range
    let y_max = snapshots.iter().flat_map(|(_, p)| p.iter().cloned()).fold(0.0_f64, f64::max).max(1.0); //largest concentration, sets the y range

    let mut chart = ChartBuilder::on(&root)
        .caption("Radial concentration profiles at successive times", ("sans-serif", 26)) //title
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(0.0..r_max, 0.0..(y_max * 1.05))?; //axes ranges

    chart.configure_mesh().x_desc("Radius r (m)   [node i at r = i*dr]").y_desc("Concentration (mol/m^3)").draw()?; //axis labels and grid

    for (k, (time, profile)) in snapshots.iter().enumerate() { //one line per saved time
        let colour = Palette99::pick(k).mix(1.0); //pick a colour for this time
        let points = profile.iter().enumerate().map(|(i, &c)| (i as f64 * delta_radius_m, c)); //node index becomes radius i*dr
        chart.draw_series(LineSeries::new(points, &colour))? //draw the profile
            .label(format!("t = {:.0} s", time)) //legend entry
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], colour)); //legend swatch
    }

    chart.configure_series_labels().background_style(WHITE.mix(0.85)).border_style(BLACK).draw()?; //draw the legend box
    root.present()?; //save the file
    Ok(())
}
