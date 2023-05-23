#![allow(unused)]
use std::{iter::FromIterator, result};
use sandbox_rust::*;
use plotly::{ImageFormat, Plot, Scatter};
use ndarray::*;
use ndarray_linalg::*;

// #[allow(dead_code)]
// #[allow(unused_mut)]
// #[allow(unused_variables)]
#[allow(non_snake_case)]


pub fn setup() {
    println!("running setup ...");
    let particles = make_particles(3, 3, 1.);
    let particles = make_mesh(1);
    // println!("{:?}", particles);
}

pub fn iterate() {
    println!("running iterate ...");
}

pub fn create_plot() {
    let mut plot = Plot::new();
    let trace = Scatter::new(vec![0, 1, 2], vec![2, 1, 0]);
    plot.add_trace(trace);

    //plot.write_html("../data/out.html");
    //plot.show();
    plot.write_image("../data/out.png", ImageFormat::PNG, 800, 600, 1.0);
}

fn mainunofficial() {
    println!("main runs ...");

    const GRIDSIZE: i32 = 2;
    setup(); 
    let mut t = 98;
    
    while t < 100 {

        iterate();

        t += 1; 
    
    }


    
    //create_plot()

}

fn mainmessy() {

    println!("{}", START);

    let mut t = 1;

    generate_mesh();

    // revamp();

    while t == 1 {
        println!("helo");
        

        t = t + 1;
    }

    println!("{}", END);
    
    // run_mpm();

    // let (x, y) = (13, 13);
    // let (mut X, mut Y, mut Z) = poisson_data(&x, &y);
    // println!("x--------------x");
    // println!("{X}");
    // println!("x--------------x");
    // println!("{Y}");
    // println!("x--------------x");
    // println!("{Z}");
}
/*/
pub fn l_step(old_local_particle_pos: Array2<f64>, node_vel: Array2<f64>) {

    let mut result: Array2<f64> = Array::zeros((0, 2));
    let vel = ArrayView::from(&node_vel);

    for particle in old_local_particle_pos.rows() {
        if particle.iter().all(|&x| x >= 0.0) {  // if lower than zero then invalid valid
            let new_local_particle_pos = particle * node_vel; 
            println!("{}", new_local_particle_pos);
            // result.push_row(ArrayView::from(&new_local_particle_pos)).unwrap();
        }
    }


    // let new_local_particle_pos = local_particle_pos * node_vel; 
    
    println!("{}", result);
    // result.push_row(ArrayView::from(&new_local_particle_pos)).unwrap();
            
    // println!("{}", result)
}
*/ 
// velocity from mesh gets transferred onto particle
pub fn f(node_pos: Array2<f64>, node_conn: Array2<usize>, node_vel: Array2<f64>, global_particle_pos: Array2<f64>) -> Array2<f64> {

    let mut result: Array2<f64> = Array::zeros((0, 2));

    // each idx contains row index values of node_pos that makes up a traingle mesh
    for idx in node_conn.rows() { //each triangle

        let pt_neg: Array1<f64> = arr1(&[node_pos.row(idx[0])[0], node_pos.row(idx[0])[1]]);

        let lhs_ul = node_pos.row(idx[1])[0] - pt_neg[0];
        let lhs_ur = node_pos.row(idx[2])[0] - pt_neg[0];
        let lhs_ll = node_pos.row(idx[1])[1] - pt_neg[1];
        let lhs_lr = node_pos.row(idx[2])[1] - pt_neg[1];

        for pt in global_particle_pos.rows() { // each point in each triangle
    

            let lhs = arr2(&[[lhs_ul, lhs_ur], 
                            [lhs_ll, lhs_lr]]);

            let rhs = arr1(&[pt[0] - pt_neg[0], 
                            pt[1] - pt_neg[1]]);

            let local_particle_pos = lhs.solve_into(rhs).unwrap(); // from global positions we get local positions
            result.push_row(ArrayView::from(&local_particle_pos)).unwrap();
            
            // applying velocity from mesh nodes to particle
            if local_particle_pos.iter().all(|&x| x >= 0.0) {  // if lower than zero then invalid valid
                println!("before giving velocity: {:?}", local_particle_pos);

                let new_local_pos = local_particle_pos * node_vel.view(); 

                println!("after giving velocity: {:?}", new_local_pos);
            }
        }

    }
        result
}

fn calc_b(node_pos: Array2<f64>, node_conn: Array2<usize>) -> Array2<f64> {

    let mut B: Array2<f64> = Array::zeros((0, 6));
    let gradient_local_N0 = arr1(&[-1, -1]);
    let gradient_local_N1 = arr1(&[1, 0]);
    let gradient_local_N2 = arr1(&[0, 1]);

    for N in node_conn.rows() {
        let J =  arr2(&[[node_pos.row(N[1])[0] - node_pos.row(N[0])[0], node_pos.row(N[1])[1] - node_pos.row(N[0])[1]], 
                        [node_pos.row(N[2])[0] - node_pos.row(N[0])[0], node_pos.row(N[3])[1] - node_pos.row(N[0])[1]]]);

        let gradient_global_N0 = J.solve_into(gradient_local_N0).unwrap();
        let gradient_global_N1 = J.solve_into(gradient_local_N1).unwrap();
        let gradient_global_N2 = J.solve_into(gradient_local_N2).unwrap();
        let intermediate = arr2(&[[gradient_global_N0[0], 0, gradient_global_N1[0], 0, gradient_global_N2[0], 0], 
                                  [0, gradient_global_N0[1], 0, gradient_global_N1[1], 0, gradient_global_N2[1]], 
                                  [gradient_global_N0[1], gradient_global_N0[0], gradient_global_N1[1], 
                                   gradient_global_N1[0], gradient_global_N2[1], gradient_global_N2[0]]]); 

        B.push_row(ArrayView::from(&intermediate)).unwrap();
        
    }

    B
}

fn mainofficial() {
    let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
    let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);
    let node_vel: Array2<f64> = arr2(&[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [2., 1.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]);
    // let pts: Array2<f64> = arr2(&[[0.25, 1.75], [1.75, 0.25], [0.75, 1.75]]);
    let global_particle_pos: Array2<f64> = arr2(&[[0.5, 0.75], [0.25, 0.5]]);

    // global particles to local particles without velocity (print statement used to display local particle pos after velocity)
    let local_particle_pos = f(node_pos, node_conn, node_vel, global_particle_pos);
    
    
}





fn main() {
    mainofficial()
}
