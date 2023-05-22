#![allow(unused)]
use std::iter::FromIterator;
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

fn mainofficial() {
let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
let node_idx: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);
let node_vel: Array2<f64> = arr2(&[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [2., 1.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]);
// let pts: Array2<f64> = arr2(&[[0.25, 1.75], [1.75, 0.25], [0.75, 1.75]]);
let pts: Array2<f64> = arr2(&[[0.5, 0.75], [0.25, 0.5]]);


// each idx contains row index values of node_pos that makes up a traingle mesh
for idx in node_idx.rows() {

    let pt_neg: Array1<f64> = arr1(&[node_pos.row(idx[0])[0], node_pos.row(idx[0])[1]]);

    let lhs_ul = node_pos.row(idx[1])[0] - pt_neg[0];
    let lhs_ur = node_pos.row(idx[2])[0] - pt_neg[0];
    let lhs_ll = node_pos.row(idx[1])[1] - pt_neg[1];
    let lhs_lr = node_pos.row(idx[2])[1] - pt_neg[1];

    for pt in pts.rows() {
    

        let lhs = arr2(&[[lhs_ul, lhs_ur], 
                         [lhs_ll, lhs_lr]]);

        let rhs = arr1(&[pt[0] - pt_neg[0], 
                         pt[1] - pt_neg[1]]);

        let mut x = lhs.solve_into(rhs).unwrap();

        if x.iter().all(|&x| x >= 0.0) {
            println!("before: {:?}", x);
            x = x * node_vel.row(4); 
            println!("after {:?}", x);
        }
    }

}

// println!("{:?}", pts_flattened);


}

fn main() {
    mainofficial()
}
