#![allow(unused)]
#[allow(non_snake_case)]
// #[allow(dead_code)]
// #[allow(unused_mut)]
// #[allow(unused_variables)]

// use std::{iter::FromIterator, result};
use sandbox_rust::*;
// use plotly::{ImageFormat, Plot, Scatter};
use ndarray::*;
mod tests;


fn mainofficial() {
    let nodes = Mesh {
        pos: arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]),
        conn: arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]),
    };
    let nodes_vel: Array2<f64> = arr2(&[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [2., 1.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]);
    // let pts: Array2<f64> = arr2(&[[0.25, 1.75], [1.75, 0.25], [0.75, 1.75]]);
    // let pt: Array1<f64> = arr1(&[0.5, 0.75]);
    let pt1: Array2<f64> = arr2(&[[0.25, 0.75], [1.50, 0.75]]);
    let pt2: Array2<f64> = arr2(&[[0.75, 0.25], [1.25, 1.50]]);
    let pt3: Array2<f64> = arr2(&[[1.75, 1.25], [1.50, 1.25]]);

    // returns result [[0.50, 0.25], [1.0, 0.5]]
    let var1 = get_particle_velocity(&nodes, &nodes_vel, &pt1); 
    let var2 = get_particle_velocity(&nodes, &nodes_vel, &pt2); 
    let var3 = get_particle_velocity(&nodes, &nodes_vel, &pt3); 
    
    
    // let var1 = g_fn(&node_pos, &node_conn, &pt); 
    // let (idx, xpos, ypos) = get_particle_position(&node_pos, &node_conn, &pt); 
    println!("{:?}", var1);
    println!("{:?}", var2);
    println!("{:?}", var3);
    
}

fn mainofficial1() {
    let f = get_mesh(3, 2);

    // [[0.0, 0.0],
    //  [1.0, 0.0],
    //  [2.0, 0.0],
    //  [0.0, 1.0],
    //  [1.0, 1.0],
    //  [2.0, 1.0]]
    println!("{:?}", f.pos);

    // [[0, 4, 3],
    //  [0, 1, 4],
    //  [1, 5, 4],
    //  [1, 2, 5]]
    println!("{:?}", f.conn);
}
fn main() {
    mainofficial1();
}
