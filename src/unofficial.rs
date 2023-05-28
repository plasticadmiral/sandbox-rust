
#![allow(unused)]
#![allow(non_snake_case)]
use ndarray::{*, Array3};
use ndarray_linalg::*;
// use plotters::prelude::*;




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




pub const START: &'static str = "running... ";
pub const END: &'static str = "complete! ";

// min max step returns node positions and connectivity matrix
pub fn generate_mesh() {
    let mut all = array![[1.,2.,3.], [4.,5.,6.]]; 
    let mut alo = Array::range(0., 10., 0.5);
    println!("{}", alo);
}



pub fn is_point_in_element() {
    

}

pub fn make_particles(ylen: usize, xlen: usize, gaps: f32) -> Array3<f32> { 

    let yvars:Array1<f32> = Array::range(0., ylen as f32, gaps);
    let xvars:Array1<f32> = Array::range(0., xlen as f32, gaps);

    let mut particle:Array3<f32> = Array::zeros((2, yvars.len(), xvars.len()));

    let mut y_layer = particle.index_axis_mut(Axis(0), 0);
    for i in 0..y_layer.shape()[1] {
        let mut col = y_layer.index_axis_mut(Axis(1), i);
        col += &yvars;
    }
    drop(y_layer);

    let mut x_layer = particle.index_axis_mut(Axis(0), 1);
    for i in 0..x_layer.shape()[0] {
        let mut row = x_layer.index_axis_mut(Axis(0), i);
        row += &xvars;
    }
    particle
}

pub fn make_mesh(factor: usize) {

    // Global coordinates for particles
    // let particle_positions = arr2(&[[0.25,0.75],
                                    // [1.75,0.25],
                                    // [1.75,1.75]]);

    // Holds mesh indices
    // let node_positions = arr2(&[[0,0],[1,0],[2,0],
    //                             [0,1],[1,1],[2,1],
    //                             [0,2],[1,2],[2,2]]);
    
    let yvars:Array1<f32> = Array::range(0., (3 + 2 * factor) as f32, 1.);
    let xvars:Array1<f32> = Array::range(0., (3 + 2 * factor) as f32, 1.);

    let mut node_positions:Array3<f32> = Array::zeros((2, yvars.len(), xvars.len()));

    let mut y_layer = node_positions.index_axis_mut(Axis(0), 0);
    for i in 0..y_layer.shape()[1] {
        let mut col = y_layer.index_axis_mut(Axis(1), i);
        col += &yvars;
    }
    drop(y_layer);

    let mut x_layer = node_positions.index_axis_mut(Axis(0), 1);
    for i in 0..x_layer.shape()[0] {
        let mut row = x_layer.index_axis_mut(Axis(0), i);
        row += &xvars;
    }

    // Holds mesh_positions indices
    // let node_connectivity = arr2(&[[0,1,4],[1,2,4],[0,3,4],
    //                                [0,3,4],[2,5,4],[3,4,6],
    //                                [4,5,8],[4,7,6],[4,8,7]]);
        
    // let mut node_connectivity = arr3(&[[[0,1,0],
    //                                     [0,2,3],
    //                                     [4,4,4]],

    //                                    [[1,2,3],
    //                                     [3,5,4],
    //                                     [5,7,8]],

    //                                    [[4,4,4],
    //                                     [4,4,6],
    //                                     [8,6,7]]]);

    let mut node_connectivity:Array3<i32>;
    if factor == 0 {
        node_connectivity = Array::zeros((3, 4, 2));
    } else {
        node_connectivity = Array::zeros((3, 4 * 2 * factor, 2 * 2 * factor));
    }

    // for mut layer in node_connectivity.axis_iter_mut(Axis(0)) {
    //     println!("{}", layer(node_position));
    // }
    // println!("{}", node_positions);
    
    let mut pt1_layer = node_connectivity.index_axis_mut(Axis(0), 0);
    for y in 0..pt1_layer.shape()[0] {
        for x in 0..pt1_layer.shape()[1] {
            if x % 2 == 0 && y % 2 == 0 {
                println!("pass")
                //pt1_layer[y][x] = node_positions[y][x]
            }
        }
        // let mut col = y_layer.index_axis_mut(Axis(1), i);
        // col += &yvars;
    }
    // drop(y_layer);

    let mut pt2_layer = node_connectivity.index_axis_mut(Axis(0), 1);
    // for i in 0..x_layer.shape()[0] {
    //     let mut row = x_layer.index_axis_mut(Axis(0), i);
    //     row += &xvars;
    // }

    let mut pt3_layer = node_connectivity.index_axis_mut(Axis(0), 2);
    // for i in 0..x_layer.shape()[0] {
    //     let mut row = x_layer.index_axis_mut(Axis(0), i);
    //     row += &xvars;
    // }

}
// MPM code
pub fn run_mpm() {

let meshx = arr2(&[[0.,1.,0.],
                   [0.,1.,1.],
                   [1.,2.,1.],
                   [2.,2.,1.],
                   [1.,2.,2.],
                   [1.,2.,1.],
                   [1.,1.,0.],
                   [0.,1.,0.]]);

let meshy = arr2(&[[0.,1.,1.],
                   [0.,0.,1.],
                   [0.,0.,1.],
                   [0.,1.,1.], 
                   [1.,1.,2.],
                   [1.,2.,2.],
                   [1.,2.,2.],
                   [1.,1.,2.]]);

let ptx = arr2(&[[0.,0.25],
                 [0.,0.], 
                 [0.,1.75], 
                 [0.,0.], 
                 [0.,0.], 
                 [0.,0.], 
                 [0.,0.75],
                 [0.,0.]]);

let pty = arr2(&[[0.,0.75],
                 [0.,0.], 
                 [0.,0.25], 
                 [0.,0.], 
                 [0.,0.], 
                 [0.,0.], 
                 [0.,1.75],
                 [0.,0.]]);

let aa = arr1(&[1, 2, 43]);
// let elems = Array::<f64, Ix2>::zeros((8, 2).f());
let ul:Array::<f64, Ix1> = &meshx.slice(s![..,1]) - &meshx.slice(s![..,0]); // LHS
let ur:Array::<f64, Ix1> = &meshx.slice(s![..,2]) - &meshx.slice(s![..,0]); // LHS
let ll:Array::<f64, Ix1> = &meshy.slice(s![..,1]) - &meshy.slice(s![..,0]); // LHS
let lr:Array::<f64, Ix1> = &meshy.slice(s![..,2]) - &meshy.slice(s![..,0]); // LHS
let uu:Array::<f64, Ix1> = &ptx.slice(s![..,1]) - &meshx.slice(s![..,0]);   // RHS 
let dd:Array::<f64, Ix1> = &pty.slice(s![..,1]) - &meshy.slice(s![..,0]);   // RHS

println!("{}", meshx.slice(s![..,2]));

// let mut n = 0;
// while n < 8 {
//     if ptx[(n, 1)] != 0. {
//         //inv matrix here
//         let A = arr2(&[[ul[n], ur[n]], 
//                          [ll[n], lr[n]]]);
//         let B = arr1(&[uu[n], dd[n]]);
        
//         let x = A.solve_into(B).unwrap();

//         println!("{}", x);
//     }
//     else {
//         println!("no particle in {}th element", n);
//     }

//     n+=1
// }


//println!("arrray is: {:}",meshx.slice(s![..,1]));
// let elex = concatenate![Axis(0), &[ele1.view(), ele2.view()]];

println!("complete");
         
}