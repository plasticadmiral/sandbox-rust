#![allow(unused)]
#![allow(non_snake_case)]
use ndarray::{*, Array3};
use ndarray_linalg::*;
// use plotters::prelude::*;


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



// poisson equation code
pub fn meshgrid(x: &Array<f64,Ix1>, y: &Array<f64,Ix1>) -> 
                    (Array<f64,Ix2>,Array<f64,Ix2>) {
    
    let X = x.broadcast((y.dim(), x.dim()))
                .unwrap().to_owned();
    let Y = y.broadcast((x.dim(), y.dim()))
                .unwrap().reversed_axes().to_owned();
    (X, Y)
}


pub fn poisson_data(nxx: &usize, nyy: &usize) -> 
        (Array<f64,Ix2>,Array<f64,Ix2>,Array<f64,Ix2>) {
    
    // lattice space values
    // let (nx, ny) = (10, 10);
    let (nx, ny) = (*nxx, *nyy);
    let n_max = (nx-1) * (ny-1);
    let (dx, dy) = (1./((nx as f64)-1.), 1./((ny as f64)-1.));
    assert_eq!(nx, ny);

    // loading tridiagonal values
    let mut block = Array::<f64,_>::eye(nx - 1)
         * -(2./f64::powi(dx,2) + 2./f64::powi(dy,2));
    block.slice_mut(s![(1 as usize).., ..]).diag_mut().fill(1./f64::powi(dy,2));
    block.slice_mut(s![.., (1 as usize)..]).diag_mut().fill(1./f64::powi(dy,2));
    
    // loading offdiagonal values
    let mut  A = linalg::kron(&Array::<f64,_>::eye(nx-1), &block);
    A.slice_mut(s![(ny-1).., ..]).diag_mut().fill(1./f64::powi(dy,2));
    A.slice_mut(s![.., (ny-1)..]).diag_mut().fill(1./f64::powi(dy,2));
    
    // calculating z = A^-1 b with b = const
    let b = Array::<f64,_>::ones(n_max) * -3.;
    let op = A.solve(&b).unwrap();

    // prepping to plot data
    let x = Array::linspace(1., (nx as f64) - 1., nx - 1) * dx;
    let y = Array::linspace(1., (ny as f64) - 1., ny - 1) * dy;
    let (X, Y) = meshgrid(&x, &y);
    let Z = op.into_shape((nx-1, ny-1)).unwrap();

    (X,Y,Z)
}
