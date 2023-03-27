#![allow(unused)]
#![allow(non_snake_case)]
use ndarray::*;
use plotly::*;
use ndarray_linalg::*;
// use plotters::prelude::*;


pub const GREETING: &'static str = "hello, Rust library here!";


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
let ul:Array::<f64, Ix1> = &meshx.slice(s![..,1]) - &meshx.slice(s![..,0]);
let ur:Array::<f64, Ix1> = &meshx.slice(s![..,2]) - &meshx.slice(s![..,0]);
let ll:Array::<f64, Ix1> = &meshy.slice(s![..,1]) - &meshy.slice(s![..,0]);
let lr:Array::<f64, Ix1> = &meshy.slice(s![..,2]) - &meshy.slice(s![..,0]);
let uu:Array::<f64, Ix1> = &ptx.slice(s![..,1]) - &meshx.slice(s![..,0]);
let dd:Array::<f64, Ix1> = &pty.slice(s![..,1]) - &meshy.slice(s![..,0]);

let mut n = 0;
while n < 8 {
    if ptx[(n, 1)] != 0. {
        //inv matrix here
        let A = arr2(&[[ul[n], ur[n]], 
                         [ll[n], lr[n]]]);
        let B = arr1(&[uu[n], dd[n]]);
        
        let x = A.solve_into(B).unwrap();

        println!("{}", x);
    }
    else {
        println!("no particle in {}th element", n);
    }

    n+=1
}


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
