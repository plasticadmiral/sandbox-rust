#![allow(unused)]
#![allow(non_snake_case)]

use std::{iter::FromIterator, result};
// use sandbox_rust::*;
use plotly::{ImageFormat, Plot, Scatter};
use ndarray::*;
use ndarray_linalg::*;
use ndarray::{*, Array3};
use ndarray_linalg::*;
// use plotters::prelude::*;

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


pub fn f_fn(node_pos: Array2<f64>, node_conn: Array2<usize>, node_vel: Array2<f64>, pts: Array2<f64>) -> Array2<f64> {

    let mut pts_vel:Array2<f64> = Array::zeros(node_vel.dim()); // allocating memory

    for pt in pts.rows() {

        let output = g_fn(&node_pos, &node_conn, &pt.to_owned()); // returns [node_conn index, local xpos, local ypos]
        let nval = node_conn.row(output[0] as usize); // getting N position index for triangle containing the point

        // V_p = N(X_p)V^e
        pts_vel.row_mut(output[0] as usize)[0] = output[1] as f64 * node_vel.row(nval[1] as usize)[0];
        pts_vel.row_mut(output[0] as usize)[1] = output[2] as f64 * node_vel.row(nval[2] as usize)[1];
    }
    pts_vel 
}

pub fn g_fn(node_pos: &Array2<f64>, node_conn: &Array2<usize>, pt: &Array1<f64>) -> Array1<f64> {

    let ne = node_conn.shape()[0]; // number of mesh elements
    let es = node_conn.shape()[1]; // number of points in mesh element

    // memory for final calculation
    let mut lhs: Array2<f64> = Array2::zeros((2, 2 * ne));
    let mut rhs: Array2<f64> = Array2::zeros((2, ne));
    let mut r: Array2<f64> = Array2::zeros((es-1, ne));

    let mut pos = 0; // position in matrix 
    while pos < lhs.len() {

        if pos < 2*ne { //upper half of the matrix

            for n in 1..es { if pos == 2*ne {break;} //for each N of mesh element
                lhs[(0, pos)] = node_pos.row(node_conn.row(pos/2)[n])[0] 
                                        - node_pos.row(node_conn.row(pos/2)[0])[0];
                pos += 1;
            }
        }

        if 2*ne <= pos && pos < 4*ne { //lower half of the matrix

            for n in 1..es { if pos == 4*ne {break;} //for each N of mesh element
                lhs[(1, pos-2*ne)] = node_pos.row(node_conn.row((pos-2*ne)/2)[n])[1] 
                                        - node_pos.row(node_conn.row((pos-2*ne)/2)[0])[1];
                pos += 1;
            }
        }
    } 

    pos = 0;
    while pos < rhs.len() {
        if pos < ne {
            rhs[(0, pos)] = pt[0] - node_pos.row(node_conn.row(pos)[0])[0];
            pos += 1;
        }
        if ne <= pos && pos < 2*ne {
            rhs[(1, pos-ne)] = pt[1] - node_pos.row(node_conn.row(pos-ne)[0])[1];
            pos += 1;
        }
    }

    let mut lhs_slice: Array2<f64> = Array::zeros((ne-1, ne-1));
    let mut rhs_slice: Array1<f64> = Array::zeros(ne-1);
    let mut r_slice: Array1<f64> = Array::zeros(ne-1);

    pos = 0;
    while pos < rhs.row(0).len() { //could be replaced by axis_chunks_iter
        lhs_slice = lhs.slice(s![0..2,(0+2*pos)..(2+2*pos)]).to_owned();
        rhs_slice = rhs.slice(s![0..2, pos]).to_owned();
        r_slice = lhs_slice.solve_into(rhs_slice).unwrap();

        r.column_mut(pos)[0] = r_slice[0]; 
        r.column_mut(pos)[1] = r_slice[1]; 

        pos += 1;
    }
    let mut idx_and_position: Array1<f64> = Array::zeros(3);
    for (i, elem) in r.axis_iter_mut(Axis((1))).enumerate() {
        if elem[0] > 0. && elem[1] > 0. && 1.-elem[0]-elem[1] > 0. {
            idx_and_position[0] = i as f64;
            idx_and_position[1] = elem[0];
            idx_and_position[2] = elem[1];
        }
    } 
    idx_and_position 
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
