#![allow(unused)]
#![allow(non_snake_case)]

use std::{iter::FromIterator, result, usize, process::id, f32::consts::E};
// use sandbox_rust::*;
use plotly::{ImageFormat, Plot, Scatter};
use ndarray::{*, linalg::Dot};
use ndarray_linalg::*;
use ndarray::{*, Array3};
use ndarray_linalg::*;
// use plotters::prelude::*;

pub fn get_internal_forces(n_pos: &Array2<f64>, n_conn: &Array2<usize>, n_vel: &Array2<f64>, moe: f64, nu: f64) -> Array2<f64> {

    let b = get_b(n_pos, n_conn);
    let stress = get_stress(n_pos, n_conn, n_vel, moe, nu);
    let mut int_forces: Array2<f64> = Array2::zeros((b.shape()[0], b.shape()[2])); //8x6

    for (idx, b_layer) in b.axis_iter(Axis(0)).enumerate() {

        int_forces.row_mut(idx).assign(&b_layer.t().dot(&stress.row(idx)));
    }

    int_forces
}


pub fn get_stress(n_pos: &Array2<f64>, n_conn: &Array2<usize>, n_vel: &Array2<f64>, moe: f64, nu: f64) -> Array2<f64> {

    let mut stress: Array2<f64> = Array2::zeros((n_conn.shape()[0], 3));
    let strain: Array2<f64> = get_strain(n_pos, n_conn, n_vel);
    let D: Array2<f64> = moe/(1.-nu.powi(2)) * arr2(&[[1., nu, 0.], [nu, 1., 0.], [0., 0., (1.-nu)/2.]]);    

    stress = D.dot(&strain); 
    stress
}

pub fn get_strain(n_pos: &Array2<f64>, n_conn: &Array2<usize>, n_vel: &Array2<f64>) -> Array2<f64> {

    let mut strain: Array2<f64> = Array2::zeros((n_conn.shape()[0], 3));
    let mut v: Array2<f64> = Array2::zeros((3, 6));
    let b: Array3<f64> = get_b(&n_pos, &n_conn);

    for (idx, elmnt) in n_conn.axis_iter(Axis(0)).enumerate() {

        let mut r = Array2::zeros((3, 2)); // holds velocities of each element

        for (pos, node_idx) in elmnt.iter().enumerate() {

            r.row_mut(pos).assign(&n_vel.row(*node_idx)); 
        }
        for (pos, val) in r.iter().enumerate() {v[[idx, pos]] = *val;} // flattens r into v      
    }
    
    for (idx, b_layer) in b.axis_iter(Axis(0)).enumerate() {

        strain.row_mut(idx).assign(&b_layer.dot(&v.row(idx))); 
    }
    
    strain
}


pub fn get_delta_e(n_pos: &Array2<f64>, n_conn: &Array2<usize>, n_vel: &Array2<f64>, pts: &Array2<f64>, delta_t: &f64) -> Array2<f64> {
    let mut delta_e: Array2<f64> = Array2::zeros((pts.shape()[0], 3));
    let mut v: Array2<f64> = Array2::zeros((0, 2));
    let mut idx_mesh: usize = 0;

    let B = get_b(&n_pos, &n_conn);

    for (idx_pt, pt) in pts.axis_iter(Axis(0)).enumerate() {

        (idx_mesh, _, _) = get_particle_position(&n_pos, &n_conn, &pt.to_owned());
        
        for idx_node in n_conn.row(idx_mesh) {v.push_row(n_vel.row(*idx_node));}

        let v_flat = Array::from_iter(v.iter().cloned());
        let B_calc = B.index_axis(Axis(0), idx_mesh);

        delta_e.row_mut(idx_pt).assign(&B_calc.dot(&v_flat));
    } 

    delta_e = delta_e * *delta_t;
    delta_e
}


pub fn get_b(n_pos: &Array2<f64>, n_conn: &Array2<usize>) -> Array3<f64> {

    let mut B: Array3<f64> = Array3::zeros((n_conn.shape()[0], 3, 6));
    let mut J: Array2<f64> = Array2::zeros((2,2));
    let mut grd_N: Array2<f64> = arr2(&[[-1., -1.], [1., 0.], [0., 1.]]);

    for (idx, N) in n_conn.axis_iter(Axis(0)).enumerate() {

        grd_N = arr2(&[[-1., -1.], [1., 0.], [0., 1.]]);

        for r in 0..2 {
            for c in 0..2 {

                J[[r,c]] = n_pos.row(N[r+1])[c] - n_pos.row(N[0])[c];
            }
        }
        for layer in 0..3 { J.solve_inplace(&mut grd_N.row_mut(layer)); }

        for c in 0..6 {
            if c % 2 == 0 { 
                B[[idx, 0, c]] = grd_N[[c/2, 0]]; 
                B[[idx, 2, c]] = grd_N[[c/2, 1]];
            } else {
                B[[idx, 1, c]] = grd_N[[c/2, 1]]; 
                B[[idx, 2, c]] = grd_N[[c/2, 0]];
            }
        }
    }
    B
}

pub fn get_particle_velocity(n_pos: &Array2<f64>, n_conn: &Array2<usize>, n_vel: &Array2<f64>, pts: &Array2<f64>) -> Array2<f64> {

    let mut pts_vel:Array2<f64> = Array2::zeros(pts.dim()); 
    let mut N = Array1::zeros(n_conn.shape()[1]);
    let mut v = Array2::zeros((n_conn.shape()[1], n_vel.shape()[1]));

    for (idx_pt,pt) in pts.axis_iter(Axis(0)).enumerate() {

        let (idx_elmnt, xpos, ypos) = get_particle_position(&n_pos, &n_conn, &pt.to_owned()); 

        N.assign(&arr1(&[1.-xpos-ypos, xpos, ypos]));
        for (idx, node) in n_conn.row(idx_elmnt).iter().enumerate() {
            v.row_mut(idx).assign(&n_vel.row(*node));
        }

        // V_p = N(X_p)V^e
        let r = N.dot(&v);
        pts_vel.row_mut(idx_pt).assign(&arr1(&[r[0], r[1]]));
    }
    pts_vel 
}

pub fn get_particle_position(n_pos: &Array2<f64>, n_conn: &Array2<usize>, pt: &Array1<f64>) -> (usize, f64, f64) {
    
    let (mut idx_elmnt, mut pos_x, mut pos_y): (usize, f64, f64) = (0, 0., 0.);
    
    let mut lhs: Array2<f64> = Array2::zeros((2,2));
    let mut rhs: Array1<f64> = Array1::zeros(2);

    for (idx, N) in n_conn.axis_iter(Axis(0)).enumerate() {
        for r in 0..2 { 
            for c in 0..2 {

                lhs[[r, c]] = n_pos.row(N[c+1])[r] - n_pos.row(N[0])[r];
            }

            rhs[r] = pt[r] - n_pos.row(N[0])[r];
        }

        lhs.solve_inplace(&mut rhs);

        if rhs.iter().all(|&x| x > 0.) && 1.-rhs.sum() > 0. {

            (idx_elmnt, pos_x, pos_y) = (idx, rhs[0], rhs[1]);
            break;
        }
    }
    (idx_elmnt, pos_x, pos_y)
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
