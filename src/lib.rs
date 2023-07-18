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

pub fn get_acceleration(n: &Mesh, n_vel: &Array2<f64>, ext_forces: &Array2<f64>, rho: f64, moe: f64, nu: f64) -> Array1<f64>{

    let int_forces: Array2<f64> = get_nodal_forces(&n, &n_vel, moe, nu);
    let masses: Array2<f64> = get_mass(&n, rho);
    
    let mut f_int: Array1<f64> = Array1::zeros(n.pos.len());
    let mut f_ext: Array1<f64> = Array1::zeros(n.pos.len());
    let mut m: Array1<f64> = Array1::zeros(n.pos.len());
    let mut a: Array1<f64> = Array1::zeros(n.pos.len());

    for (row, pos) in (0..n.pos.len()).step_by(2).enumerate() {

        (m[pos], m[pos+1]) = (masses[[row, 0]], masses[[row, 1]]);

        (f_int[pos], f_int[pos+1]) = (int_forces[[row, 0]], int_forces[[row, 1]]);

        (f_ext[pos], f_ext[pos+1]) = (ext_forces[[row, 0]], ext_forces[[row, 1]]);
    }
    a = (f_int - f_ext) / m;
    a
}

pub fn get_mass(n: &Mesh, rho: f64) -> Array2<f64> {

    let mut n_mass: Array2<f64> = Array2::zeros(n.pos.dim());
    let mut element_mass: Array2<f64> = Array2::zeros((6,6));
    let mut J: Array2<f64> = Array2::zeros((2,2)); 

    for (idx, elmnt) in n.conn.axis_iter(Axis(0)).enumerate() {

        for r in 0..2 {
            for c in 0..2 {
                J[[r,c]] = n.pos.row(elmnt[r+1])[c] - n.pos.row(elmnt[0])[c];
            }
        }
        let area: f64 = (J[(0,0)] * J[(1,1)] - J[(0,1)] * J[(1,0)]) * 0.5;

        element_mass = rho * area * 1./3. * Array::eye(6);

        for N in 0..3 {
           n_mass.row_mut(elmnt[N])[0] += element_mass.row(N*2)[N*2]; 
           n_mass.row_mut(elmnt[N])[1] += element_mass.row(N*2+1)[N*2+1]; 
        }
    }
    n_mass
}


pub fn get_nodal_forces(n: &Mesh, n_vel: &Array2<f64>, moe: f64, nu: f64) -> Array2<f64> {

    let forces: Array2<f64> = get_internal_forces(&n, &n_vel, moe, nu);
    let mut n_force: Array2<f64> = Array2::zeros(n.pos.dim());

    for (idx, elmnt) in n.conn.axis_iter(Axis(0)).enumerate() {
    
        for N in 0..3 {
            n_force.row_mut(elmnt[N])[0] += forces.row(idx)[N*2]; 
            n_force.row_mut(elmnt[N])[1] += forces.row(idx)[N*2+1]; 
        }
    }
    n_force
}

pub fn get_internal_forces(n: &Mesh, n_vel: &Array2<f64>, moe: f64, nu: f64) -> Array2<f64> {

    let b: Array3<f64> = get_b(&n);
    let stress: Array2<f64> = get_stress(&n, &n_vel, moe, nu);
    let mut int_forces: Array2<f64> = Array2::zeros((b.shape()[0], b.shape()[2])); //8x6
    let mut J: Array2<f64> = Array2::zeros((2,2)); 

    for (idx, b_layer) in b.axis_iter(Axis(0)).enumerate() {

        let N: Array1<usize> = n.conn.row(idx).to_owned();

        for r in 0..2 {
            for c in 0..2 {
                J[[r,c]] = n.pos.row(N[r+1])[c] - n.pos.row(N[0])[c];
            }
        }
        let det_J: f64 = J[(0,0)] * J[(1,1)] - J[(0,1)] * J[(1,0)];

        let var = b_layer.t().dot(&stress.row(idx)) * det_J * 0.5;
        int_forces.row_mut(idx).assign(&var);
            
    }

    int_forces
}


pub fn get_stress(n: &Mesh, n_vel: &Array2<f64>, moe: f64, nu: f64) -> Array2<f64> {

    let mut stress: Array2<f64> = Array2::zeros((n.conn.shape()[0], 3));
    let strain: Array2<f64> = get_strain(&n, n_vel);
    let D: Array2<f64> = moe/(1.-nu.powi(2)) * arr2(&[[1., nu, 0.], [nu, 1., 0.], [0., 0., (1.-nu)/2.]]);    

    let var = D.dot(&strain.t()); 
    stress = var.t().to_owned();
    stress
}

pub fn get_strain(n: &Mesh, n_vel: &Array2<f64>) -> Array2<f64> {

    let mut strain: Array2<f64> = Array2::zeros((n.conn.shape()[0], 3));
    let mut v: Array2<f64> = Array2::zeros((n.conn.shape()[0], 6));
    let b: Array3<f64> = get_b(&n);

    for (idx, elmnt) in n.conn.axis_iter(Axis(0)).enumerate() {

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

pub fn get_b(n: &Mesh) -> Array3<f64> {

    let mut B: Array3<f64> = Array3::zeros((n.conn.shape()[0], 3, 6));
    let mut J: Array2<f64> = Array2::zeros((2,2));
    let mut grd_N: Array2<f64> = arr2(&[[-1., -1.], [1., 0.], [0., 1.]]);

    for (idx, N) in n.conn.axis_iter(Axis(0)).enumerate() {

        grd_N = arr2(&[[-1., -1.], [1., 0.], [0., 1.]]);

        for r in 0..2 {
            for c in 0..2 {

                J[[r,c]] = n.pos.row(N[r+1])[c] - n.pos.row(N[0])[c];
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

pub fn get_particle_velocity(n: &Mesh, n_vel: &Array2<f64>, pts: &Array2<f64>) -> Array2<f64> {

    let mut pts_vel:Array2<f64> = Array2::zeros(pts.dim()); 
    let mut N = Array1::zeros(n.conn.shape()[1]);
    let mut v = Array2::zeros((n.conn.shape()[1], n_vel.shape()[1]));

    for (idx_pt,pt) in pts.axis_iter(Axis(0)).enumerate() {

        let (idx_elmnt, xpos, ypos) = get_particle_position(&n, &pt.to_owned()); 

        N.assign(&arr1(&[1.-xpos-ypos, xpos, ypos]));
        for (idx, node) in n.conn.row(idx_elmnt).iter().enumerate() {
            v.row_mut(idx).assign(&n_vel.row(*node));
        }

        // V_p = N(X_p)V^e
        let r = N.dot(&v);
        pts_vel.row_mut(idx_pt).assign(&arr1(&[r[0], r[1]]));
    }
    pts_vel 
}


pub fn get_particle_position(n: &Mesh,  pt: &Array1<f64>) -> (usize, f64, f64) {
    
    let (mut idx_elmnt, mut pos_x, mut pos_y): (usize, f64, f64) = (0, 0., 0.);
    
    let mut lhs: Array2<f64> = Array2::zeros((2,2));
    let mut rhs: Array1<f64> = Array1::zeros(2);

    for (idx, N) in n.conn.axis_iter(Axis(0)).enumerate() {
        for r in 0..2 { 
            for c in 0..2 {

                lhs[[r, c]] = n.pos.row(N[c+1])[r] - n.pos.row(N[0])[r];
            }

            rhs[r] = pt[r] - n.pos.row(N[0])[r];
        }

        lhs.solve_inplace(&mut rhs);

        if rhs.iter().all(|&x| x > 0.) && 1.-rhs.sum() > 0. {

            (idx_elmnt, pos_x, pos_y) = (idx, rhs[0], rhs[1]);
            break;
        }
    }
    (idx_elmnt, pos_x, pos_y)
}
pub struct Mesh {
    pub pos: Array2<f64>,
    pub conn: Array2<usize>, 
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
