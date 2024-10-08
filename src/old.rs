
// fn calc_b(node_pos: Array2<f64>, node_conn: Array2<usize>) -> Array2<f64> {

//     let mut B: Array2<f64> = Array::zeros((0, 6));
//     let gradient_local_N0 = arr1(&[-1, -1]);
//     let gradient_local_N1 = arr1(&[1, 0]);
//     let gradient_local_N2 = arr1(&[0, 1]);

//     for N in node_conn.rows() {
//         let J =  arr2(&[[node_pos.row(N[1])[0] - node_pos.row(N[0])[0], node_pos.row(N[1])[1] - node_pos.row(N[0])[1]], 
//                         [node_pos.row(N[2])[0] - node_pos.row(N[0])[0], node_pos.row(N[3])[1] - node_pos.row(N[0])[1]]]);

//         let gradient_global_N0 = J.solve_into(gradient_local_N0).unwrap();
//         let gradient_global_N1 = J.solve_into(gradient_local_N1).unwrap();
//         let gradient_global_N2 = J.solve_into(gradient_local_N2).unwrap();
//         let intermediate = arr2(&[[gradient_global_N0[0], 0, gradient_global_N1[0], 0, gradient_global_N2[0], 0], 
//                                   [0, gradient_global_N0[1], 0, gradient_global_N1[1], 0, gradient_global_N2[1]], 
//                                   [gradient_global_N0[1], gradient_global_N0[0], gradient_global_N1[1], 
//                                    gradient_global_N1[0], gradient_global_N2[1], gradient_global_N2[0]]]); 

//         B.push_row(ArrayView::from(&intermediate)).unwrap();
        
//     }

//     B
// }

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

// velocity from mesh gets transferred onto particle
pub fn f(node_pos: Array2<f64>, node_conn: Array2<usize>, node_vel: Array2<f64>, global_particle_pos: Array1<f64>) -> Array2<f64> {

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
            // println!("{:?}", local_particle_pos.shape())
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