
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