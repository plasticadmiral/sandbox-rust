

#[cfg(test)]
mod tests {
    use sandbox_rust::*;
    use ndarray::*;



    fn particle_velocity_test() {

    }

    #[test]
    fn b_fn_test() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]);
        let node_conn: Array2<usize> = arr2(&[[0,2,3], [0,1,3]]);
        let node_vel: Array2<f64> = arr2(&[[0., 0.], [1., 0.], [0., 1.], [1., 1.]]);
        let result = arr1(&[1., 1., 0.]);
        let b = get_b(&node_pos, &node_conn);
        let mut v: Array2<f64> = Array2::zeros((3, 6));

        for (idx, elmnt) in node_conn.axis_iter(Axis(0)).enumerate() {

            let mut r = Array2::zeros((3, 2));

            for (pos, node_idx) in elmnt.iter().enumerate() {

                r.row_mut(pos).assign(&node_vel.row(*node_idx)); 
            }
            for (pos, val) in r.iter().enumerate() {v[[idx, pos]] = *val;}        
        }
        
        for (idx, b_layer) in b.axis_iter(Axis(0)).enumerate() {
            
           assert_eq!(b_layer.dot(&v.row(idx)), result); 
        }
    }

    #[test]
    fn velocity_test() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);
        let node_vel: Array2<f64> = arr2(&[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [2., 1.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]);

        let pt1: Array2<f64> = arr2(&[[0.25, 0.75], [1.50, 0.75]]);
        let pt2: Array2<f64> = arr2(&[[0.75, 0.25], [1.25, 1.50]]);
        let pt3: Array2<f64> = arr2(&[[1.75, 1.25], [1.50, 1.25]]);

        // returns result [[0.50, 0.25], [1.0, 0.5]]
        let var1 = get_particle_velocity(&node_pos, &node_conn, &node_vel, &pt1); 
        let var2 = get_particle_velocity(&node_pos, &node_conn, &node_vel, &pt2); 
        let var3 = get_particle_velocity(&node_pos, &node_conn, &node_vel, &pt3); 

        for row in var1.rows() {assert_eq!(row[0]/row[1], 2.);}
        for row in var2.rows() {assert_eq!(row[0]/row[1], 2.);}
        for row in var3.rows() {assert_eq!(row[0]/row[1], 2.);}
        
    }


    #[test]
    fn single_pt_in_mesh_2() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);

        let pt = arr1(&[0.25, 0.50]);
        // let output = get_particle_position(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[2., 0.25, 0.25]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (2, 0.25, 0.25);
        assert_eq!(output, should_be);
    }

    #[test]
    fn single_pt_in_mesh_3() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);

        let pt = arr1(&[1.50, 0.75]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[3., 0.25, 0.50]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (3, 0.25, 0.50);
        assert_eq!(output, should_be);

    }

    #[test]
    fn single_pt_in_mesh_7() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);

        let pt = arr1(&[1.25, 1.75]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[7., 0.25, 0.50]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (7, 0.25, 0.50);
        assert_eq!(output, should_be);

    }

    #[test]
    fn mult_pts_in_mesh_6() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);

        let pt = arr1(&[0.50, 1.75]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[6., 0.25, 0.50]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (6, 0.25, 0.50);
        assert_eq!(output, should_be);

        let pt = arr1(&[0.75, 1.75]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[6., 0.5, 0.25]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (6, 0.5, 0.25);
        assert_eq!(output, should_be);

    }
    
    #[test]
    fn pts_in_mesh_0_and_2() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);
        let pt = arr1(&[0.60, 0.25]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[0., 0.35, 0.25]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (0, 0.35, 0.25);
        assert_eq!(output, should_be);

        let pt = arr1(&[0.50, 0.75]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[2., 0.25, 0.50]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (2, 0.25, 0.50);
        assert_eq!(output, should_be);
    }

    #[test]
    fn pts_in_mesh_2_and_7() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);

        let pt = arr1(&[0.30, 0.80]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[2., 0.5, 0.3]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (2, 0.5, 0.3);
        assert_eq!(output, should_be);

        let pt = arr1(&[1.25, 1.75]);
        // let output = g_fn(&node_pos, &node_conn, &pt);
        // let should_be = arr1(&[7., 0.25, 0.50]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (7, 0.25, 0.50);
        assert_eq!(output, should_be);
    }
    
}