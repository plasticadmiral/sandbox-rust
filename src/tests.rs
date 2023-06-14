

#[cfg(test)]
mod tests {
    use sandbox_rust::*;
    use ndarray::*;

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