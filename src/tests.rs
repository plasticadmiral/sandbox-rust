

#[cfg(test)]
mod tests {
    use sandbox_rust::*;
    use ndarray::*;

    #[test]
    fn strain_test() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]);
        let node_conn: Array2<usize> = arr2(&[[0,2,3], [0,1,3]]);
        let node_vel: Array2<f64> = arr2(&[[0., 0.], [1., 0.], [0., 1.], [1., 1.]]);
        let result = arr1(&[1., 1., 0.]);

        let output = get_strain(&node_pos, &node_conn, &node_vel);

        for row in output.axis_iter(Axis(0)) {
           assert_eq!(row, result); 
        }
    }

    #[test]
    fn velocity_test() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);
        let node_vel: Array2<f64> = arr2(&[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [2., 1.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]);

        let pts: Array2<f64> = arr2(&[[0.75, 0.25], [0.25, 0.75], // element 0 and 2
                                      [1.50, 0.75], [1.75, 1.25], // element 3 and 5
                                      [1.25, 1.50]]); // element 7

        let should_be: Array2<f64> = arr2(&[[0.50, 0.25], [0.50, 0.25], // element 0 and 2
                                            [1.00, 0.50], [0.50, 0.25], // element 3 and 5
                                            [1.00, 0.50]]); // element 7

        let vars = get_particle_velocity(&node_pos, &node_conn, &node_vel, &pts); 
        
        assert_eq!(vars, should_be);
        
    }


    #[test]
    fn global_to_local_test() {
        let node_pos: Array2<f64> = arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]);
        let node_conn: Array2<usize> = arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]);

        let pt = arr1(&[0.60, 0.25]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (0, 0.35, 0.25);
        assert_eq!(output, should_be);
        
        let pt = arr1(&[0.25, 0.50]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (2, 0.25, 0.25);
        assert_eq!(output, should_be);

        let pt = arr1(&[1.50, 0.75]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (3, 0.25, 0.50);
        assert_eq!(output, should_be);

        let pt = arr1(&[0.50, 1.75]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (6, 0.25, 0.50);
        assert_eq!(output, should_be);

        let pt = arr1(&[1.25, 1.75]);
        let output = get_particle_position(&node_pos, &node_conn, &pt);
        let should_be = (7, 0.25, 0.50);
        assert_eq!(output, should_be);
    }
    
}