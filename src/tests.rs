

#[cfg(test)]
mod tests {
    use sandbox_rust::*;
    use ndarray::*;

    #[test]
    fn acceleration_test() {
        let nodes = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]),
            conn: arr2(&[[0,3,2], [0,1,3]]),
        };
        let elmnt_width: f64 = 1.;
        let (rho, moe, nu): (f64, f64, f64) = (1.0, 130., 0.2);

        let n_vel: Array2<f64> = arr2(&[[0., 0.], [1., 0.], [0., 1.], [1., 1.]]);
        let ext_f: Array2<f64> = Array2::zeros(n_vel.dim()); 

        let accel = get_acceleration(&nodes, &n_vel, &ext_f, rho, moe, nu);

        let f: f64 = moe / (1.-nu.powi(2)) * (1.+nu) * elmnt_width * 0.5;
        let m: f64 = 1./2. * 1./3.;

        let num: Array2<f64> = arr2(&[[-f, -f], [ f, -f], [-f,  f], [ f,  f]]);
        let den: Array2<f64> = arr2(&[[m*2., m*2.], [m, m], [m, m], [m*2., m*2.]]);

        let r: Array2<f64> = num/den;
        let mut r_flat: Array1<f64> = Array1::zeros(n_vel.len());

        for (pos, val) in r.iter().enumerate() {r_flat[pos] = *val};    

        assert_eq!(accel, r_flat);
    }

    #[test]
    fn mass_test() {
        let nodes = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]),
            conn: arr2(&[[0,3,2], [0,1,3]]),
        };
        let rho: f64 = 1.0;
        let m: f64 = 1./2. * 1./3.;

        let mass: Array2<f64> = get_mass(&nodes, rho);
        let result: Array2<f64> = arr2(&[[m*2., m*2.], [m, m], [m, m], [m*2., m*2.]]);
        assert_eq!(mass, result);
    } 

    #[test]
    fn internal_forces_test() {
        let nodes = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]),
            conn: arr2(&[[0,3,2], [0,1,3]]),
        };
        let nodes_vel: Array2<f64> = arr2(&[[0., 0.], [1., 0.], [0., 1.], [1., 1.]]);
        let moe: f64 = 130.;
        let nu: f64 = 0.2;
        let elmnt_width: f64 = 1.;

        let f: f64 = moe / (1.-nu.powi(2)) * (1.+nu) * elmnt_width * 0.5;

        let forces: Array2<f64> = get_nodal_forces(&nodes, &nodes_vel, moe, nu);
        let result: Array2<f64> = arr2(&[[-f, -f], [ f, -f], [-f,  f], [ f,  f]]);

        assert_eq!(forces, result); 

    }

    #[test]
    fn strain_test() {
        let nodes = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]),
            conn: arr2(&[[0,2,3], [0,1,3]]),
        };
        let nodes_vel: Array2<f64> = arr2(&[[0., 0.], [1., 0.], [0., 1.], [1., 1.]]);
        let result = arr1(&[1., 1., 0.]);

        let output = get_strain(&nodes, &nodes_vel);

        for row in output.axis_iter(Axis(0)) {
           assert_eq!(row, result); 
        }
    }

    #[test]
    fn velocity_test() {
        let nodes_vel: Array2<f64> = arr2(&[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [2., 1.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]);

        let nodes = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]),
            conn: arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]),
        };

        let pts: Array2<f64> = arr2(&[[0.75, 0.25], [0.25, 0.75], // element 0 and 2
                                      [1.50, 0.75], [1.75, 1.25], // element 3 and 5
                                      [1.25, 1.50]]); // element 7

        let should_be: Array2<f64> = arr2(&[[0.50, 0.25], [0.50, 0.25], // element 0 and 2
                                            [1.00, 0.50], [0.50, 0.25], // element 3 and 5
                                            [1.00, 0.50]]); // element 7

        let vars = get_particle_velocity(&nodes, &nodes_vel, &pts); 
        
        assert_eq!(vars, should_be);
        
    }


    #[test]
    fn global_to_local_test() {

        let nodes = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [2.,0.], [0.,1.], [1.,1.], [2.,1.], [0.,2.], [1.,2.], [2.,2.]]),
            conn: arr2(&[[0,1,4], [1,2,4], [0,3,4], [2,5,4], [3,4,6], [4,5,8], [4,7,6], [4,8,7]]),
        };

        let pt = arr1(&[0.60, 0.25]);
        let output = get_particle_position(&nodes, &pt);
        let should_be = (0, 0.35, 0.25);
        assert_eq!(output, should_be);
        
        let pt = arr1(&[0.25, 0.50]);
        let output = get_particle_position(&nodes, &pt);
        let should_be = (2, 0.25, 0.25);
        assert_eq!(output, should_be);

        let pt = arr1(&[1.50, 0.75]);
        let output = get_particle_position(&nodes, &pt);
        let should_be = (3, 0.25, 0.50);
        assert_eq!(output, should_be);

        let pt = arr1(&[0.50, 1.75]);
        let output = get_particle_position(&nodes, &pt);
        let should_be = (6, 0.25, 0.50);
        assert_eq!(output, should_be);

        let pt = arr1(&[1.25, 1.75]);
        let output = get_particle_position(&nodes, &pt);
        let should_be = (7, 0.25, 0.50);
        assert_eq!(output, should_be);
    }
    
}