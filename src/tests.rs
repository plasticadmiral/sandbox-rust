

#[cfg(test)]
mod tests {
    use sandbox_rust::*;
    use ndarray::*;

    #[test]
    fn sim() {
        let n: Mesh = get_mesh(8, 2);
        let (rho, moe, nu, dt): (f64, f64, f64, f64) = (1., 130., 0.2, 1.); 

        let mut n_vel: Array2<f64> = Array2::zeros(n.pos.dim());
        let mut a: Array2<f64> = Array2::zeros(n.pos.dim());
        let mut disp: Array2<f64> = Array2::zeros(n.pos.dim());
        let mut ext_f: Array2<f64> = Array2::zeros(n.pos.dim());

        //row 7 and row 15 of ext_f make x to 1 and 2
        ext_f.row_mut(7)[0] = 1.;
        ext_f.row_mut(15)[0] = 2.;

        for l in 0..10 {
            verletstep1(&mut disp, &mut n_vel, &a, dt);
            a = get_acceleration(&n, &n_vel, &ext_f, rho, moe, nu);
            
            //row 0 & 8 of a make x and y to 0
            a.row_mut(0).assign(&arr1(&[0.,0.]));
            a.row_mut(8).assign(&arr1(&[0.,0.]));
            verletstep2(&mut n_vel, &a, dt);
        }
        // 5th timestep
        // [[0.0, 0.0],
        // [0.0, 0.0],
        // [0.0, 0.0],
        // [0.0, 0.0],
        // [9178005.642361114, 6317328.559027783],
        // [-80650005.42534725, -26724408.637152787],
        // [428203347.98177105, -70345816.94878478],
        // [-1310375332.2526047, 960797656.5755213],
        // [0.0, 0.0],
        // [0.0, 0.0],
        // [0.0, 0.0],
        // [0.0, 0.0],
        // [17760036.89236112, -4886990.017361114],
        // [-100185129.12326394, 39597455.5121528],
        // [147218842.77343753, -39751741.64496529],
        // [22900108.014322974, -336707568.5221356]],


        // 10th timestep
        // [[0.0, 0.0],
        // [-2.8298363877830957e20, 3.21527849046472e20],
        // [1.2074140023625346e21, -1.5118699354617996e21],
        // [2.1111280185852586e20, 1.986962450935603e21],
        // [-2.8397887017787098e22, 1.737275332340676e22],
        // [1.599423868598454e23, -1.2282805450955085e23],
        // [-5.2631194763474104e23, 4.4425388753637895e23],
        // [1.4083450825013607e24, -1.4233424646050618e24],
        // [0.0, 0.0],
        // [3.223442605923869e20, -1.696911709807679e20],
        // [-1.6862446891930752e21, 1.1316589836272437e21],
        // [3.5568891506146473e21, -3.3141355541587223e21],
        // [9.02981761759937e21, -4.050265342023486e21],
        // [-8.692373941346358e22, 7.3775565796026714e22],
        // [2.7721622655346035e23, -3.038352601226352e23],
        // [-4.1599425817342005e23, 5.5700181500277426e23]]
        println!("{:?}", disp);

    }
    #[test]
    fn verlet_test() {
        let n = Mesh {
            pos: arr2(&[[0.,0.], [1.,0.], [0.,1.], [1.,1.]]),
            conn: arr2(&[[0,3,2], [0,1,3]]),
        };

        let (rho, moe, nu, dt): (f64, f64, f64, f64) = (1., 130., 0.2, 1.); 

        let mut n_vel: Array2<f64> = Array2::zeros(n.pos.dim());
        let mut a: Array2<f64> = Array2::zeros(n.pos.dim());
        let mut disp: Array2<f64> = Array2::zeros(n.pos.dim());

        let  ext_f: Array2<f64> = arr2(&[[2.,2.], [1.,1.], [1.,1.], [2.,2.]]);

        for l in 0..2 {
            verletstep1(&mut disp, &mut n_vel, &a, dt);

            // disp + n.pos for acceleration ?
            a = get_acceleration(&n, &n_vel, &ext_f, rho, moe, nu);
            verletstep2(&mut n_vel, &a, dt);
        }

        assert_eq!(n_vel, 9. * Array2::ones((4, 2)));
        assert_eq!(disp, 6. * Array2::ones((4, 2)));
        
    }

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

        // should i make this negative ?
        let f: f64 = moe / (1.-nu.powi(2)) * (1.+nu) * elmnt_width * 0.5;
        let m: f64 = 1./2. * 1./3.;

        let num: Array2<f64> = arr2(&[[ f,  f], [-f,  f], [f, -f], [-f, -f]]);
        let den: Array2<f64> = arr2(&[[m*2., m*2.], [m, m], [m, m], [m*2., m*2.]]);

        let r: Array2<f64> = num/den;

        assert_eq!(accel, r);
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