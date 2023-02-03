#![allow(unused)]
use sandbox_rust::*;



// #[allow(dead_code)]
// #[allow(unused_mut)]
// #[allow(unused_variables)]
#[allow(non_snake_case)]

fn main() {
    println!("{}", GREETING);

    let (x, y) = (13, 13);
let (mut X, mut Y, mut Z) = poisson_data(&x, &y);
println!("x--------------x");
println!("{X}");
println!("x--------------x");
println!("{Y}");
println!("x--------------x");
println!("{Z}");
}


