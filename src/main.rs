extern crate gf256;

use gf256::Gf256;

fn main() {
    let a = Gf256::from_byte(3);
    let b = Gf256::from_byte(4);
    let c = a + b;

    let d = a / b;

    println!("c = {:?}", c);
    println!("d = {:?}", d);
}
