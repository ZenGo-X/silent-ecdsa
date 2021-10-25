#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]

mod dpf;
pub mod dspf;
pub mod fft;

#[macro_use]
mod utils;
mod keygen;
mod poly;
mod sign;
mod tests;

const LAMBDA: usize = 128;
const LAMBDA_BYTES_LEN: usize = LAMBDA / 8;
//Possible factors:
//  2 - up to power of 5.
//  3, 149, 631
const N: usize = 4 * 149;
const t: usize = 4; // todo: validate,
const c: usize = 2; // todo: validate
const n: usize = 2; // num of parties
pub const use_cyclotomic: bool = true;
pub const use_fft: bool = true;

#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;
