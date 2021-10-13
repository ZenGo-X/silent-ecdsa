#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]

mod dpf;
pub mod dspf;
pub mod fft;
mod poly;
#[macro_use]
mod utils;
mod keygen;
mod sign;
mod tests;

const LAMBDA: usize = 128;
const LAMBDA_BYTES_LEN: usize = LAMBDA / 8;
const N: usize = 32;
const t: usize = N / 16; // todo: validate
const c: usize = 2; // todo: validate
const n: usize = 3; // num of parties

#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;
