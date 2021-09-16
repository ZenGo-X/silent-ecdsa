#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]

mod dpf;
pub mod dspf;
pub mod fft;
mod poly;
mod utils;

const LAMBDA: usize = 128;
const LAMBDA_BYTES_LEN: usize = LAMBDA / 8;
const N: usize = 32;

#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;
