#![allow(non_snake_case)]

mod dpf;
mod dspf;
mod utils;

const LAMBDA: usize = 128;
const LAMBDA_BYTES_LEN: usize = LAMBDA / 8;
const N: usize = 32;
