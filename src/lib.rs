#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
mod dpf;
mod dspf;
mod utils;
mod poly;

const LAMBDA: usize = 128;
const LAMBDA_BYTES_LEN: usize = LAMBDA / 8;
const N: usize = 32;
