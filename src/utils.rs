#![allow(non_snake_case)]

use crate::LAMBDA_BYTES_LEN;
use chacha20::ChaCha20Rng;
use rand_core::{RngCore, SeedableRng};
use std::cmp::{max, min, Ordering};
use std::convert::TryInto;
use std::iter::FromIterator;
use std::ops::BitXor;

// prng from |lambda| to |2 lambda|
pub fn prng(key: &[u8]) -> Vec<u8> {
    // TODO: handle 128 bits key
    let mut new_key = key.to_vec();
    while new_key.len() < 32 {
        new_key.extend_from_slice(key);
    }
    let mut rng = ChaCha20Rng::from_seed(new_key[0..32].try_into().unwrap());
    let mut out = [0; LAMBDA_BYTES_LEN * 2];
    let _ = rng.try_fill_bytes(&mut out[..]);
    out.into()
}

pub fn xor<T>(a: &[T], b: &[T]) -> Vec<T>
where
    T: BitXor + Clone,
    Vec<T>: FromIterator<<T as BitXor>::Output>,
{
    let min = min(a.len(), b.len());
    let max = max(a.len(), b.len());
    let mut res: Vec<_> = a[0..min]
        .iter()
        .zip(b[0..min].to_vec())
        .map(|byte| byte.0.clone().bitxor(byte.1.clone()))
        .collect();
    match a.len().cmp(&b.len()) {
        Ordering::Less => res.extend_from_slice(&b[min..max]),
        Ordering::Greater => res.extend_from_slice(&a[min..max]),
        Ordering::Equal => (),
    }
    res
}

mod tests {
    use crate::utils::prng;

    #[test]
    fn test_prng_zero() {
        let key = [0u8; 32].to_vec();
        let out = prng(&key[..]);
        println!("out : {:?}", out);
    }
}
