#![allow(non_snake_case)]

use crate::LAMBDA_BYTES_LEN;
use chacha20::ChaCha20Rng;
//use curv::elliptic::curves::{Scalar, Secp256k1};
//use curv::BigInt;
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
/*
pub fn outer_mul(u: &[Scalar<Secp256k1>], v: &[Scalar<Secp256k1>]) -> Vec<Scalar<Secp256k1>> {
    let m = u.len();
    let l = v.len();
    let mut res = Vec::new();
    for i in 0..m {
        for j in 0..l {
            res.push(&u[i] * &v[j])
        }
    }
    return res;
}

pub fn outer_add(u: &[BigInt], v: &[BigInt]) -> Vec<BigInt> {
    let m = u.len();
    let l = v.len();
    let mut res = Vec::new();
    for i in 0..m {
        for j in 0..l {
            res.push(&u[i] + &v[j])
        }
    }
    return res;
}

 */

#[macro_export]
macro_rules! make_array {
    ($n:expr, $constructor:expr) => {{
        let mut items: [_; $n] = mem::MaybeUninit::uninit().assume_init();
        for (_, place) in items.iter_mut().enumerate() {
            ptr::write(place, $constructor);
        }
        items
    }};
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
