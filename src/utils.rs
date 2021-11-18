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
pub fn prng(key: [u8; LAMBDA_BYTES_LEN]) -> [u8; LAMBDA_BYTES_LEN * 2] {
    // TODO: handle 128 bits key
    let mut extended_key = [0; LAMBDA_BYTES_LEN * 2];
    extended_key[..LAMBDA_BYTES_LEN].copy_from_slice(&key);
    let mut rng = ChaCha20Rng::from_seed(extended_key);
    let mut out = [0; LAMBDA_BYTES_LEN * 2];
    rng.try_fill_bytes(&mut out[..]);
    out
}

// pub fn xor<T>(a: &[T], b: &[T]) -> Vec<T>
// where
//     T: BitXor + Clone,
//     Vec<T>: FromIterator<<T as BitXor>::Output>,
// {
//     let min = min(a.len(), b.len());
//     let max = max(a.len(), b.len());
//     let mut res: Vec<_> = a[0..min]
//         .iter()
//         .zip(b[0..min].to_vec())
//         .map(|byte| byte.0.clone().bitxor(byte.1.clone()))
//         .collect();
//     match a.len().cmp(&b.len()) {
//         Ordering::Less => res.extend_from_slice(&b[min..max]),
//         Ordering::Greater => res.extend_from_slice(&a[min..max]),
//         Ordering::Equal => (),
//     }
//     res
// }
pub fn xor_slices_in_place<const N: usize>(lhs: &mut [u8; N], rhs: &[u8; N]) {
    lhs.iter_mut().zip(rhs.iter()).for_each(|(l, r)| *l ^= *r);
}
pub fn xor_slices<const N: usize>(op1: &[u8; N], op2: &[u8; N]) -> [u8; N] {
    let mut out: [u8; N] = op1.clone();
    out.iter_mut()
        .zip(op2.iter())
        .for_each(|(out, op2_e)| *out ^= *op2_e);
    out
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

// Can be replaced with https://github.com/rust-lang/rust/issues/89379 once stabilized.
pub fn array_from_fn<F, T, const N: usize>(mut cb: F) -> [T; N]
where
    F: FnMut(usize) -> T,
{
    let mut idx = 0;
    [(); N].map(|_| {
        let res = cb(idx);
        idx += 1;
        res
    })
}

mod tests {
    use crate::utils::prng;

    #[test]
    fn test_prng_zero() {
        let key = [0u8; 16];
        let out = prng(key);
        println!("out : {:?}", out);
    }
}
