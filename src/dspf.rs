use curv::arithmetic::Zero;
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::dpf::DPF;
use crate::N;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DSPF {
    pub key: Vec<DPF>,
}

impl DSPF {
    pub fn init(l: usize) -> Self {
        let mut dspf = DSPF {
            key: Vec::with_capacity(l),
        };
        for _ in 0..l {
            dspf.key.push(DPF::init());
        }
        dspf
    }
    pub fn gen(alpha_vec: &[BigInt], beta_vec: &[Scalar<Secp256k1>]) -> (Self, Self) {
        // TODO: make sure there are no repetitions in alpha_vec?
        let alpha_vec_len = alpha_vec.len();
        assert_eq!(alpha_vec_len, beta_vec.len());
        let (key0, key1): (Vec<_>, Vec<_>) = (0..alpha_vec_len)
            .map(|i| {
                assert!(alpha_vec[i] >= BigInt::zero());
                assert!(alpha_vec[i] < BigInt::from(2) * BigInt::from(N as u32));
                DPF::gen(&alpha_vec[i], &beta_vec[i])
            })
            .unzip();
        return (DSPF { key: key0.clone() }, DSPF { key: key1.clone() });
    }

    pub fn eval(&self, b: &u8, x: &BigInt, two_n_flag: bool) -> Scalar<Secp256k1> {
        if two_n_flag {
            assert!(x < &BigInt::from(2 * N as u32));
        } else {
            assert!(x < &BigInt::from(N as u32));
        }
        assert!(x >= &BigInt::zero());
        self.key
            .iter()
            .fold(Scalar::zero(), |acc, dpf| acc + dpf.eval(*b, x))
    }

    pub fn full_eval_N(&self, b: u8) -> Vec<Scalar<Secp256k1>> {
        self.key
            .par_iter()
            .map(|dpf_key| dpf_key.full_eval_N(b))
            .reduce(
                || vec![Scalar::zero(); N],
                |mut a, b| {
                    a.iter_mut()
                        .zip(&b)
                        .for_each(|(a_i, b_i)| *a_i = &*a_i + b_i);
                    a
                },
            )
    }

    pub fn full_eval_2N(&self, b: u8) -> Vec<Scalar<Secp256k1>> {
        self.key
            .par_iter()
            .map(|dpf_key| dpf_key.full_eval_2N(b))
            .reduce(
                || vec![Scalar::zero(); 2 * N],
                |mut a, b| {
                    a.iter_mut()
                        .zip(&b)
                        .for_each(|(a_i, b_i)| *a_i = &*a_i + b_i);
                    a
                },
            )
    }
}

#[cfg(test)]
mod tests {
    use crate::dspf::DSPF;
    use crate::N;
    use curv::elliptic::curves::{Scalar, Secp256k1};
    use curv::BigInt;

    #[test]
    fn test_simple_dspf_eval_and_full_eval() {
        let alpha_vec = &[BigInt::from(10), BigInt::from(20), BigInt::from(30)].to_vec();
        let beta_vec: &[Scalar<Secp256k1>] =
            &[Scalar::random(), Scalar::random(), Scalar::random()].to_vec();
        let mut fe0_vec = Vec::new();
        let mut fe1_vec = Vec::new();

        let (key0, key1) = DSPF::gen(&alpha_vec[..], &beta_vec[..]);
        let full_eval_0 = key0.full_eval_N(0);
        let full_eval_1 = key1.full_eval_N(1);
        for x in 0..N {
            let x_bn = BigInt::from(x as u32);
            let fe0 = &full_eval_0[x];
            let fe1 = &full_eval_1[x];

            fe0_vec.push(fe0.clone());
            fe1_vec.push(fe1.clone());

            if x_bn == alpha_vec[0] {
                assert_eq!(fe0 + fe1, beta_vec[0]);
            } else if x_bn == alpha_vec[1] {
                assert_eq!(fe0 + fe1, beta_vec[1]);
            } else if x_bn == alpha_vec[2] {
                assert_eq!(fe0 + fe1, beta_vec[2]);
            } else {
                assert_eq!(
                    fe0.to_bigint(),
                    Scalar::<Secp256k1>::group_order() - fe1.to_bigint()
                );
            }
        }
        let p0_shares = key0.full_eval_N(0u8);
        let p1_shares = key1.full_eval_N(1u8);
        assert_eq!(p0_shares, fe0_vec);
        assert_eq!(p1_shares, fe1_vec);
    }
}
