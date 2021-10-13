#![allow(non_snake_case)]
// implementing https://eprint.iacr.org/2018/707.pdf figure 1

use crate::utils::{prng, xor};
use crate::{LAMBDA, LAMBDA_BYTES_LEN, N};
use curv::arithmetic::{Converter, Modulo, One, Samplable, Zero};
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;
use serde::{Deserialize, Serialize};
use std::convert::TryInto;

#[derive(Clone, Debug)]
pub struct PrngOut {
    pub s_i_L: Vec<u8>,
    pub t_i_L: u8,
    pub s_i_R: Vec<u8>,
    pub t_i_R: u8,
}

#[derive(Clone, Debug, Copy, Serialize, Deserialize)]
pub struct CWi {
    pub s_CW: [u8; LAMBDA_BYTES_LEN],
    pub t_CW_L: u8,
    pub t_CW_R: u8,
}

impl CWi {
    pub fn init() -> Self {
        CWi {
            s_CW: [0; LAMBDA_BYTES_LEN],
            t_CW_L: 0u8,
            t_CW_R: 0u8,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CWnp1(Scalar<Secp256k1>);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Key {
    pub s_i_0: [u8; LAMBDA_BYTES_LEN],
    pub CW_n_plus_1: CWnp1,
    pub CWs: [CWi; N],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DPF(pub(crate) Key);

impl DPF {
    pub fn init() -> Self {
        let key = Key {
            s_i_0: [0; LAMBDA_BYTES_LEN],
            CW_n_plus_1: CWnp1(Scalar::zero()),
            CWs: [CWi::init(); N],
        };
        DPF(key)
    }

    pub fn gen(alpha: &BigInt, beta: &Scalar<Secp256k1>) -> (Self, Self) {
        let s_0_0_vec = BigInt::sample(LAMBDA * 2).to_bytes()[0..LAMBDA_BYTES_LEN].to_vec();
        let s_1_0_vec = BigInt::sample(LAMBDA * 2).to_bytes()[0..LAMBDA_BYTES_LEN].to_vec();
        let mut s_0_0: [u8; LAMBDA_BYTES_LEN] = [0; LAMBDA_BYTES_LEN];
        let mut s_1_0: [u8; LAMBDA_BYTES_LEN] = [0; LAMBDA_BYTES_LEN];
        for i in 0..LAMBDA_BYTES_LEN {
            s_0_0[i] = s_0_0_vec[i].clone();
            s_1_0[i] = s_1_0_vec[i].clone();
        }

        let t_0_0: u8 = 0;
        let t_1_0: u8 = 1;
        s_0_0[0] = s_0_0[0] & 0xfe;
        s_1_0[0] = s_1_0[0] & 0xfe;
        let n = N;
        let mut s_0_i_minus_1 = s_0_0_vec;
        let mut s_1_i_minus_1 = s_1_0_vec;
        let mut t_0_i_minus_1 = t_0_0;
        let mut t_1_i_minus_1 = t_1_0;
        let mut CW_i_vec: Vec<CWi> = Vec::new();
        for i in 0..n {
            let prng_out_0_i = G(&s_0_i_minus_1);
            let prng_out_1_i = G(&s_1_i_minus_1);
            let alpha_i_bn = alpha >> i & BigInt::one();
            let alpha_i = if alpha_i_bn > BigInt::zero() {
                true
            } else {
                false
            };
            let (s_0_Lose, s_0_Keep) = if alpha_i {
                (prng_out_0_i.s_i_L.clone(), prng_out_0_i.s_i_R.clone())
            } else {
                (prng_out_0_i.s_i_R.clone(), prng_out_0_i.s_i_L.clone())
            };
            let (s_1_Lose, s_1_Keep) = if alpha_i {
                (prng_out_1_i.s_i_L.clone(), prng_out_1_i.s_i_R.clone())
            } else {
                (prng_out_1_i.s_i_R.clone(), prng_out_1_i.s_i_L.clone())
            };

            let (t_0_Keep, t_1_Keep) = if alpha_i {
                (prng_out_0_i.t_i_R.clone(), prng_out_1_i.t_i_R.clone())
            } else {
                (prng_out_0_i.t_i_L.clone(), prng_out_1_i.t_i_L.clone())
            };
            let s_CW: [u8; LAMBDA_BYTES_LEN] = xor(&s_0_Lose, &s_1_Lose).try_into().unwrap();
            let t_CW_L = prng_out_0_i.t_i_L ^ prng_out_1_i.t_i_L ^ alpha_i_bn.to_bytes()[0] ^ 1u8;

            let t_CW_R = prng_out_0_i.t_i_R ^ prng_out_1_i.t_i_R ^ alpha_i_bn.to_bytes()[0];

            CW_i_vec.push(CWi {
                s_CW: s_CW.clone(),
                t_CW_L,
                t_CW_R,
            });
            let t_CW_Keep = if alpha_i { t_CW_R } else { t_CW_L };
            s_0_i_minus_1 = if t_0_i_minus_1 > 0 {
                xor(&s_0_Keep, &s_CW)
            } else {
                s_0_Keep.clone()
            };
            s_1_i_minus_1 = if t_1_i_minus_1 > 0 {
                xor(&s_1_Keep, &s_CW)
            } else {
                s_1_Keep.clone()
            };

            t_0_i_minus_1 = t_0_Keep ^ (t_0_i_minus_1 * t_CW_Keep);
            t_1_i_minus_1 = t_1_Keep ^ (t_1_i_minus_1 * t_CW_Keep);
        }

        let s_0_n = s_0_i_minus_1;
        let s_1_n = s_1_i_minus_1;

        let t_1_n = t_1_i_minus_1;

        let s_0_n_fe = convert(s_0_n);
        let s_1_n_fe = convert(s_1_n);

        let mut CW_n_plus_1 = s_1_n_fe + beta - &s_0_n_fe;

        CW_n_plus_1 = if t_1_n == 1u8 {
            let neg_cw_bn = Scalar::<Secp256k1>::group_order() - CW_n_plus_1.to_bigint();
            Scalar::from(&neg_cw_bn)
        } else {
            CW_n_plus_1
        };
        let key0 = Key {
            s_i_0: s_0_0.clone(),
            CW_n_plus_1: CWnp1(CW_n_plus_1.clone()),
            CWs: CW_i_vec.clone().try_into().unwrap(),
        };
        let key1 = Key {
            s_i_0: s_1_0.clone(),
            CW_n_plus_1: CWnp1(CW_n_plus_1.clone()),
            CWs: CW_i_vec.try_into().unwrap(),
        };
        (DPF(key0), DPF(key1))
    }

    pub fn eval(&self, b: &u8, x: &BigInt) -> Scalar<Secp256k1> {
        assert!(b <= &1u8);
        let mut s_i_minus_1 = self.0.s_i_0.to_vec();
        let mut t_i_minus_1 = b.clone();
        for i in 0..N {
            let stst_i_concat = concat_s_t_values(
                &self.0.CWs[i].s_CW.to_vec(),
                &self.0.CWs[i].t_CW_L,
                &self.0.CWs[i].s_CW.to_vec(),
                &self.0.CWs[i].t_CW_R,
            );
            let prng_out_i = G(&s_i_minus_1);
            let prng_i_concat = concat_s_t_values(
                &prng_out_i.s_i_L,
                &prng_out_i.t_i_L,
                &prng_out_i.s_i_R,
                &prng_out_i.t_i_R,
            );

            let tau_i_vec = if t_i_minus_1 == 0u8 {
                prng_i_concat
            } else {
                xor(&prng_i_concat, &stst_i_concat)
            };

            let tau_i = split_s_t_values(&tau_i_vec);
            let x_i_bn = x >> i & BigInt::one();
            let x_i = if x_i_bn > BigInt::zero() { true } else { false };
            s_i_minus_1 = if x_i { tau_i.s_i_R } else { tau_i.s_i_L };
            t_i_minus_1 = if x_i { tau_i.t_i_R } else { tau_i.t_i_L };
        }
        let s_n = s_i_minus_1.clone();
        let t_n = t_i_minus_1.clone();
        let t_n_bn = BigInt::from_bytes(&[t_n]);
        let t_n_cw_n1 = BigInt::mod_mul(
            &self.0.CW_n_plus_1.0.to_bigint(),
            &t_n_bn,
            &Scalar::<Secp256k1>::group_order(),
        );

        let mut output = BigInt::mod_add(
            &convert(s_n).to_bigint(),
            &t_n_cw_n1,
            &Scalar::<Secp256k1>::group_order(),
        );

        if b == &1u8 {
            output = Scalar::<Secp256k1>::group_order() - output;
        }
        Scalar::from(&output)
    }

    pub fn full_eval_N(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
        (0..N)
            .map(|i| self.eval(b, &BigInt::from(i as u32)))
            .collect()
    }

    pub fn full_eval_2N(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
        (0..2 * N)
            .map(|i| self.eval(b, &BigInt::from(i as u32)))
            .collect()
    }
}

fn G(prng_in: &Vec<u8>) -> PrngOut {
    let mut prng_in_with_set_bit = prng_in.clone();

    if prng_in_with_set_bit.len() < LAMBDA_BYTES_LEN {
        let mut template = vec![0; LAMBDA_BYTES_LEN - prng_in_with_set_bit.len()];
        template.extend_from_slice(&prng_in_with_set_bit);
        prng_in_with_set_bit = template;
    }
    prng_in_with_set_bit[0] = prng_in_with_set_bit[0] & 0xfe;
    let output = prng(&prng_in_with_set_bit);
    split_s_t_values(&output)
}

// figure 3 in the paper. takes s and converts to Fq
fn convert(s: Vec<u8>) -> Scalar<Secp256k1> {
    let fe_vec = prng(&s);
    let bn = BigInt::from_bytes(&fe_vec[..]);
    Scalar::from(&bn)
}

// takes  [s1||t1||s2||t2] and outputs s1,t1,s2,t2
fn split_s_t_values(input: &Vec<u8>) -> PrngOut {
    assert_eq!(input.len(), LAMBDA_BYTES_LEN * 2);

    let mut stst = input.clone();
    let mut s_R = stst.split_off(LAMBDA_BYTES_LEN);
    let mut s_L = stst;

    let t_R = s_R[0] & 1;
    let t_L = s_L[0] & 1;
    s_R[0] = s_R[0] & 0xfe;
    s_L[0] = s_L[0] & 0xfe;

    PrngOut {
        s_i_L: s_L,
        t_i_L: t_L,
        s_i_R: s_R,
        t_i_R: t_R,
    }
}
// takes s1,t1,s2,t2 and outputs [s1||t1||s2||t2]
fn concat_s_t_values(s1: &Vec<u8>, t1: &u8, s2: &Vec<u8>, t2: &u8) -> Vec<u8> {
    let mut s_L = s1.clone();
    s_L[0] = s_L[0] & 0xfe;
    s_L[0] = s_L[0] + t1;
    let mut s_R = s2.clone();
    s_R[0] = s_R[0] & 0xfe;
    s_R[0] = s_R[0] + t2;
    [s_L, s_R].concat()
}

mod tests {
    use crate::dpf::DPF;
    use crate::N;
    use curv::elliptic::curves::{Scalar, Secp256k1};
    use curv::BigInt;

    #[test]
    fn test_simple_dpf() {
        let alpha = BigInt::from(10);
        let beta = Scalar::random();
        let (dpf0, dpf1) = DPF::gen(&alpha, &beta);
        for x in 0..N {
            let x_bn = BigInt::from(x as u32);
            let fe0 = dpf0.eval(&0u8, &x_bn);
            let fe1 = dpf1.eval(&1u8, &x_bn);
            if x_bn == alpha {
                assert_eq!(fe0 + fe1, beta);
            } else {
                assert_eq!(
                    fe0.to_bigint(),
                    Scalar::<Secp256k1>::group_order() - fe1.to_bigint()
                );
            }
        }
    }
}
