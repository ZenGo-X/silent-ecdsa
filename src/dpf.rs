#![allow(non_snake_case)]
// implementing https://eprint.iacr.org/2018/707.pdf figure 1

use crate::utils::{prng, xor};
use crate::{LAMBDA, LAMBDA_BYTES_LEN, N};
use curv::arithmetic::{Converter, Modulo, One, Samplable, Zero};
use curv::elliptic::curves::secp256_k1::FE;
use curv::elliptic::curves::traits::ECScalar;
use curv::BigInt;

#[derive(Clone, Debug)]
pub struct PrngOut {
    pub s_i_L: Vec<u8>,
    pub t_i_L: u8,
    pub s_i_R: Vec<u8>,
    pub t_i_R: u8,
}

#[derive(Clone, Debug)]
pub struct CWi {
    pub s_CW: Vec<u8>,
    pub t_CW_L: u8,
    pub t_CW_R: u8,
}

pub struct CWnp1(FE);

pub struct Key {
    pub s_i_0: Vec<u8>,
    pub CW_n_plus_1: CWnp1,
    pub CWs: Vec<CWi>,
}

pub struct DPF;

impl DPF {
    pub fn gen(alpha: &BigInt, beta: &FE) -> (Key, Key) {
        let mut s_0_0 = BigInt::sample(LAMBDA * 2).to_bytes()[0..LAMBDA_BYTES_LEN].to_vec();
        let mut s_1_0 = BigInt::sample(LAMBDA * 2).to_bytes()[0..LAMBDA_BYTES_LEN].to_vec();

        let t_0_0: u8 = 0;
        let t_1_0: u8 = 1;
        s_0_0[0] = s_0_0[0] & 0xfe;
        s_1_0[0] = s_1_0[0] & 0xfe;
        let n = N;
        let mut s_0_i_minus_1 = s_0_0.clone();
        let mut s_1_i_minus_1 = s_1_0.clone();
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
            let s_CW = xor(&s_0_Lose, &s_1_Lose);
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

        let mut CW_n_plus_1 = s_1_n_fe + beta.sub(&s_0_n_fe.get_element());

        CW_n_plus_1 = if t_1_n == 1u8 {
            let neg_cw_bn = FE::q() - CW_n_plus_1.to_big_int();
            ECScalar::from(&neg_cw_bn)
        } else {
            CW_n_plus_1
        };
        let key0 = Key {
            s_i_0: s_0_0.clone(),
            CW_n_plus_1: CWnp1(CW_n_plus_1.clone()),
            CWs: CW_i_vec.clone(),
        };
        let key1 = Key {
            s_i_0: s_1_0.clone(),
            CW_n_plus_1: CWnp1(CW_n_plus_1.clone()),
            CWs: CW_i_vec.clone(),
        };
        (key0, key1)
    }

    pub fn eval(b: &u8, key_b: &Key, x: &BigInt) -> FE {
        assert!(b <= &1u8);
        let mut s_i_minus_1 = key_b.s_i_0.clone();
        let mut t_i_minus_1 = b.clone();
        for i in 0..N {
            let stst_i_concat = concat_s_t_values(
                &key_b.CWs[i].s_CW,
                &key_b.CWs[i].t_CW_L,
                &key_b.CWs[i].s_CW,
                &key_b.CWs[i].t_CW_R,
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
        let t_n_cw_n1 = BigInt::mod_mul(&key_b.CW_n_plus_1.0.to_big_int(), &t_n_bn, &FE::q());

        let mut output = BigInt::mod_add(&convert(s_n).to_big_int(), &t_n_cw_n1, &FE::q());

        if b == &1u8 {
            output = FE::q() - output;
        }
        ECScalar::from(&output)
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
fn convert(s: Vec<u8>) -> FE {
    let fe_vec = prng(&s);
    let bn = BigInt::from_bytes(&fe_vec[..]);
    ECScalar::from(&bn)
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
    use curv::elliptic::curves::secp256_k1::FE;
    use curv::elliptic::curves::traits::ECScalar;
    use curv::BigInt;

    #[test]
    fn test_simple_dpf() {
        let alpha = BigInt::from(10);
        let beta = ECScalar::new_random();
        let (key0, key1) = DPF::gen(&alpha, &beta);
        for x in 0..N {
            let x_bn = BigInt::from(x as u32);
            let fe0 = DPF::eval(&0u8, &key0, &x_bn);
            let fe1 = DPF::eval(&1u8, &key1, &x_bn);
            if x_bn == alpha {
                assert_eq!(fe0 + fe1, beta);
            } else {
                assert_eq!(fe0.to_big_int(), FE::q() - fe1.to_big_int());
            }
        }
    }
}
