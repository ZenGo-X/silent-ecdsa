#![allow(non_snake_case)]
// implementing https://eprint.iacr.org/2018/707.pdf figure 1

use crate::utils::{prng, xor_slices, xor_slices_in_place};
use crate::{LAMBDA, LAMBDA_BYTES_LEN, N};
use curv::arithmetic::{BitManipulation, Converter, Modulo, One, Samplable, Zero};
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;
use serde::{Deserialize, Serialize};
use std::convert::TryInto;

#[derive(Debug)]
pub struct PrngOut {
    pub s_i_L: [u8; LAMBDA_BYTES_LEN],
    pub t_i_L: u8,
    pub s_i_R: [u8; LAMBDA_BYTES_LEN],
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
    pub CWs: Vec<CWi>,
    pub tree_depth: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DPF(pub(crate) Key);

impl DPF {
    pub fn init() -> Self {
        let key = Key {
            s_i_0: [0; LAMBDA_BYTES_LEN],
            CW_n_plus_1: CWnp1(Scalar::zero()),
            CWs: vec![CWi::init(); u32::BITS as usize],
            tree_depth: 0,
        };
        DPF(key)
    }

    pub fn gen(alpha: &BigInt, beta: &Scalar<Secp256k1>) -> (Self, Self) {
        let s_0_0_vec: [u8; LAMBDA_BYTES_LEN] = BigInt::sample(LAMBDA * 2).to_bytes()
            [0..LAMBDA_BYTES_LEN]
            .try_into()
            .expect("Shouldn't happen, lengths are correct!");
        let s_1_0_vec: [u8; LAMBDA_BYTES_LEN] = BigInt::sample(LAMBDA * 2).to_bytes()
            [0..LAMBDA_BYTES_LEN]
            .try_into()
            .expect("Shouldn't happen, length are correct!");
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
        let n = BigInt::from(2 * N as u32).bit_length();
        let mut s_0_i_minus_1 = s_0_0_vec;
        let mut s_1_i_minus_1 = s_1_0_vec;
        let mut t_0_i_minus_1 = t_0_0;
        let mut t_1_i_minus_1 = t_1_0;
        let mut CW_i_vec: Vec<CWi> = Vec::new();
        let bit_size = n;
        for i in 0..bit_size {
            let prng_out_0_i = G(s_0_i_minus_1);
            let prng_out_1_i = G(s_1_i_minus_1);
            let alpha_i_bn = alpha >> (bit_size - 1 - i) & BigInt::one();
            let alpha_i = if alpha_i_bn > BigInt::zero() {
                true
            } else {
                false
            };
            let (s_0_Lose, s_0_Keep) = if alpha_i {
                (prng_out_0_i.s_i_L, prng_out_0_i.s_i_R)
            } else {
                (prng_out_0_i.s_i_R, prng_out_0_i.s_i_L)
            };
            let (s_1_Lose, s_1_Keep) = if alpha_i {
                (prng_out_1_i.s_i_L, prng_out_1_i.s_i_R)
            } else {
                (prng_out_1_i.s_i_R, prng_out_1_i.s_i_L)
            };

            let (t_0_Keep, t_1_Keep) = if alpha_i {
                (prng_out_0_i.t_i_R, prng_out_1_i.t_i_R)
            } else {
                (prng_out_0_i.t_i_L, prng_out_1_i.t_i_L)
            };
            let s_CW: [u8; LAMBDA_BYTES_LEN];

            s_CW = xor_slices(&s_0_Lose, &s_1_Lose);
            let t_CW_L = prng_out_0_i.t_i_L ^ prng_out_1_i.t_i_L ^ alpha_i_bn.to_bytes()[0] ^ 1u8;

            let t_CW_R = prng_out_0_i.t_i_R ^ prng_out_1_i.t_i_R ^ alpha_i_bn.to_bytes()[0];

            CW_i_vec.push(CWi {
                s_CW: s_CW.clone(),
                t_CW_L,
                t_CW_R,
            });
            let t_CW_Keep = if alpha_i { t_CW_R } else { t_CW_L };
            // s_0_i_minus_1 =
            s_0_i_minus_1 = if t_0_i_minus_1 > 0 {
                xor_slices(&s_0_Keep, &s_CW)
            } else {
                s_0_Keep
            };
            s_1_i_minus_1 = if t_1_i_minus_1 > 0 {
                xor_slices(&s_1_Keep, &s_CW)
            } else {
                s_1_Keep
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
            tree_depth: n,
        };
        let key1 = Key {
            s_i_0: s_1_0.clone(),
            CW_n_plus_1: CWnp1(CW_n_plus_1.clone()),
            CWs: CW_i_vec.try_into().unwrap(),
            tree_depth: n,
        };
        (DPF(key0), DPF(key1))
    }

    pub fn eval(&self, b: u8, x: &BigInt) -> Scalar<Secp256k1> {
        assert!(b <= 1u8);
        let mut s_i_minus_1 = self.0.s_i_0;
        let mut t_i_minus_1 = b;
        let input_bit_length = BigInt::from(2 * N as u32).bit_length();
        for i in 0..input_bit_length {
            let stst_i_concat = concat_s_t_values(
                self.0.CWs[i].s_CW,
                self.0.CWs[i].t_CW_L,
                self.0.CWs[i].s_CW,
                self.0.CWs[i].t_CW_R,
            );
            let prng_out_i = G(s_i_minus_1);
            let prng_i_concat = concat_s_t_values(
                prng_out_i.s_i_L,
                prng_out_i.t_i_L,
                prng_out_i.s_i_R,
                prng_out_i.t_i_R,
            );

            let tau_i_vec = if t_i_minus_1 == 0u8 {
                prng_i_concat
            } else {
                xor_slices(&prng_i_concat, &stst_i_concat)
            };

            let tau_i = split_s_t_values(tau_i_vec);
            let x_i_bn = x >> (input_bit_length - 1 - i) & BigInt::one();
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

        if b == 1u8 {
            output = Scalar::<Secp256k1>::group_order() - output;
        }
        Scalar::from(&output)
    }

    // This is an optimization that instead of iteratively evaluate the DpfKey
    // in each address, we save time by computing the DPF over a range.
    // This optimization is not explicitly stated in the original paper.
    // But figure 4 in the paper gives part of the general sense
    // behind the optimization.

    fn full_eval_optimized_recursive(
        &self,
        id: &bool,
        all_ones_path_bits: &Vec<bool>,
        upper_bound_path_bits: &Vec<bool>,
        control_bit: bool,
        seed: [u8; LAMBDA_BYTES_LEN],
        current_depth: usize,
    ) -> Vec<Scalar<Secp256k1>> {
        if current_depth == self.0.tree_depth {
            let mut output_val = convert(seed);
            if control_bit {
                output_val = output_val + &self.0.CW_n_plus_1.0;
            }
            if *id {
                output_val = -output_val;
            }
            return vec![output_val];
        }
        let expanded_seed = G(seed);
        let (mut seed_left, mut control_bit_left, mut seed_right, mut control_bit_right) = (
            expanded_seed.s_i_L,
            expanded_seed.t_i_L,
            expanded_seed.s_i_R,
            expanded_seed.t_i_R,
        );
        if control_bit {
            xor_slices_in_place(&mut seed_left, &self.0.CWs[current_depth].s_CW);
            xor_slices_in_place(&mut seed_right, &self.0.CWs[current_depth].s_CW);
            control_bit_left ^= self.0.CWs[current_depth].t_CW_L;
            control_bit_right ^= self.0.CWs[current_depth].t_CW_R;
        }

        if !upper_bound_path_bits[current_depth] {
            self.full_eval_optimized_recursive(
                id,
                all_ones_path_bits,
                upper_bound_path_bits,
                control_bit_left != 0,
                seed_left,
                current_depth + 1,
            )
        } else {
            // We will set here all ones to that the left subtree will be evaluated fully.
            let mut left_eval = self.full_eval_optimized_recursive(
                id,
                all_ones_path_bits,
                all_ones_path_bits,
                control_bit_left != 0,
                seed_left,
                current_depth + 1,
            );
            let mut right_eval = self.full_eval_optimized_recursive(
                id,
                all_ones_path_bits,
                upper_bound_path_bits,
                control_bit_right != 0,
                seed_right,
                current_depth + 1,
            );
            left_eval.append(&mut right_eval);
            left_eval
        }
    }
    // Evals for all i such that: 0 <= i < upper_bound
    pub fn full_eval_optimized(&self, b: &bool, upper_bound: &BigInt) -> Vec<Scalar<Secp256k1>> {
        let upper_bound_path_bits: Vec<bool> = (0..self.0.tree_depth)
            .rev()
            .map(|idx| upper_bound.test_bit(idx))
            .collect();
        let all_ones_path_bits: Vec<bool> = upper_bound_path_bits.iter().map(|_| true).collect();
        let seed = self.0.s_i_0.clone();
        let control_bit = *b;
        self.full_eval_optimized_recursive(
            b,
            &all_ones_path_bits,
            &upper_bound_path_bits,
            control_bit,
            seed,
            0,
        )
    }
    pub fn full_eval_N(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
        self.full_eval_optimized(&(*b != 0), &BigInt::from((N - 1) as u64))
    }
    // pub fn full_eval_N_naive(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
    //     (0..N)
    //         .map(|i| self.eval(b, &BigInt::from(i as u32)))
    //         .collect()
    // }

    pub fn full_eval_2N(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
        self.full_eval_optimized(&(*b != 0), &BigInt::from((2 * N - 1) as u64))
    }

    // pub fn full_eval_2N_naive(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
    //     (0..2 * N)
    //         .map(|i| self.eval(b, &BigInt::from(i as u32)))
    //         .collect()
    // }
}

fn G(mut prng_in: [u8; LAMBDA_BYTES_LEN]) -> PrngOut {
    // if prng_in_with_set_bit.len() < LAMBDA_BYTES_LEN {
    //     let mut template = vec![0; LAMBDA_BYTES_LEN - prng_in_with_set_bit.len()];
    //     template.extend_from_slice(&prng_in_with_set_bit);
    //     prng_in_with_set_bit = template;
    // }
    // let prev_byte = prng_in[0];
    prng_in[0] &= 0xfe;
    let output = prng(prng_in);
    // prng_in[0] = prev_byte;
    split_s_t_values(output)
}

// figure 3 in the paper. takes s and converts to Fq
fn convert(s: [u8; LAMBDA_BYTES_LEN]) -> Scalar<Secp256k1> {
    let fe_vec = prng(s);
    let bn = BigInt::from_bytes(&fe_vec[..]);
    Scalar::from(&bn)
}

// takes  [s1||t1||s2||t2] and outputs s1,t1,s2,t2
fn split_s_t_values(input: [u8; LAMBDA_BYTES_LEN * 2]) -> PrngOut {
    PrngOut {
        s_i_L: input[..LAMBDA_BYTES_LEN]
            .try_into()
            .expect("Length is correct, shouldn't happen!"),
        t_i_L: input[0] & 1,
        s_i_R: input[LAMBDA_BYTES_LEN..LAMBDA_BYTES_LEN * 2]
            .try_into()
            .expect("Length is correct, shouldn't happen!"),
        t_i_R: input[LAMBDA_BYTES_LEN] & 1,
    }

    // let s_R = &mut input[..LAMBDA_BYTES_LEN];
    // let s_L = &mut input[LAMBDA_BYTES_LEN..LAMBDA_BYTES_LEN * 2];

    // let t_R = s_R[0] & 1;
    // let t_L = s_L[0] & 1;
    // s_R[0] = s_R[0] & 0xfe;
    // s_L[0] = s_L[0] & 0xfe;

    // PrngOut {
    //     s_i_L: s_L,
    //     t_i_L: t_L,
    //     s_i_R: s_R,
    //     t_i_R: t_R,
    // }
}
// takes s1,t1,s2,t2 and outputs [s1||t1||s2||t2]
fn concat_s_t_values(
    s1: [u8; LAMBDA_BYTES_LEN],
    t1: u8,
    s2: [u8; LAMBDA_BYTES_LEN],
    t2: u8,
) -> [u8; LAMBDA_BYTES_LEN * 2] {
    let mut out: [u8; LAMBDA_BYTES_LEN * 2] = [s1, s2]
        .concat()
        .try_into()
        .expect("Should be of correct length here!");
    out[0] = out[0] & 0xfe + t1;
    out[LAMBDA_BYTES_LEN] = out[LAMBDA_BYTES_LEN] & 0xfe + t2;
    out
}

#[cfg(test)]
mod tests {
    use crate::dpf::DPF;
    use crate::N;
    use curv::elliptic::curves::{Scalar, Secp256k1};
    use curv::BigInt;

    #[test]
    fn test_simple_dpf() {
        let alpha = 10;
        let beta = Scalar::random();
        let (dpf0, dpf1) = DPF::gen(&BigInt::from(alpha as u64), &beta);
        let full_eval_0 = dpf0.full_eval_2N(&0u8);
        let full_eval_1 = dpf1.full_eval_2N(&1u8);
        for x in 0..2 * N {
            let fe0 = &full_eval_0[x];
            let fe1 = &full_eval_1[x];
            if x == alpha {
                assert_eq!(fe0 + fe1, beta);
            } else {
                assert_eq!(
                    fe0.to_bigint(),
                    Scalar::<Secp256k1>::group_order() - fe1.to_bigint()
                );
            }
        }
    }
    // #[test]
    // fn test_full_eval() {
    //     let alpha = BigInt::from(10);
    //     let beta = Scalar::from(1);
    //     let (dpf0, _) = DPF::gen(&alpha, &beta);
    //     assert_eq!(dpf0.full_eval_2N_naive(&0u8), dpf0.full_eval_2N(&0u8));
    // }
}
