use crate::dspf::DSPF;
use crate::{c, n, t, N};
use curv::arithmetic::{Samplable, Zero};
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;
use std::convert::TryInto;
use std::mem;
use std::ptr;

#[derive(Debug)]
struct LongTermKey {
    pub alpha_i: Scalar<Secp256k1>,
    pub sk_i: Scalar<Secp256k1>,
    pub w_i: [BigInt; c * t],
    pub eta_i: [BigInt; c * t],
    pub beta_i: [Scalar<Secp256k1>; c * t],
    pub gamma_i: [Scalar<Secp256k1>; c * t],
    pub u_i_0: [DSPF; c * (n - 1)],
    pub u_i_1: [DSPF; c * (n - 1)],
    pub v_i_0: [DSPF; c * (n - 1)],
    pub v_i_1: [DSPF; c * (n - 1)],
    pub c_i_0: [DSPF; c * c * (n - 1)],
    pub c_i_1: [DSPF; c * c * (n - 1)],
}

impl LongTermKey {
    pub fn init() -> Self {
        LongTermKey {
            alpha_i: Scalar::zero(),
            sk_i: Scalar::zero(),
            w_i: unsafe { make_array!(c * t, BigInt::zero()) },
            eta_i: unsafe { make_array!(c * t, BigInt::zero()) },
            beta_i: unsafe { make_array!(c * t, Scalar::zero()) },
            gamma_i: unsafe { make_array!(c * t, Scalar::zero()) },
            u_i_0: unsafe { make_array!(c * (n - 1), DSPF::init(t*c * (n - 1))) },
            u_i_1: unsafe { make_array!(c * (n - 1), DSPF::init(t*c * (n - 1))) },
            v_i_0: unsafe { make_array!(c * (n - 1), DSPF::init(t*c * (n - 1))) },
            v_i_1: unsafe { make_array!(c * (n - 1), DSPF::init(t*c * (n - 1))) },
            c_i_0: unsafe { make_array!(c * c * (n - 1), DSPF::init(t*t*c * c * (n - 1))) },
            c_i_1: unsafe { make_array!(c * c * (n - 1), DSPF::init(t*t*c * c * (n - 1))) },
        }
    }

    pub fn trusted_key_gen() -> [Self; n] {
        let mut long_term_keys = unsafe { make_array!(n, LongTermKey::init()) };
        for i in 0..n {
            long_term_keys[i].alpha_i = Scalar::<Secp256k1>::random();
            long_term_keys[i].sk_i = Scalar::<Secp256k1>::random();
        }

        for i in 0..n {
            for r in 0..c {
                let beta_i_vec = pick_Fq_t();
                let gamma_i_vec = pick_Fq_t();
                let w_i_vec = pick_N_t().to_vec();
                let eta_i_vec = pick_N_t().to_vec();
                for z in 0..t {
                    long_term_keys[i].beta_i[t * r + z] = beta_i_vec[z].clone();
                    long_term_keys[i].gamma_i[t * r + z] = gamma_i_vec[z].clone();
                    long_term_keys[i].w_i[t * r + z] = w_i_vec[z].clone();
                    long_term_keys[i].eta_i[t * r + z] = eta_i_vec[z].clone();
                }
            }
        }

        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                for r in 0..c {
                    let alpha_j_beta_i_r: Vec<_> = (0..t)
                        .map(|z| &long_term_keys[j].alpha_i * &long_term_keys[i].beta_i[t * r + z])
                        .collect();
                    let w_i_r = &long_term_keys[i].w_i[t * r..t * r + t];
                    let (u_ij_0, u_ij_1) = DSPF::gen(&w_i_r[..], &alpha_j_beta_i_r[..], false);
                    let sk_j_gamma_i_r: Vec<_> = (0..t)
                        .map(|z| &long_term_keys[j].sk_i * &long_term_keys[i].gamma_i[t * r + z])
                        .collect();
                    let eta_i_r = &long_term_keys[i].eta_i[t * r..t * r + t];
                    let (v_ij_0, v_ij_1) = DSPF::gen(&eta_i_r[..], &sk_j_gamma_i_r[..], false);
                    let index = if j < i { j } else { j - 1 };
                    long_term_keys[i].u_i_0[c * index + r] = u_ij_0;
                    long_term_keys[i].u_i_1[c * index + r] = u_ij_1;
                    long_term_keys[i].v_i_0[c * index + r] = v_ij_0;
                    long_term_keys[i].v_i_1[c * index + r] = v_ij_1;
                }
            }
        }
        for i in 0..n {
            for j in 0..n {
                for r in 0..c {
                    for s in 0..c {
                        if j == i {
                            continue;
                        }
                        let w_i_r_outer_sum_eta_j_s = outer_sum(
                            &long_term_keys[i].w_i[t * r..t * r + t],
                            &long_term_keys[j].eta_i[t * s..t * s + t],
                        );
                        let beta_i_r_outer_mul_gamma_j_s = outer_product(
                            &long_term_keys[i].beta_i[t * r..t * r + t],
                            &long_term_keys[j].gamma_i[t * s..t * s + t],
                        );
                        let (c_ij_rs_0, c_ij_rs_1) = DSPF::gen(
                            &w_i_r_outer_sum_eta_j_s[..],
                            &beta_i_r_outer_mul_gamma_j_s[..],
                            true,
                        );
                        let index = if j < i { j } else { j - 1 };
                        long_term_keys[i].c_i_0[r * index + s] = c_ij_rs_0;
                        long_term_keys[i].c_i_1[r * index + s] = c_ij_rs_1;
                    }
                }
            }
        }

        return long_term_keys;
    }
}

fn pick_N_t() -> [BigInt; t] {
    (0..t)
        .map(|_| BigInt::sample_below(&BigInt::from(N as u32)))
        .collect::<Vec<BigInt>>()
        .try_into()
        .unwrap()
}

fn pick_Fq_t() -> [Scalar<Secp256k1>; t] {
    (0..t)
        .map(|_| Scalar::<Secp256k1>::random())
        .collect::<Vec<Scalar<Secp256k1>>>()
        .try_into()
        .unwrap()
}

// We recall that if u and v have dimensions m and l, the outer product and the outer sum
// are the ml-dimensional vectors whose (im + j)-th entry is ui · vj and ui + vj respectively.
fn outer_product(u: &[Scalar<Secp256k1>], v: &[Scalar<Secp256k1>]) -> Vec<Scalar<Secp256k1>> {
    let mut output: [Scalar<Secp256k1>; t * t] = unsafe { make_array!(t * t, Scalar::zero()) };
    for i in 0..t {
        for j in 0..t {
            output[i * t + j] = &u[i] * &v[j];
        }
    }
    output.to_vec()
}

fn outer_sum(u: &[BigInt], v: &[BigInt]) -> Vec<BigInt> {
    let mut output: [BigInt; t * t] = unsafe { make_array!(t * t, BigInt::zero()) };
    for i in 0..t {
        for j in 0..t {
            output[i * t + j] = &u[i] + &v[j];
        }
    }
    output.to_vec()
}

mod tests {
    use crate::keygen::LongTermKey;

    #[test]
    fn test_run_keygen() {
        let key = LongTermKey::trusted_key_gen();
       // println!("key : {:?}", key);
        println!("c_i_0 : {:?}", key[0].c_i_0.clone());
        println!("c_i_1 : {:?}", key[0].c_i_1.clone());

        println!("{:?}", key[0].u_i_0.len());
        println!("{:?}", key[0].u_i_0[0].key.len());
        println!("{:?}", key[0].u_i_0[0].key[0].0.CWs.len());
        assert!(false);

    }
}
