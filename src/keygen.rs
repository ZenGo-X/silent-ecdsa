use crate::dspf::DSPF;
use crate::poly::{poly_add_f, poly_mod, poly_mul_f, poly_mul_f_naive};
use crate::utils::array_from_fn;
use crate::{c, n, t, N};
use curv::arithmetic::{BasicOps, Converter, Modulo, One, Samplable, Zero};
use curv::cryptographic_primitives::secret_sharing::{ffts::PowerIterator, Polynomial};
use curv::elliptic::curves::{ECScalar, Point, Scalar, Secp256k1};
use curv::BigInt;
use rand::{random, Rng};
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use std::convert::TryFrom;
use std::convert::TryInto;
use std::fs;
use std::ops::{Add, Mul, Neg};

#[serde_as]
#[derive(Debug, Serialize, Deserialize)]
pub struct LongTermKey {
    pub alpha_i: Scalar<Secp256k1>,
    pub sk_i: Scalar<Secp256k1>,
    #[serde_as(as = "[[_; t]; c]")]
    pub w_i: [[u32; t]; c],
    #[serde_as(as = "[[_; t]; c]")]
    pub eta_i: [[u32; t]; c],
    #[serde_as(as = "[[_; t]; c]")]
    pub beta_i: [[Scalar<Secp256k1>; t]; c],
    #[serde_as(as = "[[_; t]; c]")]
    pub gamma_i: [[Scalar<Secp256k1>; t]; c],
    pub u_i_0: [[DSPF; (n - 1)]; c],
    pub u_i_1: [[DSPF; (n - 1)]; c],
    pub v_i_0: [[DSPF; (n - 1)]; c],
    pub v_i_1: [[DSPF; (n - 1)]; c],
    #[serde_as(as = "[[_; (n - 1)]; c * c]")]
    pub c_i_0: [[DSPF; (n - 1)]; c * c],
    #[serde_as(as = "[[_; (n - 1)]; c * c]")]
    pub c_i_1: [[DSPF; (n - 1)]; c * c],
    pub pk: Point<Secp256k1>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tuple {
    pub x_i: Vec<Scalar<Secp256k1>>,
    pub y_i: Vec<Scalar<Secp256k1>>,
    pub z_i: Vec<Scalar<Secp256k1>>,
    pub d_i: Vec<Scalar<Secp256k1>>,
    pub M_i_j: [Vec<Scalar<Secp256k1>>; n - 1],
    pub K_j_i: [Vec<Scalar<Secp256k1>>; n - 1],
}

impl LongTermKey {
    fn init() -> Self {
        LongTermKey {
            alpha_i: Scalar::zero(),
            sk_i: Scalar::zero(),
            w_i: array_from_fn(|_| array_from_fn(|_| 0u32)),
            eta_i: array_from_fn(|_| array_from_fn(|_| 0u32)),
            beta_i: array_from_fn(|_| array_from_fn(|_| Scalar::zero())),
            gamma_i: array_from_fn(|_| array_from_fn(|_| Scalar::zero())),
            u_i_0: array_from_fn(|_| array_from_fn(|_| DSPF::init(t * (n - 1)))),
            u_i_1: array_from_fn(|_| array_from_fn(|_| DSPF::init(t * (n - 1)))),
            v_i_0: array_from_fn(|_| array_from_fn(|_| DSPF::init(t * (n - 1)))),
            v_i_1: array_from_fn(|_| array_from_fn(|_| DSPF::init(t * (n - 1)))),
            c_i_0: array_from_fn(|_| array_from_fn(|_| DSPF::init(t * t * (n - 1)))),
            c_i_1: array_from_fn(|_| array_from_fn(|_| DSPF::init(t * t * (n - 1)))),
            pk: Point::<Secp256k1>::zero(),
        }
    }

    // todo: compute from seed, should be a random oracle
    pub fn sample_a() -> Vec<Vec<Scalar<Secp256k1>>> {
        let mut a = vec![vec![Scalar::<Secp256k1>::zero(); N]; c];

        for i in 0..c - 1 {
            a[i] = pick_R();
        }
        a[c - 1][0] = Scalar::from_bigint(&BigInt::one());
        a
    }

    pub fn trusted_key_gen() -> [Self; n] {
        let mut long_term_keys = array_from_fn(|_| LongTermKey::init());
        for i in 0..n {
            long_term_keys[i].alpha_i = Scalar::<Secp256k1>::random();
            long_term_keys[i].sk_i = Scalar::<Secp256k1>::random();
            //  long_term_keys[i].sk_i = Scalar::<Secp256k1>::from_bigint(&BigInt::one());
        }
        let pk = long_term_keys
            .iter()
            .fold(Point::<Secp256k1>::zero(), |acc, x| {
                acc + &x.sk_i * Point::<Secp256k1>::generator()
            });
        for i in 0..n {
            long_term_keys[i].pk = pk.clone();
        }
        for i in 0..n {
            for r in 0..c {
                let beta_i = pick_Fq_t();
                let gamma_i = pick_Fq_t();
                let w_i = pick_N_t();
                let eta_i = pick_N_t();
                long_term_keys[i].beta_i[r] = beta_i;
                long_term_keys[i].gamma_i[r] = gamma_i;
                long_term_keys[i].w_i[r] = w_i;
                long_term_keys[i].eta_i[r] = eta_i;
            }
        }

        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                for r in 0..c {
                    let alpha_j_beta_i_r: Vec<_> = (0..t)
                        .map(|z| &long_term_keys[j].alpha_i * &long_term_keys[i].beta_i[r][z])
                        .collect();
                    let w_i_r = &long_term_keys[i].w_i[r];
                    let (u_ij_0, u_ij_1) = DSPF::gen(&w_i_r[..], &alpha_j_beta_i_r[..]);
                    let sk_j_gamma_i_r: Vec<_> = (0..t)
                        .map(|z| &long_term_keys[j].sk_i * &long_term_keys[i].gamma_i[r][z])
                        .collect();
                    let eta_i_r = &long_term_keys[i].eta_i[r];
                    let (v_ij_0, v_ij_1) = DSPF::gen(&eta_i_r[..], &sk_j_gamma_i_r[..]);
                    let index_j = if j < i { j } else { j - 1 };
                    let index_i = if i < j { i } else { i - 1 };
                    long_term_keys[i].u_i_0[r][index_j] = u_ij_0;
                    long_term_keys[j].u_i_1[r][index_i] = u_ij_1;
                    long_term_keys[i].v_i_0[r][index_j] = v_ij_0;
                    long_term_keys[j].v_i_1[r][index_i] = v_ij_1;
                }
            }
        }

        for i in 0..n {
            for j in 0..n {
                if j == i {
                    continue;
                }
                for r in 0..c {
                    for s in 0..c {
                        let w_i_r_outer_sum_eta_j_s =
                            outer_sum(&long_term_keys[i].w_i[r], &long_term_keys[j].eta_i[s]);

                        let beta_i_r_outer_mul_gamma_j_s = outer_product_t(
                            &long_term_keys[i].beta_i[r],
                            &long_term_keys[j].gamma_i[s],
                        );
                        let (c_ij_rs_0, c_ij_rs_1) = DSPF::gen(
                            &w_i_r_outer_sum_eta_j_s[..],
                            &beta_i_r_outer_mul_gamma_j_s[..],
                        );
                        // let index = if j < i { j } else { j - 1 };
                        let index_j = if j < i { j } else { j - 1 };
                        let index_i = if i < j { i } else { i - 1 };
                        long_term_keys[i].c_i_0[r * c + s][index_j] = c_ij_rs_0;
                        long_term_keys[j].c_i_1[r * c + s][index_i] = c_ij_rs_1;
                    }
                }
            }
        }

        // write down keys to files
        for id in 0..n {
            let long_term_key_json = serde_json::to_string(&long_term_keys[id]).unwrap();
            fs::write(
                "long_term_key_".to_string().add(&id.to_string()),
                long_term_key_json,
            )
            .expect("Unable to save !");
        }

        return long_term_keys;
    }

    pub fn get_tuple(
        &self,
        a: Vec<Vec<Scalar<Secp256k1>>>,
        f_x: &Polynomial<Secp256k1>,
        id: usize,
    ) {
        let mut M_i_j = array_from_fn(|_| vec![Scalar::zero(); N]);
        let mut K_j_i = array_from_fn(|_| vec![Scalar::zero(); N]);
        let mut d_i = vec![Scalar::zero(); N];
        let mut x_i = vec![Scalar::zero(); N];
        let mut y_i = vec![Scalar::zero(); N];
        let mut z_i = vec![Scalar::zero(); N];

        for r in 0..c {
            let u_i_r = set_poly(&self.beta_i[r], &self.w_i[r]);
            let v_i_r = set_poly(&self.gamma_i[r], &self.eta_i[r]);
            let mut v_i_r_tilde: Vec<_> = (0..N).map(|k| &self.sk_i * &v_i_r[k]).collect();

            let a_r_mul_u_i_r = poly_mod(&poly_mul_f(&a[r], &u_i_r[..]), f_x.coefficients()).1;

            let a_r_mul_v_i_r = poly_mod(&poly_mul_f(&a[r], &v_i_r[..]), f_x.coefficients()).1;

            x_i = poly_add_f(&x_i[..], &a_r_mul_u_i_r[..]).try_into().unwrap();

            y_i = poly_add_f(&y_i[..], &a_r_mul_v_i_r[..]).try_into().unwrap();

            for j in 0..n - 1 {
                let M_i_j_r_tilde = &self.u_i_0[r][j].full_eval_N(0u8);
                let K_i_j_r_tilde = &self.u_i_1[r][j].full_eval_N(1u8);

                v_i_r_tilde = poly_add_f(&v_i_r_tilde, &self.v_i_0[r][j].full_eval_N(0u8));
                v_i_r_tilde = poly_mod(&v_i_r_tilde, f_x.coefficients()).1;
                v_i_r_tilde = poly_add_f(&v_i_r_tilde, &self.v_i_1[r][j].full_eval_N(1u8));
                v_i_r_tilde = poly_mod(&v_i_r_tilde, f_x.coefficients()).1;

                let a_r_mul_M_i_j_r_tilde =
                    poly_mod(&poly_mul_f(&a[r], &M_i_j_r_tilde[..]), f_x.coefficients()).1;
                let a_r_mul_K_i_j_r_tilde =
                    poly_mod(&poly_mul_f(&a[r], &K_i_j_r_tilde[..]), f_x.coefficients()).1;

                M_i_j[j] = poly_add_f(&M_i_j[j], &a_r_mul_M_i_j_r_tilde[..])
                    .try_into()
                    .unwrap();

                K_j_i[j] = poly_add_f(&K_j_i[j], &a_r_mul_K_i_j_r_tilde[..])
                    .try_into()
                    .unwrap();
            }

            let a_r_mul_v_i_r_tilde =
                poly_mod(&poly_mul_f(&a[r], &v_i_r_tilde[..]), f_x.coefficients()).1;

            d_i = poly_mod(
                &poly_add_f(&d_i[..], &a_r_mul_v_i_r_tilde[..]),
                f_x.coefficients(),
            )
            .1
            .try_into()
            .unwrap();
        }

        // -K_j_i
        for j in 0..n - 1 {
            for f in 0..N {
                K_j_i[j][f] = Scalar::from_bigint(
                    &(Scalar::<Secp256k1>::group_order() - K_j_i[j][f].to_bigint()),
                );
            }
        }

        let a_mul_a = outer_product_c(&a, &a);
        for r in 0..c {
            for s in 0..c {
                let u_i_r = set_poly(&self.beta_i[r], &self.w_i[r]);
                let v_i_s = set_poly(&self.gamma_i[s], &self.eta_i[s]);
                let mut w_i_r_s =
                    poly_mod(&poly_mul_f(&u_i_r[..], &v_i_s[..]), f_x.coefficients()).1;
                for j in 0..n - 1 {
                    w_i_r_s = poly_mod(
                        &poly_add_f(
                            &w_i_r_s[..],
                            &poly_mod(
                                &self.c_i_0[r * c + s][j].full_eval_2N(0u8),
                                f_x.coefficients(),
                            )
                            .1[..],
                        ),
                        f_x.coefficients(),
                    )
                    .1;
                    w_i_r_s = poly_mod(
                        &poly_add_f(
                            &w_i_r_s[..],
                            &poly_mod(
                                &self.c_i_1[r * c + s][j].full_eval_2N(1u8),
                                f_x.coefficients(),
                            )
                            .1[..],
                        ),
                        f_x.coefficients(),
                    )
                    .1;
                }
                let mut a_r_mul_a_s_mod = poly_mod(&a_mul_a[r * c + s], &f_x.coefficients()).1;
                for _ in 0..N - a_r_mul_a_s_mod.len() {
                    a_r_mul_a_s_mod.push(Scalar::<Secp256k1>::zero());
                }
                let a_r_a_s_mul_w_i_r_s = poly_mod(
                    &poly_mul_f(&w_i_r_s[..], &a_r_mul_a_s_mod),
                    f_x.coefficients(),
                )
                .1;
                z_i = poly_mod(
                    &poly_add_f(&z_i[..], &a_r_a_s_mul_w_i_r_s[..]),
                    f_x.coefficients(),
                )
                .1
                .try_into()
                .unwrap();
            }
        }

        let tuple = Tuple {
            x_i: x_i.to_vec(),
            y_i: y_i.to_vec(),
            z_i: z_i.to_vec(),
            d_i: d_i.to_vec(),
            M_i_j,
            K_j_i,
        };
        let tuple_json =
            serde_json::to_string(&(tuple, self.alpha_i.clone(), self.sk_i.clone())).unwrap();
        fs::write("tuple_".to_string().add(&id.to_string()), tuple_json).expect("Unable to save !");
    }
}

// define the Ring we work in, right now the function is deterministic
// todo: move to cyclotomic rings
pub fn pick_f_x() -> (Polynomial<Secp256k1>, Vec<Scalar<Secp256k1>>) {
    let mut roots = Vec::new();
    let mut a = vec![Scalar::<Secp256k1>::zero(); N + 1];
    let mut b = vec![Scalar::<Secp256k1>::zero(); N + 1];
    let small_factor_threshold = BigInt::from(1000);
    let small_factors: Vec<&(BigInt, u32)> =
        Scalar::<Secp256k1>::multiplicative_group_order_factorization()
            .iter()
            .filter(|(factor, _)| factor < &small_factor_threshold)
            .collect();
    let smooth_order = small_factors
        .iter()
        .fold(BigInt::one(), |acc, (factor, exp)| acc * (factor.pow(*exp)));
    let group_order = Scalar::<Secp256k1>::group_order();
    let multiplicative_smooth_group_generator = BigInt::mod_pow(
        &Scalar::<Secp256k1>::primitive_root_of_unity().to_bigint(),
        &(&(group_order - BigInt::one()) / &smooth_order),
        group_order,
    );
    let u32_N = u32::try_from(N).expect("N should fit into u32");
    if !&smooth_order.modulus(&BigInt::from(u32_N)).is_zero() {
        panic!("Order must devide multiplicative group order of Secp256k1");
    }
    a[0] = Scalar::from_bigint(&BigInt::one());
    PowerIterator::new(
        Scalar::<Secp256k1>::from_bigint(&BigInt::mod_pow(
            &multiplicative_smooth_group_generator,
            &(&smooth_order / &BigInt::from(u32_N)),
            &Scalar::<Secp256k1>::group_order(),
        )),
        N,
    )
    .enumerate()
    .for_each(|(_i, root)| {
        // pick a root

        let mut root = root;
        if !crate::use_cyclotomic {
            root = Scalar::<Secp256k1>::random();
            b[0] = Scalar::from(&(Scalar::<Secp256k1>::group_order() - root.to_bigint()));
            b[1] = Scalar::from_bigint(&BigInt::one());
            a = poly_mul_f_naive(&a[..], &b[..])[0..N + 1].to_vec();
        }
        roots.push(root);
    });
    if crate::use_cyclotomic {
        let one = Scalar::from(1);
        a[0] = one.clone().neg();
        a[N] = one;
    }
    return (Polynomial::from_coefficients(a), roots);
}

fn pick_N_t() -> [u32; t] {
    let mut trng = rand::thread_rng();
    let u32_N = u32::try_from(N).expect("N should fit into u32");
    assert!(t < N);
    let mut candidate_vec = [0u32; t];
    let mut i = 0;
    while i < t {
        let candidate: u32 = trng.gen_range(0..u32_N);
        if candidate_vec[..i].contains(&candidate) {
            continue;
        }
        candidate_vec[i] = candidate;
        i += 1;
    }
    candidate_vec.sort();
    candidate_vec
}

fn pick_Fq_t() -> [Scalar<Secp256k1>; t] {
    (0..t)
        .map(|_| Scalar::<Secp256k1>::random())
        .collect::<Vec<Scalar<Secp256k1>>>()
        .try_into()
        .unwrap()
}

fn pick_R() -> Vec<Scalar<Secp256k1>> {
    (0..N)
        .map(|_| Scalar::<Secp256k1>::random())
        .collect::<Vec<Scalar<Secp256k1>>>()
        .try_into()
        .unwrap()
}

// set up a ploynomial of degree N from t<N coeffs at locations locs
fn set_poly(coeffs: &[Scalar<Secp256k1>; t], locs: &[u32; t]) -> Vec<Scalar<Secp256k1>> {
    let locs_usize: Vec<_> = (0..locs.len())
        .map(|i| {
            usize::try_from(locs[i])
                .expect("usize should be large enough to represent locs elements")
            // let mut bytes_arr = [0u8; 4];
            // // let mut bytes = locs[i].to_bytes();
            // let mut bytes = locs[i].to_be_bytes();
            // bytes.reverse();

            // for j in 0..bytes.len() {
            //     bytes_arr[j] = bytes[j].clone();
            // }
            // u32::from_le_bytes(bytes_arr) as usize
        })
        .collect();
    for i in 0..locs_usize.len() {
        assert!(locs_usize[i] < N)
    }
    // we assume correctness of input
    let mut poly_new = vec![Scalar::<Secp256k1>::zero(); N];
    for i in 0..t {
        poly_new[locs_usize[i].clone()] = coeffs[i].clone();
    }
    return poly_new;
}

// We recall that if u and v have dimensions m and l, the outer product and the outer sum
// are the ml-dimensional vectors whose (im + j)-th entry is ui · vj and ui + vj respectively.
fn outer_product_t(u: &[Scalar<Secp256k1>], v: &[Scalar<Secp256k1>]) -> Vec<Scalar<Secp256k1>> {
    let mut output: [Scalar<Secp256k1>; t * t] = array_from_fn(|_| Scalar::zero());
    for i in 0..t {
        for j in 0..t {
            output[i * t + j] = &u[i] * &v[j];
        }
    }
    output.to_vec()
}

fn outer_product_c(
    u: &Vec<Vec<Scalar<Secp256k1>>>,
    v: &Vec<Vec<Scalar<Secp256k1>>>,
) -> [Vec<Scalar<Secp256k1>>; c * c] {
    let mut output: [Vec<Scalar<Secp256k1>>; c * c] =
        array_from_fn(|_| vec![Scalar::zero(); 2 * N]);
    for i in 0..c {
        for j in 0..c {
            let mul_ij = poly_mul_f(&u[i], &v[j]);
            for k in 0..mul_ij.len() {
                output[i * c + j][k] = mul_ij[k].clone();
            }
        }
    }
    output
}

fn outer_sum(u: &[u32; t], v: &[u32; t]) -> [u32; t * t] {
    let mut output = [0u32; t * t];
    for i in 0..t {
        for j in 0..t {
            output[i * t + j] = &u[i] + &v[j];
        }
    }
    output
}
