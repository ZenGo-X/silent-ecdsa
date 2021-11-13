use crate::fft::{bit_rev_radix_2_intt, bit_rev_radix_2_ntt};
use curv::arithmetic::{Converter, Modulo, Zero};
use curv::cryptographic_primitives::secret_sharing::PolynomialDegree;
use curv::cryptographic_primitives::secret_sharing::{ffts::multiply_polynomials, Polynomial};
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;
use std::convert::TryInto;
use std::ops::{Add, Neg};

// a mod (b) = (d,r)  such that a = b*d + r
pub fn poly_mod(
    a: &[Scalar<Secp256k1>],
    b: &[Scalar<Secp256k1>],
) -> (Vec<Scalar<Secp256k1>>, Vec<Scalar<Secp256k1>>) {
    if is_cyclotomic(b) {
        let mut b_deg = 0;
        for x in (0..b.len()).rev() {
            if b[x] != Scalar::<Secp256k1>::zero() {
                b_deg = x;
                break;
            }
        }
        return poly_mod_cyclotomic(a, b_deg);
    }
    let a_poly = Polynomial::from_coefficients(a.to_vec());
    let b_poly = Polynomial::from_coefficients(b.to_vec());
    let b_deg = b_poly.degree();
    assert!(b_poly.degree() > PolynomialDegree::Finite(0));

    let mut a_trim = a.clone().to_vec();
    a_trim = trim(&a_trim);
    let a_trim_poly = Polynomial::from_coefficients(a_trim.clone());

    if a_poly.degree() < b_poly.degree() {
        return ([Scalar::<Secp256k1>::zero()].to_vec(), a_trim);
    }
    let a_deg = match a_poly.degree() {
        PolynomialDegree::Finite(n) => n,
        PolynomialDegree::Infinity => panic!("b_poly shouldn't be the zero polynomial"),
    };
    let b_deg = match b_poly.degree() {
        PolynomialDegree::Finite(n) => n,
        PolynomialDegree::Infinity => panic!("b_poly shouldn't be the zero polynomial"),
    };
    let lc_inv = b_poly.coefficients()[b_deg].clone().invert().unwrap();
    let mut coef = vec![Scalar::<Secp256k1>::zero(); (a_deg - b_deg + 1) as usize];
    let mut reminder_poly = a_trim_poly.clone();
    let mut reminder_coef = reminder_poly.coefficients().to_vec();

    while reminder_poly.degree() >= b_poly.degree() {
        let d = match reminder_poly.degree() {
            PolynomialDegree::Finite(n) => n,
            PolynomialDegree::Infinity => break,
        };
        let d = d - b_deg;
        let c = &lc(&reminder_poly) * &lc_inv;
        for i in 0..b_deg as usize {
            reminder_coef[i + d] = &reminder_coef[i + d] - &c * &b_poly.coefficients().to_vec()[i];
        }

        reminder_coef.pop();
        reminder_coef = trim(&reminder_coef);

        reminder_poly = Polynomial::from_coefficients(reminder_coef.clone());
        coef[d] = c;
    }
    return (coef, reminder_coef);
}

fn is_cyclotomic(a: &[Scalar<Secp256k1>]) -> bool {
    let one = Scalar::<Secp256k1>::from_bigint(&BigInt::from(1));
    let minus_one = one.clone().neg();
    let mut has_seen_one = false;
    if a[0] != minus_one {
        return false;
    }
    for i in 1..a.len() {
        if a[i] != Scalar::<Secp256k1>::zero() {
            if a[i] != one {
                return false;
            }
            if has_seen_one {
                return false;
            }
            has_seen_one = true;
        }
    }
    true
}
// computes a mod (x^deg - 1).
pub fn poly_mod_cyclotomic(
    a: &[Scalar<Secp256k1>],
    deg: usize,
) -> (Vec<Scalar<Secp256k1>>, Vec<Scalar<Secp256k1>>) {
    if a.len() <= deg {
        return (vec![Scalar::<Secp256k1>::zero()], a.to_vec());
    }

    let mut rem = a.to_vec();
    let mut div_result: Vec<Scalar<Secp256k1>> = Vec::with_capacity(a.len() - deg);
    (deg..a.len()).rev().for_each(|i| {
        rem[i - deg] = &rem[i - deg] + &a[i];
        div_result.push(a[i].clone());
    });
    (
        div_result.into_iter().rev().collect(),
        rem.into_iter().take(deg).collect(),
    )
}

//fn trim(poly: &Vec<Scalar<Secp256k1>>) -> Vec<Scalar<Secp256k1>>{
//    poly[0..Polynomial::from_coefficients(poly.clone()).degree() as usize].to_vec()
//}

fn trim(poly: &Vec<Scalar<Secp256k1>>) -> Vec<Scalar<Secp256k1>> {
    let poly_len = poly.len();
    let mut counter = 0;
    for j in 0..poly_len {
        if poly[poly_len - 1 - j] == Scalar::zero() {
            counter += 1;
        } else {
            break;
        }
    }
    return poly[0..poly_len - counter].to_vec();
}

// largest coefficient
fn lc(poly: &Polynomial<Secp256k1>) -> Scalar<Secp256k1> {
    match poly.degree() {
        PolynomialDegree::Infinity => panic!("Shouldn't get zero polynomial here"),
        PolynomialDegree::Finite(n) => poly.coefficients()[n].clone(),
    }
}

// wrapper around poly_mul
pub fn poly_mul_f(a: &[Scalar<Secp256k1>], b: &[Scalar<Secp256k1>]) -> Vec<Scalar<Secp256k1>> {
    if !crate::use_fft {
        return poly_mul_f_naive(a, b);
    }
    let va = Vec::from(a);
    let pa = Polynomial::from_coefficients(va);

    let vb = Vec::from(b);
    let pb = Polynomial::from_coefficients(vb);
    multiply_polynomials(pa, pb)
        .into_coefficients()
        .into_iter()
        .take(a.len() + b.len() - 1)
        .collect()
}

pub fn poly_add_f(a: &[Scalar<Secp256k1>], b: &[Scalar<Secp256k1>]) -> Vec<Scalar<Secp256k1>> {
    let a_poly = Polynomial::from_coefficients(a.to_vec());
    let b_poly = Polynomial::from_coefficients(b.to_vec());
    let c_poly = a_poly.add(&b_poly);
    c_poly.coefficients().to_vec()
}

// wrapper around poly_mul
pub fn poly_mul_f_naive(
    a: &[Scalar<Secp256k1>],
    b: &[Scalar<Secp256k1>],
) -> Vec<Scalar<Secp256k1>> {
    let mut a_bn: Vec<_> = a.to_vec();
    let mut b_bn: Vec<_> = b.to_vec();

    if a_bn.len() < b_bn.len() {
        for _ in 0..(b_bn.len() - a_bn.len()) {
            a_bn.push(Scalar::zero());
        }
    }
    if b_bn.len() < a_bn.len() {
        for _ in 0..(a_bn.len() - b_bn.len()) {
            b_bn.push(Scalar::zero());
        }
    }

    // todo: propagate errors
    let c = poly_mul(&a_bn[..], &b_bn[..]);
    c
}
// c = a * b
pub fn poly_mul(a: &[Scalar<Secp256k1>], b: &[Scalar<Secp256k1>]) -> Vec<Scalar<Secp256k1>> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    let c: Vec<_> = (0..2 * n - 1)
        .map(|i| {
            let mut acc = Scalar::zero();
            let mut tmp;
            for j in 0..=i {
                tmp = Scalar::zero();
                if j < n && i - j < n {
                    tmp = &a[j] * &b[i - j];
                }
                acc = &acc + &tmp;
            }
            acc
        })
        .collect();
    c
}

pub fn poly_mul_fft(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let p = Scalar::<Secp256k1>::group_order();
    let w = BigInt::from_str_radix(
        "8538e806a1ff6810b4181dbeca3ee5d2095a19ff2489c3a80693427f238d0ba4",
        16,
    )
    .unwrap();
    let s = 6; // two-arity of secp256k1
    let n = 6; //bits
    let n_elements = 2_i32.pow(n) as usize;
    assert!(n <= s);

    let mut F_points_vec = a.to_vec();
    let mut G_points_vec = b.to_vec();
    assert_eq!(F_points_vec.len(), G_points_vec.len());
    for _ in 0..n_elements - F_points_vec.len() {
        F_points_vec.push(BigInt::zero());
        G_points_vec.push(BigInt::zero());
    }
    bit_rev_radix_2_ntt(&mut F_points_vec, &p, &w);
    bit_rev_radix_2_ntt(&mut G_points_vec, &p, &w);

    // mul points :
    let FG_points_vec: Vec<_> = (0..64)
        .map(|i| &F_points_vec[i] * &G_points_vec[i])
        .collect();

    // go back to coeffs
    let mut FG_coeffs_vec = FG_points_vec.clone();
    bit_rev_radix_2_intt(&mut FG_coeffs_vec, &p, &w);
    FG_coeffs_vec
}

mod tests {

    use crate::fft::{bit_rev_radix_2_intt, bit_rev_radix_2_ntt};
    use crate::poly::{poly_mod, poly_mod_cyclotomic, poly_mul, poly_mul_f, poly_mul_fft};
    use crate::utils::array_from_fn;
    use curv::arithmetic::{Converter, Modulo, One, Samplable, Zero};
    use curv::cryptographic_primitives::secret_sharing::Polynomial;
    use curv::elliptic::curves::{Scalar, Secp256k1};
    use curv::BigInt;
    use std::time::Instant;

    // 1) https://crypto.stackexchange.com/questions/63614/finding-the-n-th-root-of-unity-in-a-finite-field
    // 2) the factorization of q -1 in secp256k1 is 2^6 * 3 * 149 * 631 * 107361793816595537 * 174723607534414371449 * 341948486974166000522343609283189
    // 3) candidate root: "8538e806a1ff6810b4181dbeca3ee5d2095a19ff2489c3a80693427f238d0ba4"
    #[test]
    fn test_find_primitive_root_of_unity() {
        let mut test = BigInt::one();
        let mut candidate_root = BigInt::one();
        while test == BigInt::one() {
            let x = BigInt::sample_below(&Scalar::<Secp256k1>::group_order());
            let n_inv =
                BigInt::mod_inv(&BigInt::from(64), &Scalar::<Secp256k1>::group_order()).unwrap();
            let q_min_1_n_inv = BigInt::mod_mul(
                &(Scalar::<Secp256k1>::group_order() - BigInt::one()),
                &n_inv,
                &Scalar::<Secp256k1>::group_order(),
            );
            candidate_root =
                BigInt::mod_pow(&x, &q_min_1_n_inv, &Scalar::<Secp256k1>::group_order());
            test = BigInt::mod_pow(
                &candidate_root,
                &BigInt::from(32),
                &Scalar::<Secp256k1>::group_order(),
            );
        }
        // check
        let mut g = candidate_root.clone();
        assert_ne!(g, BigInt::one());
        for _ in 1..6 {
            g = BigInt::mod_mul(&g, &g, &Scalar::<Secp256k1>::group_order());
            assert_ne!(g, BigInt::one());
        }
        g = BigInt::mod_mul(&g, &g, &Scalar::<Secp256k1>::group_order());
        assert_eq!(g, BigInt::one());
        println!(
            "candidate root: {:?}",
            candidate_root.clone().to_str_radix(16)
        );
    }

    #[test]
    fn test_poly_mul_over_secp256k1() {
        let p = Scalar::<Secp256k1>::group_order();
        let w = BigInt::from_str_radix(
            "8538e806a1ff6810b4181dbeca3ee5d2095a19ff2489c3a80693427f238d0ba4",
            16,
        )
        .unwrap();
        let s = 6; // two-arity of secp256k1
        let n = 6; //bits
        let n_elements = 2_i32.pow(n) as usize;
        assert!(n <= s);
        // TODO: use polynomial
        let F_coeffs_vec: Vec<_> = (0..n_elements)
            .map(|i| {
                // we fill only half the slots with actual numbers to allow for the multiplication
                if i < n_elements / 2 {
                    Scalar::<Secp256k1>::random().to_bigint()
                } else {
                    BigInt::zero()
                }
            })
            .collect();
        let mut F_points_vec = F_coeffs_vec.clone();
        bit_rev_radix_2_ntt(&mut F_points_vec, &p, &w);

        let G_coeffs_vec: Vec<_> = (0..n_elements)
            .map(|i| {
                if i < n_elements / 2 {
                    Scalar::<Secp256k1>::random().to_bigint()
                } else {
                    BigInt::zero()
                }
            })
            .collect();
        let mut G_points_vec = G_coeffs_vec.clone();
        bit_rev_radix_2_ntt(&mut G_points_vec, &p, &w);

        // mul points :
        let FG_points_vec: Vec<_> = (0..64)
            .map(|i| &F_points_vec[i] * &G_points_vec[i])
            .collect();

        // go back to coeffs
        let mut FG_coeffs_vec = FG_points_vec.clone();
        bit_rev_radix_2_intt(&mut FG_coeffs_vec, &p, &w);

        // manual mul for comparison
        let coeff0 = BigInt::mod_mul(
            &G_coeffs_vec[0],
            &F_coeffs_vec[0],
            &Scalar::<Secp256k1>::group_order(),
        );
        let coeff1_1 = BigInt::mod_mul(
            &G_coeffs_vec[0],
            &F_coeffs_vec[1],
            &Scalar::<Secp256k1>::group_order(),
        );
        let coeff1_2 = BigInt::mod_mul(
            &G_coeffs_vec[1],
            &F_coeffs_vec[0],
            &Scalar::<Secp256k1>::group_order(),
        );
        let coeff1 = BigInt::mod_add(&coeff1_1, &coeff1_2, &Scalar::<Secp256k1>::group_order());
        let coeff62 = BigInt::mod_mul(
            &G_coeffs_vec[31],
            &F_coeffs_vec[31],
            &Scalar::<Secp256k1>::group_order(),
        );
        let coeff61_1 = BigInt::mod_mul(
            &G_coeffs_vec[31],
            &F_coeffs_vec[30],
            &Scalar::<Secp256k1>::group_order(),
        );
        let coeff61_2 = BigInt::mod_mul(
            &G_coeffs_vec[30],
            &F_coeffs_vec[31],
            &Scalar::<Secp256k1>::group_order(),
        );
        let coeff61 = BigInt::mod_add(&coeff61_1, &coeff61_2, &Scalar::<Secp256k1>::group_order());
        assert_eq!(coeff0, FG_coeffs_vec[0]);
        assert_eq!(coeff1, FG_coeffs_vec[1]);
        assert_eq!(coeff61, FG_coeffs_vec[61]);
        assert_eq!(coeff62, FG_coeffs_vec[62]);
    }

    #[test]
    fn test_simple_poly_mul() {
        let a = [1, 2, 3, 4].map(|i| Scalar::from(i));
        let b = [5, 6, 7, 8].map(|i| Scalar::from(i));
        let c = poly_mul(&a, &b);
        println!("c: {:?}", c);
    }

    #[test]
    fn test_compare_simple_poly_and_fft_poly() {
        let mut dur_av_simple = 0;
        let mut dur_av_fft = 0;
        let rep = 100 as usize;
        for _ in 0..rep {
            let a: Vec<_> = (0..32).map(|_| Scalar::<Secp256k1>::random()).collect();
            let b: Vec<_> = (0..32).map(|_| Scalar::<Secp256k1>::random()).collect();

            let start_simple = Instant::now();
            let mut c_simple = poly_mul(&a[..], &b[..]);
            let duration_simple = start_simple.elapsed();
            dur_av_simple = dur_av_simple + duration_simple.as_millis();

            let a_bn: Vec<_> = a.iter().map(Scalar::to_bigint).collect();
            let b_bn: Vec<_> = b.iter().map(Scalar::to_bigint).collect();
            let start_fft = Instant::now();
            let c_fft = poly_mul_fft(&a_bn[..], &b_bn[..]);
            let duration_fft = start_fft.elapsed();
            dur_av_fft = dur_av_fft + duration_fft.as_millis();

            c_simple.push(Scalar::zero());
            let c_fft_scalar: Vec<_> = c_fft.iter().map(Scalar::from_bigint).collect();
            assert_eq!(c_simple, c_fft_scalar);
        }
        println!("duration simple: {:?}", dur_av_simple);
        println!("duration fft: {:?}", dur_av_fft);
    }

    use std::mem;
    use std::ops::Neg;
    use std::ptr;

    #[test]
    fn test_poly_mod() {
        let a: Vec<_> = (0..63).map(|_| Scalar::<Secp256k1>::random()).collect();
        let b: Vec<_> = (0..33).map(|_| Scalar::<Secp256k1>::random()).collect();
        let (mut d, r) = poly_mod(&a, &b);
        for _ in 0..(b.len() - d.len()) {
            d.push(Scalar::zero());
        }
        let dq = poly_mul_f(&d[..], &b[..]);
        let dq_poly = Polynomial::from_coefficients(dq);
        let r_poly = Polynomial::from_coefficients(r);
        let dq_plus_r = &dq_poly + &r_poly;
        let mut dq_plus_r_coef = dq_plus_r.coefficients().to_vec();
        for _ in 0..dq_plus_r_coef.len() - a.len() {
            dq_plus_r_coef.pop();
        }

        assert_eq!(dq_plus_r_coef, a);
    }

    #[test]
    fn test_poyl_mod_cyclotomic() {
        let one_scalar = Scalar::<Secp256k1>::from_bigint(&BigInt::from(1));
        let deg = 320;
        let a: Vec<_> = (0..(2 * deg - 1))
            .map(|_| Scalar::<Secp256k1>::random())
            .collect();
        let mut b: Vec<_> = (0..(deg + 1))
            .map(|_| Scalar::<Secp256k1>::zero())
            .collect();
        b[0] = one_scalar.clone().neg();
        b[deg] = one_scalar.clone();
        let (mut d, r) = poly_mod_cyclotomic(&a, deg);
        for _ in 0..(b.len() - d.len()) {
            d.push(Scalar::zero());
        }
        let dq = poly_mul_f(&d[..], &b[..]);
        let dq_poly = Polynomial::from_coefficients(dq);
        let r_poly = Polynomial::from_coefficients(r);
        let dq_plus_r = &dq_poly + &r_poly;
        let mut dq_plus_r_coef = dq_plus_r.coefficients().to_vec();
        for _ in 0..dq_plus_r_coef.len() - a.len() {
            dq_plus_r_coef.pop();
        }

        assert_eq!(dq_plus_r_coef, a);
    }

    #[test]
    fn test_poly_mod_edge_case() {
        // let a: Vec<_> = (0..63).map(|_| Scalar::<Secp256k1>::random()).collect();
        let mut a: [Scalar<Secp256k1>; 34] = array_from_fn(|_| Scalar::zero());
        a[33] = Scalar::<Secp256k1>::random();

        let b: Vec<_> = (0..33).map(|_| Scalar::<Secp256k1>::random()).collect();
        let (mut d, r) = poly_mod(&a, &b);
        for _ in 0..(b.len() - d.len()) {
            d.push(Scalar::zero());
        }
        let dq = poly_mul_f(&d[..], &b[..]);
        let dq_poly = Polynomial::from_coefficients(dq);
        let r_poly = Polynomial::from_coefficients(r);
        let dq_plus_r = &dq_poly + &r_poly;
        let mut dq_plus_r_coef = dq_plus_r.coefficients().to_vec();
        for _ in 0..dq_plus_r_coef.len() - a.len() {
            dq_plus_r_coef.pop();
        }

        assert_eq!(dq_plus_r_coef, a);
    }
}
