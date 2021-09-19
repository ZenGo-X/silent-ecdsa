use crate::fft::{bit_rev_radix_2_intt, bit_rev_radix_2_ntt};
use curv::arithmetic::{Converter, Modulo, Zero};
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;

// c = a * b
pub fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    let c: Vec<_> = (0..2 * n - 1)
        .map(|i| {
            let mut acc = BigInt::zero();
            let mut tmp = BigInt::zero();
            for j in 0..=i {
                tmp = BigInt::zero();
                if j < n && i - j < n {
                    tmp = BigInt::mod_mul(&a[j], &b[i - j], &Scalar::<Secp256k1>::group_order());
                }
                acc = BigInt::mod_add(&acc, &tmp, &Scalar::<Secp256k1>::group_order());
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
    for i in 0..n_elements - F_points_vec.len() {
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
    use crate::poly::{poly_mul, poly_mul_fft};
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
        let a = [1, 2, 3, 4].map(|i| BigInt::from(i));
        let b = [5, 6, 7, 8].map(|i| BigInt::from(i));
        let c = poly_mul(&a, &b);
        println!("c: {:?}", c);
    }

    #[test]
    fn test_compare_simple_poly_and_fft_poly() {
        let mut dur_av_simple = 0;
        let mut dur_av_fft = 0;
        let rep = 100 as usize;
        for _ in 0..rep {
            let a: Vec<_> = (0..32)
                .map(|_| Scalar::<Secp256k1>::random().to_bigint())
                .collect();
            let b: Vec<_> = (0..32)
                .map(|_| Scalar::<Secp256k1>::random().to_bigint())
                .collect();

            let start_simple = Instant::now();
            let mut c_simple = poly_mul(&a[..], &b[..]);
            let duration_simple = start_simple.elapsed();
            dur_av_simple = dur_av_simple + duration_simple.as_millis();

            let start_fft = Instant::now();
            let c_fft = poly_mul_fft(&a[..], &b[..]);
            let duration_fft = start_fft.elapsed();
            dur_av_fft = dur_av_fft + duration_fft.as_millis();

            c_simple.push(BigInt::zero());
            assert_eq!(c_simple, c_fft);
        }
        println!("duration simple: {:?}", dur_av_simple);
        println!("duration fft: {:?}", dur_av_fft);
    }
}