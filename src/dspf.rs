use crate::dpf::{Key, DPF};
use crate::N;
use curv::arithmetic::Zero;
use curv::elliptic::curves::secp256_k1::FE;
use curv::elliptic::curves::traits::ECScalar;
use curv::BigInt;

pub struct DSPF {
    key: Vec<Key>,
}

impl DSPF {
    pub fn gen(alpha_vec: &[BigInt], beta_vec: &[FE]) -> (Self, Self) {
        // TODO: make sure there are no repetitions in alpha_vec?
        let t = alpha_vec.len();
        assert_eq!(t, beta_vec.len());

        let (key0, key1): (Vec<_>, Vec<_>) = (0..t)
            .map(|i| {
                assert!(alpha_vec[i] < BigInt::from(N as u32));
                assert!(alpha_vec[i] >= BigInt::zero());
                DPF::gen(&alpha_vec[i], &beta_vec[i])
            })
            .unzip();
        (DSPF { key: key0 }, DSPF { key: key1 })
    }

    pub fn eval(&self, b: &u8, x: &BigInt) -> FE {
        assert!(x < &BigInt::from(N as u32));
        assert!(x >= &BigInt::zero());
        self.key
            .iter()
            .fold(FE::zero(), |acc, dpf_key| acc + DPF::eval(b, dpf_key, x))
    }
}

mod tests {
    use crate::dspf::DSPF;
    use crate::N;
    use curv::elliptic::curves::secp256_k1::FE;
    use curv::elliptic::curves::traits::ECScalar;
    use curv::BigInt;

    #[test]
    fn test_simple_dspf() {
        let alpha_vec = &[BigInt::from(10), BigInt::from(20), BigInt::from(30)].to_vec();
        let beta_vec: &[FE] = &[
            ECScalar::new_random(),
            ECScalar::new_random(),
            ECScalar::new_random(),
        ]
        .to_vec();
        let (key0, key1) = DSPF::gen(&alpha_vec[..], &beta_vec[..]);
        for x in 0..N {
            let x_bn = BigInt::from(x as u32);
            let fe0 = key0.eval(&0u8, &x_bn);
            let fe1 = key1.eval(&1u8, &x_bn);

            if x_bn == alpha_vec[0] {
                assert_eq!(fe0 + fe1, beta_vec[0]);
            } else if x_bn == alpha_vec[1] {
                assert_eq!(fe0 + fe1, beta_vec[1]);
            } else if x_bn == alpha_vec[2] {
                assert_eq!(fe0 + fe1, beta_vec[2]);
            } else {
                assert_eq!(fe0.to_big_int(), FE::q() - fe1.to_big_int());
            }
        }
    }
}
