use crate::dpf::DPF;
use crate::N;
use curv::arithmetic::Zero;
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;

#[derive(Debug, Clone)]
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
    pub fn gen(
        alpha_vec: &[BigInt],
        beta_vec: &[Scalar<Secp256k1>],
        two_n_flag: bool,
    ) -> (Self, Self) {
        // TODO: make sure there are no repetitions in alpha_vec?
        let t = alpha_vec.len();
        assert_eq!(t, beta_vec.len());

        let (key0, key1): (Vec<_>, Vec<_>) = (0..t)
            .map(|i| {
                if two_n_flag {
                    assert!(alpha_vec[i] < BigInt::from(2 * N as u32));
                } else {
                    assert!(alpha_vec[i] < BigInt::from(N as u32));
                }
                assert!(alpha_vec[i] >= BigInt::zero());
                DPF::gen(&alpha_vec[i], &beta_vec[i])
            })
            .unzip();
        (DSPF { key: key0 }, DSPF { key: key1 })
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
            .fold(Scalar::zero(), |acc, dpf| acc + dpf.eval(b, x))
    }

    pub fn full_eval_N(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
        let zero_vec = vec![Scalar::zero(); N];
        self.key.iter().fold(zero_vec, |acc, dpf_key| {
            let dpf_i_full_eval = dpf_key.full_eval_N(b);
            (0..N).map(|i| &acc[i] + &dpf_i_full_eval[i]).collect()
        })
    }

    pub fn full_eval_2N(&self, b: &u8) -> Vec<Scalar<Secp256k1>> {
        let zero_vec = vec![Scalar::zero(); 2 * N];
        self.key.iter().fold(zero_vec, |acc, dpf_key| {
            let dpf_i_full_eval = dpf_key.full_eval_2N(b);
            (0..2 * N).map(|i| &acc[i] + &dpf_i_full_eval[i]).collect()
        })
    }
}

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

        let (key0, key1) = DSPF::gen(&alpha_vec[..], &beta_vec[..], false);
        for x in 0..N {
            let x_bn = BigInt::from(x as u32);
            let fe0 = key0.eval(&0u8, &x_bn, false);
            let fe1 = key1.eval(&1u8, &x_bn, false);
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
        let p0_shares = key0.full_eval_N(&0u8);
        let p1_shares = key1.full_eval_N(&1u8);
        assert_eq!(p0_shares, fe0_vec);
        assert_eq!(p1_shares, fe1_vec);
    }
}
