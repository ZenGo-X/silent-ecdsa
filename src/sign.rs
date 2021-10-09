
use crate::keygen::Tuple;
use std::fs;
use std::ops::Add;
use curv::elliptic::curves::{Scalar, Secp256k1, Point};
use curv::cryptographic_primitives::secret_sharing::Polynomial;
use crate::n;
use curv::BigInt;


#[derive(Debug)]
pub struct PreSigningKey{
    pub alpha_i: Scalar<Secp256k1>,
    pub sk_i: Scalar<Secp256k1>,
    pub x_i : Scalar<Secp256k1>,
    pub y_i : Scalar<Secp256k1>,
    pub z_i : Scalar<Secp256k1>,
    pub d_i : Scalar<Secp256k1>,
    pub M_i_j : Vec<Scalar<Secp256k1>>,
    pub K_j_i : Vec<Scalar<Secp256k1>>,
}

pub struct PreSignMessage{
    pub x_i_G: Point<Secp256k1>,
    pub M_i_j_G_vec: Vec<Point<Secp256k1>>,
}



impl PreSigningKey{
    pub fn get_key_from_tuple(id: usize)-> Self{
        let data = fs::read_to_string("tuple_".to_string().add(&id.to_string()))
            .expect("Unable to load Tuple, did you run get tuple first? ");
        let (tuple, alpha_i, sk_i): (Tuple , Scalar<Secp256k1>, Scalar<Secp256k1>) = serde_json::from_str(&data).unwrap();
        let rand_point = Scalar::<Secp256k1>::random();
        let x_i = Polynomial::from_coefficients(tuple.x_i.clone().to_vec()).evaluate(&rand_point);
        let y_i = Polynomial::from_coefficients(tuple.y_i.clone().to_vec()).evaluate(&rand_point);
        let z_i = Polynomial::from_coefficients(tuple.z_i.clone().to_vec()).evaluate(&rand_point);
        let d_i = Polynomial::from_coefficients(tuple.d_i.clone().to_vec()).evaluate(&rand_point);
        let mut M_i_j : Vec<Scalar<Secp256k1>>  = Vec::new();
        let mut K_j_i : Vec<Scalar<Secp256k1>>  = Vec::new();
        for i in 0..n-1{
            M_i_j.push(
                Polynomial::from_coefficients(tuple.M_i_j[i].clone().to_vec()).evaluate(&rand_point)
            );
            K_j_i.push(
                Polynomial::from_coefficients(tuple.K_j_i[i].clone().to_vec()).evaluate(&rand_point)
            );
        }
        return PreSigningKey{
            alpha_i,
            sk_i,
            x_i,
            y_i,
            z_i,
            d_i,
            M_i_j,
            K_j_i,
        }

    }

    pub fn presign_message(&self) -> PreSignMessage {
        let mut M_i_j_G_vec : Vec<Point<Secp256k1>>  = Vec::new();
        for val in &self.M_i_j{
            M_i_j_G_vec.push(val * Point::<Secp256k1>::generator())
        }
        PreSignMessage{
            x_i_G: &self.x_i * Point::<Secp256k1>::generator(),
            M_i_j_G_vec,
        }
    }

    // helper function to collect presign messages intended for party "id"
    pub fn pre_sign_message_collect(presign_messages : &[PreSignMessage], id: usize) -> (Vec<Point<Secp256k1>>,Vec<Point<Secp256k1>>){
        //todo: check id is valid
        let mut M_j_i_G_vec : Vec<Point<Secp256k1>>  = Vec::new();
        let mut x_i_G_vec : Vec<Point<Secp256k1>>  = Vec::new();
        for i in 0..n{
            if i != id{
                x_i_G_vec.push(presign_messages[i].x_i_G.clone());
                let ind = if i<id {id-1} else {id};
                M_j_i_G_vec.push(presign_messages[i].M_i_j_G_vec[ind].clone())
            }
        }

        (x_i_G_vec, M_j_i_G_vec)
    }

    // if valid, output  r
    pub fn presign_message_verify(&self, x_i_G_vec : &[Point<Secp256k1>], M_j_i_vec: &[Point<Secp256k1>] ) -> Result<BigInt, ()>{
        if M_j_i_vec.len() != n - 1 {return Err(())};
        if x_i_G_vec.len() != n - 1 {return Err(())};

        for j in 0..n-1{
            let a_i_x_j_G = &self.alpha_i * &x_i_G_vec[j];
            let K_j_i_G_plus_a_i_x_j_G = &self.K_j_i[j] * Point::<Secp256k1>::generator() + &a_i_x_j_G;
            if M_j_i_vec[j] != K_j_i_G_plus_a_i_x_j_G{
                return Err(())
            }

        }
        let R = x_i_G_vec.iter()
            .fold(&self.x_i * Point::<Secp256k1>::generator(), |acc, x|&acc + x);

        return Ok(R.x_coord().unwrap());

    }

}

