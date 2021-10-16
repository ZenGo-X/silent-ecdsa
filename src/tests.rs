use crate::keygen::{pick_f_x, LongTermKey};
use crate::n;
use crate::sign::{PreSignMessage, PreSigningKey, Signature};
use curv::elliptic::curves::Scalar;
use curv::BigInt;
use std::fs;
use std::ops::Add;

#[test]
fn test_run_keygen() {
    let key = LongTermKey::trusted_key_gen();
    let a = LongTermKey::sample_a();
    let f_x = pick_f_x();
    key[0].get_tuple(a.clone(), &f_x.0, 0);
    key[1].get_tuple(a.clone(), &f_x.0, 1);
}

#[test]
fn test_presign() {
    // note: we assume keygen run already, otherwise an error.
    let mut presign_keys: Vec<PreSigningKey> = Vec::new();
    let mut presign_messages: Vec<PreSignMessage> = Vec::new();
    let sample_at_point = Scalar::from_bigint(&BigInt::from(4)); // our ring has zeros in 1,2,3,4,...N
    for id in 0..n {
        let presign_key_id = PreSigningKey::get_key_from_tuple(&sample_at_point, id);
        presign_messages.push(presign_key_id.presign_message());
        presign_keys.push(presign_key_id);
    }
    // for two parties
    /*
    let x = presign_keys[0].x_i.clone() + presign_keys[1].x_i.clone();
    let y = presign_keys[0].y_i.clone() + presign_keys[1].y_i.clone();
    let z = presign_keys[0].z_i.clone() + presign_keys[1].z_i.clone();
    let d = presign_keys[0].d_i.clone() + presign_keys[1].d_i.clone();
    let sk = presign_keys[0].sk_i.clone() + presign_keys[1].sk_i.clone();
    //assert_eq!(sk * &y , d);
    //assert_eq!(x * &y , z);

     */

    for id in 0..n {
        let (x_i_G_vec, M_j_i_G_vec) =
            PreSigningKey::pre_sign_message_collect(&presign_messages[..], id);
        presign_keys[id]
            .presign_message_verify(&x_i_G_vec, &M_j_i_G_vec)
            .expect("error presign message verify party");
    }
}

#[test]
fn test_presign_and_sign() {
    // note: we assume keygen run already, otherwise an error.
    let mut presign_keys: Vec<PreSigningKey> = Vec::new();
    let mut presign_messages: Vec<PreSignMessage> = Vec::new();
    let sample_at_point = Scalar::from_bigint(&BigInt::from(4)); // our ring has zeros in 1,2,3,4,...N
    for id in 0..n {
        let presign_key_id = PreSigningKey::get_key_from_tuple(&sample_at_point, id);
        presign_messages.push(presign_key_id.presign_message());
        presign_keys.push(presign_key_id);
    }
    let mut r_vec = Vec::new();
    for id in 0..n {
        let (x_i_G_vec, M_j_i_G_vec) =
            PreSigningKey::pre_sign_message_collect(&presign_messages[..], id);
        r_vec.push(
            presign_keys[id]
                .presign_message_verify(&x_i_G_vec, &M_j_i_G_vec)
                .unwrap(),
        );
    }

    let message = BigInt::from(1234);
    let mut local_sig_vec = Vec::new();
    for id in 0..n {
        local_sig_vec.push(Signature::local_signature(
            &message,
            &presign_keys[id],
            &r_vec[id],
        ));
    }

    let mut sig_vec = Vec::new();
    for id in 0..n {
        sig_vec.push(Signature::output(&local_sig_vec[..], &r_vec[id]));
    }

    for id in 0..n {
        let data = fs::read_to_string("long_term_key_".to_string().add(&id.to_string()))
            .expect("Unable to load long term key, did you run get keygen first? ");
        let key: LongTermKey = serde_json::from_str(&data).unwrap();
        sig_vec[id]
            .verify(&key.pk, &message)
            .expect("failed verification");
    }
}
