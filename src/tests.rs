use crate::keygen::{pick_f_x, LongTermKey};
use crate::n;
use crate::sign::{PreSignMessage, PreSigningKey};
use curv::elliptic::curves::Scalar;
use curv::BigInt;

#[test]
fn test_run_keygen() {
    let key = LongTermKey::trusted_key_gen();
    let a = LongTermKey::sample_a();
    let f_x = pick_f_x();
    key[0].get_tuple(a.clone(), &f_x.0, 0);
    key[1].get_tuple(a.clone(), &f_x.0, 1);
    key[2].get_tuple(a.clone(), &f_x.0, 2);
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
    for id in 0..n {
        let (x_i_G_vec, M_j_i_G_vec) =
            PreSigningKey::pre_sign_message_collect(&presign_messages[..], id);
        presign_keys[id]
            .presign_message_verify(&x_i_G_vec, &M_j_i_G_vec)
            .expect("error presign message verify party");
    }
}
