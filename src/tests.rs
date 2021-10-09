use crate::keygen::{LongTermKey, pick_f_x};
use crate::sign::{PreSigningKey, PreSignMessage};

#[test]
fn test_run_keygen() {
    let key = LongTermKey::trusted_key_gen();
    let a = LongTermKey::sample_a();
    let f_x = pick_f_x();
    let _tuple1 = key[0].get_tuple(a, &f_x, 0);

}


#[test]
fn test_presign() {
    // note: we assume keygen run already, otherwise an error.
    let mut presign_messages : Vec<PreSignMessage> = Vec::new();
    let presignkey_0 = PreSigningKey::get_key_from_tuple(0);
    let presignkey_1 = PreSigningKey::get_key_from_tuple(1);
    presign_messages.push(presignkey_0.presign_message());
    presign_messages.push(presignkey_1.presign_message());

    let (x_i_G_vec_0, M_j_i_G_vec_0) =  PreSigningKey::pre_sign_message_collect(&presign_messages[..], 0);
    let (x_i_G_vec_1, M_j_i_G_vec_1) =  PreSigningKey::pre_sign_message_collect(&presign_messages[..], 1);

    presignkey_0.presign_message_verify(&x_i_G_vec_0[..], &M_j_i_G_vec_0[..]).expect("error presign message verify party 0");
    presignkey_1.presign_message_verify(&x_i_G_vec_1[..], &M_j_i_G_vec_1[..]).expect("error presign message verify party 0");


    //  assert!(false);
}
