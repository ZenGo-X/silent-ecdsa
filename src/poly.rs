

mod tests {

    #![allow(non_upper_case_globals)]
    #![allow(non_camel_case_types)]
    #![allow(non_snake_case)]
    #![allow(dead_code)]

    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));



    use std::ffi::CStr;
    use std::mem::swap;
    use std::ops::Neg;
    use std::{str, ptr};




    #[test]
    fn test_simple_binding() {
        let bytes: [u8;2] = [0x02, 0x03];
        let bytes_ptr : *const u8 = &bytes[0];



            //let res2 : *mut NTL_ZZ =  *res;
           // let z : *mut NTL_ZZ = &mut res;
          //  et my_speed_ptr: *mut i32 = &mut my_speed;


        unsafe {

            let res = NTL_ZZFromBytes1( bytes_ptr, 3);
        }

    }
}
