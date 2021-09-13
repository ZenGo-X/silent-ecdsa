use std::collections::HashSet;
extern crate bindgen;

use std::env;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn main(){


    let ignored_macros = IgnoreMacros(
        vec![
            "FP_INFINITE".into(),
            "FP_NAN".into(),
            "FP_NORMAL".into(),
            "FP_SUBNORMAL".into(),
            "FP_ZERO".into(),
            "IPPORT_RESERVED".into(),
        ]
            .into_iter()
            .collect(),
    );

    // Tell cargo to tell rustc to link the system bzip2
    // shared library.
    println!("cargo:rustc-link-lib=ntl");
    println!("cargo:rustc-link-lib=c++");

    println!("cargo:rustc-link-lib=m");
   // println!("cargo:rustc-link-lib=gmp");
    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=wrapper.hpp");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.hpp")
        .clang_arg("-std=c++14")


        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.

        //  .generate_inline_functions(false)
        // Finish the builder and generate the bindings.
           .allowlist_recursively(true)

        .allowlist_type("*.ZZ_p.*")
        .allowlist_type("*.ZZ_pX.*")
        .allowlist_function("*.ZZFromBytes.*")
        .allowlist_type("*.ZZ.*")
        .allowlist_function("*.ZZ_p_init.*")
        .allowlist_function("*.to_ZZ_p.*")
          .allowlist_function("*.mul.*")
        .blocklist_type("*.Mat_value_type.*")
        .blocklist_type("*.Mat_reference.*")
        .blocklist_type("*.Mat_const_reference.*")
        .blocklist_type("*.OptionalVal.*")
        .blocklist_type("*.ZZ_pXModulus.*")
        .blocklist_type("*.ZZ_pX_modulus.*")
        .blocklist_function("*.ZZ_pXModulus.*")
        .blocklist_type("*ZZ_pXMultiplier.*")
        .blocklist_type("*ZZ_pX_multiplier.*")
        .blocklist_function("*ZZ_pXMultiplier.*")
        .parse_callbacks(Box::new(ignored_macros))
        .rustfmt_bindings(true)
        .generate_inline_functions(true)
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

/*
fn main() {
    let ntl_dir = "./depend/ntl-11.5.1/src";
    let ntl_install = format!("{}/ntl", env::var("OUT_DIR").unwrap());
    {
        let path = fs::canonicalize(format!("{}/Configure", ntl_dir)).unwrap();
        let output = Command::new(path.to_str().unwrap())
        //    .arg(format!("--prefix={}", &ntl_install))
            .current_dir(ntl_dir)
            .output()
            .expect("failed to execute process");

        if !output.status.success() {
            std::io::stderr().write_all(&output.stderr).unwrap();
            panic!("NTL: failed to configure");
        }
    }

    {
        let output = Command::new("make")
         //   .arg("install-nodata")
            .current_dir(ntl_dir)
            .output()
            .expect("failed to make");

        if !output.status.success() {
            std::io::stderr().write_all(&output.stderr).unwrap();
            panic!("NTL: failed to ‘make install-nodata’");
        }
    }

    {
        let output = Command::new("sudo")
            .arg("make")
            .arg("install")
            .current_dir(ntl_dir)
            .output()
            .expect("failed to make");

        if !output.status.success() {
            std::io::stderr().write_all(&output.stderr).unwrap();
            panic!("NTL: failed to ‘make install’");
        }
    }
/*
    let ignored_macros = IgnoreMacros(
        vec![
            "FP_INFINITE".into(),
            "FP_NAN".into(),
            "FP_NORMAL".into(),
            "FP_SUBNORMAL".into(),
            "FP_ZERO".into(),
            "IPPORT_RESERVED".into(),
        ]
            .into_iter()
            .collect(),
    );
*/
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.h")
        .clang_arg(format!("-I{}/include", &ntl_install))
        /*
        .whitelist_type("GEN")
        .whitelist_function("mkintn")
        .whitelist_function("GENtostr")
        .whitelist_function("compo")
        .whitelist_function("qfi")
        .whitelist_function("nupow")
        .whitelist_function("qfbcompraw")
        .whitelist_function("primeform")
        .whitelist_function("pari_init")
        .whitelist_function("gneg")
        .whitelist_function("gadd")
        .whitelist_function("shifti")
        .whitelist_function("isprime")
        .parse_callbacks(Box::new(ignored_macros))

         */
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");

    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rustc-link-search=native={}/lib", ntl_install);
    println!("cargo:rustc-link-lib=static=libntl");


}
*/

#[derive(Debug)]
struct IgnoreMacros(HashSet<String>);

impl bindgen::callbacks::ParseCallbacks for IgnoreMacros {
    fn will_parse_macro(&self, name: &str) -> bindgen::callbacks::MacroParsingBehavior {
        if self.0.contains(name) {
            bindgen::callbacks::MacroParsingBehavior::Ignore
        } else {
            bindgen::callbacks::MacroParsingBehavior::Default
        }
    }
}
