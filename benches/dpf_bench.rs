use criterion::PlotConfiguration;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use curv::elliptic::curves::{Scalar, Secp256k1};
use curv::BigInt;
use silent_t_ecdsa::dpf::DPF;

fn criterion_benchmark(c: &mut Criterion) {
    let bench_sizes: Vec<u32> = (0..15).map(|i| 1 << i).collect();
    let mut group = c.benchmark_group("dpf");
    let plot_config = PlotConfiguration::default().summary_scale(criterion::AxisScale::Logarithmic);
    group.plot_config(plot_config);
    let alpha = &BigInt::from(10);
    let beta = &Scalar::random();
    let (dpf0, _) = DPF::gen(alpha, beta);
    for size in bench_sizes {
        group.bench_with_input(BenchmarkId::new("Naive", &size), &size, |b, size| {
            b.iter(|| {
                let x: Vec<Scalar<Secp256k1>> = (0..*size)
                    .map(|i| dpf0.eval(0u8, &BigInt::from(i as u32)))
                    .collect();
                black_box(x);
            })
        });
        group.bench_with_input(BenchmarkId::new("Optimized", &size), &size, |b, size| {
            b.iter(|| {
                let x: Vec<Scalar<Secp256k1>> =
                    dpf0.full_eval_optimized(&false, &BigInt::from(size - 1));
                black_box(x);
            })
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
