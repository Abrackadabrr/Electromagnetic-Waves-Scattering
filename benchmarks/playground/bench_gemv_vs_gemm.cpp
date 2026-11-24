#include <benchmark/benchmark.h>
#include <eigen3/Eigen/Dense>

constexpr int N = 512;

// K последовательных GEMV
static void BM_Matvec(benchmark::State &state) {
    const int K = static_cast<int>(state.range(0));

    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
    const Eigen::VectorXd B = Eigen::VectorXd::Random(K * N);

    Eigen::VectorXd res(K * N);

    for (auto _ : state) {
        for (int i = 0; i < K; ++i) {
            res.segment(i * N, N).noalias() = A * B.segment(i * N, N);
        }

        benchmark::DoNotOptimize(res.data());
        benchmark::ClobberMemory();
    }
}

// Одна GEMM через Eigen::Map
static void BM_Map(benchmark::State &state) {
    const int K = static_cast<int>(state.range(0));

    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
    const Eigen::VectorXd B = Eigen::VectorXd::Random(K * N);

    Eigen::VectorXd res(K * N);

    const Eigen::Map<const Eigen::MatrixXd> B_reshaped(B.data(), N, K);

    Eigen::Map<Eigen::MatrixXd> res_reshaped(res.data(), N, K);

    for (auto _ : state) {
        res_reshaped.noalias() = A * B_reshaped;

        benchmark::DoNotOptimize(res.data());
        benchmark::ClobberMemory();
    }
}

static void BenchmarkArguments(::benchmark::Benchmark *b) {
    for (int K : {1, 2, 3, 4, 5, 7, 8, 9, 10, 20, 30, 40, 50}) {
        b->Arg(K);
    }
}

BENCHMARK(BM_Matvec)->Apply(BenchmarkArguments)->Unit(benchmark::kMicrosecond);

BENCHMARK(BM_Map)->Apply(BenchmarkArguments)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();