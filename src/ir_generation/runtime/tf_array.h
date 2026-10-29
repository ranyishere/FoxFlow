// FoxFlow runtime: fixed-size stack arrays (TFArray<N>).
//
// `FixedList<<N, Float>>` in FoxFlow source compiles to TFArray<N>: N doubles
// stored inline (no heap allocation, no dispatch), with elementwise arithmetic
// and the handful of element/reduction accessors the generated code uses.
// Use `TFTensor` instead when you need a torch::Tensor (deep learning models,
// multi-dimensional data, autograd). The two are never mixed implicitly; use
// tf::to_tensor / tf::to_array to convert.
#ifndef FOXFLOW_TF_ARRAY_H
#define FOXFLOW_TF_ARRAY_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <torch/torch.h>

namespace tf {

// Read-only element value that also offers the tensor-style item<T>() accessor.
struct Scalar {
    double v;
    constexpr Scalar(double x = 0.0) : v(x) {}
    constexpr operator double() const { return v; }
    template <typename T = double>
    constexpr T item() const { return static_cast<T>(v); }
};

// Writable element reference: x[i] = value, x[i] += value, x[i].item<double>().
struct Elem {
    double* p;
    explicit Elem(double* ptr) : p(ptr) {}
    operator double() const { return *p; }
    Elem& operator=(double x) { *p = x; return *this; }
    Elem& operator=(const Elem& o) { *p = *o.p; return *this; }
    Elem& operator+=(double x) { *p += x; return *this; }
    Elem& operator-=(double x) { *p -= x; return *this; }
    Elem& operator*=(double x) { *p *= x; return *this; }
    Elem& operator/=(double x) { *p /= x; return *this; }
    template <typename T = double>
    T item() const { return static_cast<T>(*p); }
};

template <std::size_t N>
struct Array {
    double v[N] = {};

    static constexpr std::size_t size_v = N;
    constexpr std::size_t numel() const { return N; }
    constexpr std::size_t size(int = 0) const { return N; }

    Elem operator[](std::size_t i) { return Elem(&v[i]); }
    Scalar operator[](std::size_t i) const { return Scalar(v[i]); }

    // Value semantics: clone/detach are copies, kept for source compatibility
    // with tensor-style generated code.
    Array clone() const { return *this; }
    Array detach() const { return *this; }

    Scalar sum() const { double s = 0.0; for (std::size_t i = 0; i < N; ++i) s += v[i]; return Scalar(s); }
    Scalar mean() const { return Scalar(sum().v / static_cast<double>(N)); }
    Scalar norm() const { double s = 0.0; for (std::size_t i = 0; i < N; ++i) s += v[i] * v[i]; return Scalar(std::sqrt(s)); }
    Scalar min() const { return Scalar(*std::min_element(v, v + N)); }
    Scalar max() const { return Scalar(*std::max_element(v, v + N)); }

    Array& operator+=(const Array& o) { for (std::size_t i = 0; i < N; ++i) v[i] += o.v[i]; return *this; }
    Array& operator-=(const Array& o) { for (std::size_t i = 0; i < N; ++i) v[i] -= o.v[i]; return *this; }
    Array& operator*=(const Array& o) { for (std::size_t i = 0; i < N; ++i) v[i] *= o.v[i]; return *this; }
    Array& operator/=(const Array& o) { for (std::size_t i = 0; i < N; ++i) v[i] /= o.v[i]; return *this; }
    Array& operator+=(double s) { for (std::size_t i = 0; i < N; ++i) v[i] += s; return *this; }
    Array& operator-=(double s) { for (std::size_t i = 0; i < N; ++i) v[i] -= s; return *this; }
    Array& operator*=(double s) { for (std::size_t i = 0; i < N; ++i) v[i] *= s; return *this; }
    Array& operator/=(double s) { for (std::size_t i = 0; i < N; ++i) v[i] /= s; return *this; }

    template <class Archive>
    void serialize(Archive& archive) {
        for (std::size_t i = 0; i < N; ++i) archive(v[i]);
    }
};

#define TF_ARRAY_BINOP(OP)                                                              \
    template <std::size_t N>                                                            \
    inline Array<N> operator OP(const Array<N>& a, const Array<N>& b) {                 \
        Array<N> r;                                                                     \
        for (std::size_t i = 0; i < N; ++i) r.v[i] = a.v[i] OP b.v[i];                  \
        return r;                                                                       \
    }                                                                                   \
    template <std::size_t N>                                                            \
    inline Array<N> operator OP(const Array<N>& a, double b) {                          \
        Array<N> r;                                                                     \
        for (std::size_t i = 0; i < N; ++i) r.v[i] = a.v[i] OP b;                       \
        return r;                                                                       \
    }                                                                                   \
    template <std::size_t N>                                                            \
    inline Array<N> operator OP(double a, const Array<N>& b) {                          \
        Array<N> r;                                                                     \
        for (std::size_t i = 0; i < N; ++i) r.v[i] = a OP b.v[i];                       \
        return r;                                                                       \
    }

TF_ARRAY_BINOP(+)
TF_ARRAY_BINOP(-)
TF_ARRAY_BINOP(*)
TF_ARRAY_BINOP(/)
#undef TF_ARRAY_BINOP

template <std::size_t N>
inline Array<N> operator-(const Array<N>& a) {
    Array<N> r;
    for (std::size_t i = 0; i < N; ++i) r.v[i] = -a.v[i];
    return r;
}

#define TF_ARRAY_UNARY(NAME, EXPR)                                                      \
    template <std::size_t N>                                                            \
    inline Array<N> NAME(const Array<N>& a) {                                           \
        Array<N> r;                                                                     \
        for (std::size_t i = 0; i < N; ++i) r.v[i] = EXPR(a.v[i]);                      \
        return r;                                                                       \
    }

TF_ARRAY_UNARY(exp, std::exp)
TF_ARRAY_UNARY(log, std::log)
TF_ARRAY_UNARY(sqrt, std::sqrt)
TF_ARRAY_UNARY(abs, std::fabs)
TF_ARRAY_UNARY(sin, std::sin)
TF_ARRAY_UNARY(cos, std::cos)
TF_ARRAY_UNARY(tanh, std::tanh)
#undef TF_ARRAY_UNARY

template <std::size_t N>
inline Scalar dot(const Array<N>& a, const Array<N>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < N; ++i) s += a.v[i] * b.v[i];
    return Scalar(s);
}

template <std::size_t N>
inline Array<N> ones() {
    Array<N> r;
    for (std::size_t i = 0; i < N; ++i) r.v[i] = 1.0;
    return r;
}

template <std::size_t N>
inline Array<N> zeros() { return Array<N>{}; }

// Explicit conversions between arrays and tensors (copies).
template <std::size_t N>
inline torch::Tensor to_tensor(const Array<N>& a) {
    return torch::from_blob(const_cast<double*>(a.v), {static_cast<int64_t>(N)},
                            torch::kFloat64).clone();
}

template <std::size_t N>
inline Array<N> to_array(const torch::Tensor& t) {
    torch::Tensor c = t.to(torch::kFloat64).contiguous();
    TORCH_CHECK(static_cast<std::size_t>(c.numel()) == N,
                "to_array: tensor has ", c.numel(), " elements, expected ", N);
    Array<N> r;
    const double* src = c.data_ptr<double>();
    std::copy(src, src + N, r.v);
    return r;
}

}  // namespace tf

template <std::size_t N>
using TFArray = tf::Array<N>;

#endif
