#ifndef SYMENGINE_MP_INT64_H
#define SYMENGINE_MP_INT64_H

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <climits>
#include <numeric>
#include <symengine/symengine_rcp.h>

namespace SymEngine
{

class int64_wrapper
{
    int64_t val_;

public:
    int64_wrapper() : val_(0) {}
    int64_wrapper(int i) : val_(i) {}
    int64_wrapper(long i) : val_(i) {}
    int64_wrapper(long long i) : val_(i) {}
    int64_wrapper(unsigned i) : val_(static_cast<int64_t>(i)) {}
    int64_wrapper(unsigned long i) : val_(static_cast<int64_t>(i)) {}
    int64_wrapper(unsigned long long i) : val_(static_cast<int64_t>(i)) {}
    int64_wrapper(double d) : val_(static_cast<int64_t>(d)) {}

    explicit int64_wrapper(const std::string &s)
    {
        val_ = std::stoll(s, nullptr, 0);
    }

    explicit int64_wrapper(const char *s)
    {
        val_ = std::stoll(s, nullptr, 0);
    }

    int64_t get() const
    {
        return val_;
    }

    double get_d() const
    {
        return static_cast<double>(val_);
    }

    unsigned long get_ui() const
    {
        return static_cast<unsigned long>(val_ < 0 ? -val_ : val_);
    }

    long get_si() const
    {
        return static_cast<long>(val_);
    }

    bool fits_ulong_p() const
    {
        return val_ >= 0
               && static_cast<uint64_t>(val_) <= static_cast<uint64_t>(ULONG_MAX);
    }

    bool fits_slong_p() const
    {
        return val_ >= static_cast<int64_t>(LONG_MIN)
               && val_ <= static_cast<int64_t>(LONG_MAX);
    }

    std::string str(int /*unused*/ = 0, std::ios_base::fmtflags /*unused*/ = {}) const
    {
        return std::to_string(val_);
    }

    int sign() const
    {
        return (val_ > 0) - (val_ < 0);
    }

    template <typename T>
    T convert_to() const
    {
        return static_cast<T>(val_);
    }

    // Arithmetic operators
    int64_wrapper operator-() const { return int64_wrapper(-val_); }

    int64_wrapper &operator+=(const int64_wrapper &o) { val_ += o.val_; return *this; }
    int64_wrapper &operator-=(const int64_wrapper &o) { val_ -= o.val_; return *this; }
    int64_wrapper &operator*=(const int64_wrapper &o) { val_ *= o.val_; return *this; }
    int64_wrapper &operator/=(const int64_wrapper &o) { val_ /= o.val_; return *this; }
    int64_wrapper &operator%=(const int64_wrapper &o) { val_ %= o.val_; return *this; }

    int64_wrapper &operator+=(long o) { val_ += o; return *this; }
    int64_wrapper &operator-=(long o) { val_ -= o; return *this; }
    int64_wrapper &operator*=(long o) { val_ *= o; return *this; }
    int64_wrapper &operator/=(long o) { val_ /= o; return *this; }
    int64_wrapper &operator%=(long o) { val_ %= o; return *this; }

    int64_wrapper &operator++() { ++val_; return *this; }
    int64_wrapper operator++(int) { int64_wrapper t(*this); ++val_; return t; }
    int64_wrapper &operator--() { --val_; return *this; }
    int64_wrapper operator--(int) { int64_wrapper t(*this); --val_; return t; }

    // Shift operators
    int64_wrapper operator<<(unsigned long n) const { return int64_wrapper(val_ << n); }
    int64_wrapper operator>>(unsigned long n) const { return int64_wrapper(val_ >> n); }
    int64_wrapper &operator<<=(unsigned long n) { val_ <<= n; return *this; }
    int64_wrapper &operator>>=(unsigned long n) { val_ >>= n; return *this; }

    // Comparison
    friend bool operator==(const int64_wrapper &a, const int64_wrapper &b) { return a.val_ == b.val_; }
    friend bool operator!=(const int64_wrapper &a, const int64_wrapper &b) { return a.val_ != b.val_; }
    friend bool operator<(const int64_wrapper &a, const int64_wrapper &b) { return a.val_ < b.val_; }
    friend bool operator>(const int64_wrapper &a, const int64_wrapper &b) { return a.val_ > b.val_; }
    friend bool operator<=(const int64_wrapper &a, const int64_wrapper &b) { return a.val_ <= b.val_; }
    friend bool operator>=(const int64_wrapper &a, const int64_wrapper &b) { return a.val_ >= b.val_; }

    // Arithmetic binary operators
    friend int64_wrapper operator+(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ + b.val_); }
    friend int64_wrapper operator-(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ - b.val_); }
    friend int64_wrapper operator*(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ * b.val_); }
    friend int64_wrapper operator/(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ / b.val_); }
    friend int64_wrapper operator%(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ % b.val_); }

    // Bitwise
    friend int64_wrapper operator&(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ & b.val_); }
    friend int64_wrapper operator|(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ | b.val_); }
    friend int64_wrapper operator^(const int64_wrapper &a, const int64_wrapper &b) { return int64_wrapper(a.val_ ^ b.val_); }

    // Stream
    friend std::ostream &operator<<(std::ostream &os, const int64_wrapper &i)
    {
        return os << i.val_;
    }

    // Hash support
    friend struct std::hash<int64_wrapper>;
};

class rational_int64
{
public:
    int64_wrapper num_;
    int64_wrapper den_;

private:
    void reduce()
    {
        int64_t g = num_.get();
        int64_t d = den_.get();
        int64_t a = g < 0 ? -g : g;
        int64_t b = d < 0 ? -d : d;
        while (b) {
            int64_t t = b;
            b = a % b;
            a = t;
        }
        if (a > 1) {
            num_ = int64_wrapper(g / a);
            den_ = int64_wrapper(d / a);
        }
        if (den_.get() < 0) {
            num_ = -num_;
            den_ = -den_;
        }
    }

public:
    rational_int64() : num_(0), den_(1) {}
    rational_int64(int64_wrapper n) : num_(n), den_(1) {}
    rational_int64(int64_wrapper n, int64_wrapper d) : num_(n), den_(d)
    {
        reduce();
    }

    // Enable construction from any pair of integer types
    template <typename T, typename U,
              typename = typename std::enable_if<
                  std::is_integral<T>::value && std::is_integral<U>::value>::type>
    rational_int64(T n, U d)
        : num_(static_cast<int64_t>(n)), den_(static_cast<int64_t>(d))
    {
        reduce();
    }

    double get_d() const
    {
        return static_cast<double>(num_.get()) / static_cast<double>(den_.get());
    }

    int sign() const
    {
        return num_.sign();
    }

    // Arithmetic
    rational_int64 &operator+=(const rational_int64 &o)
    {
        num_ = num_ * o.den_ + o.num_ * den_;
        den_ = den_ * o.den_;
        reduce();
        return *this;
    }

    rational_int64 &operator+=(const int64_wrapper &o)
    {
        num_ = num_ + o * den_;
        reduce();
        return *this;
    }

    rational_int64 &operator-=(const rational_int64 &o)
    {
        num_ = num_ * o.den_ - o.num_ * den_;
        den_ = den_ * o.den_;
        reduce();
        return *this;
    }

    rational_int64 &operator*=(const rational_int64 &o)
    {
        num_ = num_ * o.num_;
        den_ = den_ * o.den_;
        reduce();
        return *this;
    }

    rational_int64 &operator*=(const int64_wrapper &o)
    {
        num_ = num_ * o;
        reduce();
        return *this;
    }

    rational_int64 &operator/=(const rational_int64 &o)
    {
        num_ = num_ * o.den_;
        den_ = den_ * o.num_;
        reduce();
        return *this;
    }

    rational_int64 &operator/=(const int64_wrapper &o)
    {
        den_ = den_ * o;
        reduce();
        return *this;
    }

    friend rational_int64 operator+(const rational_int64 &a, const rational_int64 &b)
    {
        return rational_int64(a.num_ * b.den_ + b.num_ * a.den_, a.den_ * b.den_);
    }

    friend rational_int64 operator+(const rational_int64 &a, const int64_wrapper &b)
    {
        return rational_int64(a.num_ + b * a.den_, a.den_);
    }

    friend rational_int64 operator+(const int64_wrapper &a, const rational_int64 &b)
    {
        return b + a;
    }

    friend rational_int64 operator-(const rational_int64 &a, const rational_int64 &b)
    {
        return rational_int64(a.num_ * b.den_ - b.num_ * a.den_, a.den_ * b.den_);
    }

    friend rational_int64 operator-(const rational_int64 &a, const int64_wrapper &b)
    {
        return rational_int64(a.num_ - b * a.den_, a.den_);
    }

    friend rational_int64 operator*(const rational_int64 &a, const rational_int64 &b)
    {
        return rational_int64(a.num_ * b.num_, a.den_ * b.den_);
    }

    friend rational_int64 operator*(const rational_int64 &a, const int64_wrapper &b)
    {
        return rational_int64(a.num_ * b, a.den_);
    }

    friend rational_int64 operator*(const int64_wrapper &a, const rational_int64 &b)
    {
        return b * a;
    }

    friend rational_int64 operator/(const rational_int64 &a, const rational_int64 &b)
    {
        return rational_int64(a.num_ * b.den_, a.den_ * b.num_);
    }

    friend rational_int64 operator/(const rational_int64 &a, const int64_wrapper &b)
    {
        return rational_int64(a.num_, a.den_ * b);
    }

    friend rational_int64 operator/(const int64_wrapper &a, const rational_int64 &b)
    {
        return rational_int64(a * b.den_, b.num_);
    }

    rational_int64 operator-() const
    {
        return rational_int64(-num_, den_);
    }

    // Comparison (rational vs rational)
    friend bool operator==(const rational_int64 &a, const rational_int64 &b)
    {
        return a.num_ * b.den_ == b.num_ * a.den_;
    }
    friend bool operator!=(const rational_int64 &a, const rational_int64 &b)
    {
        return !(a == b);
    }
    friend bool operator<(const rational_int64 &a, const rational_int64 &b)
    {
        return a.num_ * b.den_ < b.num_ * a.den_;
    }
    friend bool operator>(const rational_int64 &a, const rational_int64 &b)
    {
        return b < a;
    }
    friend bool operator<=(const rational_int64 &a, const rational_int64 &b)
    {
        return !(b < a);
    }
    friend bool operator>=(const rational_int64 &a, const rational_int64 &b)
    {
        return !(a < b);
    }

    // Comparison (rational vs int64_wrapper, handles all integer types)
    friend bool operator==(const rational_int64 &a, const int64_wrapper &b)
    {
        return a.den_.get() == 1 && a.num_ == b;
    }
    friend bool operator==(const int64_wrapper &a, const rational_int64 &b)
    {
        return b == a;
    }
    friend bool operator!=(const rational_int64 &a, const int64_wrapper &b)
    {
        return !(a == b);
    }
    friend bool operator!=(const int64_wrapper &a, const rational_int64 &b)
    {
        return !(b == a);
    }
    friend bool operator<(const rational_int64 &a, const int64_wrapper &b)
    {
        return a.num_ < b * a.den_;
    }
    friend bool operator<(const int64_wrapper &a, const rational_int64 &b)
    {
        return a * b.den_ < b.num_;
    }
    friend bool operator>(const rational_int64 &a, const int64_wrapper &b)
    {
        return b < a;
    }
    friend bool operator>(const int64_wrapper &a, const rational_int64 &b)
    {
        return b < a;
    }
    friend bool operator<=(const rational_int64 &a, const int64_wrapper &b)
    {
        return !(b < a);
    }
    friend bool operator<=(const int64_wrapper &a, const rational_int64 &b)
    {
        return !(b < a);
    }
    friend bool operator>=(const rational_int64 &a, const int64_wrapper &b)
    {
        return !(a < b);
    }
    friend bool operator>=(const int64_wrapper &a, const rational_int64 &b)
    {
        return !(a < b);
    }

    // Comparison with double
    friend bool operator==(const rational_int64 &a, double b)
    {
        return a.get_d() == b;
    }
    friend bool operator==(double a, const rational_int64 &b)
    {
        return b.get_d() == a;
    }
    friend bool operator!=(const rational_int64 &a, double b)
    {
        return !(a == b);
    }

    // Stream
    friend std::ostream &operator<<(std::ostream &os, const rational_int64 &r)
    {
        if (r.den_ == 1)
            return os << r.num_;
        return os << r.num_ << "/" << r.den_;
    }
};

} // namespace SymEngine

namespace std
{
template <>
struct hash<SymEngine::int64_wrapper> {
    size_t operator()(const SymEngine::int64_wrapper &i) const
    {
        return std::hash<int64_t>()(i.val_);
    }
};
} // namespace std

#endif // SYMENGINE_MP_INT64_H
