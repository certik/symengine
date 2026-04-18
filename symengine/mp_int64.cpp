#include "mp_class.h"
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_exception.h>
#include <symengine/prime_sieve.h>

#if SYMENGINE_INTEGER_CLASS == SYMENGINE_INT64

namespace SymEngine
{

static int64_t i64_gcd(int64_t a, int64_t b)
{
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    while (b) {
        int64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Floor division: quotient rounded toward -infinity
static int64_t i64_fdiv_q(int64_t a, int64_t b)
{
    int64_t q = a / b;
    int64_t r = a % b;
    // If remainder is nonzero and signs differ, subtract 1
    if (r != 0 && ((r ^ b) < 0))
        q -= 1;
    return q;
}

// Floor division remainder
static int64_t i64_fdiv_r(int64_t a, int64_t b)
{
    int64_t r = a % b;
    if (r != 0 && ((r ^ b) < 0))
        r += b;
    return r;
}

// Ceiling division: quotient rounded toward +infinity
static int64_t i64_cdiv_q(int64_t a, int64_t b)
{
    int64_t q = a / b;
    int64_t r = a % b;
    if (r != 0 && ((r ^ b) >= 0))
        q += 1;
    return q;
}

static int64_t i64_pow(int64_t base, unsigned long exp)
{
    int64_t result = 1;
    while (exp > 0) {
        if (exp & 1)
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

void mp_fdiv_qr(integer_class &q, integer_class &r, const integer_class &a,
                const integer_class &b)
{
    q = i64_fdiv_q(a.get(), b.get());
    r = i64_fdiv_r(a.get(), b.get());
}

void mp_cdiv_qr(integer_class &q, integer_class &r, const integer_class &a,
                const integer_class &b)
{
    q = i64_cdiv_q(a.get(), b.get());
    r = a.get() - q.get() * b.get();
}

void mp_pow_ui(rational_class &res, const rational_class &i, unsigned long n)
{
    int64_t num = i64_pow(i.num_.get(), n);
    int64_t den = i64_pow(i.den_.get(), n);
    res = rational_class(num, den);
}

void mp_powm(integer_class &res, const integer_class &base,
             const integer_class &exp, const integer_class &m)
{
    int64_t mod = m.get() < 0 ? -m.get() : m.get();
    int64_t e = exp.get();

    if (e < 0) {
        integer_class base_inv;
        if (!mp_invert(base_inv, base, m)) {
            throw SymEngineException("negative exponent undefined "
                                     "in powm if base is not "
                                     "invertible mod m");
        }
        e = -e;
        int64_t b = i64_fdiv_r(base_inv.get(), mod);
        int64_t result = 1;
        b = i64_fdiv_r(b, mod);
        while (e > 0) {
            if (e & 1)
                result = i64_fdiv_r(result * b, mod);
            b = i64_fdiv_r(b * b, mod);
            e >>= 1;
        }
        res = i64_fdiv_r(result, mod);
        return;
    }

    int64_t b = i64_fdiv_r(base.get(), mod);
    int64_t result = 1;
    while (e > 0) {
        if (e & 1)
            result = i64_fdiv_r(result * b, mod);
        b = i64_fdiv_r(b * b, mod);
        e >>= 1;
    }
    res = i64_fdiv_r(result, mod);
}

void mp_gcdext(integer_class &gcd, integer_class &s, integer_class &t,
               const integer_class &a, const integer_class &b)
{
    int64_t old_r = a.get(), r = b.get();
    int64_t old_s = 1, ss = 0;
    int64_t old_t = 0, tt = 1;

    while (r != 0) {
        int64_t q = old_r / r;
        int64_t tmp = r;
        r = old_r - q * r;
        old_r = tmp;
        tmp = ss;
        ss = old_s - q * ss;
        old_s = tmp;
        tmp = tt;
        tt = old_t - q * tt;
        old_t = tmp;
    }
    if (old_r < 0) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
    }
    gcd = old_r;
    s = old_s;
    t = old_t;
}

bool mp_invert(integer_class &res, const integer_class &a,
               const integer_class &m)
{
    integer_class gcd, s, t;
    mp_gcdext(gcd, s, t, a, m);
    if (gcd.get() != 1) {
        res = 0;
        return false;
    }
    int64_t r = i64_fdiv_r(s.get(), m.get());
    if (r < 0)
        r += (m.get() < 0 ? -m.get() : m.get());
    res = r;
    return true;
}

bool mp_root(integer_class &res, const integer_class &i, unsigned long n)
{
    if (n == 0)
        throw std::runtime_error("0th root is undefined");
    if (n == 1) {
        res = i;
        return true;
    }
    if (i.get() == 0) {
        res = 0;
        return true;
    }

    int64_t val = i.get();
    bool neg = val < 0;
    if (neg && n % 2 == 0)
        throw std::runtime_error("even root of a negative is non-real");
    if (neg)
        val = -val;

    // Newton's method for integer nth root
    int64_t x = static_cast<int64_t>(std::pow(static_cast<double>(val), 1.0 / n));
    // Adjust for floating point imprecision
    if (x > 0) {
        while (i64_pow(x + 1, n) <= val)
            x++;
        while (i64_pow(x, n) > val)
            x--;
    }
    res = neg ? -x : x;
    return i64_pow(x, n) == val;
}

integer_class mp_sqrt(const integer_class &i)
{
    integer_class res;
    mp_root(res, i, 2);
    return res;
}

void mp_rootrem(integer_class &a, integer_class &b, const integer_class &i,
                unsigned long n)
{
    mp_root(a, i, n);
    b = i.get() - i64_pow(a.get(), n);
}

void mp_sqrtrem(integer_class &a, integer_class &b, const integer_class &i)
{
    a = mp_sqrt(i);
    b = i.get() - a.get() * a.get();
}

int mp_probab_prime_p(const integer_class &i, unsigned retries)
{
    int64_t n = i.get();
    if (n < 2)
        return 0;
    if (n < 4)
        return 2; // definitely prime
    if (n % 2 == 0)
        return 0;

    // Trial division for small primes
    static const int small_primes[]
        = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    for (int p : small_primes) {
        if (n == p)
            return 2;
        if (n % p == 0)
            return 0;
    }

    // Miller-Rabin with deterministic witnesses for 64-bit
    // These witnesses are sufficient to determine primality for all n < 3,317,044,064,679,887,385,961,981
    int64_t d = n - 1;
    int r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }

    // For n < 3.3*10^24, witnesses {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} suffice
    static const int64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int64_t a : witnesses) {
        if (a >= n)
            continue;

        // Compute a^d mod n using modular exponentiation
        // Use __int128 for intermediate products to avoid overflow
        __int128 x = 1;
        __int128 base = a;
        int64_t exp = d;
        __int128 mod = n;
        base %= mod;
        while (exp > 0) {
            if (exp & 1)
                x = (x * base) % mod;
            base = (base * base) % mod;
            exp >>= 1;
        }

        int64_t result = static_cast<int64_t>(x);
        if (result == 1 || result == n - 1)
            continue;

        bool composite = true;
        for (int j = 0; j < r - 1; j++) {
            x = (static_cast<__int128>(result) * result) % mod;
            result = static_cast<int64_t>(x);
            if (result == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite)
            return 0;
    }
    return 1; // probably prime
}

void mp_nextprime(integer_class &res, const integer_class &i)
{
    int64_t n = i.get();
    if (n < 2) {
        res = 2;
        return;
    }
    int64_t candidate = (n % 2 == 0) ? n + 1 : n + 2;
    while (!mp_probab_prime_p(integer_class(candidate), 25)) {
        candidate += 2;
    }
    res = candidate;
}

unsigned long mp_scan1(const integer_class &i)
{
    if (i.get() == 0)
        return ULONG_MAX;
    uint64_t v = static_cast<uint64_t>(i.get());
    unsigned long pos = 0;
    while ((v & 1) == 0) {
        v >>= 1;
        pos++;
    }
    return pos;
}

void mp_fib_ui(integer_class &res, unsigned long n)
{
    if (n == 0) {
        res = 0;
        return;
    }
    int64_t a = 0, b = 1;
    for (unsigned long i = 2; i <= n; i++) {
        int64_t t = a + b;
        a = b;
        b = t;
    }
    res = b;
}

void mp_fib2_ui(integer_class &a, integer_class &b, unsigned long n)
{
    if (n == 0) {
        a = 0;
        b = 1; // F(-1) = 1 by convention
        return;
    }
    int64_t prev = 0, curr = 1;
    for (unsigned long i = 2; i <= n; i++) {
        int64_t t = prev + curr;
        prev = curr;
        curr = t;
    }
    a = curr;
    b = prev;
}

void mp_lucnum_ui(integer_class &res, unsigned long n)
{
    if (n == 0) {
        res = 2;
        return;
    }
    if (n == 1) {
        res = 1;
        return;
    }
    int64_t a = 2, b = 1;
    for (unsigned long i = 2; i <= n; i++) {
        int64_t t = a + b;
        a = b;
        b = t;
    }
    res = b;
}

void mp_lucnum2_ui(integer_class &a, integer_class &b, unsigned long n)
{
    if (n == 0)
        throw std::runtime_error("index of lucas number cannot be negative");
    int64_t prev = 2, curr = 1;
    for (unsigned long i = 2; i <= n; i++) {
        int64_t t = prev + curr;
        prev = curr;
        curr = t;
    }
    a = curr;
    b = prev;
}

void mp_fac_ui(integer_class &res, unsigned long n)
{
    int64_t r = 1;
    for (unsigned long i = 2; i <= n; i++)
        r *= static_cast<int64_t>(i);
    res = r;
}

void mp_bin_ui(integer_class &res, const integer_class &n, unsigned long r)
{
    int64_t x = n.get();
    int64_t result = 1;
    for (unsigned long i = 1; i <= r; i++) {
        result *= (x - static_cast<int64_t>(r) + static_cast<int64_t>(i));
        result /= static_cast<int64_t>(i);
    }
    res = result;
}

bool mp_perfect_power_p(const integer_class &i)
{
    int64_t n = i.get();
    if (n == 0 || n == 1 || n == -1)
        return true;

    int64_t absn = n < 0 ? -n : n;
    unsigned long max_exp
        = static_cast<unsigned long>(std::log2(static_cast<double>(absn)));

    // Check if perfect square
    if (n > 0) {
        int64_t s = static_cast<int64_t>(std::sqrt(static_cast<double>(absn)));
        for (int64_t c = std::max(s - 1, (int64_t)1); c <= s + 1; c++) {
            if (c * c == absn)
                return true;
        }
    }

    // Check prime exponents
    integer_class p(2);
    integer_class root;
    while (true) {
        mp_nextprime(p, p);
        if (static_cast<unsigned long>(p.get()) > max_exp)
            return false;
        if (n < 0 && p.get() % 2 == 0)
            continue;
        try {
            if (mp_root(root, i, static_cast<unsigned long>(p.get())))
                return true;
        } catch (...) {
            continue;
        }
    }
}

bool mp_perfect_square_p(const integer_class &i)
{
    if (i.get() < 0)
        return false;
    integer_class root;
    return mp_root(root, i, 2);
}

integer_class mp_primorial(unsigned long n)
{
    int64_t res = 1;
    Sieve::iterator pi(static_cast<unsigned>(n));
    unsigned int p;
    while ((p = pi.next_prime()) <= n) {
        res *= static_cast<int64_t>(p);
    }
    return integer_class(res);
}

// Floored modulus helper
static int64_t fmod_int64(int64_t a, int64_t mod)
{
    int64_t r = a % mod;
    if (r < 0)
        r += mod;
    return r;
}

int mp_legendre(const integer_class &a, const integer_class &n)
{
    integer_class res;
    mp_powm(res, a, integer_class((n.get() - 1) / 2), n);
    return res.get() <= 1 ? static_cast<int>(res.get()) : -1;
}

static int unchecked_jacobi(int64_t a, int64_t n)
{
    if (a == 1)
        return 1;

    a = fmod_int64(a, n);

    unsigned long factors_of_two = 0;
    while (a % 2 == 0 && a != 0) {
        a /= 2;
        ++factors_of_two;
    }

    int product_of_twos = 1;
    int n_mod_8 = static_cast<int>(fmod_int64(n, 8));
    if ((factors_of_two % 2 == 1) && (n_mod_8 == 3 || n_mod_8 == 5))
        product_of_twos = -1;

    if (a == 1)
        return product_of_twos;
    if (i64_gcd(a, n) != 1)
        return 0;

    int qr_factor = 1;
    if ((fmod_int64(a, 4) == 3) && (fmod_int64(n, 4) == 3))
        qr_factor = -1;

    return product_of_twos * qr_factor * unchecked_jacobi(n, a);
}

int mp_jacobi(const integer_class &a, const integer_class &n)
{
    if (n.get() < 0)
        throw std::runtime_error("jacobi denominator must be positive");
    if (n.get() % 2 == 0)
        throw std::runtime_error("jacobi denominator must be odd");
    return unchecked_jacobi(a.get(), n.get());
}

int mp_kronecker(const integer_class &a, const integer_class &n)
{
    int64_t nv = n.get();
    if (nv == 0) {
        int64_t av = a.get();
        return (av == 1 || av == -1) ? 1 : 0;
    }

    int kr_a_u = 1;
    if (nv < 0 && a.get() < 0)
        kr_a_u = -1;

    int64_t m = nv < 0 ? -nv : nv;
    unsigned long j = 0;
    while (m % 2 == 0 && m != 0) {
        m /= 2;
        ++j;
    }

    int kr_a_2_to_j = 1;
    if (j > 0) {
        if (a.get() % 2 == 0) {
            kr_a_2_to_j = 0;
        } else {
            int a_mod_8 = static_cast<int>(fmod_int64(a.get(), 8));
            int kr_a_2 = (a_mod_8 == 1 || a_mod_8 == 7) ? 1 : -1;
            kr_a_2_to_j = (kr_a_2 == -1 && (j % 2 != 0)) ? -1 : 1;
        }
    }

    if (nv % 2 == 0)
        return kr_a_u * kr_a_2_to_j * unchecked_jacobi(a.get(), m);
    else
        return kr_a_u * unchecked_jacobi(a.get(), m);
}

} // namespace SymEngine

#endif // SYMENGINE_INTEGER_CLASS == SYMENGINE_INT64
