template <typename T>
T inverse(T a, T m) {
    T u = 0, v = 1;
    while (a != 0) {
        T t = m / a;
        m -= t * a; swap(a, m);
        u -= t * v; swap(u, v);
    }
    assert(m == 1);
    return u;
}
 
template<int mod>
class Modular {
private:
    int a;
public:
    constexpr Modular() : a(0) {}
    constexpr Modular(long long b) : a(b >= 0 ? b % mod : ((b % mod) + mod) % mod) {}
    constexpr Modular& operator+=(const Modular& rhs) noexcept { if ((a += rhs.a) >= mod) a -= mod; return *this; }
    constexpr Modular& operator-=(const Modular& rhs) noexcept { if ((a -= rhs.a) < 0) a += mod; return *this; }
    constexpr Modular& operator*=(const Modular& rhs) noexcept { a = (int) (1LL * a * rhs.a % mod); return *this; }
    constexpr Modular& operator/=(const Modular& rhs) noexcept { return *this *= Modular(inverse(rhs.a, mod)); }
    constexpr Modular operator+(const Modular& rhs) const noexcept { return Modular(*this) += rhs; }
    constexpr Modular operator-(const Modular& rhs) const noexcept { return Modular(*this) -= rhs; }
    constexpr Modular operator*(const Modular& rhs) const noexcept { return Modular(*this) *= rhs; }
    constexpr Modular operator/(const Modular& rhs) const noexcept { return Modular(*this) /= rhs; }
    friend std::ostream& operator<<(std::ostream& out, const Modular& rhs) { return out << rhs.a; }
};
