long long inv(long long a, long long m) { return a == 1 ? a : m + (-m * inv(m % a, a) + 1) / a; }

template<int mod>
class Modular {
private:
    int a;
public:
    Modular() : a(0) {}
    Modular(long long b) : a(b >= 0 ? b % mod : ((b % mod) + mod) % mod) {}
    Modular& operator+=(Modular& rhs) { if ((a += rhs.a) >= mod) a -= mod; return *this; }
    Modular& operator-=(Modular& rhs) { if ((a -= rhs.a) < 0) a += mod; return *this; }
    Modular& operator*=(Modular& rhs) { a = (int) (1LL * a * rhs.a % mod); return *this; }
    Modular& operator/=(Modular& rhs) { return *this *= Modular(inv(rhs.a, mod)); }
    Modular operator+(Modular& rhs) { return Modular(*this) += rhs; }
    Modular operator-(Modular& rhs) { return Modular(*this) -= rhs; }
    Modular operator*(Modular& rhs) { return Modular(*this) *= rhs; }
    Modular operator/(Modular& rhs) { return Modular(*this) /= rhs; }
    friend std::ostream& operator<<(std::ostream& out, const Modular& rhs) { return out << rhs.a; }
};
