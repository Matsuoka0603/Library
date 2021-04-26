template<int mod>
class Modular {
private:
    int x;
public:
    Modular() : x(0) {}
    Modular(long long y) : x(y >= 0 ? y % mod : (mod - (-y) % mod) % mod) {}
    Modular &operator+=(const Modular& p) { if ((x += p.x) >= mod) x -= mod; return *this; }
    Modular &operator-=(const Modular& p) { if ((x -= p.x) < 0) x += mod; return *this; }
    Modular &operator*=(const Modular& p) { x = (int) (1LL * x * p.x % mod); return *this; }
    Modular &operator/=(const Modular& p) { *this *= p.Inverse(); return *this; }
    Modular operator-() const { return Modular(-x); }
    Modular operator+(const Modular& p) { return Modular(*this) += p; }
    Modular operator-(const Modular& p) { return Modular(*this) -= p; }
    Modular operator*(const Modular& p) { return Modular(*this) *= p; }
    Modular operator/(const Modular& p) { return Modular(*this) /= p; }
    bool operator==(const Modular& p) { return x == p.x; }
    bool operator!=(const Modular& p) { return x != p.x; }
    friend std::ostream& operator<<(std::ostream& out, const Modular& p) { return out << p.x; }
    friend std::istream& operator>>(std::istream& in, Modular& a) {
        long long t;
        in >> t;
        a = Modular<mod>(t);
        return in;
    }
    Modular Inverse() const {
        int a = x, m = mod, u = 1, v = 0, t;
        while (m > 0) {
            t = a / m;
            std::swap(a -= t * m, m);
            std::swap(u -= t * v, v);
        }
        return Modular(u);
    }
};
