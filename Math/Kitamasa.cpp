#include <vector>

// 0-index
template<typename T>
class Kitamasa {
    std::vector<T> a, d;
    int k;
public:
    Kitamasa(std::vector<T> _a, std::vector<T> _d, int _k) :a(_a), d(_d), k(_k) {}
private:
    // f_{n} -> f_{n+1}
    std::vector<T> Calc1(std::vector<T> x) {
        std::vector<T> res(k);
        for (int i = 0; i < k; ++i) {
            res[i] = x[k - 1] * d[i] + (i > 0 ? x[i - 1] : 0);
        }
        return res;
    }
    // f_{n} -> f_{2*n}
    std::vector<T> Calc2(std::vector<T> x) {
        std::vector<T> res(k), y = x;
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                res[j] += y[j] * x[i];
            }
            y = Calc1(y);
        }
        return res;
    }
public:
    // return a[n]
    T Solve(long long n) {
        if (n < k) {
            return a[n];
        }
        std::vector<T> f = Rec(n);
        T res = 0;
        for (int i = 0; i < k; ++i) {
            res += f[i] * a[i];
        }
        return res;
    }
private:
    std::vector<T> Rec(long long n) {
        if (n == k) {
            return d;
        } else if (!(n & 1) && n >= (k << 1)) {
            return Calc2(Rec(n >> 1));
        } else {
            return Calc1(Rec(n - 1));
        }
    }
};

/*
 * @see https://scrapbox.io/daikimatsuoka/Kitamasa_Method
 */
