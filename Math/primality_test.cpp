#include <random>
#include <chrono>

std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

// 合成数を素数と判定する確率は k を試行回数として 4^{-k}
bool Miller_Rabin(long long n) {
    if (n == 2) {
        return true;
    }
    if (~n & 1) {
        return false;
    }
    long long d = n - 1;
    while (d & 1) {
        d /= 2;
    }
    int k = 20;
    while (k--) {
        long long a = rng() % (n - 1) + 1;
        if (std::gcd(a, n) != 1) {
            return false;
        }
        long long f = mod_pow(a, d, n);
        if (f == 1 || f == n - 1) {
            continue;
        }
        bool flag = false;
        while (d < n - 1) {
            d *= 2;
            f *= 2;
            f %= n;
            if (f == n - 1) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            return false;
        }
    }
    return true;
}
