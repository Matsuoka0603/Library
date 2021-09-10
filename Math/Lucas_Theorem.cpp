/**
 * binom(n, r) mod p を Lucas's Theorem により計算する
 * 空間計算量 : O(p^2)
 * 時間計算量 : 前計算 O(p^2), クエリ O(log n)
 */

#include <vector>
#include <cassert>

class Lucas_Theorem {
    int p;
    std::vector<std::vector<int>> c;
public:
    Lucas_Theorem(int _p) : p(_p) {
        assert(2 <= p && p <= 10000);
        // p 以下の二項係数テーブルを作る
        c.resize(p + 1, std::vector<int> (p + 1));
        c[0][0] = 1;
        for (int i = 0; i <= p; ++i) {
            for (int j = 0; j <= i; ++j) {
                if (i - 1 >= 0 && j - 1 >= 0) {
                    c[i][j] += c[i - 1][j - 1];
                }
                if (i - 1 >= 0) {
                    c[i][j] += c[i - 1][j];
                }
            }
        }
    }
    // binom(n, r) mod p を計算
    int solve(long long n, long long r) {
        if (n < r) {
            return 0;
        }
        int res = 1;
        while (n > 0 && r > 0) {
            res *= c[n % p][r % p];
            res %= p;
            n /= p;
            r /= p;
        }
        return res;
    }
};
