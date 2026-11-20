
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <thread>
#include <mutex>
#include <atomic>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

using Clock = chrono::high_resolution_clock;

cpp_int mod(cpp_int a, const cpp_int& m) {
    a %= m;
    if (a < 0) a += m;
    return a;
}

cpp_int gcdBig(cpp_int a, cpp_int b) {
    while (b != 0) {
        cpp_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

cpp_int modPow(cpp_int base, cpp_int exp, cpp_int m) {
    base = mod(base, m);
    cpp_int result = 1;

    while (exp > 0) {
        if ((exp & 1) != 0)
            result = mod(result * base, m);

        base = mod(base * base, m);
        exp >>= 1;
    }

    return result;
}

cpp_int modInverse(cpp_int a, cpp_int m) {
    cpp_int t = 0, newT = 1;
    cpp_int r = m, newR = mod(a, m);

    while (newR != 0) {
        cpp_int q = r / newR;

        cpp_int tempT = t - q * newT;
        t = newT;
        newT = tempT;

        cpp_int tempR = r - q * newR;
        r = newR;
        newR = tempR;
    }

    if (r != 1)
        throw runtime_error("Modular inverse does not exist.");

    return mod(t, m);
}

bool isPrime(long long n) {
    if (n < 2) return false;

    for (long long i = 2; i * i <= n; ++i)
        if (n % i == 0)
            return false;

    return true;
}

vector<long long> factorBase(long long B) {
    vector<long long> primes;

    for (long long i = 2; i <= B; ++i)
        if (isPrime(i))
            primes.push_back(i);

    return primes;
}

bool isGenerator(cpp_int alpha, cpp_int p) {
    cpp_int order = p - 1;
    cpp_int n = order;

    vector<cpp_int> factors;

    for (cpp_int q = 2; q * q <= n; ++q) {
        if (n % q == 0) {
            factors.push_back(q);

            while (n % q == 0)
                n /= q;
        }
    }

    if (n > 1)
        factors.push_back(n);

    for (const cpp_int& q : factors) {
        if (modPow(alpha, order / q, p) == 1)
            return false;
    }

    return modPow(alpha, order, p) == 1;
}

bool factorSmooth(cpp_int value, const vector<long long>& base,
                 vector<int>& exponents) {
    exponents.assign(base.size(), 0);

    for (size_t i = 0; i < base.size(); ++i) {
        while (value % base[i] == 0) {
            value /= base[i];
            exponents[i]++;
        }
    }

    return value == 1;
}

struct Relation {
    cpp_int exponent;
    vector<int> exponents;
};

bool solvePrimeField(const vector<vector<cpp_int>>& A,
                     const vector<cpp_int>& b,
                     cpp_int q,
                     vector<cpp_int>& solution) {
    if (A.empty())
        return false;

    int rows = (int)A.size();
    int cols = (int)A[0].size();

    vector<vector<cpp_int>> mat(rows, vector<cpp_int>(cols + 1));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j)
            mat[i][j] = mod(A[i][j], q);

        mat[i][cols] = mod(b[i], q);
    }

    int pivotRow = 0;
    vector<int> pivotCol;

    for (int col = 0; col < cols && pivotRow < rows; ++col) {
        int selected = -1;

        for (int r = pivotRow; r < rows; ++r) {
            if (mat[r][col] != 0) {
                selected = r;
                break;
            }
        }

        if (selected == -1)
            continue;

        swap(mat[pivotRow], mat[selected]);

        cpp_int inv = modInverse(mat[pivotRow][col], q);

        for (int j = col; j <= cols; ++j)
            mat[pivotRow][j] = mod(mat[pivotRow][j] * inv, q);

        for (int r = 0; r < rows; ++r) {
            if (r == pivotRow)
                continue;

            cpp_int factor = mat[r][col];

            for (int j = col; j <= cols; ++j)
                mat[r][j] = mod(mat[r][j] -
                                  factor * mat[pivotRow][j], q);
        }

        pivotCol.push_back(col);
        pivotRow++;
    }

    for (int r = pivotRow; r < rows; ++r) {
        bool allZero = true;

        for (int c = 0; c < cols; ++c)
            if (mat[r][c] != 0)
                allZero = false;

        if (allZero && mat[r][cols] != 0)
            return false;
    }

    solution.assign(cols, 0);

    for (int r = 0; r < (int)pivotCol.size(); ++r)
        solution[pivotCol[r]] = mat[r][cols];

    return true;
}

vector<cpp_int> primeFactorsOfModulus(cpp_int m) {
    vector<cpp_int> factors;

    for (cpp_int i = 2; i * i <= m; ++i) {
        if (m % i == 0) {
            factors.push_back(i);

            while (m % i == 0)
                m /= i;
        }
    }

    if (m > 1)
        factors.push_back(m);

    return factors;
}

vector<pair<cpp_int, int>> primePowerFactors(cpp_int m) {
    vector<pair<cpp_int, int>> result;

    for (cpp_int q = 2; q * q <= m; ++q) {
        if (m % q == 0) {
            int exponent = 0;

            while (m % q == 0) {
                m /= q;
                exponent++;
            }

            result.push_back({q, exponent});
        }
    }

    if (m > 1)
        result.push_back({m, 1});

    return result;
}

bool solvePrimePower(const vector<vector<cpp_int>>& A,
                     const vector<cpp_int>& b,
                     cpp_int q,
                     int power,
                     vector<cpp_int>& solution) {
    cpp_int modulus = 1;

    for (int i = 0; i < power; ++i)
        modulus *= q;

    vector<cpp_int> x;

    if (!solvePrimeField(A, b, q, x))
        return false;

    cpp_int currentModulus = q;

    for (int k = 1; k < power; ++k) {
        vector<cpp_int> correctionRHS(b.size());

        for (size_t i = 0; i < A.size(); ++i) {
            cpp_int residual = b[i];

            for (size_t j = 0; j < x.size(); ++j)
                residual -= A[i][j] * x[j];

            residual = mod(residual, currentModulus * q);
            residual /= currentModulus;

            correctionRHS[i] = mod(residual, q);
        }

        vector<cpp_int> delta;

        if (!solvePrimeField(A, correctionRHS, q, delta))
            return false;

        for (size_t j = 0; j < x.size(); ++j)
            x[j] += currentModulus * delta[j];

        currentModulus *= q;
    }

    for (size_t i = 0; i < A.size(); ++i) {
        cpp_int value = 0;

        for (size_t j = 0; j < x.size(); ++j)
            value += A[i][j] * x[j];

        if (mod(value - b[i], modulus) != 0)
            return false;
    }

    solution.resize(x.size());

    for (size_t j = 0; j < x.size(); ++j)
        solution[j] = mod(x[j], modulus);

    return true;
}

bool solveRelationSystem(const vector<Relation>& relations,
    int baseSize,
    cpp_int modulus,
    vector<cpp_int>& logs) {
if ((int)relations.size() < baseSize)
return false;

vector<vector<cpp_int>> A;
vector<cpp_int> b;

for (const auto& relation : relations) {
vector<cpp_int> row(baseSize);

for (int j = 0; j < baseSize; ++j)
row[j] = relation.exponents[j];

A.push_back(row);
b.push_back(relation.exponent);
}

auto factors = primePowerFactors(modulus);

vector<cpp_int> result(baseSize, 0);
cpp_int accumulatedModulus = 1;

for (const auto& factor : factors) {
cpp_int q = factor.first;
int power = factor.second;

vector<cpp_int> partial;

if (!solvePrimePower(A, b, q, power, partial))
return false;

cpp_int primePower = 1;

for (int i = 0; i < power; ++i)
primePower *= q;

cpp_int inverse = modInverse(accumulatedModulus, primePower);

for (int j = 0; j < baseSize; ++j) {
cpp_int delta = mod(
(partial[j] - result[j]) * inverse,
primePower
);

result[j] += accumulatedModulus * delta;
result[j] = mod(result[j], accumulatedModulus * primePower);
}

accumulatedModulus *= primePower;
}

logs = result;

return true;
}
void collectRelationsParallel(
    cpp_int p,
    cpp_int alpha,
    const vector<long long>& base,
    int requiredRelations,
    vector<Relation>& relations,
    long long maxAttempts,
    int threadCount
) {
    mutex relationMutex;
    atomic<long long> attempts(0);
    atomic<bool> stop(false);

    cpp_int order = p - 1;

    vector<thread> workers;

    for (int t = 0; t < threadCount; ++t) {
        workers.emplace_back([&, t]() {
            random_device rd;
            mt19937_64 rng(rd() + t);

            while (!stop.load()) {
                long long attempt = attempts.fetch_add(1);

                if (attempt >= maxAttempts) {
                    stop.store(true);
                    break;
                }

                cpp_int exponent = rng() % order;
                cpp_int value = modPow(alpha, exponent, p);

                vector<int> exponents;

                if (!factorSmooth(value, base, exponents))
                    continue;

                Relation relation;
                relation.exponent = exponent;
                relation.exponents = exponents;

                lock_guard<mutex> lock(relationMutex);

                if (stop.load())
                    continue;

                relations.push_back(relation);

                if ((int)relations.size() >= requiredRelations)
                    stop.store(true);
            }
        });
    }

    for (auto& worker : workers)
        worker.join();
}

bool findTargetRelation(cpp_int p,
                       cpp_int alpha,
                       cpp_int beta,
                       const vector<long long>& base,
                       const vector<cpp_int>& logs,
                       cpp_int& result,
                       long long maxAttempts) {
    cpp_int order = p - 1;

    random_device rd;
    mt19937_64 rng(rd());

    for (long long attempt = 0; attempt < maxAttempts; ++attempt) {
        cpp_int r = rng() % order;
        cpp_int value = mod(beta * modPow(alpha, r, p), p);

        vector<int> exponents;

        if (!factorSmooth(value, base, exponents))
            continue;

            cpp_int x = -r;

            for (size_t i = 0; i < logs.size(); ++i)
                x += cpp_int(exponents[i]) * logs[i];
            
            result = mod(x, order);

        if (modPow(alpha, result, p) == beta)
            return true;
    }

    return false;
}

int main() {
    try {
        cpp_int p, alpha, beta;
        long long B;
        int threadCount;

        cout << "Enter prime p: ";
        cin >> p;

        cout << "Enter generator alpha: ";
        cin >> alpha;

        cout << "Enter beta: ";
        cin >> beta;

        cout << "Enter factor base bound B: ";
        cin >> B;

        cout << "Enter number of threads: ";
        cin >> threadCount;

        if (p <= 2 || p > 1000000000LL) {
            cout << "Unsupported p. Use 2 < p <= 10^9.\n";
            return 1;
        }

        if (!isPrime(p.convert_to<long long>())) {
            cout << "p must be prime.\n";
            return 1;
        }

        alpha = mod(alpha, p);
        beta = mod(beta, p);

        if (alpha == 0 || beta == 0) {
            cout << "alpha and beta must be nonzero modulo p.\n";
            return 1;
        }

        if (B < 2 || B >= p) {
            cout << "Invalid factor base bound B.\n";
            return 1;
        }

        if (threadCount < 1) {
            cout << "Number of threads must be positive.\n";
            return 1;
        }

        if (!isGenerator(alpha, p)) {
            cout << "alpha must be a generator of Z_p*.\n";
            return 1;
        }

        vector<long long> base = factorBase(B);

        if (base.empty()) {
            cout << "Factor base is empty.\n";
            return 1;
        }

        cout << "\nFactor base: ";

        for (auto prime : base)
            cout << prime << " ";

        cout << "\n";

        if (threadCount > (int)thread::hardware_concurrency() && thread::hardware_concurrency() > 0) {
            threadCount = (int)thread::hardware_concurrency();
            cout << "Thread count limited to available hardware threads: " << threadCount << "\n";
        }

        int requiredRelations = (int)base.size() + 20;
        long long maxAttempts = 2000000;
        vector<Relation> relations;
        vector<cpp_int> logs;
        bool solved = false;

        auto totalStart = Clock::now();
        auto relationStart = Clock::now();
        double relationAccumulated = 0.0;
        double solvingAccumulated = 0.0;

        cout << "\nCollecting relations in parallel...\n";

        for (int round = 0; round < 5 && !solved; ++round) {
            int targetRelations = requiredRelations + round * 20;
            cout << "Collection round " << (round + 1)
                 << ": target " << targetRelations << " relations.\n";

            auto collectRoundStart = Clock::now();
            collectRelationsParallel(
                p, alpha, base, targetRelations,
                relations, maxAttempts, threadCount
            );
            relationAccumulated += chrono::duration<double>(Clock::now() - collectRoundStart).count();

            cout << "Collected " << relations.size() << " relations.\n";

            if ((int)relations.size() < (int)base.size())
                continue;

            cout << "Solving relation system...\n";
            auto solveRoundStart = Clock::now();
            bool solveResult = solveRelationSystem(relations, (int)base.size(), p - 1, logs);
            solvingAccumulated += chrono::duration<double>(Clock::now() - solveRoundStart).count();
            if (!solveResult)
                continue;

            solved = true;

            for (size_t i = 0; i < base.size(); ++i) {
                if (modPow(alpha, logs[i], p) != base[i]) {
                    solved = false;
                    break;
                }
            }
        }

        auto relationEnd = Clock::now();
        double relationTimeMeasured = relationAccumulated;
        double solvingTimeMeasured = solvingAccumulated;

        if (!solved) {
            cout << "Unable to obtain a verified solution after multiple collection rounds.\n";
            return 1;
        }

        cout << "\nFactor base logarithms:\n";

        for (size_t i = 0; i < base.size(); ++i)
            cout << "log_alpha(" << base[i] << ") = "
                 << logs[i] << "\n";

        cout << "\nSearching for target relation...\n";

        auto targetStart = Clock::now();

        cpp_int result;

        bool found = findTargetRelation(
            p, alpha, beta, base, logs,
            result, maxAttempts
        );

        auto targetEnd = Clock::now();

        if (!found) {
            cout << "Target relation was not found.\n";
            return 1;
        }

        cout << "\nDiscrete logarithm:\n";
        cout << "x = " << result << "\n";
        cout << "Verification: alpha^x mod p = beta.\n";

        auto totalEnd = Clock::now();

        auto relationTime = relationTimeMeasured;
        auto solvingTime = solvingTimeMeasured;

        auto targetTime =
            chrono::duration<double>(targetEnd - targetStart).count();

        auto totalTime =
            chrono::duration<double>(totalEnd - totalStart).count();

        cout << "\nExecution time:\n";
        cout << "Relation collection: " << relationTime << " s\n";
        cout << "System solving: " << solvingTime << " s\n";
        cout << "Target relation: " << targetTime << " s\n";
        cout << "Total: " << totalTime << " s\n";

    } catch (const exception& e) {
        cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
 
