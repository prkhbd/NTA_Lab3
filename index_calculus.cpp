
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <set>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

cpp_int mod(cpp_int a, const cpp_int& m) {
    a %= m;
    if (a < 0)
        a += m;
    return a;
}

cpp_int modPow(cpp_int a, cpp_int e, const cpp_int& m) {
    cpp_int result = 1;
    a = mod(a, m);

    while (e > 0) {
        if ((e & 1) != 0)
            result = mod(result * a, m);
        a = mod(a * a, m);
        e >>= 1;
    }

    return result;
}

long long modLL(long long a, long long m) {
    a %= m;
    if (a < 0)
        a += m;
    return a;
}

long long extendedGCD(long long a, long long b, long long& x, long long& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }

    long long x1, y1;
    long long gcd = extendedGCD(b, a % b, x1, y1);

    x = y1;
    y = x1 - (a / b) * y1;

    return gcd;
}

bool inverseMod(long long a, long long m, long long& inverse) {
    long long x, y;
    long long gcd = extendedGCD(a, m, x, y);

    if (gcd != 1)
        return false;

    inverse = modLL(x, m);
    return true;
}

bool isPrime(long long n) {
    if (n < 2)
        return false;

    if (n % 2 == 0)
        return n == 2;

    for (long long i = 3; i <= n / i; i += 2) {
        if (n % i == 0)
            return false;
    }

    return true;
}

vector<long long> primeFactors(long long n) {
    vector<long long> factors;

    for (long long i = 2; i <= n / i; i++) {
        if (n % i == 0) {
            factors.push_back(i);
            while (n % i == 0)
                n /= i;
        }
    }

    if (n > 1)
        factors.push_back(n);

    return factors;
}

vector<pair<long long, long long>> primePowerFactors(long long n) {
    vector<pair<long long, long long>> result;

    for (long long q = 2; q <= n / q; q++) {
        if (n % q == 0) {
            long long power = 1;

            while (n % q == 0) {
                n /= q;
                power *= q;
            }

            result.push_back({q, power});
        }
    }

    if (n > 1)
        result.push_back({n, n});

    return result;
}

bool isPrimitiveRoot(long long alpha, long long p) {
    if (!isPrime(p))
        return false;

    long long n = p - 1;
    alpha = modLL(alpha, p);

    if (alpha == 0)
        return false;

    vector<long long> factors = primeFactors(n);

    for (long long q : factors) {
        if (modPow(alpha, n / q, p) == 1)
            return false;
    }

    return true;
}

vector<long long> generateFactorBase(long long B, long long p) {
    vector<long long> factorBase;

    for (long long i = 2; i <= B; i++) {
        if (i < p && isPrime(i))
            factorBase.push_back(i);
    }

    return factorBase;
}

struct SmoothResult {
    bool smooth;
    vector<long long> exponents;
};

SmoothResult checkSmoothness(
    cpp_int value,
    const vector<long long>& factorBase
) {
    SmoothResult result;
    result.smooth = false;
    result.exponents.assign(factorBase.size(), 0);

    if (value <= 0)
        return result;

    for (size_t i = 0; i < factorBase.size(); i++) {
        long long prime = factorBase[i];

        while (value % prime == 0) {
            value /= prime;
            result.exponents[i]++;
        }
    }

    result.smooth = (value == 1);
    return result;
}

struct Relation {
    long long k;
    vector<long long> exponents;
};

vector<Relation> collectRelations(
    long long p,
    long long alpha,
    long long n,
    const vector<long long>& factorBase,
    int requiredRelations,
    long long maxAttempts
) {
    vector<Relation> relations;
    set<long long> usedK;

    mt19937_64 generator(
        chrono::steady_clock::now().time_since_epoch().count()
    );

    uniform_int_distribution<long long> distribution(0, n - 1);
    long long attempts = 0;

    while ((int)relations.size() < requiredRelations &&
           attempts < maxAttempts) {
        attempts++;

        long long k = distribution(generator);

        if (!usedK.insert(k).second)
            continue;

        cpp_int value = modPow(alpha, k, p);
        SmoothResult result = checkSmoothness(value, factorBase);

        if (result.smooth) {
            Relation relation;
            relation.k = k;
            relation.exponents = result.exponents;
            relations.push_back(relation);
        }
    }

    return relations;
}

bool solvePrimeField(
    const vector<vector<long long>>& A,
    const vector<long long>& b,
    long long q,
    int variables,
    vector<long long>& solution
) {
    int rows = (int)A.size();

    vector<vector<long long>> matrix(
        rows,
        vector<long long>(variables + 1, 0)
    );

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < variables; j++)
            matrix[i][j] = modLL(A[i][j], q);

        matrix[i][variables] = modLL(b[i], q);
    }

    int pivotRow = 0;
    vector<int> pivotColumns;

    for (int col = 0; col < variables; col++) {
        int selected = -1;

        for (int row = pivotRow; row < rows; row++) {
            if (matrix[row][col] != 0) {
                selected = row;
                break;
            }
        }

        if (selected == -1)
            continue;

        swap(matrix[pivotRow], matrix[selected]);

        long long inverse;

        if (!inverseMod(matrix[pivotRow][col], q, inverse))
            return false;

        for (int j = col; j <= variables; j++)
            matrix[pivotRow][j] =
                modLL(matrix[pivotRow][j] * inverse, q);

        for (int row = 0; row < rows; row++) {
            if (row == pivotRow)
                continue;

            long long coefficient = matrix[row][col];

            for (int j = col; j <= variables; j++) {
                matrix[row][j] = modLL(
                    matrix[row][j] -
                    coefficient * matrix[pivotRow][j],
                    q
                );
            }
        }

        pivotColumns.push_back(col);
        pivotRow++;
    }

    for (int row = pivotRow; row < rows; row++) {
        bool zeroRow = true;

        for (int col = 0; col < variables; col++) {
            if (matrix[row][col] != 0) {
                zeroRow = false;
                break;
            }
        }

        if (zeroRow && matrix[row][variables] != 0)
            return false;
    }

    if ((int)pivotColumns.size() != variables)
        return false;

    solution.assign(variables, 0);

    for (int row = 0; row < variables; row++)
        solution[pivotColumns[row]] = matrix[row][variables];

    return true;
}

bool solvePrimePower(
    const vector<vector<long long>>& A,
    const vector<long long>& b,
    long long prime,
    long long modulus,
    int variables,
    vector<long long>& solution
) {
    if (!solvePrimeField(A, b, prime, variables, solution))
        return false;

    long long currentModulus = prime;

    while (currentModulus < modulus) {
        vector<long long> correctionRHS(b.size());

        for (size_t i = 0; i < A.size(); i++) {
            cpp_int sum = 0;

            for (int j = 0; j < variables; j++)
                sum += cpp_int(A[i][j]) * solution[j];

            cpp_int residual = cpp_int(b[i]) - sum;

            if (residual % currentModulus != 0)
                return false;

            cpp_int quotient = residual / currentModulus;
            correctionRHS[i] = mod(quotient, prime).convert_to<long long>();
        }

        vector<long long> correction;

        if (!solvePrimeField(
                A,
                correctionRHS,
                prime,
                variables,
                correction
            ))
            return false;

        for (int j = 0; j < variables; j++) {
            cpp_int updated =
                cpp_int(solution[j]) +
                cpp_int(currentModulus) * correction[j];

            solution[j] =
                mod(updated, modulus).convert_to<long long>();
        }

        currentModulus *= prime;
    }

    return true;
}

bool solveRelationSystem(
    const vector<Relation>& relations,
    int variables,
    long long n,
    vector<long long>& logarithms
) {
    vector<vector<long long>> A;
    vector<long long> b;

    for (const Relation& relation : relations) {
        A.push_back(relation.exponents);
        b.push_back(relation.k);
    }

    vector<pair<long long, long long>> components =
        primePowerFactors(n);

    vector<vector<long long>> componentSolutions;

    for (auto component : components) {
        long long prime = component.first;
        long long modulus = component.second;

        vector<long long> solution;

        if (!solvePrimePower(
                A,
                b,
                prime,
                modulus,
                variables,
                solution
            ))
            return false;

        componentSolutions.push_back(solution);
    }

    logarithms.assign(variables, 0);

    for (int i = 0; i < variables; i++) {
        cpp_int x = 0;
        cpp_int M = 1;

        for (size_t j = 0; j < components.size(); j++) {
            long long q = components[j].second;
            long long residue = componentSolutions[j][i];

            long long inverse;

            if (!inverseMod(
                    mod(M, q).convert_to<long long>(),
                    q,
                    inverse
                ))
                return false;

            long long adjustment = mod(
                cpp_int(residue) - x,
                q
            ).convert_to<long long>();

            adjustment = modLL(adjustment * inverse, q);

            x += M * adjustment;
            M *= q;
            x = mod(x, M);
        }

        logarithms[i] =
            mod(x, n).convert_to<long long>();
    }

    return true;
}

bool verifyRelationSystem(
    const vector<Relation>& relations,
    const vector<long long>& logarithms,
    long long n
) {
    for (const Relation& relation : relations) {
        cpp_int lhs = 0;

        for (size_t j = 0; j < logarithms.size(); j++)
            lhs += cpp_int(relation.exponents[j]) * logarithms[j];

        if (mod(lhs, n) != mod(cpp_int(relation.k), n))
            return false;
    }

    return true;
}

bool findTargetRelation(
    long long p,
    long long alpha,
    long long beta,
    long long n,
    const vector<long long>& factorBase,
    const vector<long long>& logarithms,
    long long maxAttempts,
    long long& answer
) {
    mt19937_64 generator(
        chrono::steady_clock::now().time_since_epoch().count()
    );

    uniform_int_distribution<long long> distribution(0, n - 1);
    set<long long> usedK;

    for (long long attempts = 0; attempts < maxAttempts; attempts++) {
        long long k = distribution(generator);

        if (!usedK.insert(k).second)
            continue;

        cpp_int value = mod(
            cpp_int(beta) * modPow(alpha, k, p),
            p
        );

        SmoothResult result = checkSmoothness(value, factorBase);

        if (!result.smooth)
            continue;

        cpp_int x = 0;

        for (size_t i = 0; i < logarithms.size(); i++)
            x += cpp_int(result.exponents[i]) * logarithms[i];

        answer = mod(x - k, n).convert_to<long long>();

        if (modPow(alpha, answer, p) == beta)
            return true;
    }

    return false;
}

int main() {
    long long p, alpha, beta, B;

    cout << "Enter prime p: ";
    cin >> p;

    if (!isPrime(p)) {
        cout << "Error: p must be prime.\n";
        return 1;
    }

    if (p > 1000000000LL) {
        cout << "Error: this implementation is limited to p <= 10^9.\n";
        return 1;
    }

    long long n = p - 1;

    cout << "Enter generator alpha: ";
    cin >> alpha;

    alpha = modLL(alpha, p);

    if (!isPrimitiveRoot(alpha, p)) {
        cout << "Error: alpha must be a primitive root modulo p.\n";
        return 1;
    }

    cout << "Enter beta: ";
    cin >> beta;

    beta = modLL(beta, p);

    if (beta == 0) {
        cout << "Error: beta must be nonzero modulo p.\n";
        return 1;
    }

    cout << "Enter factor base bound B: ";
    cin >> B;

    if (B < 2 || B >= p) {
        cout << "Error: invalid factor base bound.\n";
        return 1;
    }

    vector<long long> factorBase = generateFactorBase(B, p);

    if (factorBase.empty()) {
        cout << "Error: factor base is empty.\n";
        return 1;
    }

    cout << "\nFactor base: ";

    for (long long prime : factorBase)
        cout << prime << " ";

    cout << "\n";

    int variables = (int)factorBase.size();
    long long maxAttempts = 2000000;
    int requiredRelations = variables + 20;
    vector<Relation> relations;
    vector<long long> logarithms;
    bool solved = false;
    double relationTime = 0.0;
    double solveTime = 0.0;

    cout << "\nCollecting relations and solving system...\n";

    for (int round = 0; round < 5 && !solved; round++) {
        auto relationStart = chrono::high_resolution_clock::now();

        relations = collectRelations(
            p,
            alpha,
            n,
            factorBase,
            requiredRelations,
            maxAttempts
        );

        auto relationEnd = chrono::high_resolution_clock::now();
        relationTime += chrono::duration<double>(
            relationEnd - relationStart
        ).count();

        cout << "Collected " << relations.size()
             << " relations (attempt " << round + 1 << ").\n";

        if ((int)relations.size() < variables) {
            requiredRelations += variables;
            continue;
        }

        auto solveStart = chrono::high_resolution_clock::now();
        solved = solveRelationSystem(
            relations,
            variables,
            n,
            logarithms
        );
        auto solveEnd = chrono::high_resolution_clock::now();
        solveTime += chrono::duration<double>(
            solveEnd - solveStart
        ).count();

        if (!solved)
            requiredRelations += variables + 20;
    }

    if (!solved) {
        cout << "Error: relation system could not be solved after retries.\n";
        cout << "Try a different factor base bound B.\n";
        return 1;
    }

    if (!verifyRelationSystem(relations, logarithms, n)) {
        cout << "Error: relation system verification failed.\n";
        return 1;
    }

    cout << "\nFactor base logarithms:\n";

    for (size_t i = 0; i < factorBase.size(); i++) {
        cout << "log_alpha(" << factorBase[i] << ") = "
             << logarithms[i] << "\n";
    }

    long long answer = -1;

    cout << "\nSearching for target relation...\n";

    auto targetStart = chrono::high_resolution_clock::now();

    bool found = findTargetRelation(
        p,
        alpha,
        beta,
        n,
        factorBase,
        logarithms,
        maxAttempts,
        answer
    );

    auto targetEnd = chrono::high_resolution_clock::now();

    double targetTime = chrono::duration<double>(
        targetEnd - targetStart
    ).count();

    if (!found) {
        cout << "Error: target relation was not found.\n";
        return 1;
    }

    cout << "\nDiscrete logarithm:\n";
    cout << "x = " << answer << "\n";

    if (modPow(alpha, answer, p) == beta)
        cout << "Verification: alpha^x mod p = beta.\n";
    else
        cout << "Verification failed.\n";

    cout << "\nExecution time:\n";
    cout << "Relation collection: " << relationTime << " s\n";
    cout << "System solving: " << solveTime << " s\n";
    cout << "Target relation: " << targetTime << " s\n";
    cout << "Total: "
         << relationTime + solveTime + targetTime
         << " s\n";

    return 0;
}
