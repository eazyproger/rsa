import random
from sympy import isprime
import math
import decimal


def nums_generation():
    prime_num = []
    while True:
        a = random.randrange(10 ** 18, 10 ** 19)
        if isprime(a):
            prime_num.append(a)
            if len(prime_num) == 2:
                if prime_num[0] == prime_num[1]:
                    prime_num.remove(prime_num[1])
                else:
                    break
    return prime_num


def n_bit_random(n):
    return random.randrange(2 ** (n - 1) + 1, 2 ** n - 1)


first_primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                     53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
                     109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                     173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                     233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
                     293, 307, 311, 313, 317, 331, 337, 347, 349]


def check_low_level_prime(n):
    while True:
        prime_candidate = n_bit_random(n)
        for divisor in first_primes_list:
            if prime_candidate % divisor == 0 and divisor ** 2 <= prime_candidate:
                break
            return prime_candidate


def check_miller_rabin_primality(mrc):
    max_divisions_by_two = 0
    ec = mrc - 1
    while ec % 2 == 0:
        ec >>= 1
        max_divisions_by_two += 1
    assert (2 ** max_divisions_by_two * ec == mrc - 1)  # lame

    def trial_composite(round_tester):
        if pow(round_tester, ec, mrc) == 1:
            return False
        for i in range(max_divisions_by_two):
            if pow(round_tester, 2 ** i * ec, mrc) == mrc - 1:
                return False
        return True

    trials_amount = 20
    for i in range(trials_amount):
        round_tester = random.randrange(2, mrc)
        if trial_composite(round_tester):
            return False
    return True


def n_bit_random_prime(n):
    while True:
        possible_prime = check_low_level_prime(n)
        if not check_miller_rabin_primality(possible_prime):
            continue
        return possible_prime


decimal.getcontext().prec=10000
_d = lambda x: decimal.Decimal(x)


def modular_pow(base, exponent, modulus):
    base, modulus = _d(base), _d(modulus)
    base %= modulus
    result = 1
    while exponent > 0:
        if int(exponent) & 1:
            result = (result * base) % modulus
        base = (base * base) % modulus
        exponent //= 2
    return result


def encrypt(message, exponent, modulus):
    out = []
    for c in message:
        encv = modular_pow(ord(c), exponent, modulus)
        out.append(encv)
    return out


def decrypt(encrypted, pkey, modulus):
    out = []
    for num in encrypted:
        decv = modular_pow(num, pkey, modulus)
        out.append(chr(decv))
    return "".join(out)


def lcm(a, b):
    return abs(a*b) // math.gcd(a, b)


def choose_e(lcmv):
    while 1:
        e = random.randrange(3, 2**16+2)
        if e < lcmv and math.gcd(e, lcmv) == 1:
            return e


def extended_gcd(a, b):
    x, old_x = 0, 1
    y, old_y = 1, 0

    while (b != 0):
        quotient = a // b
        a, b = b, a - quotient * b
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y

    return a, old_x, old_y


P, Q = nums_generation()

n = P * Q
lcmv = lcm(P-1, Q-1)
e = choose_e(lcmv)

gcd, x, y = extended_gcd(e, lcmv)
if x < 0:
    d = x + lcmv
else:
    d = x

print(e, d, n, len(str(n)), sep="\n")

testing_list = ["test1", "1234567890", "big mess"*20]
for mes in testing_list:
    print(mes == decrypt(encrypt(mes, e, n), d, n))
