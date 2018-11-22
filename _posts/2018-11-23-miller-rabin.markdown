---
title: "Miller Rabin"
layout: post
date: 2018-11-23 03:18
image: /assets/images/markdown.jpg
headerImage: false
tag:
- algorithm
- crypto
category: blog
# jemoji: '<img class="emoji" title=":ramen:" alt=":ramen:" src="https://assets.github.com/images/icons/emoji/unicode/1f35c.png" height="20" width="20" align="absmiddle">'
---

## Problem

Test whether a given positive integer $$ N $$ is prime. 

## Candidate Solutions

1. [Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes) ($$ O(n lg lg n) $$ preprocessing, $$ O(1) $$ query)
2. [Linear Sieve](https://cp-algorithms.com/algebra/prime-sieve-linear.html) ($$ O(n) $$)
3. [Sieve of Atkin](https://en.wikipedia.org/wiki/Sieve_of_Atkin) ($$ O(n / lg lg n) $$ preprocessing, $$ O(1) $$ query)
4. [Miller-Rabin](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test) ($$ O(k lg n) $$ non-deterministic)

Overall, as for general purposes considering that $$ n $$ is stored in less that 64bits, Miller-Rabin is the best option as it is fast enough and simple to write. 

## Miller-Rabin Primality Test

Let us start with [Fermat's Little Theorem](https://en.wikipedia.org/wiki/Fermat%27s_little_theorem), which states that if $$ p $$ is a prime number, then for any integer $$ a $$ not divisible by $$ p $$, $$ a^{p-1} \equiv 1 \pmod p $$.

The contrapositive gives that, if there exists some $$ a $$ not divisible by $$ p $$ such that $$ a^{p-1}\not\equiv 1\pmod p $$, $$ p $$ is not prime. From here we could reasonably derive a probabilistic primality test by randomly choosing some number of different $$ a $$ and check whether there is a witness $$ a $$ that finds $$ p $$ having a behaviour not conforming the prime number property. This is called the [Fermat Primality Test](https://en.wikipedia.org/wiki/Fermat_primality_test). However, there are some numbers that behaves like prime for every possible $$ a $$,
i.e. $$ 561 $$, which are called [Carmichael Numbers](https://en.wikipedia.org/wiki/Carmichael_number).

In fact, Miller-Rabin is an algorithm that improves Fermats primality test by being able to test such numbers without losing efficiency.

Let us first consider this lemma:

> For any $$ x^2 \equiv 1\pmod p $$ where $$ p $$ is prime, $$ x \equiv 1\pmod p $$ or $$ x \equiv -1\pmod p $$ must hold.\\
> This is for $$ x^2 - 1 = (x + 1)(x - 1) = pk $$ for some integer $$ k $$, where $$ x + 1 $$ and $$ x - 1 $$ cannot both be divisible by $$ p $$ (unless $$ p = 2 $$ which is trivial). So at least one of $$ x + 1 $$ and $$ x - 1 $$ is divislbe by $$ p $$, which implies the lemma.

Assume that the number we are to test is $$ n $$. We know that if $$ n $$ is prime, then for any $$ a $$ not divisible by $$ n $$, $$ a^{n-1} \equiv 1 \pmod n $$. Let $$ n - 1 = 2^d s $$, $$ 2\nmid s $$. Then by the lemma, as we decrement $$ d $$, we are taking square root of $$ 1 $$, and we should only have either $$ -1 $$ or $$ 1 $$ as the result. If we get $$ -1 $$, then it could be a prime and we cannot get anything more from the $$ a $$ being tested; if we get $$ 1 $$, then we should
continue doing square root (unless $$ d = 0 $$ which again implies that $$ a $$ cannot ensure us anything); otherwise any other result implies that $$ n $$ is composite.

It is said that by uniformly and independently choosing $$ a $$ to test $$ w $$ times, the probability of getting false positive is $$ 4^{-w} $$. However, in practice, the set of candidate witnesses $$ a $$ is often fixed. For example, it has been [verified](http://miller-rabin.appspot.com/) that the set {2, 325, 9375, 28178, 450775, 9780504, 1795265022} can cover all composite numbers upto $$ 2^{64} $$, which is sufficient for common usage. 

## Snippet

{% highlight cpp %}
struct PrimeFactor {
  static long long mul(long long x, long long y, long long mod) {
    long long m = x, s = 0;
    for (; y; y >>= 1, m <<= 1, m = m >= mod ? m - mod : m)
      if (y & 1) s += m, s = s >= mod ? s - mod : s;
    return s;
  }
  static long long powmod(long long x, long long p, long long mod) {
    long long s = 1, m = x % mod;
    for (; p; m = mul(m, m, mod), p >>= 1)
      if (p & 1) s = mul(s, m, mod);
    return s;
  }
  static bool miller_rabin(long long n, int s = 7) {
    static const long long wits[7] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    auto witness = [&](long long a, long long n, long long u, int t) {
      long long x = powmod(a, u, n), nx;
      for (int i = 0; i < t; ++i, x = nx) {
        nx = mul(x, x, n);
        if (nx == 1 && x != 1 && x != n - 1) return true;
      }
      return x != 1;
    };
    if (n < 2) return 0;
    if (~n & 1) return n == 2;
    long long u = n - 1, t = 0, a;
    for (; ~u & 1; u >>= 1) ++t;
    while (s--) if ((a = wits[s] % n) && witness(a, n, u, t)) return false;
    return true;
  }
};
{% endhighlight %}
