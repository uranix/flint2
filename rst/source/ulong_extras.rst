.. _ulong_extras:

**ulong_extras** -- machine word arithmetic
:::::::::::::::::::::::::::::::::::::::::::

The *ulong_extras* module provides fast primitives for unsigned machine word
arithmetic:

- bit manipulation
- division/modular arithmetic with precomputed inverses
- greatest common divisor
- integer square/cube/n-th roots
- integer logarithms
- jacobi symbols
- probable prime testing
- primality proving
- integer factorisation
- number theoretic functions (Jacobi, Euler phi, Mobius)
- random generation

Also see *longlong.h* for assembly optimised low level word operations.

Types, macros and constants
---------------------------

The following types, macros and constants relevant to *ulong_extras* are
defined in *flint.h* and available throughout the whole of Flint.

.. type:: ulong

    This is the basic type that *ulong_extras* functions operate on. It is
    defined to be the same as GMP's *mp_limb_t*, which is to say a 32 or 64 bit
    unsigned integer depending on the ABI of the machine. It's usually
    implemented as either an *unsigned long* or *unsigned long long*.

    This type is intended as a portable replacement for *unsigned long* which
    is not portable, especially to Windows and MIPS.

.. type:: slong

    As per *ulong*, but signed. This type is more often used for loop counters
    and lengths of arrays than for arithmetic. It is defined to be the same as
    GMP's *mp_limb_signed_t* and is usually implemented as either a *long* or
    *long long*.

    This type is intended as a portable replacement for *long* which is not
    portable, especially to Windows and MIPS.

.. macro:: UWORD_MAX

    The largest integer able to be stored in a *ulong*. Usually either
    `2^{32} - 1` or `2^{64} - 1`.

.. macro:: UWORD_MIN

    The smallest integer able to be stored in a *ulong*. Usually `0`.

.. macro:: WORD_MAX

    The largest integer able to be stored in an *slong*. Usually either
    `2^{31} - 1` or `2^{63} - 1`.

.. macro:: WORD_MIN

    The smallest integer able to be stored in an *slong*. Usually either
    `-2^{31}` or `-2^{63}`.

Input/output
------------

Flint provides its own printing/reading functions which can deal with the
*ulong* and *slong* types. This makes printing, reading, file and stream
I/O portable between unixes and other operating systems such as Windows.

.. function:: int flint_printf(char * fmt, ...)

    As per the standard C *printf* function, but supports the following
    additional format specifiers:

    - *%wd* : print an *slong* in decimal format
    - *%wu* : print a *ulong* in decimal format
    - *%wx* : print a *ulong* in hexadecimal format
    - *%\*wd* : print an *slong*, space padded in a field of the given width

    Returns the number of characters printed.

.. function:: int flint_sprintf(char * str, const char * fmt, ...)

    As per the standard C *sprintf* function, but with the additional format
    specifiers provided by *flint_printf*.

    Returns the number of characters printed.

.. function:: int flint_fprintf(FILE * stream, const char * fmt, ...)

    As per the standard C *fprintf* function, but with the additional format
    specifiers provided by *flint_printf*.

    Returns the number of characters printed.

.. function:: int flint_scanf(const char * fmt, ...)

    As per the standard C *scanf* function, but supports the following
    additional format specifiers:

    - *%wd* : read an *slong* in decimal format
    - *%wu* : read a *ulong* in decimal format
    - *%wx* : read a *ulong* in hexadecimal format

    Returns the number of items in the argument list successfully filled.

.. function:: int flint_sscanf(const char * str, const char * fmt, ...)

    As per the standard C *sscanf* function, but supports the additional format
    specifiers provided by *flint_scanf*.

    Returns the number of items in the argument list successfully filled.

.. function:: int flint_fscanf(FILE * stream, const char * fmt, ...)

    As per the standard C *fscanf* function, but supports the additional format
    specifiers provided by *flint_scanf*.

    Returns the number of items in the argument list successfully filled.

Bit manipulation
----------------

.. function:: ulong n_revbin(ulong n, ulong b)

    Considering `n` to be a binary number with `b` bits, return the number with
    the reverse binary bit representation. For example *n_revbin(3, 4)*
    would consider `3` as the `4` bit binary number `0011`, which it would
    reverse to give the return value of `12`, with binary representation `1100`.

    **Conditions:** We require `b \leq` *FLINT_BITS*. Only the lower `b` bits of
    `n` are read and the remaining bits can be arbitrary.

    **Algorithm:** If `b \leq 8` we swap the bits in the least significant byte
    using a lookup table, then shift to the required number of bits, otherwise
    we swap all bits in the word and then shift to the required number of bits.

    To swap all the bits in a word we mask and swap alternate bits, then
    alternate pairs of bits, then alternate nibbles and finally swap all bytes
    in the word using the *byte_swap* macro from *longlong.h*.

Arithmetic with precomputed inverse
-----------------------------------

.. function:: ulong n_preinvert_limb(ulong n)

    Return a Moller-Granlund precomputed inverse of the word `n`.

    **Conditions:** We require `n > 0`.

    **Algorithm:** This function normalises `n` (shifts it so that its most
    significant bit is set) and then computes the precomputed inverse
    `(\beta^2 - 1)/n - \beta` where `\beta = 2^w` with `w =` *FLINT_BITS*,
    using the *invert_limb* macro from longlong.h.

    For details of the precomputed inverse see [MolGra2011]_. 

.. function:: ulong n_div2_preinv(ulong a, ulong n, ulong ninv)

    Return the Euclidean quotient of `a` by `n` given a precomputed inverse
    *ninv* provided by *n_preinvert_limb*.

    **Conditions:** We require `n > 0`.

    **Algorithm:** Both `a` and `n` are shifted left by the same number of
    bits `b` so that `n` is normalised. In general `a` now occupies two limbs.
    We then compute the quotient using Algorithm 4 of [MolGra2011]_.

.. function:: ulong n_mod2_preinv(ulong a, ulong n, ulong ninv)

    Return the Euclidean remainder of `a` divided by `n` given a precomputed
    inverse *ninv* provided by *n_preinvert_limb*.

    **Conditions:** We require `n > 0`.

    **Algorithm:** Both `a` and `n` are shifted left by the same number of
    bits `b` so that `n` is normalised. In general `a` now occupies two limbs.
    We then reduce `a` modulo `n` using Algorithm 4 of [MolGra2011]_.
    Finally the resulting remainder is shifted right by `b` bits.

.. function:: ulong n_divrem2_preinv(ulong * q, ulong a, ulong n, ulong ninv)

    Return the Euclidean remainder and set `q` to the Euclidean quotient of `a`
    by `n`, given a precomputed inverse *ninv* provided by *n_preinvert_limb*.

    **Conditions:** We require `n > 0`.

    **Algorithm:** Both `a` and `n` are shifted left by the same number of
    bits `b` so that `n` is normalised. In general `a` now occupies two limbs.
    We then reduce `a` modulo `n` using Algorithm 4 of [MolGra2011]_.
    Finally the resulting remainder is shifted right by `b` bits.

.. function:: ulong n_ll_mod_preinv(ulong a1, ulong a0, ulong n, ulong ninv)

    Given a two word input `a = \langle a_1, a_0 \rangle`, return the
    remainder upon division of `a` by `n`, given a precomputed inverse *ninv*
    provided by *n_preinvert_limb*.

    This function is useful for delayed reduction and reduction of products
    that have accumulated in no more than two words.

    **Conditions:** We require `n > 0`.

    **Algorithm:** If the word `a_1` is not reduced modulo `n` we reduce it
    as per *n_mod2_preinv*. The remainder of `a` divided by `n` is then
    computed using Algorithm 4 of [MolGra2011]_.

.. function:: ulong n_lll_mod_preinv(ulong a2, ulong a1, ulong a0, ulong n, ulong ninv)

    Given a three word input `a = \langle a_2, a_1, a_0 \rangle`, return the
    remainder upon division of `a` by `n`, given a precomputed inverse *ninv*
    provided by *n_preinvert_limb*.

    This function is useful for delayed reduction and reduction of products
    that have accumulated in three words.

    **Conditions:** We require `n > 0` and `a_2 < n`.

    **Algorithm:** As the word `a_2` is reduced modulo `n` we reduce
    `\langle a_2, a_1 \rangle` modulo `n`. Now as `a_2` is reduced modulo `n`
    we reduce `\langle a_1, a_0 \rangle` modulo `n`. The remainders modulo `n`
    are computed using Algorithm 4 of [MolGra2011]_.

For additional functions which take a precomputed inverse, see the modular
arithmetic section below.

Greatest common divisor
-----------------------

.. function:: ulong n_gcd(ulong x, ulong y)

    Return the greatest common divisor of `x` and `y`. If `x = 0` we define
    `\gcd(x, y) = y` and if `y = 0` we define `\gcd(x, y) = x`.

    **Conditions:** None.

    **Algorithm:** Two algorithms are used, the first on machines with a fast
    *count_trailing_zeros* function (currently *x86* and *x86_64*), the other
    as a fallback on other architectures.

    Algorithm 1: First deal with the special cases where either `x = 0` or
    `y = 0`. Now determine the greatest power of `2` dividing both inputs. Call
    this value `2^k`. This value must be the power of `2` dividing the greatest
    common divisor.

    From this point on, any powers of two dividing the two values can be
    divided out, since they do not contribute to the result. In particular we
    begin with the two values shifted right until they are both odd. Call these
    positive, odd values `r_0` and `r_1`.

    At each iteration we start with two unequal, odd values. We subtract the
    smaller from the larger, which doesn't change their GCD, but it makes the
    larger number even. We shift it to the right again so that it is odd, and
    repeat.

    The loop terminates when both of the values are the same. This must happen
    eventually since the sum of the two values is always decreasing and both
    numbers are always positive. The final common value must be
    `\gcd(r_0, r_1)`. We multiply this by `2^k` to get `\gcd(x, y)`.

    Algorithm 2: This is a variant of the ordinary Euclidean algorithm. We
    begin with `r_0 = x` and `r_1 = y` and keep applying the division algorithm
    in order to obtain a remainder sequence `\{r_i\}`. The last nonzero
    remainder is the greatest common divisor (if both inputs are zero all tests
    fall through and zero is returned).

    To minimize the number of divisions performed, the algorithm deals
    specially with the cases were `r_i < 4r_{i+1}`, i.e. where the quotient is
    `1`, `2` or `3`.

    We first compute `s = r_i - r_{i+1}`. If `s < r_{i+1}`, i.e.
    `r_i < 2r_{i+1}`, we know the quotient is `1`, else if `s < 2r_{i+1}`, i.e.
    `r_i < 3r_{i+1}` we know the quotient is `2`. In the remaining cases, the
    quotient must be `3`.

    When the quotient is `4` or above, we use division. However this happens
    rarely for generic inputs.

    To prevent overflows in the arithmetic the values are first reordered so
    that `r_0 \geq r_1`. The special case where both have top bit set is then
    dealt with, followed by the case where the second value has its second most
    significant bit set. It is then safe to multiply the second value by `4` as
    required by the algorithm, without overflow.

.. function:: ulong n_gcdinv(ulong * a, ulong x, ulong y)

    Return the greatest common divisor of `x` and `y` and set `a` to a value in
    the range `[0, y)` such that `ax \equiv \gcd(x, y) \pmod{y}`.

    If `y = 1` then the greatest common divisor is `1` and `a` is set to `0`.

    **Conditions:** We require `x < y`. In particular `y \neq 0`. Aliasing of
    `a` with any of the inputs is permitted.

    **Algorithm:** The algorithm to compute the greatest common divisor is as
    per Algorithm 2 of *n_gcd*.

    In order to compute the cofactor, we start with `v_1 = 0` and `v_2 = 1`.
    Each time we compute `x = qy + r` in the Euclidean algorithm, we set 

    .. math::
        \left(\begin{array}{c}v_1\\ v_2\end{array}\right) = 
        \left(\begin{array}{cc}0 & 1\\ 1 & -q\end{array}\right)
        \left(\begin{array}{c}v_1\\ v_2\end{array}\right).

    Upon termination of the Euclidean algorithm, `v_1` is a cofactor in the
    range `[-y/2, y/2]`. If it is negative we add `y` to it so that it is in
    the range `[0, y)`.

    For a proof that the cofactors never overflow, see *n_xgcd*.

.. function:: ulong n_xgcd(ulong * s, ulong * t, ulong x, ulong y)

    Return the greatest common divisor of `x` and `y` and set `s` and `t` to
    non-negative values such that `\gcd(x, y) = sx - ty`.

    If `y \neq 0` we will have `s \leq y` and `t \leq x`.

    In the case that `y = 0` we will have `s = 1` and `t = 0`.

    **Conditions:** We require `x \geq y`. Aliasing of `s` and `t` with any
    of the inputs is permitted.

    **Algorithm:** The algorithm to compute the greatest common divisor is as
    per Algorithm 2 of *n_gcd*.

    We compute the cofactors by starting with a matrix with signed entries

    .. math::
        M = \left(\begin{array}{cc}u_1 & v_1\\ u_2 & v_2\end{array}\right)
        = \left(\begin{array}{cc}1 & 0\\ 0 & 1\end{array}\right)

    At each iteration of the algorithm we compute 
    `r_{i - 1} = q_ir_i + r_{i + 1}`.

    We multiply the matrix `M` on the left by

    .. math::
        \left(\begin{array}{cc}0 & 1\\ 1 & -q_i\end{array}\right)
    
    After each step of the algorithm we will have

    .. math::
        \left(\begin{array}{c}r_i\\ r_{i + 1}\end{array}\right) = 
        \left(\begin{array}{cc}u_1 & v_1\\ u_2 & v_2\end{array}\right)
        \left(\begin{array}{c}x\\ y\end{array}\right)
          
    We claim that if the greatest common divisor is computed via the Euclidean
    algorithm, starting with `x \geq y > 0` and `x` not a multiple of `y` then
    we always have `|s| \leq y/2` and `|t| < x/2`.

    Recall that the cofactors are obtained by backsubstituting the steps of the
    Euclidean algorithm. We first prove the result for the case
    `\gcd(x, y) = 1`. We proceed by induction.

    We will show that at each step in the backsubstitution we have
    `1 = \pm s r_{i-1} \mp t r_i` with `s \leq r_i/2, t < r_{i-1}/2`.

    At the final step of the Euclidean algorithm we have
    `r_{n-1} = q_nr_n + 1`. We can rewrite this as `1 = r_{n-1} - q_nr_n`.
    As `r_n > 1` we must have `q_n < r_{n-1}/2`. Similarly, as
    `r_n > 1` we have `1 \leq r_n/2`. Thus the claim holds at
    the first step of the backsubstitution.

    Suppose the claim is true at some point in the backsubstitution. The
    next equation to backsubstitute is `r_{i-2} = q_{i-1}r_{i-1} + r_i`.
    This we rewrite as `r_i = r_{i-2} - q_{i-1}r_{i-1}`.

    Making the substitution yields
    `1 = s r_{i-1} - t (r_{i-2} - q_{i-1}r_{i-1})
    = (s + t q_{i-1}) r_{i-1} - t r_{i-2}`.

    It suffices to show that `s + t q_{i-1} < r_{i-2}/2` since
    `t < r_{i-1}/2` which would complete the induction.

    But `s + t q_{i-1} < r_i/2 + r_{i-1}/2 q_{i-1} = r_{i-2}/2`, so the claim
    is proved.

    In the case where the greatest common divisor is greater than `1`, all the
    equations are simply multiplied through by the GCD and the result holds
    there too.

    The important consequence of this theorem is that the cofactors can never
    overflow a signed word and comparison of the cofactors with zero is always
    permitted.

    This means that at the end of the algorithm, if we have `1 = -s x + t y`
    for `s, t > 0` we can replace `s` with `s + y` and `t` with `t - x`. Then
    the first cofactor is guaranteed to be positive.

    In the case where `y = 0` the algorithm terminates immediately with
    cofactors `s = 1` and `t = 0`.

    In the case where `x` is a multiple of `y` the algorithm terminates
    after one step with `s = y` and `t = 1 - x`.

Modular arithmetic
------------------

.. function:: ulong n_addmod(ulong a, ulong b, ulong n)

    Returns `a + b \pmod{n}`.

    **Conditions:** Requires that `a` and `b` are reduced modulo `n` and that
    `n \neq 0`.

    **Algorithm:** We subtract `y` from `n` and if the result is greater than
    `x` we know that the result of `x + y` will be less than `n`. But this
    comparison has the advantage of not overflowing the word.

    If `x + y` will be less than `n` we return that result, otherwise we return
    `x + y - n`. It does not matter if overflow occurs during this computation
    since the result is computed modulo `2^B` where `B` is the number of bits
    in a word, and the result fits in a word.
    
.. function:: ulong n_submod(ulong a, ulong b, ulong n)

    Returns `a - b \pmod{n}`.

    **Conditions:** Requires that `a` and `b` are reduced modulo `n` and that
    `n \neq 0`.

    **Algorithm:** If `y > x` we compute `x - y + n`, otherwise the we compute
    `x - y`. It doesn't matter if overflow occurs during the computation of
    `x - y + n` since the result is computed modulo `2^B` where `B` is the
    number of bits in a word, and the result fits in a word.

.. function:: ulong n_negmod(ulong a, ulong n)

    Returns `-a \pmod{n}`.

    **Conditions:** Requires that `a` is reduced modulo `n` and that
    `n \neq 0`.

    **Algorithm:** This function calls *n_submod* with first argument `0`. The
    call is inlined.

.. function:: ulong n_mulmod_preinv(ulong a, ulong b, ulong n, ulong ninv, ulong norm)

    Returns `ab/2^m \pmod{n}` where `m =` *norm*, given a precomputed inverse
    *ninv* provided by *n_preinvert_limb*. This function is intended to be used
    to compute `ab \pmod{n}` as described below.

    This is the fastest but least convenient mulmod function. It requires `n`
    to be normalised (most significant bit set). In this case it can be used
    with *norm* equal to `0` and the function will then compute `ab \pmod{n}`.

    However, the function is designed to be used with other values of `n` as
    follows. Let *norm* be the number of leading zeroes of `n`. Before using
    this function we shift `a`, `b` and `n` to the left by *norm* bits. Then
    after using the function, we shift the result to the right by *norm* bits.

    Using the function in this way will result in `ab \pmod{n}` for any `n`.

    Note that the function performs an additional shift right by *norm* bits
    internally so that the result is meaningful after the user also performs
    such a shift (two such shifts are required in total since both `a` and
    `b` will have been shifted left by *norm* bits).
    
    The function is generally intended to be used in cases where a long
    computation is carried out in shifted representation. This function is
    designed to leave the result in such representation so that no additional
    shifts are required before the next operation.

    **Conditions**: We require `n \neq 0`, `a, b < n` and `n` to be normalised,
    i.e. most significant bit set. However, see the description for how to use
    this function with other values of `n`.

    **Algorithm**: The function first shifts `a` right by *norm* bits, then
    computes the product `ab`. The result is then reduced modulo `n` using 
    Algorithm 4 of [MolGra2011]_.

.. function:: ulong n_mulmod2_preinv(ulong a, ulong b, ulong n, ulong ninv)

    Returns `ab \pmod{n}` given a precomputed *ninv* provided by 
    *n_preinvert_limb*.

    **Conditions:** We require `n \neq 0`.

    **Algorithm:** The product `ab` is reduced modulo `n` using
    *n_ll_mod_preinv*.

.. function:: ulong n_mulmod2(ulong a, ulong b, ulong n)

    Returns `ab \pmod{n}`. This is the most convenient, but least efficient
    mulmod function. It is provided for convenience only. As it doesn't take a
    precomputed inverse it is only useful in code that is not time critical.

    **Conditions:** We require `n \neq 0`.

    **Algorithm:** A precomputed inverse is computed after which the
    implementation is the same as *n_mulmod2_preinv*.

.. function:: ulong n_powmod_ui_preinv(ulong a, ulong m, ulong n, ulong ninv, ulong norm)

    Returns `a^m \pmod{n}`. For convenience we define everything
    modulo `1` to be `0` and otherwise `a^0 = 1 \pmod{n}` for all `n`.

    This is the fastest but least convenient powmod function. It requires `n`
    to be normalised. However it can be used with other values of `n` with
    inputs in shifted representation as per *n_mulmod_preinv*.

    **Conditions:** We require `n \neq 0`, `a < n` and `n` to be normalised,
    i.e. most significant bit set. However, see the description for how to use
    this function with other values of `n`. Note that `m` is unsigned.

    **Algorithm:** This uses binary exponentiation with an optimisation to
    save a multiplication after the first `1` is encountered in the binary
    representation of `m`.

.. function:: ulong n_powmod2_ui_preinv(ulong a, ulong m, ulong n, ulong ninv)

    Returns `a^m \pmod{n}`. For convenience we define everything
    modulo `1` to be `0` and otherwise `a^0 = 1 \pmod{n}` for all `n`.

    **Conditions:** We require `n \neq 0`. Note that `m` is unsigned.

    **Algorithm:** If `a \geq n` we first reduce it modulo `n`. We then use
    binary exponentiation with an optimisation to save a multiplication after
    the first `1` is encountered in the binary representation of `m`.

.. function:: ulong n_powmod2_preinv(ulong a, slong m, ulong n, ulong ninv)

    Returns `a^m \pmod{n}`. For convenience we define everything
    modulo `1` to be `0` and otherwise `a^0 = 1 \pmod{n}` for all `n`.

    **Conditions:** We require `n \neq 0`. Note that `m` is signed. If `a`
    is not invertible modulo `n` and exception is raised.

    **Algorithm:** If `m < 0` we first invert `a` modulo `n`. If `a \geq n` we
    reduce it modulo `n`. We then call *n_powmod_ui_preinv*.
