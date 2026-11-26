"""Finite-block calculations for Gisin--Wolf Example 3.

This implements the parity-code filters in Supplemental Material Sec. S9.6
of "Bipartite Bound Information Exists".

Choose m coordinates, repeat each coordinate r times, and let

    L_m = {u in F_3^m : sum_i u_i = 0}.

For c in {1,2}, the four GGK filtering sets are the expanded copies of

    A1=L_m, A2=B1=L_m+c-shift, B2=L_m-c-shift.

Only the affine coset C_c={u: sum u_i=c} matters.  A word in C_c has a
disjoint-support partner in C_c iff it contains at least one zero.

The weighted residue sums are computed by a three-state dynamic program, so
the calculation uses O(m) recurrence steps rather than enumerating O(3^m) words.
Rational alpha values are handled with exact integers.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction


Triple = tuple[int, int, int]


def convolve_residues(state: Triple, weights: Triple) -> Triple:
    out = [0, 0, 0]
    for residue, value in enumerate(state):
        for symbol, weight in enumerate(weights):
            out[(residue + symbol) % 3] += value * weight
    return tuple(out)  # type: ignore[return-value]


@dataclass(frozen=True)
class ExactCertificate:
    alpha: Fraction
    m: int
    r: int
    c: int
    good_overlap: int
    inverse_coset_mass: int
    code_mass: int

    @property
    def squared_margin(self) -> int:
        return (
            self.good_overlap * self.good_overlap
            - self.inverse_coset_mass * self.code_mass
        )

    @property
    def ratio_squared(self) -> Fraction:
        return Fraction(
            self.good_overlap * self.good_overlap,
            self.inverse_coset_mass * self.code_mass,
        )

    @property
    def blocklength(self) -> int:
        return self.m * self.r


def exact_certificate(alpha: Fraction, m: int, r: int, c: int | None = None) -> ExactCertificate:
    if alpha < 2 or alpha > 5:
        raise ValueError("alpha must satisfy 2 <= alpha <= 5")
    if m < 2 or r < 1:
        raise ValueError("m >= 2 and r >= 1 are required")
    if c is None:
        if m % 3 == 0:
            c = 2
        elif m % 3 == 2:
            c = 1
        else:
            raise ValueError("the asymptotic construction skips m == 1 mod 3")
    if c not in (1, 2):
        raise ValueError("c must be 1 or 2")

    # If alpha=p/q, multiplying each repeated-coordinate weight by q^r
    # produces the exact integers below.  The common factor cancels.
    p, q = alpha.numerator, alpha.denominator
    A = (2 * q) ** r
    B = p**r
    C = (5 * q - p) ** r

    all_residues: Triple = (1, 0, 0)
    full_support: Triple = (1, 0, 0)
    for _ in range(m):
        all_residues = convolve_residues(all_residues, (A, B, C))
        full_support = convolve_residues(full_support, (0, B, C))

    good = all_residues[c] - full_support[c]
    inverse = all_residues[(-c) % 3]
    code = all_residues[0]
    return ExactCertificate(alpha, m, r, c, good, inverse, code)


def first_exact_witness(
    alpha: Fraction,
    max_m: int,
    max_r: int,
) -> ExactCertificate | None:
    best: ExactCertificate | None = None
    for r in range(1, max_r + 1):
        for m in range(2, max_m + 1):
            if m % 3 == 1:
                continue
            certificate = exact_certificate(alpha, m, r)
            if certificate.squared_margin <= 0:
                continue
            if best is None or certificate.blocklength < best.blocklength:
                best = certificate
    return best


def decimal_squared_margin(
    alpha: Decimal,
    m: int,
    r: int,
    c: int | None = None,
    precision: int = 80,
) -> Decimal:
    """High-precision normalized margin, useful for threshold bisection."""
    if c is None:
        if m % 3 == 0:
            c = 2
        elif m % 3 == 2:
            c = 1
        else:
            raise ValueError("the asymptotic construction skips m == 1 mod 3")
    with localcontext() as context:
        context.prec = precision
        x = (Decimal(2) / alpha) ** r
        y = ((Decimal(5) - alpha) / alpha) ** r
        all_residues = (Decimal(1), Decimal(0), Decimal(0))
        full_support = (Decimal(1), Decimal(0), Decimal(0))
        for _ in range(m):
            all_residues = convolve_residues(all_residues, (x, Decimal(1), y))
            full_support = convolve_residues(
                full_support, (Decimal(0), Decimal(1), y)
            )
        good = all_residues[c] - full_support[c]
        return +(good * good - all_residues[(-c) % 3] * all_residues[0])


def threshold_decimal(m: int, r: int, precision: int = 60) -> Decimal | None:
    """Estimate a crossing of the criterion in [3,4] by bisection.

    This numerical estimate does not establish uniqueness of the crossing
    or the full parameter range satisfying the criterion.
    """
    with localcontext() as context:
        context.prec = precision
        lo, hi = Decimal(3), Decimal(4)
        if decimal_squared_margin(hi, m, r, precision=precision) <= 0:
            return None
        for _ in range(4 * precision):
            mid = (lo + hi) / 2
            if decimal_squared_margin(mid, m, r, precision=precision) > 0:
                hi = mid
            else:
                lo = mid
        return +(lo + hi) / 2


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def main() -> None:
    # Exact certificates can contain more digits than Python's default limit.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=parse_fraction, default=Fraction(4))
    parser.add_argument("--m", type=int)
    parser.add_argument("--r", type=int)
    parser.add_argument("--max-m", type=int, default=200)
    parser.add_argument("--max-r", type=int, default=15)
    parser.add_argument("--print-exact-margin", action="store_true")
    args = parser.parse_args()

    if args.m is not None or args.r is not None:
        if args.m is None or args.r is None:
            parser.error("--m and --r must be supplied together")
        certificate = exact_certificate(args.alpha, args.m, args.r)
    else:
        certificate = first_exact_witness(args.alpha, args.max_m, args.max_r)
        if certificate is None:
            print("No positive-rate example found within the requested search bounds.")
            return

    print("alpha                 =", certificate.alpha)
    print("m, r, n               =", certificate.m, certificate.r, certificate.blocklength)
    print("coset residue c       =", certificate.c)
    print("B^2-PQ > 0            =", certificate.squared_margin > 0)
    print("B^2/PQ                =", float(certificate.ratio_squared))
    if args.print_exact_margin:
        print("exact squared margin  =", certificate.squared_margin)
    else:
        print(
            "margin decimal digits =",
            len(str(abs(certificate.squared_margin))),
            "(use --print-exact-margin to display it)",
        )


if __name__ == "__main__":
    main()
