# Bipartite Bound Information Exists

Finite-block calculations accompanying **Bipartite Bound Information Exists**
by Jef Pauwels, Nicolas Gisin, and Renato Renner.

`example3_parity_family.py` implements the recurrences in Supplemental Material
Sec. S9.6 for the Gisin–Wolf Example 3 source, with `2 <= alpha <= 5`.
It requires **Python 3.10 or later** and uses only the standard library.

## Reproduce the certificates

Run these commands from this directory:

```sh
python3 example3_parity_family.py --alpha 18/5 --m 5 --r 3
python3 example3_parity_family.py --alpha 7/2 --m 6 --r 4
python3 example3_parity_family.py --alpha 16/5 --m 119 --r 10
python3 example3_parity_family.py --alpha 4 --m 2 --r 1
```

The first three commands reproduce the examples in Sec. S9.6 with blocklengths
15, 24, and 1190. The fourth reproduces the two-copy certificate at `alpha = 4`
in Sec. S9.4. Each prints `True` for `B^2-PQ > 0`.

The filter repeats each of `m` coordinates `r` times, giving blocklength
`n = m*r`. The coset residue is chosen as in Sec. S9.6: `c = 1` when
`m = 2 mod 3`, and `c = 2` when `m = 0 mod 3`.

For rational `alpha`, the criterion is checked with exact integer arithmetic.
A positive result certifies a strictly positive secret-key rate; it does not
calculate that rate. The displayed `B^2/PQ` is a floating-point approximation.
Add `--print-exact-margin` to display the exact integer margin after clearing
denominators.

## Bounded search and numerical threshold

Omit `--m` and `--r` to search within specified bounds:

```sh
python3 example3_parity_family.py --alpha 18/5 --max-m 20 --max-r 5
```

This returns the shortest blocklength found within those bounds and the residue
choices above. It is not a search over all possible filters. An unsuccessful
finite search does not imply a zero secret-key rate.

The numerical helper estimates a crossing of the criterion in `[3,4]`:

```sh
python3 -c "from example3_parity_family import threshold_decimal; print(threshold_decimal(2, 1))"
```

For the two-copy filter this gives approximately `3.917737700`. The proof in
Sec. S9.4 establishes that its criterion holds above this threshold, up to and
including `alpha = 5`. For general `m` and `r`, the helper does not prove
uniqueness of the crossing or the full working range.

## Version and license

This is version **v1.0.0** of the accompanying code, distributed under the
[MIT License](LICENSE).
