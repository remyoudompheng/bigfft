package bigfft

import (
	. "math/big"
)

// Arithmetic modulo 2^n+1.

// A fermat of length w+1 represents a number modulo 2^(w*_W) + 1. The last
// word is zero or one. A number has at most two representatives satisfying the
// 0-1 last word constraint.
type fermat nat

func (n fermat) String() string { return nat(n).String() }

// Shift computes (x << k) mod (2^n+1).
func (z fermat) Shift(x fermat, k int) {
	if len(z) != len(x) {
		println(len(z), len(x))
		panic("len(z) != len(x) in Shift")
	}
	n := len(x) - 1
	// Shift by n*_W is Neg.
	k %= 2 * n * _W
	if k < 0 {
		k += 2 * n * _W
	}
	neg := false
	if k >= n*_W {
		k -= n * _W
		neg = true
	}
	// Shift left by kw words.
	// x = a·2^(n-k) + b
	// x<<k = (b<<k) - a
	kw, kb := k/_W, k%_W
	copy(z[kw:], x[:n-kw])
	z[n] = 1 // Add (-1)
	b := subVV(z[:kw+1], z[:kw+1], x[n-kw:])
	subVW(z[kw+1:], z[kw+1:], b)
	addVW(z, z, 1) // Add back 1.
	// Shift left by kb bits
	shlVU(z, z, uint(kb))
	c := z[n]
	if z[0] >= c {
		z[n] = 0
		z[0] -= c
	} else {
		z[n] = 1
		subVW(z, z, c-1)
	}
	if neg {
		z.Neg()
	}
}

// Neg computes the opposite of x mod 2^n+1.
func (x fermat) Neg() {
	n := len(x) - 1
	c := x[n]
	// 2^n + 1 - x = ^x + 2.
	for i, a := range x {
		x[i] = ^a
	}
	x[n] = 0
	c = addVW(x[:n], x[:n], c+2)
	if x[0] >= c {
		x[0] -= c
	} else {
		x[n] = 1
		subVW(x, x, c-1)
	}
	return
}

// Add computes addition mod 2^n+1.
func (z fermat) Add(x, y fermat) fermat {
	z = x
	addVV(z, x, y) // there cannot be a carry here.
	c := z[len(z)-1]
	if z[0] >= c {
		z[len(z)-1] = 0
		z[0] -= c
	} else {
		z[len(z)-1] = 1
		subVW(z, z, c-1)
	}
	return z
}

func (z fermat) Mul(x, y fermat) fermat {
	var xi, yi, zi Int
	xi.SetBits(x)
	yi.SetBits(y)
	zi.SetBits(z)
	z = zi.Mul(&xi, &yi).Bits()
	n := len(x) - 1
	// len(z) is at most 2n+1.
	if len(z) > 2*n+1 {
		panic("len(z) > 2n+1")
	}
	i := len(z) - (n + 1)
	c := subVV(z[1:i+1], z[1:i+1], z[n+1:])
	z = z[:n+1]
	z[n]++ // Add -1.
	subVW(z[i+1:], z[i+1:], c)
	// Normalize.
	c = z[n]
	if z[0] >= c {
		z[n] = 0
		z[0] -= c
	} else {
		z[n] = 1
		subVW(z, z, c-1)
	}
	// Add 1.
	if z[n] == 1 {
		z[n] = 0
	} else {
		addVW(z, z, 1)
	}
	return z
}
