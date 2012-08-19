package bigfft

import (
	. "math/big"
	"unsafe"
)

const _W = int(unsafe.Sizeof(Word(0)) * 8)

type nat []Word

func (n nat) String() string {
	v := new(Int)
	v.SetBits(n)
	return v.String()
}

// returns the FFT length k, m the number of words per chunk
func fftSize(x, y nat) (k uint, m int) {
	words := len(x) + len(y)
	bits := int64(words) * int64(_W)
	switch {
	case bits < 1<<12:
		k = 6 // chunks have less than 64 bits.
	case bits < 1<<18:
		k = 8 // chunks have less than 1024 bits.
	case bits < 1<<24:
		k = 10 // chunks have less than 64*N>>k < 8192 bits.
	case bits < 1<<30:
		k = 14 // chunks have less than 64*(1<<10) < 64kbits
	default:
		k = 16
	}
	m = words>>k + 1
	return
}

// poly represents an integer via a polynomial in Z[x]/(x^K+1)
// where K is the FFT length and b is the computation basis 1<<(m*_W).
// If P = a[0] + a[1] x + ... a[n] x^(K-1), the associated natural number
// is P(b^m).
type poly struct {
	k uint  // k is such that K = 1<<k.
	m int   // the m such that P(b^m) is the original number.
	a []nat // a slice of at most K m-word coefficients.
}

// polyFromNat slices the number x into a polynomial
// with 1<<k coefficients made of m words.
func polyFromNat(x nat, k uint, m int) poly {
	p := poly{k: k, m: m}
	length := len(x)/m + 1
	p.a = make([]nat, length)
	for i := range p.a {
		if len(x) < m {
			p.a[i] = make(nat, m)
			copy(p.a[i], x)
			break
		}
		p.a[i] = x[:m]
		x = x[m:]
	}
	return p
}

// Int evaluates back a poly to its integer value.
func (p *poly) Int() nat {
	length := len(p.a)*p.m + 1
	if na := len(p.a); na > 0 {
		length += len(p.a[na-1])
	}
	n := make(nat, length)
	m := p.m
	np := n
	for i := range p.a {
		l := len(p.a[i])
		c := addVV(np[:l], np[:l], p.a[i])
		if np[l] < ^Word(0) {
			np[l] += c
		} else {
			addVW(np[l:], np[l:], c)
		}
		np = np[m:]
	}
	n = trim(n)
	return n
}

func trim(n nat) nat {
	for i := range n {
		if n[len(n)-1-i] != 0 {
			return n[:len(n)-i]
		}
	}
	return nil
}

func (p *poly) valueSize() int {
	n := 2*p.m*_W + int(p.k)
	n = ((n >> p.k) + 1) << p.k
	return n / _W
}

// A polValues represents the value of a poly at the odd powers of a
// (2K)-th root of unity θ=2^l in Z/(b^n+1)Z, where b^n = 2^Kl.
type polValues struct {
	k      uint     // k is such that K = 1<<k.
	n      int      // the length of coefficients, n*_W a multiple of 1<<k.
	values []fermat // a slice of K (n+1)-word values
}

// Transform evaluates p at θω^i for i = 0...K-1, where
// θ is a (2K)-th primitive root of unity in Z/(b^n+1)Z
// and ω = θ².
func (p *poly) Transform(n int) polValues {
	k := p.k
	// θ is represented as a shift.
	θshift := (n * _W) >> k
	// p(x) = a_0 + a_1 x + ... + a_{K-1} x^(K-1)
	// p(θx) = q(x) where
	// q(x) = a_0 + θa_1 x + ... + θ^(K-1) a_{K-1} x^(K-1)
	//
	// Twist p by θ to obtain q.
	tbits := make([]Word, (n+1)<<k)
	twisted := make([]fermat, 1<<k)
	src := make(fermat, n+1)
	for i := range twisted {
		twisted[i] = fermat(tbits[i*(n+1) : (i+1)*(n+1)])
		if i < len(p.a) {
			for i := range src {
				src[i] = 0
			}
			copy(src, p.a[i])
			twisted[i].Shift(src, θshift*i)
		}
	}

	// Now computed q(ω^i) for i = 0 ... K-1
	valbits := make([]Word, (n+1)<<k)
	values := make([]fermat, 1<<k)
	for i := range values {
		values[i] = fermat(valbits[i*(n+1) : (i+1)*(n+1)])
	}
	fourier(values, twisted, false, n, k)
	return polValues{k, n, values}
}

// InvTransform reconstructs a polynomial from its values at
// roots of x^K+1. The m field of the returned polynomial
// is unspecified.
func (v *polValues) InvTransform() poly {
	k := v.k
	n := v.n
	θshift := (n * _W) >> k

	// Perform an inverse Fourier transform to recover q.
	qbits := make([]Word, (n+1)<<k)
	q := make([]fermat, 1<<k)
	for i := range q {
		q[i] = fermat(qbits[i*(n+1) : (i+1)*(n+1)])
	}
	fourier(q, v.values, true, n, k)

	// Divide by K, and untwist q to recover p.
	u := make(fermat, n+1)
	a := make([]nat, 1<<k)
	for i := range q {
		u.Shift(q[i], -int(k)-i*θshift)
		copy(q[i], u)
		a[i] = nat(q[i])
	}
	return poly{k: k, m: 0, a: a}
}

// fourier performs an unnormalized Fourier transform
// of src, a length 1<<k vector of numbers modulo b^n+1
// where b = 1<<_W.
//
// The root of unity used in the transform is ω=1<<ωshift.
// The source array may use shifted indices (i.e. the i-th
// element is src[i << idxShift]).
func fourier(dst []fermat, src []fermat, backward bool, n int, k uint) {
	var rec func(dst, src []fermat, size uint)
	u, v := make(fermat, n+1), make(fermat, n+1) // pre-allocate temporary variables.
	rec = func(dst, src []fermat, size uint) {
		idxShift := k - size
		ωshift := (2 * n * _W) >> size
		if backward {
			ωshift = -ωshift
		}

		// Easy cases.
		if len(src[0]) != n+1 || len(dst[0]) != n+1 {
			panic("len(src[0]) != n+1 || len(dst[0]) != n+1")
		}
		switch size {
		case 0:
			copy(dst[0], src[0])
			return
		case 1:
			dst[0].Add(src[0], src[1<<idxShift]) // dst[0] = src[0] + src[1]
			copy(dst[1], src[1<<idxShift])
			dst[1].Neg()
			dst[1].Add(src[0], dst[1]) // dst[1] = src[0] - src[1]
			return
		}

		// Let P(x) = src[0] + src[1<<idxShift] * x + ... + src[K-1 << idxShift] * x^(K-1)
		// The P(x) = Q1(x²) + x*Q2(x²)
		// where Q1's coefficients are src with indices shifted by 1
		// where Q2's coefficients are src[1<<idxShift:] with indices shifted by 1

		// Split destination vectors in halves.
		dst1 := dst[:1<<(size-1)]
		dst2 := dst[1<<(size-1):]
		// Transform Q1 and Q2 in the halves.
		rec(dst1, src, size-1)
		rec(dst2, src[1<<idxShift:], size-1)

		// Reconstruct P's transform from transforms of Q1 and Q2.
		// dst[i]            is dst1[i] + ω^i * dst2[i]
		// dst[i + 1<<(k-1)] is dst1[i] + ω^(i+K/2) * dst2[i]
		//
		for i := range dst1 {
			copy(u, dst1[i])
			v.Shift(dst2[i], i*ωshift) // ω^i * dst2[i]
			dst1[i].Add(u, v)
			v.Neg()
			dst2[i].Add(u, v)
		}
	}
	rec(dst, src, k)
}

