package bigfft

import (
	. "math/big"
	"unsafe"
)

const _W = int(unsafe.Sizeof(Word(0)) * 8)

type nat []Word

// returns the FFT length k, m the number of words per chunk
func fftSize(x, y nat) (k uint, m int) {
	words := len(x) + len(y)
	bits := int64(words) * int64(_W)
	switch {
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

// poly represents an integer via a polynomial in Z[x]/(x^K+1)
// where K is the FFT length and b is the computation basis 1<<(m*_W).
// If P = a[0] + a[1] x + ... a[n] x^(K-1), the associated natural number
// is P(b^m).
type poly struct {
	k uint  // k is such that K = 1<<k.
	m int   // the m such that P(b^m) is the original number.
	a []nat // a slice of at most K m-word coefficients.
}

// Int evaluates back a poly to its integer value.
func (p *poly) Int() nat {
	n := make(nat, len(p.a)*p.m+1)
	m := p.m
	np := n
	for i := range p.a {
		c := addVV(np[:m], np[:m], p.a[i])
		np = np[m:]
		if np[0] > 1 {
			panic("impossible")
		}
		np[0] += c
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

