package bigfft

import (
	. "math/big"
	"math/rand"
	"testing"
)

func cmpnat(t *testing.T, x, y nat) int {
	var a, b Int
	a.SetBits(x)
	b.SetBits(y)
	c := a.Cmp(&b)
	if c != 0 {
		t.Logf("a.len=%d, b.len=%d", a.BitLen(), b.BitLen())
		for i := 0; i < len(x) || i < len(y); i++ {
			var u, v Word
			if i < len(x) {
				u = x[i]
			}
			if i < len(y) {
				v = y[i]
			}
			if diff := u ^ v; diff != 0 {
				t.Logf("diff at word %d: %x", i, diff)
			}
		}
	}
	return c
}

func TestRoundTripIntPoly(t *testing.T) {
	N := 10
	step := 500
      // Sizes 12800 and 34300 may cause problems.
	for size := 300; size < 50000; size += step {
            n := make(nat, size)
		for i := 0; i < N; i++ {
			for p := range n {
				n[p] = Word(rand.Int63())
			}
			k, m := fftSize(n, nil)
			pol := polyFromNat(n, k, m)
			n2 := pol.Int()
			if cmpnat(t, n, n2) != 0 {
				t.Errorf("different n and n2, size=%d, iter=%d", size, i)
			}
		}
	}
}
