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
	if testing.Short() {
		N = 2
	}
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

func testFourier(t *testing.T, N int, k uint) {
	t.Logf("testFourier(t, %d, %d)", N, k)
	ωshift := (2 * N * _W) >> k
	src := make([]fermat, 1<<k)
	dst1 := make([]fermat, 1<<k)
	dst2 := make([]fermat, 1<<k)
	// random inputs.
	for i := range src {
		src[i] = make(fermat, N+1)
		dst1[i] = make(fermat, N+1)
		dst2[i] = make(fermat, N+1)
		for p := 0; p < N; p++ {
			src[i][p] = Word(rand.Int63())
		}
	}

	// naive transform
	tmp := make(fermat, N+1)
	for i := range src {
		for j := range dst1 {
			tmp.Shift(src[i], i*j*ωshift)
			dst1[j].Add(dst1[j], tmp)
		}
	}

	// fast transform
	fourier(dst2, src, false, N, k)

	for i := range src {
		if cmpnat(t, nat(dst1[i]), nat(dst2[i])) != 0 {
			var x, y Int
			x.SetBits(dst1[i])
			y.SetBits(dst2[i])
			t.Errorf("difference in dst[%d]: %x %x", i, &x, &y)
		}
	}
}

func TestFourier(t *testing.T) {
	// 1-word transforms.
	testFourier(t, 1, 2)
	testFourier(t, 1, 3)
	testFourier(t, 1, 4)

	// 2-word transforms
	testFourier(t, 2, 2)
	testFourier(t, 2, 3)
	testFourier(t, 2, 4)

	testFourier(t, 4, 4)
	testFourier(t, 4, 5)
	testFourier(t, 4, 6)
}

// Tests Fourier transform and its reciprocal.
func TestRoundTripPolyValues(t *testing.T) {
	Size := 100000
	if testing.Short() {
		Size = 50
	}
	// Build a polynomial from an integer.
	n := make(nat, Size)
	for p := range n {
		n[p] = Word(rand.Int63())
	}
	k, m := fftSize(n, nil)
	pol := polyFromNat(n, k, m)

	// Transform it.
	f := 2*pol.m*_W + int(pol.k)
	f = ((f >> pol.k) + 1) << pol.k
	f /= _W
	values := pol.Transform(f)

	// Inverse transform.
	pol2 := values.InvTransform()
	pol2.m = m

	t.Logf("k=%d, m=%d", k, m)

	// Evaluate and compare.
	n2 := pol2.Int()
	if cmpnat(t, n, n2) != 0 {
		t.Errorf("different n and n2")
	}
}

