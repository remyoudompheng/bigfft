package bigfft

import (
	"fmt"
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
	N := 4
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

func TestFourierSizes(t *testing.T) {
	sizes := []int{1e3, 2e3, 5e3, 10e3, 20e3, 100e3,
		1e6, 2e6, 5e6, 10e6, 20e6, 50e6, 100e6}
	for _, s := range sizes {
		k, m := fftSize(make(nat, s/_W), make(nat, s/_W))
		v := valueSize(k, m)
		t.Logf("bits=%d => FFT size %d, chunk size = %d, value size = %d",
			s, 1<<k, m, v)
		if v > 3*m {
			t.Errorf("FFT word size %d >> input word size %d", v, m)
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
	f := valueSize(k, m)
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

var rnd = rand.New(rand.NewSource(0x43de683f473542af))

func rndNat(n int) nat {
	x := make(nat, n)
	for i := 0; i < n; i++ {
		x[i] = Word(rnd.Int63()<<1 + rnd.Int63n(2))
	}
	return x
}

func TestMul(t *testing.T) {
	sizes := []int{1e3, 5e3, 15e3, 25e3, 70e3, 200e3, 500e3}
	iters := 10
	if testing.Short() {
		iters = 1
	}

	var x, y Int
	for i := 0; i < iters; i++ {
		for _, size1 := range sizes {
			for _, size2 := range sizes {
				x.SetBits(rndNat(size1 / _W))
				y.SetBits(rndNat(size2 / _W))
				z := new(Int).Mul(&x, &y)
				z2 := Mul(&x, &y)
				if z.Cmp(z2) != 0 {
					t.Errorf("z (%d bits) != z2 (%d bits)", z.BitLen(), z2.BitLen())
					logbig(t, new(Int).Xor(z, z2))
				}
			}
		}
	}
}

func logbig(t *testing.T, n *Int) {
	s := fmt.Sprintf("%x", n)
	for len(s) > 64 {
		t.Log(s[:64])
		s = s[64:]
	}
	t.Log(s)
}

func benchmarkMulBig(b *testing.B, sizex, sizey int) {
	mulx := rndNat(sizex / _W)
	muly := rndNat(sizey / _W)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var x, y, z Int
		x.SetBits(mulx)
		y.SetBits(muly)
		z.Mul(&x, &y)
	}
}

func benchmarkMulFFT(b *testing.B, sizex, sizey int) {
	mulx := rndNat(sizex / _W)
	muly := rndNat(sizey / _W)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var x, y Int
		x.SetBits(mulx)
		y.SetBits(muly)
		_ = Mul(&x, &y)
	}
}

func BenchmarkMulBig_1kb(b *testing.B)   { benchmarkMulBig(b, 1e3, 1e3) }
func BenchmarkMulBig_10kb(b *testing.B)  { benchmarkMulBig(b, 1e4, 1e4) }
func BenchmarkMulBig_100kb(b *testing.B) { benchmarkMulBig(b, 1e5, 1e5) }
func BenchmarkMulBig_1Mb(b *testing.B)   { benchmarkMulBig(b, 1e6, 1e6) }
func BenchmarkMulBig_5Mb(b *testing.B)   { benchmarkMulBig(b, 5e6, 5e6) }
func BenchmarkMulBig_10Mb(b *testing.B)  { benchmarkMulBig(b, 10e6, 10e6) }
func BenchmarkMulBig_20Mb(b *testing.B)  { benchmarkMulBig(b, 20e6, 20e6) }

// 50Mb multiplication takes about 1 minute.
//func BenchmarkMulBig_50Mb(b *testing.B)  { benchmarkMulBig(b, 50e6, 50e6) }

func BenchmarkMulFFT_1kb(b *testing.B)   { benchmarkMulFFT(b, 1e3, 1e3) }
func BenchmarkMulFFT_10kb(b *testing.B)  { benchmarkMulFFT(b, 1e4, 1e4) }
func BenchmarkMulFFT_100kb(b *testing.B) { benchmarkMulFFT(b, 1e5, 1e5) }

func BenchmarkMulFFT_1Mb(b *testing.B)   { benchmarkMulFFT(b, 1e6, 1e6) }
func BenchmarkMulFFT_5Mb(b *testing.B)   { benchmarkMulFFT(b, 5e6, 5e6) }
func BenchmarkMulFFT_10Mb(b *testing.B)  { benchmarkMulFFT(b, 10e6, 10e6) }
func BenchmarkMulFFT_20Mb(b *testing.B)  { benchmarkMulFFT(b, 20e6, 20e6) }
func BenchmarkMulFFT_50Mb(b *testing.B)  { benchmarkMulFFT(b, 50e6, 50e6) }
func BenchmarkMulFFT_100Mb(b *testing.B) { benchmarkMulFFT(b, 100e6, 100e6) }
