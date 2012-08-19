// Usage: go test -run=TestCalibrate -calibrate

package bigfft

import (
	"flag"
	"fmt"
	"testing"
	"time"
)

var calibrate = flag.Bool("calibrate", false, "run calibration test")

func measureMul(th int) (tBig, tFFT time.Duration) {
	bigLoad := func(b *testing.B) { benchmarkMulBig(b, th, th) }
	fftLoad := func(b *testing.B) { benchmarkMulFFT(b, th, th) }

	res1 := testing.Benchmark(bigLoad)
	res2 := testing.Benchmark(fftLoad)
	tBig = time.Duration(res1.NsPerOp())
	tFFT = time.Duration(res2.NsPerOp())
	return
}

func TestCalibrate(t *testing.T) {
	if !*calibrate {
		t.Log("not calibrating, use -calibrate to do so.")
		return
	}

	lower := int(1e3) // math/big is faster at this size.
	upper := int(1e6) // FFT is really faster at this size.

	big, fft := measureMul(lower)
	lowerX := float64(big) / float64(fft)
	fmt.Printf("speedup at size %d: %.2f\n", lower, lowerX)
	big, fft = measureMul(upper)
	upperX := float64(big) / float64(fft)
	fmt.Printf("speedup at size %d: %.2f\n", upper, upperX)
	for {
		mid := (lower + upper) / 2
		big, fft := measureMul(mid)
		X := float64(big) / float64(fft)
		fmt.Printf("speedup at size %d: %.2f\n", mid, X)
		switch {
		case X < 0.9:
			lower = mid
			lowerX = X
		case X > 1.1:
			upper = mid
			upperX = X
		default:
			fmt.Printf("speedup at size %d: %.2f\n", lower, lowerX)
			fmt.Printf("speedup at size %d: %.2f\n", upper, upperX)
			return
		}
	}
}
