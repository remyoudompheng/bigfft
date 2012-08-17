package bigfft

import (
	"math/big"
	"testing"
)

func parseHex(s string) fermat {
	z := new(big.Int)
	z, ok := z.SetString(s, 0)
	if !ok {
		panic(s)
	}
	return append(fermat(z.Bits()), 0)
}

func compare(t *testing.T, a, b fermat) {
	var x, y big.Int
	x.SetBits(a)
	y.SetBits(b)
	if x.Cmp(&y) != 0 {
		t.Errorf("%x != %x (%x)", &x, &y, new(big.Int).Xor(&x, &y))
	}
}

var (
	f1 = parseHex("0x01223344556677889001223344556778")
	f2 = parseHex("0x677889001223344556777feddccbbaaa") // f1 << 44 mod 2^128+1
)

func TestFermatShift(t *testing.T) {
	z := make(fermat, len(f1))
	z.Shift(f1, 44)
	compare(t, z, f2)
}

type test struct{ a, b, c fermat }

// addTests is a series of mod 2^256+1 tests.
var addTests = []test{
	{
		parseHex("0x5555555555555555555555555555555555555555555555555555555555555555"),
		parseHex("0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaab"),
		parseHex("0x10000000000000000000000000000000000000000000000000000000000000000"),
	},
	{
		parseHex("0x5555555555555555555555555555555555555555555555555555555555555555"),
		parseHex("0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
		parseHex("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"),
	},
	{
		parseHex("0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
		parseHex("0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
		parseHex("0x5555555555555555555555555555555555555555555555555555555555555553"),
	},
}

func TestFermatAdd(t *testing.T) {
	for _, item := range addTests {
		z := make(fermat, len(item.a))
		z = z.Add(item.a, item.b)
		compare(t, z, item.c)
	}
}

var mulTests = []test{
	{ // 3^400 = 3^200 * 3^200
		parseHex("0xc21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ae223c2578e26c486a03"),
		parseHex("0xc21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ae223c2578e26c486a03"),
		parseHex("0x0e65f4d3508036eaca8faa2b8194ace009c863e44bdc040c459a7127bf8bcc62"),
	},
}

func TestFermatMul(t *testing.T) {
	for _, item := range mulTests {
		z := make(fermat, len(item.a))
		z = z.Mul(item.a, item.b)
		compare(t, z, item.c)
	}
}
