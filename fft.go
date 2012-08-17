package bigfft

import (
	. "math/big"
	"unsafe"
)

const _W = int(unsafe.Sizeof(Word(0)) * 8)

type nat []Word

