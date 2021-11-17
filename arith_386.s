// Trampolines to math/big assembly implementations.

#include "textflag.h"

// func addVV(z, x, y []Word) (c Word)
TEXT ·addVV(SB),NOSPLIT,$0
	FUNCDATA $0, ·addVV(SB)
	JMP	math∕big·addVV(SB)

// func subVV(z, x, y []Word) (c Word)
TEXT ·subVV(SB),NOSPLIT,$0
	FUNCDATA $0, ·subVV(SB)
	JMP	math∕big·subVV(SB)

// func addVW(z, x []Word, y Word) (c Word)
TEXT ·addVW(SB),NOSPLIT,$0
	FUNCDATA $0, ·addVW(SB)
	JMP	math∕big·addVW(SB)

// func subVW(z, x []Word, y Word) (c Word)
TEXT ·subVW(SB),NOSPLIT,$0
	FUNCDATA $0, ·subVW(SB)
	JMP	math∕big·subVW(SB)

// func shlVU(z, x []Word, s uint) (c Word)
TEXT ·shlVU(SB),NOSPLIT,$0
	FUNCDATA $0, ·shlVU(SB)
	JMP	math∕big·shlVU(SB)

// func shrVU(z, x []Word, s uint) (c Word)
TEXT ·shrVU(SB),NOSPLIT,$0
	FUNCDATA $0, ·shrVU(SB)
	JMP	math∕big·shrVU(SB)

// func mulAddVWW(z, x []Word, y, r Word) (c Word)
TEXT ·mulAddVWW(SB),NOSPLIT,$0
	FUNCDATA $0, ·mulAddVWW(SB)
	JMP	math∕big·mulAddVWW(SB)

// func addMulVVW(z, x []Word, y Word) (c Word)
TEXT ·addMulVVW(SB),NOSPLIT,$0
	FUNCDATA $0, ·addMulVVW(SB)
	JMP	math∕big·addMulVVW(SB)

