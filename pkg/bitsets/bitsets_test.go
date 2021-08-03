package bitsets

import (
	"fmt"
	"testing"
)

// func tempGetByteArraysUnion(aa, ba []byte) ([]byte, error) {
// 	if len(aa) != len(ba) {
// 		return []byte{}, errors.New("different length arrays")
// 	}
// 	ca := make([]byte, len(aa))
// 	for i := range aa {
// 		ca[i] = aa[i] | ba[i]
// 	}
// 	return ca, nil
// }

// func operateOnSlice0(ba []byte, c chan int) {
// 	oa := []byte{8, 9}
// 	newValue, _ := tempGetByteArraysUnion(ba, oa)

// 	sum := 0
// 	for i := range newValue {
// 		sum += int(newValue[i])
// 	}

// 	c <- sum
// }

// func operateOnSlice1(ba []byte, c chan int) {
// 	oa := []byte{1, 7, 84}
// 	newValue, _ := tempGetByteArraysUnion(ba, oa)

// 	sum := 0
// 	for i := range newValue {
// 		sum += int(newValue[i])
// 	}

// 	c <- sum
// }

// const n = 1000

// func TestSliceThreadSafety(t *testing.T) {

// 	var chans [2]chan int
// 	for i := range chans {
// 		chans[i] = make(chan int, n)
// 	}

// 	var s [][]byte
// 	s = append(s, []byte{1, 2})
// 	s = append(s, []byte{99, 0})

// 	var results0 [n]int
// 	var results1 [n]int

// 	for j := 0; j < n; j++ {
// 		go operateOnSlice0(s[0], chans[0])
// 		go operateOnSlice1(s[1], chans[1])
// 	}

// 	for i := 0; i < n; i++ {
// 		results0[i] = <-chans[0]
// 	}

// 	for i := 0; i < n; i++ {
// 		results1[i] = <-chans[1]
// 	}

// 	var test0 int
// 	for i := 0; i < n; i++ {
// 		if i == 0 {
// 			test0 = results0[i]
// 		} else {
// 			if test0 != results0[i] {
// 				t.Errorf("%d", results0[i])
// 			}
// 		}
// 	}

// 	var test1 int
// 	for i := 0; i < n; i++ {
// 		if i == 0 {
// 			test1 = results1[i]
// 		} else {
// 			if test1 != results1[i] {
// 				t.Errorf("%d", results1[i])
// 			}
// 		}
// 	}
// }

// func BenchmarkSetByteArrayBit(b *testing.B) {
// 	b.ReportAllocs()
// 	d := []byte{0, 0}
// 	for i := 0; i < b.N; i++ {
// 		SetBit(d, 7)
// 	}
// }

// func BenchmarkIntersection(b *testing.B) {
// 	b.ReportAllocs()
// 	d := []byte{0, 0}

// 	SetBit(d, 1)
// 	SetBit(d, 2)
// 	SetBit(d, 3)
// 	SetBit(d, 9)

// 	e := []byte{0, 0}

// 	SetBit(e, 1)
// 	SetBit(e, 2)
// 	SetBit(e, 3)
// 	SetBit(e, 10)

// 	for i := 0; i < b.N; i++ {
// 		_ = Intersection(d, e)
// 	}
// }

// func BenchmarkUnion(b *testing.B) {
// 	b.ReportAllocs()
// 	d := []byte{0, 0}

// 	SetBit(d, 1)
// 	SetBit(d, 2)
// 	SetBit(d, 3)
// 	SetBit(d, 9)

// 	e := []byte{0, 0}

// 	SetBit(e, 1)
// 	SetBit(e, 2)
// 	SetBit(e, 3)
// 	SetBit(e, 10)

// 	h := []byte{0, 0}

// 	for i := 0; i < b.N; i++ {
// 		// _, _ = GetByteArraysUnion(d, e)
// 		_ = VarUnion([][]byte{d, e, h})
// 	}
// }

// func BenchmarkMapIntersection(b *testing.B) {
// 	d := make(map[int]bool)
// 	d[1] = true
// 	d[2] = true
// 	d[3] = true
// 	d[9] = true

// 	e := make(map[int]bool)
// 	e[1] = true
// 	e[2] = true
// 	e[3] = true
// 	e[10] = true

// 	f := make(map[int]bool)

// 	for i := 0; i < b.N; i++ {
// 		for k := range d {
// 			if e[k] {
// 				f[k] = true
// 			}
// 		}
// 	}
// }

// func BenchmarkMapUnion(b *testing.B) {
// 	d := make(map[int]bool)
// 	d[1] = true
// 	d[2] = true
// 	d[3] = true
// 	d[9] = true

// 	e := make(map[int]bool)
// 	e[1] = true
// 	e[2] = true
// 	e[3] = true
// 	e[10] = true

// 	f := make(map[int]bool)

// 	for i := 0; i < b.N; i++ {
// 		for k := range d {
// 			f[k] = true
// 		}
// 		for k := range e {
// 			f[k] = true
// 		}
// 	}
// }

// func BenchmarkArrayUnion(b *testing.B) {
// 	a := []int{1, 2, 3, 9}
// 	s := 8

// 	for i := 0; i < b.N; i++ {
// 		_ = intInArray(a, s)
// 	}
// }

func Test_SetBits(t *testing.T) {
	ba := []byte{0}
	SetBit(ba, 1)
	if ba[0] != []byte{128}[0] {
		t.Errorf("error in Test_SetBits")
	}
	ba = []byte{0, 0}
	SetBit(ba, 8)
	SetBit(ba, 9)
	if ba[1] != []byte{128}[0] {
		t.Errorf("error in Test_SetBits")
	}
	if ba[0] != []byte{1}[0] {
		t.Errorf("error in Test_SetBits")
	}
	ba = []byte{0}
	SetBit(ba, 8)
	if ba[0] != []byte{1}[0] {
		t.Errorf("error in Test_SetBits")
	}

	setbits := GetSetBits(ba)
	if len(setbits) != 1 {
		t.Errorf("error in Test_SetBits")
	}
	if setbits[0] != 8 {
		t.Errorf("error in Test_SetBits")
	}

	anysetbit := IsAnyBitSet(ba)
	if !anysetbit {
		t.Errorf("error in Test_SetBits")
	}
	ba = []byte{0}
	anysetbit = IsAnyBitSet(ba)
	if anysetbit {
		t.Errorf("error in Test_SetBits")
	}

	aa := []byte{8}
	if !Different(aa, ba) {
		t.Errorf("error in Test_SetBits")
	}
	if Different(aa, aa) {
		t.Errorf("error in Test_SetBits")
	}
	if Different(ba, ba) {
		t.Errorf("error in Test_SetBits")
	}
}

func Test_SetDiff(t *testing.T) {
	aa := []byte{128}
	ba := []byte{0}
	ca := SetDiff(aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_SetDiff")
	}
	if ca[0] != []byte{128}[0] {
		t.Errorf("error in Test_SetDiff")
	}

	InPlaceSetDiff(ca, aa, ba)
	if ca[0] != []byte{128}[0] {
		t.Errorf("error in Test_SetDiff")
	}
}

func Test_Union(t *testing.T) {
	aa := []byte{128}
	ba := []byte{16}
	ca := Union(aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_Union")
	}
	if ca[0] != []byte{144}[0] {
		t.Errorf("error in Test_Union")
	}
}

func Test_InPlaceUnion(t *testing.T) {
	aa := []byte{128}
	ba := []byte{16}
	ca := []byte{0}
	_ = InPlaceUnion(ca, aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceUnion")
	}
	if ca[0] != []byte{144}[0] {
		t.Errorf("error in Test_InPlaceUnion")
	}
}

func Test_UnionInt(t *testing.T) {
	aa := []byte{64}
	i := 1
	ca := UnionInt(aa, i)
	if len(ca) != 1 {
		t.Errorf("error in Test_Union_Int")
	}
	if ca[0] != []byte{192}[0] {
		t.Errorf("error in Test_Union_Int")
	}
}

func Test_VarUnion(t *testing.T) {
	aaa := make([][]byte, 0)
	aaa = append(aaa, []byte{1})
	aaa = append(aaa, []byte{2})
	aaa = append(aaa, []byte{4})
	aaa = append(aaa, []byte{128})
	ca := VarUnion(aaa)
	if len(ca) != 1 {
		t.Errorf("error in Test_VarUnion")
	}
	if ca[0] != []byte{135}[0] {
		t.Errorf("error in Test_VarUnion")
	}
}

func Test_InPlaceVarUnion(t *testing.T) {
	aaa := make([][]byte, 0)
	aaa = append(aaa, []byte{1})
	aaa = append(aaa, []byte{2})
	aaa = append(aaa, []byte{4})
	aaa = append(aaa, []byte{128})
	ca := []byte{0}
	InPlaceVarUnion(ca, aaa)
	if len(ca) != 1 {
		t.Errorf("error in Test_VarUnion")
	}
	if ca[0] != []byte{135}[0] {
		t.Errorf("error in Test_VarUnion")
	}
}

func Test_Intersection(t *testing.T) {
	aa := []byte{128}
	ba := []byte{16}
	ca := Intersection(aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_Intersection")
	}
	if ca[0] != []byte{0}[0] {
		t.Errorf("error in Test_Intersection")
	}

	aa = []byte{16}
	ba = []byte{144}
	ca = Intersection(aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_Intersection")
	}
	if ca[0] != []byte{16}[0] {
		t.Errorf("error in Test_Intersection")
	}
}

func Test_IntersectionInt(t *testing.T) {
	aa := []byte{192}
	i := 2
	ca := IntersectionInt(aa, i)
	if len(ca) != 1 {
		t.Errorf("error in Test_Intersection_Int")
	}
	if ca[0] != []byte{64}[0] {
		t.Errorf("error in Test_Intersection_Int")
	}
}

func Test_InPlaceIntersection(t *testing.T) {
	aa := []byte{128}
	ba := []byte{16}
	ca := []byte{0}
	InPlaceIntersection(ca, aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceIntersection")
	}
	if ca[0] != []byte{0}[0] {
		t.Errorf("error in Test_InPlaceIntersection")
	}

	aa = []byte{16}
	ba = []byte{144}
	ca = []byte{0}
	InPlaceIntersection(ca, aa, ba)
	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceIntersection")
	}
	if ca[0] != []byte{16}[0] {
		t.Errorf("error in Test_InPlaceIntersection")
	}
}

func Test_SymDiff(t *testing.T) {
	aa := []byte{64}
	ba := []byte{128}
	ca, err := SymDiff(aa, ba)
	if err != nil {
		t.Error(err)
	}
	if len(ca) != 1 {
		t.Errorf("error in Test_SymDiff")
	}
	if ca[0] != []byte{192}[0] {
		t.Errorf("error in Test_SymDiff")
	}
}

func Test_RandomlyChooseSetBit(t *testing.T) {
	aa := []byte{7}
	ca := RandomlyChooseSetBit(aa)
	setbits := GetSetBits(ca)
	if len(setbits) != 1 {
		t.Errorf("error in Test_RandomlyChooseSetBit")
	}
	if setbits[0] != 8 {
		if setbits[0] != 7 {
			if setbits[0] != 6 {
				t.Errorf("error in Test_RandomlyChooseSetBit")
			}
		}
	}
	fmt.Println("err")
}

func Test_InPlaceVarMax(t *testing.T) {
	aaa := make([][]byte, 0)
	aaa = append(aaa, []byte{1})
	aaa = append(aaa, []byte{2})
	aaa = append(aaa, []byte{4})
	aaa = append(aaa, []byte{2})
	aaa = append(aaa, []byte{4})
	ca := []byte{0}

	InPlaceVarMax(ca, aaa)

	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceVarMax")
	}
	if ca[0] != []byte{6}[0] {
		t.Errorf("error in Test_InPlaceVarMax")
	}

	aaa = make([][]byte, 0)
	aaa = append(aaa, []byte{1})
	aaa = append(aaa, []byte{2})
	aaa = append(aaa, []byte{0})
	aaa = append(aaa, []byte{2})
	aaa = append(aaa, []byte{4})
	ca = []byte{6}

	InPlaceVarMax(ca, aaa)

	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceVarMax")
	}
	if ca[0] != []byte{2}[0] {
		t.Errorf("error in Test_InPlaceVarMax")
	}
}

func Test_InPlaceVarCover(t *testing.T) {
	// see: Madison (1989), https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1096-0031.1989.tb00569.x
	aaa := make([][]byte, 6)
	aaa[0] = []byte{0}
	SetBit(aaa[0], 1)
	SetBit(aaa[0], 2)
	aaa[1] = []byte{0}
	SetBit(aaa[1], 4)
	aaa[2] = []byte{0}
	SetBit(aaa[2], 2)
	SetBit(aaa[2], 3)
	SetBit(aaa[2], 4)
	aaa[3] = []byte{0}
	SetBit(aaa[3], 6)
	aaa[4] = []byte{0}
	SetBit(aaa[4], 1)
	SetBit(aaa[4], 2)
	SetBit(aaa[4], 3)
	aaa[5] = []byte{0}
	SetBit(aaa[5], 5)
	SetBit(aaa[5], 6)
	SetBit(aaa[5], 7)
	SetBit(aaa[5], 8)

	ca := []byte{0}

	_ = InPlaceVarCover(ca, aaa)

	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceVarCover")
	}
	if len(GetSetBits(ca)) != len([]int{1, 2, 4, 6}) {
		t.Errorf("error in Test_InPlaceVarCover")
	}
	for i, n := range GetSetBits(ca) {
		if n != []int{1, 2, 4, 6}[i] {
			t.Errorf("error in Test_InPlaceVarCover")
		}
	}
}
