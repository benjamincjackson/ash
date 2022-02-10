package bitsets

import (
	"errors"
	"testing"
)

func tempGetByteArraysUnion(aa, ba []byte) ([]byte, error) {
	if len(aa) != len(ba) {
		return []byte{}, errors.New("different length arrays")
	}
	ca := make([]byte, len(aa))
	for i := range aa {
		ca[i] = aa[i] | ba[i]
	}
	return ca, nil
}

func operateOnSlice0(ba []byte, c chan int) {
	oa := []byte{8, 9}
	newValue, _ := tempGetByteArraysUnion(ba, oa)

	sum := 0
	for i := range newValue {
		sum += int(newValue[i])
	}

	c <- sum
}

func operateOnSlice1(ba []byte, c chan int) {
	oa := []byte{1, 7, 84}
	newValue, _ := tempGetByteArraysUnion(ba, oa)

	sum := 0
	for i := range newValue {
		sum += int(newValue[i])
	}

	c <- sum
}

const n = 1000

func TestSliceThreadSafety(t *testing.T) {

	var chans [2]chan int
	for i := range chans {
		chans[i] = make(chan int, n)
	}

	var s [][]byte
	s = append(s, []byte{1, 2})
	s = append(s, []byte{99, 0})

	var results0 [n]int
	var results1 [n]int

	for j := 0; j < n; j++ {
		go operateOnSlice0(s[0], chans[0])
		go operateOnSlice1(s[1], chans[1])
	}

	for i := 0; i < n; i++ {
		results0[i] = <-chans[0]
	}

	for i := 0; i < n; i++ {
		results1[i] = <-chans[1]
	}

	var test0 int
	for i := 0; i < n; i++ {
		if i == 0 {
			test0 = results0[i]
		} else {
			if test0 != results0[i] {
				t.Errorf("%d", results0[i])
			}
		}
	}

	var test1 int
	for i := 0; i < n; i++ {
		if i == 0 {
			test1 = results1[i]
		} else {
			if test1 != results1[i] {
				t.Errorf("%d", results1[i])
			}
		}
	}
}

func BenchmarkSetByteArrayBit(b *testing.B) {
	b.ReportAllocs()
	d := []byte{0, 0}
	for i := 0; i < b.N; i++ {
		SetBit(d, 7)
	}
}

func BenchmarkIntersection(b *testing.B) {
	b.ReportAllocs()
	d := []byte{0, 0}

	SetBit(d, 1)
	SetBit(d, 2)
	SetBit(d, 3)
	SetBit(d, 9)

	e := []byte{0, 0}

	SetBit(e, 1)
	SetBit(e, 2)
	SetBit(e, 3)
	SetBit(e, 10)

	for i := 0; i < b.N; i++ {
		_ = Intersection(d, e)
	}
}

func BenchmarkUnion(b *testing.B) {
	b.ReportAllocs()
	d := []byte{0, 0}

	SetBit(d, 1)
	SetBit(d, 2)
	SetBit(d, 3)
	SetBit(d, 9)

	e := []byte{0, 0}

	SetBit(e, 1)
	SetBit(e, 2)
	SetBit(e, 3)
	SetBit(e, 10)

	h := []byte{0, 0}

	for i := 0; i < b.N; i++ {
		// _, _ = GetByteArraysUnion(d, e)
		_ = VarUnion([][]byte{d, e, h})
	}
}

func BenchmarkMapIntersection(b *testing.B) {
	d := make(map[int]bool)
	d[1] = true
	d[2] = true
	d[3] = true
	d[9] = true

	e := make(map[int]bool)
	e[1] = true
	e[2] = true
	e[3] = true
	e[10] = true

	f := make(map[int]bool)

	for i := 0; i < b.N; i++ {
		for k := range d {
			if e[k] {
				f[k] = true
			}
		}
	}
}

func BenchmarkMapUnion(b *testing.B) {
	d := make(map[int]bool)
	d[1] = true
	d[2] = true
	d[3] = true
	d[9] = true

	e := make(map[int]bool)
	e[1] = true
	e[2] = true
	e[3] = true
	e[10] = true

	f := make(map[int]bool)

	for i := 0; i < b.N; i++ {
		for k := range d {
			f[k] = true
		}
		for k := range e {
			f[k] = true
		}
	}
}

func BenchmarkArrayUnion(b *testing.B) {
	a := []int{1, 2, 3, 9}
	s := 8

	for i := 0; i < b.N; i++ {
		_ = intInArray(a, s)
	}
}
