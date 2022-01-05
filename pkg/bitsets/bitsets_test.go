package bitsets

import (
	"fmt"
	"reflect"
	"testing"
)

func Test_SetBit(t *testing.T) {
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

func Test_IsAnyBitSet(t *testing.T) {
	aa := []byte{0, 128}
	ba := []byte{0, 0}

	if !IsAnyBitSet(aa) {
		t.Errorf("error in Test_IsAnyBitSet")
	}

	if IsAnyBitSet(ba) {
		t.Errorf("error in Test_IsAnyBitSet")
	}
}

func Test_GetSetBits(t *testing.T) {
	aa := []byte{128}
	ba := []byte{0}
	ca := []byte{128, 128}

	if !reflect.DeepEqual(GetSetBits(aa), []int{1}) {
		t.Errorf("error in Test_GetSetBits")
	}

	if !reflect.DeepEqual(GetSetBits(ba), []int{}) {
		t.Errorf("error in Test_GetSetBits")
	}

	if !reflect.DeepEqual(GetSetBits(ca), []int{1, 9}) {
		t.Errorf("error in Test_GetSetBits")
	}
}

func Test_Different(t *testing.T) {
	aa := []byte{128}
	ba := []byte{0}
	ca := []byte{128, 128}
	da := []byte{128, 128}

	if !Different(aa, ba) {
		t.Errorf("error in Test_Different")
	}

	if Different(ca, da) {
		t.Errorf("error in Test_Different")
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
}

func Test_InPlaceSetDiff(t *testing.T) {
	aa := []byte{128}
	ba := []byte{0}
	ca := make([]byte, 1)

	InPlaceSetDiff(ca, aa, ba)

	if len(ca) != 1 {
		t.Errorf("error in Test_InPlaceSetDiff")
	}
	if ca[0] != []byte{128}[0] {
		t.Errorf("error in Test_InPlaceSetDiff")
	}
}

func Test_IsSubset(t *testing.T) {
	aa := []byte{255}
	ba := []byte{3}
	ca := []byte{4}

	if !IsSubset(ba, aa) {
		t.Errorf("error in Test_IsSubset")
	}

	if IsSubset(ca, ba) {
		t.Errorf("error in Test_IsSubset")
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

func Test_VarIntersection(t *testing.T) {
	aaa := make([][]byte, 3)
	aaa[0] = []byte{1}
	aaa[1] = []byte{2}
	aaa[2] = []byte{4}

	ba := VarIntersection(aaa)

	if ba[0] != []byte{0}[0] {
		t.Errorf("error in Test_VarIntersection")
	}

	aaa = make([][]byte, 3)
	aaa[0] = []byte{255}
	aaa[1] = []byte{1}
	aaa[2] = []byte{7}

	ba = VarIntersection(aaa)

	if ba[0] != []byte{1}[0] {
		t.Errorf("error in Test_VarIntersection")
	}

	aaa = make([][]byte, 3)
	aaa[0] = []byte{255, 255}
	aaa[1] = []byte{1, 1}
	aaa[2] = []byte{0, 1}

	ba = VarIntersection(aaa)

	if len(ba) != 2 {
		t.Errorf("error in Test_VarIntersection")
	}
	if ba[0] != 0 {
		t.Errorf("error in Test_VarIntersection")
	}
	if ba[1] != 1 {
		t.Errorf("error in Test_VarIntersection")
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
	InPlaceUnion(ca, aa, ba)
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

func Test_SymDiff(t *testing.T) {
	aa := []byte{64}
	ba := []byte{128}
	ca := SymDiff(aa, ba)
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

// if the intersection is not an empty set, take it, else take the union
func Test_fitch(t *testing.T) {
	aa := []byte{1, 0}
	ba := []byte{1, 255}
	ca := fitch(aa, ba)
	if len(ca) != 2 {
		t.Errorf("error in Test_fitch")
	}
	if ca[0] != 1 {
		t.Errorf("error in Test_fitch")
	}
	if ca[1] != 0 {
		t.Errorf("error in Test_fitch")
	}

	aa = []byte{1, 0}
	ba = []byte{0, 255}
	ca = fitch(aa, ba)
	if len(ca) != 2 {
		t.Errorf("error in Test_fitch")
	}
	if ca[0] != 1 {
		t.Errorf("error in Test_fitch")
	}
	if ca[1] != 255 {
		t.Errorf("error in Test_fitch")
	}
}

// A = [(aa ⊗ ba) ⊗ ca] ∩ [(aa ⊗ ca) ⊗ ba] ∩ [(ba ⊗ ca) ⊗ aa]
// "⊗" is the fitch operation
func Test_ThreeSetMPR(t *testing.T) {
	aa := []byte{1}
	ba := []byte{2}
	ca := []byte{1}

	if !reflect.DeepEqual(ThreeSetMPR(aa, ba, ca), []byte{1}) {
		t.Errorf("error in Test_ThreeSetMPR")
	}
}

func Test_InPlaceThreeSetMPR(t *testing.T) {
	aa := []byte{1}
	ba := []byte{2}
	ca := []byte{1}

	da := make([]byte, 1)
	InPlaceThreeSetMPR(da, aa, ba, ca)

	if !reflect.DeepEqual(da, []byte{1}) {
		t.Errorf("error in Test_ThreeSetMPR")
	}
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

	aaa = make([][]byte, 0)
	aaa = append(aaa, []byte{1, 1})
	aaa = append(aaa, []byte{2, 1})
	aaa = append(aaa, []byte{1, 1})
	aaa = append(aaa, []byte{255, 255})
	aaa = append(aaa, []byte{2, 1})
	ca = []byte{0, 0}

	InPlaceVarMax(ca, aaa)

	if len(ca) != 2 {
		t.Errorf("error in Test_InPlaceVarMax")
	}
	if !reflect.DeepEqual(ca, []byte{0, 1}) {
		t.Errorf("error in Test_InPlaceVarMax")
	}
}

func Test_intInArray(t *testing.T) {
	if !intInArray([]int{1, 2, 3}, 1) {
		t.Errorf("error in Test_intInArray")
	}
	if intInArray([]int{1, 2, 3}, 4) {
		t.Errorf("error in Test_intInArray")
	}
}

func Test_combinations(t *testing.T) {
	set := []int{1, 2, 3}
	n := 2

	if !reflect.DeepEqual(combinations(set, n), [][]int{{1, 2}, {1, 3}, {2, 3}}) {
		t.Errorf("error in Test_combinations")
	}
}

// func Test_checkCoverage(t *testing.T) {

// }

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

	InPlaceVarCover(ca, aaa)

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
