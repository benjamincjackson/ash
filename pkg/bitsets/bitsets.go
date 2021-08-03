package bitsets

import (
	"errors"
	"math/bits"
	"math/rand"
	"time"
)

// Set the 8 - kth bit of a byte to 1. k ∈ {1,2,3,4,5,6,7,8}.
// E.g. if k = 0, then returned byte = byte & 10000000
func setBit(b byte, k int) byte {
	b = b | (1 << (7 - k))
	return b
}

// Set the (1-based) kth bit (from the LHS) in an array of bytes to 1.
// len(ba) * 8 must be >= k
func SetBit(ba []byte, k int) {

	byteindex := (k - 1) / 8
	bitindex := (k - 1) % 8

	ba[byteindex] = setBit(ba[byteindex], bitindex)
}

// is any bit in an array of bytes set, true/false?
func IsAnyBitSet(ba []byte) bool {
	for i := range ba {
		if ba[i] > 0 {
			return true
		}
	}
	return false
}

// which bits in an array of bytes are set? left-most bit of the left-most byte is bit # 1
func GetSetBits(ba []byte) []int {

	states := make([]int, 0)
	counter := 1

	for i := range ba {
		for j := 7; j > -1; j-- {
			if (ba[i]>>j)&1 == 1 {
				states = append(states, counter)
			}
			counter++
		}
	}

	return states
}

// Are the sets different
func Different(aa, ba []byte) bool {
	for i := range aa {
		if aa[i] != ba[i] {
			return true
		}
	}
	return false
}

// Set difference aa - ba
// is the same as aa AND the COMPLEMENT OF ba
func SetDiff(aa, ba []byte) []byte {
	ca := make([]byte, len(aa), len(aa))
	for i := range aa {
		ca[i] = aa[i] & ^ba[i]
	}
	return ca
}

// Set difference aa - ba
// is the same as aa AND the COMPLEMENT OF ba
func InPlaceSetDiff(ca, aa, ba []byte) {
	for i := range aa {
		ca[i] = aa[i] & ^ba[i]
	}
}

// is aa a subset of ba
// TO DO- test this
func IsSubset(aa, ba []byte) bool {
	if IsAnyBitSet(SetDiff(aa, ba)) {
		return false
	}
	return true
}

// Get the intersection of set bits from two equal length arrays of bytes.
// Allocates and returns a new byte array which contains the intersection of set bits
func Intersection(aa, ba []byte) []byte {
	ca := make([]byte, len(aa), len(aa))
	for i := range aa {
		ca[i] = aa[i] & ba[i]
	}
	return ca
}

// Get the intersection of set bits from two equal length arrays of bytes.
// Stores the result in ca, which has already been allocated
func InPlaceIntersection(ca, aa, ba []byte) error {
	if len(aa) != len(ba) {
		return errors.New("different length bitsets")
	}
	for i := range aa {
		ca[i] = aa[i] & ba[i]
	}
	return nil
}

// Get the intersection from an integer (representing a set bit) and a byte array.
// Allocates and returns a new byte array which contains the intersection
func IntersectionInt(ba []byte, k int) []byte {
	aa := make([]byte, len(ba), len(ba))
	SetBit(aa, k)

	ca := make([]byte, len(ba))
	for i := range ba {
		ca[i] = aa[i] & ba[i]
	}
	return ca
}

func VarIntersection(aaa [][]byte) []byte {
	ba := aaa[0]
	for _, aa := range aaa {
		for j, a := range aa {
			ba[j] = ba[j] & a
		}
	}
	return ba
}

// Get the union of set bits from two equal length arrays of bytes.
// Allocates and returns a new byte array which contains the union of set bits
func Union(aa, ba []byte) []byte {
	ca := make([]byte, len(aa), len(aa))
	for i := range aa {
		ca[i] = aa[i] | ba[i]
	}
	return ca
}

// Get the union of set bits from two equal length arrays of bytes.
// Modifies the third array of bytes in place which contains the union of set bits
func InPlaceUnion(ca, aa, ba []byte) error {
	if len(aa) != len(ba) {
		return errors.New("different length bitsets")
	}
	for i := range aa {
		ca[i] = aa[i] | ba[i]
	}
	return nil
}

// Get the union from an integer (representing a set bit) and a byte array.
// Allocates and returns a new byte array which contains the union
func UnionInt(ba []byte, k int) []byte {
	aa := make([]byte, len(ba), len(ba))
	SetBit(aa, k)

	ca := make([]byte, len(ba), len(ba))
	for i := range ba {
		ca[i] = aa[i] | ba[i]
	}
	return ca
}

// Get the Union from an arbitrary number of bitsets (which must be of the same length).
// Allocates and returns a new byte array which is the Union.
func VarUnion(aaa [][]byte) []byte {
	ba := make([]byte, len(aaa[0]), len(aaa[0]))
	for _, aa := range aaa {
		for j, a := range aa {
			ba[j] = ba[j] | a
		}
	}
	return ba
}

// Get the Union from an arbitrary number of bitsets (which must be of the same length).
// Modifies an array of bytes in place which contains the union of set bits
func InPlaceVarUnion(ba []byte, aaa [][]byte) {
	for _, aa := range aaa {
		for j, b := range aa {
			ba[j] = ba[j] | b
		}
	}
}

// Get all the bits that are set in either A or B but not both (using XOR)
// Allocates and returns a new byte array which contains the symmetric difference
func SymDiff(aa, ba []byte) ([]byte, error) {
	if len(aa) != len(ba) {
		return []byte{}, errors.New("different length bitsets")
	}
	ca := make([]byte, len(ba), len(ba))
	for i := range aa {
		ca[i] = aa[i] ^ ba[i]
	}
	return ca, nil
}

// Randomly choose one bit out of all the set bits in ba
// Allocates and returns a byte array which only has the chosen bit, set
func RandomlyChooseSetBit(ba []byte) []byte {
	ca := make([]byte, len(ba), len(ba))
	setbits := GetSetBits(ba)
	if len(setbits) == 0 {
		return ca
	}
	rand.Seed(time.Now().UnixNano())
	chosenOne := setbits[rand.Intn(len(setbits))]
	SetBit(ca, chosenOne)
	return ca
}

// if the intersection is not an empty set, take it, else take the union
// TO DO - test this
func fitch(aa, ba []byte) []byte {
	if IsAnyBitSet(Intersection(aa, ba)) {
		return Intersection(aa, ba)
	} else {
		return Union(aa, ba)
	}
}

// ThreeSetMPR
// A = [(aa ⊗ ba) ⊗ ca] ∩ [(aa ⊗ ca) ⊗ ba] ∩ [(ba ⊗ ca) ⊗ aa]
// TO DO - test this
func ThreeSetMPR(aa, ba, ca []byte) []byte {
	d := fitch(fitch(aa, ba), ca)
	e := fitch(fitch(aa, ca), ba)
	f := fitch(fitch(ba, ca), aa)

	A := VarIntersection([][]byte{d, e, f})

	return A
}

func InPlaceThreeSetMPR(da, aa, ba, ca []byte) {
	result := ThreeSetMPR(aa, ba, ca)
	for i := range result {
		da[i] = result[i]
	}
}

// Get the most common set bit(s) from a variable number of byte arrays
// This is for treating polytomies as hard - Maddison 1989
// or for dealing with the downpass when there are polytomies?
// TO DO: try to optimise this (don't use a map?)
// func InPlaceVarMax(toSet []byte, args [][]byte) error {
// 	// first, check that the bitsets are the same length
// 	l := 0
// 	for i, a := range args {
// 		if i == 0 {
// 			l = len(a)
// 			continue
// 		}
// 		if len(a) != l {
// 			return errors.New("different length bitsets")
// 		}
// 	}

// 	var states []int
// 	stateCounts := make(map[int]int)

// 	// populate a map with the counts of the states
// 	for _, a := range args {
// 		states = GetSetBits(a)
// 		for _, s := range states {
// 			if _, ok := stateCounts[s]; ok {
// 				stateCounts[s]++
// 			} else {
// 				stateCounts[s] = 1
// 			}
// 		}
// 	}

// 	// IF THE STATES ARE EMPTY HERE, WHAT DO WE DO?
// 	// if len(stateCounts) == 0 {

// 	// }

// 	// find the most common state(s)
// 	maxes := make([]int, 0)
// 	max := -1
// 	for k, v := range stateCounts {
// 		if v > max {
// 			max = v
// 			maxes = []int{k}
// 		} else if v == max {
// 			maxes = append(maxes, k)
// 		}
// 	}

// 	for i := range toSet {
// 		toSet[i] = 0
// 	}
// 	for _, k := range maxes {
// 		SetBit(toSet, k)
// 	}

// 	return nil
// }

func InPlaceVarMax(toSet []byte, args [][]byte) {
	statecounts := make([]int, len(args[0])*8, len(args[0])*8)
	for _, a := range args {
		states := GetSetBits(a)
		for _, s := range states {
			statecounts[s-1]++
		}
	}

	maxes := make([]int, 0)
	max := -1
	for i, v := range statecounts {
		if v > max {
			max = v
			maxes = []int{i + 1}
		} else if v == max {
			maxes = append(maxes, i+1)
		}
	}

	for i := range toSet {
		toSet[i] = 0
	}
	for _, k := range maxes {
		SetBit(toSet, k)
	}
}

func VarMax2(args [][]byte) []byte {
	statecounts := make([]int, len(args[0])*8, len(args[0])*8)
	for _, a := range args {
		states := GetSetBits(a)
		for _, s := range states {
			statecounts[s-1]++
		}
	}

	maxes := make([]int, 0)
	max := -1
	for i, v := range statecounts {
		if v > max {
			max = v
			maxes = []int{i + 1}
		} else if v == max {
			maxes = append(maxes, i+1)
		}
	}

	ba := make([]byte, len(args[0]), len(args[0]))
	for _, k := range maxes {
		SetBit(ba, k)
	}

	return ba
}

// helper function for variadicByteArraysCover()
func intInArray(ia []int, s int) bool {
	for _, k := range ia {
		if k == s {
			return true
		}
	}
	return false
}

//// helper function for variadicByteArraysCover()
//// modified from: https://github.com/mxschmitt/golang-combinations/blob/v1.1.0/combinations.go#L32
//// distributed under MIT licence ////
// Combinations returns combinations of n elements for a given string (Ben: changed to int) array.
// For n < 1, it equals to All and returns all combinations.
func combinations(set []int, n int) (subsets [][]int) {
	length := uint(len(set))

	if n > len(set) {
		n = len(set)
	}

	// Go through all possible combinations of objects
	// from 1 (only first object in subset) to 2^length (all objects in subset)
	for subsetBits := 1; subsetBits < (1 << length); subsetBits++ {
		if n > 0 && bits.OnesCount(uint(subsetBits)) != n {
			continue
		}

		var subset []int

		for object := uint(0); object < length; object++ {
			// checks if object is contained in subset
			// by checking if bit 'object' is set in subsetBits
			if (subsetBits>>object)&1 == 1 {
				// add object to subset
				subset = append(subset, set[object])
			}
		}
		// add subset to subsets
		subsets = append(subsets, subset)
	}
	return subsets
}

func checkCoverage(coverer []int, coverees [][]byte) (bool, error) {
	// first, check that the bitsets are the same length
	l := 0
	for i, a := range coverees {
		if i == 0 {
			l = len(a)
			continue
		}
		if len(a) != l {
			return false, errors.New("different length bitsets")
		}
	}

	testArray := make([]byte, len(coverees[0]))
	for _, k := range coverer {
		SetBit(testArray, k)
	}

	for _, a := range coverees {
		testArray = Intersection(testArray, a)
		test := false
		for _, ta := range testArray {
			if ta != 0 {
				test = true
			}
		}
		if !test {
			return false, nil
		}
	}

	return true, nil
}

// get the union of the smallest covering sets
// see: Madison (1989), https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1096-0031.1989.tb00569.x
func InPlaceVarCover(toSet []byte, args [][]byte) error {
	// first, check that the bitsets are the same length
	l := 0
	for i, a := range args {
		if i == 0 {
			l = len(a)
			continue
		}
		if len(a) != l {
			return errors.New("different length bitsets")
		}
	}

	// first, we check for singletons, because these must be in the final set
	var states []int
	var singletons []int

	// for every bitset, if the bitset only has one bit set, add this state to the array of singletons
	for _, a := range args {
		states = GetSetBits(a)
		if len(states) == 1 {
			if !intInArray(singletons, states[0]) {
				singletons = append(singletons, states[0])
			}
		}
	}

	// next, we don't further consider any of the bitsets (args) that contain any of the singletons
	newargs := make([][]byte, 0)
	var keep bool

	// for every bitset, if it contains a singleton state, don't add it to the set of bitsets we need to
	// consider for the next step
	for _, a := range args {
		keep = true
		states = GetSetBits(a)
		for _, s := range singletons {
			if intInArray(states, s) {
				keep = false
				break
			}
		}
		if keep {
			newargs = append(newargs, a)
		}
	}

	// for i := range newargs {
	// 	fmt.Println(GetSetBits(newargs[i]))
	// }
	// fmt.Println()

	// then, find the smallest set that covers the remaining bitsets
	//  * firstly, we need all the remaining states.
	//    this the union of all the remaining bitsets
	remainingstates := make([]int, 0)
	for _, a := range newargs {
		states = GetSetBits(a)
		for _, s := range states {
			if !intInArray(remainingstates, s) {
				remainingstates = append(remainingstates, s)
			}
		}
	}

	// fmt.Println(remainingstates)
	// fmt.Println()

	//  * secondly, we try combinations of these (remainingstates) to see if they
	//    cover the non-singleton bitsets (newargs).
	//    we start with all singles, then all pairs, then all triplets, etc... of states
	var doescover []int
	var stilltotest [][]int
	for n := 1; n < len(remainingstates)+1; n++ {
		// get the combinatations for of states for this value of n:
		combs := combinations(remainingstates, n)
		// then test each combination in turn
		for i, c := range combs {
			covers, err := checkCoverage(c, newargs)
			if err != nil {
				return err
			}
			// if it covers the sets, then we know that this number of states == the smallest set,
			// and we know that we only have to test the remaining combinations of this size,
			// so we keep those (stilltotest) and break out of both loops
			if covers {
				doescover = c
				if (i + 1) < len(combs) {
					stilltotest = combs[(i + 1):]
				}
				// break from the outer loop when we get there:
				n = len(remainingstates) + 1
				// and break from this loop
				break
			}
		}
	}

	// fmt.Println()
	// fmt.Println(doescover)
	// fmt.Println(stilltotest)

	finalSets := make([][]byte, 0)
	finalSets = append(finalSets, make([]byte, len(args[0])))
	// set the bits of the comb that certainly covers newargs
	for _, k := range doescover {
		SetBit(finalSets[0], k)
	}

	//  * thirdly, if there is anything in stilltotest, check these as coverers too,
	//    and if they are, take the union of the whole lot
	for _, totest := range stilltotest {
		covers, err := checkCoverage(totest, newargs)
		if err != nil {
			return err
		}
		// and set the bits of the other combs of the same size if they pass the coverage test
		if covers {
			// fmt.Println(totest)
			finalSets = append(finalSets, make([]byte, len(args[0])))
			for _, k := range totest {
				SetBit(finalSets[len(finalSets)-1], k)
			}
		}
	}

	// now take the union of all these guys
	InPlaceVarUnion(toSet, finalSets)

	// and add the singletons by setting their bits:
	for _, k := range singletons {
		SetBit(toSet, k)
	}

	// and we're done, jesus.
	return nil
}

// func main() {

// 	d1 := []byte{0}
// 	SetByteArrayBit(d1, 1)
// 	SetByteArrayBit(d1, 2)

// 	d2 := []byte{0}
// 	SetByteArrayBit(d2, 1)
// 	SetByteArrayBit(d2, 2)
// 	SetByteArrayBit(d2, 3)

// 	d3 := []byte{0}
// 	SetByteArrayBit(d3, 2)
// 	SetByteArrayBit(d3, 3)
// 	SetByteArrayBit(d3, 4)

// 	d4 := []byte{0}
// 	SetByteArrayBit(d4, 4)

// 	d5 := []byte{0}
// 	SetByteArrayBit(d5, 5)
// 	SetByteArrayBit(d5, 6)
// 	SetByteArrayBit(d5, 7)

// 	d6 := []byte{0}
// 	SetByteArrayBit(d6, 6)

// 	// d := []byte{0, 0}

// 	// SetByteArrayBit(d, 1)
// 	// SetByteArrayBit(d, 2)
// 	// SetByteArrayBit(d, 3)
// 	// SetByteArrayBit(d, 9)

// 	// e := []byte{0, 0}

// 	// SetByteArrayBit(e, 1)
// 	// SetByteArrayBit(e, 2)
// 	// SetByteArrayBit(e, 3)
// 	// SetByteArrayBit(e, 9)

// 	// f := []byte{0, 0}

// 	// SetByteArrayBit(f, 1)

// 	// v, _ := variadicByteArraysMax(d, e, f)

// 	// fmt.Println(v)

// 	// f, _ := GetByteArraysIntersection(d, e)

// 	// if !isAnyBitSet(f) {
// 	// 	GetByteArraysUnion(d, e, f)
// 	// }

// 	// fmt.Println(f)
// 	// fmt.Println(getSetBits(f))

// 	// // g, _ := GetByteArraysUnion(d, e)
// 	// h := []byte{0, 0}
// 	// _ = GetByteArraysUnion(d, e, h)

// 	// fmt.Println(f)
// 	// fmt.Println(h)

// 	// var b byte = 1

// 	// fmt.Println(b)

// 	// c := setbit(b, 4)

// 	// fmt.Println(c)

// 	// fmt.Println(b & c)
// }
