package characterio

import (
	"bufio"
	"errors"
	"os"
	"strings"

	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/gotree/tree"
)

// type CharacterStruct struct {
// 	Name     string   // name of the variant
// 	V        variant  // detailed info about this variant
// 	StateKey []string // the inner slice might be ["A", "C", "G", "T"] (if all four nucs are present at this site in the alignment)
// }

func countVariantsCSV(annoFile string) ([]CharacterStruct, error) {

	csa := make([]CharacterStruct, 0)

	f, err := os.Open(annoFile)
	defer f.Close()
	if err != nil {
		return []CharacterStruct{}, err
	}

	header := true

	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		fields := strings.Split(line, ",")
		if len(fields) < 2 {
			return []CharacterStruct{}, errors.New("badly formatted tipfile: fewer than two columns in the csv")
		}

		if header {
			for _, f := range fields[1:] {
				csa = append(csa, CharacterStruct{Name: f, StateKey: make([]string, 0)})
			}
			header = false
			continue
		}

		if len(fields[1:]) != len(csa) {
			return []CharacterStruct{}, errors.New("badly formatted tipfile: number of character columns doesn't match the length of the header")
		}

		for i, f := range fields[1:] {
			if stringInArray(f, csa[i].StateKey) || f == "" { //don't type empty strings/missing data
				continue
			} else {
				csa[i].StateKey = append(csa[i].StateKey, f)
			}
		}
	}

	err = s.Err()
	if err != nil {
		return []CharacterStruct{}, err
	}

	return csa, nil
}

// type NodeStates struct {
// 	ID     string   // name of this node (could be a tip)
// 	States [][]byte // bit-encoded character states. Hopefully will parallelise over the first dimension of the array
// }

// for every row in an csv of tip -> character relationships, return the information
// in a map from tip (node) name -> array of bit-encoded character states
func typeVariantsCSV(idx []StartStop, l int, annoFile string, csa []CharacterStruct) ([]NodeStates, error) {

	nsa := make([]NodeStates, 0)

	f, err := os.Open(annoFile)
	defer f.Close()
	if err != nil {
		return []NodeStates{}, err
	}

	header := true

	var bitToSet int
	var start, stop int

	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		fields := strings.Split(line, ",")
		if len(fields) < 2 {
			return []NodeStates{}, errors.New("badly formatted tipfile: fewer than two columns in the csv")
		}

		if header {
			header = false
			continue
		}

		if len(fields[1:]) != len(csa) {
			return []NodeStates{}, errors.New("badly formatted tipfile: number of character columns doesn't match the length of the header")
		}

		tip := fields[0]
		NS := NodeStates{ID: tip}
		NS.States = make([]byte, l, l)

		for i, f := range fields[1:] {
			start = idx[i].Start
			stop = idx[i].Stop
			if f == "" { // don't set bits for empty strings/missing data - handle them later
				continue
			}
			bitToSet, err = stringIndexInArray(f, csa[i].StateKey)
			if err != nil {
				return []NodeStates{}, err
			}
			bitsets.SetBit(NS.States[start:stop], bitToSet)
		}
		nsa = append(nsa, NS)
	}
	err = s.Err()
	if err != nil {
		return make([]NodeStates, 0), err
	}

	return nsa, nil
}

// For a CSV of character( state)s, convert the information to a
// map from tip name -> array of bit-encoded character states
func TypeTipfile(t *tree.Tree, tipfile string) ([]CharacterStruct, []StartStop, [][]byte, error) {

	characterStates, err := countVariantsCSV(tipfile)
	if err != nil {
		return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
	}

	idx, length := getIndex(characterStates)

	nodeStates, err := typeVariantsCSV(idx, length, tipfile, characterStates)
	if err != nil {
		return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
	}

	states := make([][]byte, len(t.Nodes()), len(t.Nodes()))

	for i := range states {
		states[i] = make([]byte, length, length)
	}

	for _, ns := range nodeStates {
		id, err := t.TipId(ns.ID)
		if err != nil {
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		}
		states[id] = ns.States
	}

	return characterStates, idx, states, nil
}
