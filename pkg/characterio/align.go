package characterio

import (
	"bufio"
	"errors"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/cov-ert/gofasta/pkg/alphabet"
	"github.com/cov-ert/gofasta/pkg/fastaio"
	"github.com/cov-ert/gofasta/pkg/genbank"

	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/ash/pkg/tree"
)

// CharacterStruct is the info for one character (e.g. one column in the alignment)
type CharacterStruct struct {
	Name     string   // name of the variant
	V        variant  // detailed info about this variant
	StateKey []string // the slice might be ["A", "C", "G", "T"] (if all four nucs are present at this site in the alignment)
}

// variant contains information about a character used for typing it in an alignment
type variant struct {
	vtype   string // "snp" or "aa" or "del"
	vpos    int    // 1-based start in alignment
	vgene   string // for amino acids only
	vres    int    // for amino acids only (residue number)
	vlength int    // for deletions only
}

type NodeStates struct {
	ID     string // name of this node (could be a tip)
	States []byte // bit-encoded character states. Hopefully will parallelise over the first dimension of the array
}

func getCDSPosFromAnnotation(gb genbank.Genbank) map[string]string {

	m := make(map[string]string)

	for _, F := range gb.FEATURES {
		if F.Feature == "CDS" {
			m[strings.ToLower(F.Info["gene"])] = F.Pos
		}
	}

	return m
}

// get the 1-based start position of the amino acid residue in question
func getAAStartPos(residue_pos int, nuc_pos_string string) (int, error) {

	p := 0
	A := make([]int, 0)

	if nuc_pos_string[0:4] == "join" {
		nuc_pos_string = strings.TrimLeft(nuc_pos_string, "join(")
		nuc_pos_string = strings.TrimRight(nuc_pos_string, ")")
		ranges := strings.Split(nuc_pos_string, ",")
		for _, x := range ranges {
			y := strings.Split(x, "..")
			for _, z := range y {
				temp, err := strconv.Atoi(z)
				if err != nil {
					return 0, err
				}
				A = append(A, temp)
			}
		}
	} else {
		y := strings.Split(nuc_pos_string, "..")
		for _, z := range y {
			temp, err := strconv.Atoi(z)
			if err != nil {
				return 0, err
			}
			A = append(A, temp)
		}
	}

	// if the length of A is not a non-zero multiple of 2, then something has gone wrong
	if len(A)%2 != 0 || len(A) == 0 {
		return 0, errors.New("Error parsing CDS positions")
	}

	if len(A)/2 > 1 {
		for i := 0; i < len(A); i += 2 {
			if A[i]+((residue_pos-1)*3) > A[i+1] {
				residue_pos = residue_pos - ((A[i+1] - A[i] + 1) / 3)
				continue
			} else {
				p = A[i] + ((residue_pos - 1) * 3)
				break
			}
		}
	} else {
		p = A[0] + ((residue_pos - 1) * 3)
	}

	// if the residue pos is either still 0 or out of range, the something has gone wrong
	if p == 0 || p > A[len(A)-1]-2 {
		return 0, errors.New("Error parsing CDS positions2")
	}

	return p, nil
}

// read a config file of variants to type and take also as input information about the positions
// of CDSes, and return an array of variant information structs that will be used to type the alignment
func readConfig(configFile string, cdspos map[string]string) ([]CharacterStruct, error) {

	csa := make([]CharacterStruct, 0)

	f, err := os.Open(configFile)
	if err != nil {
		return []CharacterStruct{}, err
	}
	defer f.Close()

	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		fields := strings.Split(line, ":")

		if len(fields) < 2 {
			return []CharacterStruct{}, errors.New("could not parse config line (too few fields): " + line)
		}

		switch fields[0] {
		case "aa":
			gene := strings.ToLower(fields[1])
			residuepos, err := strconv.Atoi(fields[2])
			if err != nil {
				return []CharacterStruct{}, err
			}
			cds_pos_string, ok := cdspos[gene]
			if !ok {
				return []CharacterStruct{}, errors.New("could not parse config line: " + line)
			}
			pos, err := getAAStartPos(residuepos, cds_pos_string)
			if err != nil {
				return []CharacterStruct{}, err
			}
			csa = append(csa, CharacterStruct{V: variant{vtype: "aa", vgene: gene, vpos: pos, vres: residuepos}})
		case "nuc":
			pos, err := strconv.Atoi(fields[1])
			if err != nil {
				return []CharacterStruct{}, err
			}
			csa = append(csa, CharacterStruct{V: variant{vtype: "nuc", vpos: pos}})
		case "del":
			pos, err := strconv.Atoi(fields[1])
			if err != nil {
				return []CharacterStruct{}, err
			}
			length, err := strconv.Atoi(fields[2])
			if err != nil {
				return []CharacterStruct{}, err
			}
			csa = append(csa, CharacterStruct{V: variant{vtype: "del", vpos: pos, vlength: length}})
		default:
			return []CharacterStruct{}, errors.New("could not parse config line (couldn't determine what sort of variant this is): " + line)
		}
	}

	err = s.Err()
	if err != nil {
		return make([]CharacterStruct, 0), err
	}

	return csa, nil
}

// like readConfig() but with no amino acid positions
func readConfigNoAA(configFile string) ([]CharacterStruct, error) {

	csa := make([]CharacterStruct, 0)

	f, err := os.Open(configFile)
	if err != nil {
		return []CharacterStruct{}, err
	}
	defer f.Close()

	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		fields := strings.Split(line, ":")

		if len(fields) < 2 {
			return []CharacterStruct{}, errors.New("could not parse config line (too few fields): " + line)
		}

		switch fields[0] {
		case "nuc":
			pos, err := strconv.Atoi(fields[1])
			if err != nil {
				return []CharacterStruct{}, err
			}
			csa = append(csa, CharacterStruct{V: variant{vtype: "nuc", vpos: pos}})
		case "del":
			pos, err := strconv.Atoi(fields[1])
			if err != nil {
				return []CharacterStruct{}, err
			}
			length, err := strconv.Atoi(fields[2])
			if err != nil {
				return []CharacterStruct{}, err
			}
			csa = append(csa, CharacterStruct{V: variant{vtype: "nuc", vpos: pos, vlength: length}})
		default:
			return []CharacterStruct{}, errors.New("could not parse config line (couldn't determine what sort of variant this is): " + line)
		}
	}

	err = s.Err()
	if err != nil {
		return make([]CharacterStruct, 0), err
	}

	return csa, nil
}

// make a deletion of length l
func makeDeletion(l int) string {
	s := make([]byte, l)
	for i := range s {
		s[i] = '-'
	}
	return string(s)
}

func getVariantName(v variant) (string, error) {

	var s string

	switch v.vtype {
	case "aa":
		s = v.vtype + ":" + v.vgene + ":" + strconv.Itoa(v.vres)
	case "nuc":
		s = v.vtype + ":" + strconv.Itoa(v.vpos)
	case "del":
		s = v.vtype + ":" + strconv.Itoa(v.vpos) + ":" + strconv.Itoa(v.vlength)
	default:
		return "", errors.New("unknown variant type")
	}

	return s, nil
}

func makeNucLookupArray() [][]string {

	nucArr := make([][]string, 256)
	for i := range nucArr {
		nucArr[i] = make([]string, 0)
	}

	nucArr['A'] = []string{"A"}
	nucArr['C'] = []string{"C"}
	nucArr['G'] = []string{"G"}
	nucArr['T'] = []string{"T"}
	nucArr['R'] = []string{"A", "G"}
	nucArr['Y'] = []string{"C", "T"}
	nucArr['S'] = []string{"G", "C"}
	nucArr['W'] = []string{"A", "T"}
	nucArr['K'] = []string{"G", "T"}
	nucArr['M'] = []string{"A", "C"}
	nucArr['B'] = []string{"C", "G", "T"}
	nucArr['D'] = []string{"A", "G", "T"}
	nucArr['H'] = []string{"A", "C", "T"}
	nucArr['V'] = []string{"A", "C", "G"}
	// nucArr['N'] = []string{"A", "C", "G", "T"}
	nucArr['N'] = []string{}

	return nucArr
}

// a function for getting the index of a string in an array of strings
func stringIndexInArray(s string, sa []string) (int, error) {
	for i := range sa {
		if s == sa[i] {
			return i + 1, nil
		}
	}
	return -1, errors.New("state not found in array of character states")
}

// a function for checking if a string is present in an array of strings
func stringInArray(s string, sa []string) bool {
	for i := range sa {
		if s == sa[i] {
			return true
		}
	}
	return false
}

// First pass over the alignment to get the number of states for each character.
// With this information, we will be able to type it in a more memory efficient way - and hopefully straight to the tips
func countVariantsFastaInner(variantsIn []CharacterStruct, cFR chan fastaio.FastaRecord, cVOut chan []CharacterStruct, cErr chan error) {

	var ok bool
	var err error
	var newvar string
	var rawnuc byte
	var nucs []string
	codonMap := alphabet.MakeCodonDict()
	nucArr := makeNucLookupArray()

	// Make the data structures for keeping information about the characters
	variantsOut := make([]CharacterStruct, len(variantsIn))

	// What are the characters called
	for i := range variantsOut {
		variantsOut[i].V = variantsIn[i].V
		variantsOut[i].Name, err = getVariantName(variantsOut[i].V)
		if err != nil {
			cErr <- err
		}
	}

	// StateKey might be ["A", "C", "G", "T"] (if all four nucs are present at this site in the alignment)
	for i := range variantsOut {
		variantsOut[i].StateKey = make([]string, 0)
	}

	for record := range cFR {
		for i, CS := range variantsIn {
			switch CS.V.vtype {
			case "aa":
				newvar, ok = codonMap[record.Seq[CS.V.vpos-1:CS.V.vpos+2]]
				if ok {
					if stringInArray(newvar, variantsOut[i].StateKey) {
						continue
					} else {
						variantsOut[i].StateKey = append(variantsOut[i].StateKey, newvar)
					}
				} else {
					// TO DO: error check the amino acids here
					continue
				}
			case "nuc":
				rawnuc = []byte{record.Seq[CS.V.vpos-1]}[0]
				nucs = nucArr[rawnuc]
				if len(nucs) > 0 {
					for _, nuc := range nucs {
						if stringInArray(nuc, variantsOut[i].StateKey) {
							continue
						} else {
							variantsOut[i].StateKey = append(variantsOut[i].StateKey, nuc)
						}
					}
				} else {
					// TO DO: error check the nucleotides here
					continue
				}
			case "del":
				newvar = record.Seq[CS.V.vpos-1 : CS.V.vpos-1+CS.V.vlength]
				switch newvar {
				case makeDeletion(CS.V.vlength):
					if stringInArray("del", variantsOut[i].StateKey) {
						continue
					} else {
						variantsOut[i].StateKey = append(variantsOut[i].StateKey, "del")
					}
				default:
					if stringInArray("oth", variantsOut[i].StateKey) {
						continue
					} else {
						variantsOut[i].StateKey = append(variantsOut[i].StateKey, "oth")
					}
				}
			}
		}
	}

	cVOut <- variantsOut
}

func countVariantsFasta(alignmentFile string, variants []CharacterStruct) ([]CharacterStruct, error) {

	cFR := make(chan fastaio.FastaRecord)
	cErr := make(chan error)
	cFRDone := make(chan bool)

	cCS := make(chan []CharacterStruct)

	go fastaio.ReadAlignment(alignmentFile, cFR, cErr, cFRDone)

	go countVariantsFastaInner(variants, cFR, cCS, cErr)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return []CharacterStruct{}, err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	CSA := make([]CharacterStruct, 0)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return []CharacterStruct{}, err
		case temp := <-cCS:
			CSA = temp
			n--
		}
	}

	return CSA, nil
}

type StartStop struct {
	Start int
	Stop  int
}

// get the index of where individual characters are going to occur in the slice of bytes
func getIndex(csa []CharacterStruct) ([]StartStop, int) {

	// make an array of where each character starts and stops in the byte array
	idx := make([]StartStop, len(csa))

	var start int = 0
	for i, cs := range csa {
		length := (len(cs.StateKey) / 8) + 1 // how many bytes are needed to store this character
		idx[i] = StartStop{Start: start, Stop: start + length}
		start = start + length
	}

	return idx, start
}

// Spin up a few instances of this:
func typeVariants(idx []StartStop, l int, variantsIn []CharacterStruct, cFR chan fastaio.FastaRecord, cNS chan NodeStates, cErr chan error) {

	var ok bool
	var err error
	var newvar string
	var rawnuc byte
	var nucs []string
	var bitToSet int
	var start, stop int
	codonMap := alphabet.MakeCodonDict()
	nucArr := makeNucLookupArray()

	for record := range cFR {

		NS := NodeStates{}
		NS.ID = record.ID
		NS.States = make([]byte, l, l)

		for i, CS := range variantsIn {

			start = idx[i].Start
			stop = idx[i].Stop

			switch CS.V.vtype {
			case "aa":
				newvar, ok = codonMap[record.Seq[CS.V.vpos-1:CS.V.vpos+2]]
				if ok {
					bitToSet, err = stringIndexInArray(newvar, CS.StateKey)
					if err != nil {
						cErr <- err
					}
					bitsets.SetBit(NS.States[start:stop], bitToSet)
				} else {
					// If there is missing data, we don't set any bits?
					continue
				}
			case "nuc":
				rawnuc = []byte{record.Seq[CS.V.vpos-1]}[0]
				nucs = nucArr[rawnuc]
				if len(nucs) > 0 {
					for _, nuc := range nucs {
						bitToSet, err = stringIndexInArray(nuc, CS.StateKey)
						if err != nil {
							cErr <- err
						}
						bitsets.SetBit(NS.States[start:stop], bitToSet)
					}
				} else {
					// If there is missing data, we don't set any bits?
					continue
				}
			case "del":
				newvar = record.Seq[CS.V.vpos-1 : CS.V.vpos-1+CS.V.vlength]
				switch newvar {
				case makeDeletion(CS.V.vlength):
					bitToSet, err = stringIndexInArray("del", CS.StateKey)
					if err != nil {
						cErr <- err
					}
					bitsets.SetBit(NS.States[start:stop], bitToSet)
				default:
					bitToSet, err = stringIndexInArray("oth", CS.StateKey)
					if err != nil {
						cErr <- err
					}
					bitsets.SetBit(NS.States[start:stop], bitToSet)
				}
			}
		}

		cNS <- NS
	}
	return
}

func assignNodeStatesToStatesArray(t *tree.Tree, l int, cNS chan NodeStates, cErr chan error, cResults chan [][]byte) {

	states := make([][]byte, len(t.Nodes()), len(t.Nodes()))
	for i := range states {
		states[i] = make([]byte, l, l)
	}

	for ns := range cNS {
		id, err := t.TipId(ns.ID)
		if err != nil {
			cErr <- err
		}
		states[id] = ns.States
	}

	cResults <- states
}

// for every record in an alignment, type it at each variant in a variant config file, and return
// the information in map from tip name -> array of bit-encoded character states
func TypeAlignment(t *tree.Tree, alignmentFile string, configFile string, genbankFile string) ([]CharacterStruct, []StartStop, [][]byte, error) {

	var err error
	var gb genbank.Genbank
	var cdspos map[string]string
	var config []CharacterStruct

	if len(genbankFile) > 0 {
		gb, err = genbank.ReadGenBank(genbankFile)
		if err != nil {
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		}
		cdspos = getCDSPosFromAnnotation(gb)

		config, err = readConfig(configFile, cdspos)
		if err != nil {
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		}
	} else {
		config, err = readConfigNoAA(configFile)
		if err != nil {
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		}
	}

	characterStates, err := countVariantsFasta(alignmentFile, config)
	if err != nil {
		return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
	}

	// the array of start/stop positions and the total length of the byte slice of characters
	// idx is []StartStop, length is int
	idx, length := getIndex(characterStates)

	cErr := make(chan error)
	cFR := make(chan fastaio.FastaRecord)
	cFRDone := make(chan bool)
	cTypeVariantsDone := make(chan bool)
	cNS := make(chan NodeStates)
	cNSResults := make(chan [][]byte)

	var wgTypeVariants sync.WaitGroup
	wgTypeVariants.Add(runtime.NumCPU())

	go fastaio.ReadAlignment(alignmentFile, cFR, cErr, cFRDone)

	for n := 0; n < runtime.NumCPU(); n++ {
		go func() {
			typeVariants(idx, length, characterStates, cFR, cNS, cErr)
			wgTypeVariants.Done()
		}()
	}

	go func() {
		wgTypeVariants.Wait()
		cTypeVariantsDone <- true
	}()

	go assignNodeStatesToStatesArray(t, length, cNS, cErr, cNSResults)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		case <-cTypeVariantsDone:
			close(cNS)
			n--
		}
	}

	var states [][]byte
	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		case states = <-cNSResults:
			n--
		}
	}

	return characterStates, idx, states, nil
}

func getAlignmentDims(infile string) (int, int, error) {
	n := 0
	l := 0

	f, err := os.Open(infile)
	if err != nil {
		return 0, 0, err
	}
	defer f.Close()

	s := bufio.NewScanner(f)

	for s.Scan() {
		line := s.Text()

		if string(line[0]) == ">" {
			n++
		}

		if n == 1 && string(line[0]) != ">" {
			l += len(line)
		}
	}

	err = s.Err()

	if err != nil {
		return 0, 0, err
	}

	return n, l, err
}

// for every record in an alignment, type it at each nucleotide
func TypeAlignmentNuc(t *tree.Tree, alignmentFile string) ([]CharacterStruct, []StartStop, [][]byte, error) {

	var err error
	config := make([]CharacterStruct, 0)

	_, l, err := getAlignmentDims(alignmentFile)
	if err != nil {
		return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
	}

	for i := 0; i < l; i++ {
		config = append(config, CharacterStruct{V: variant{vtype: "nuc", vpos: i + 1}})
	}

	characterStates, err := countVariantsFasta(alignmentFile, config)
	if err != nil {
		return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
	}

	// the array of start/stop positions and the total length of the byte slice of characters
	// idx is []StartStop, length is int
	idx, length := getIndex(characterStates)

	cErr := make(chan error)
	cFR := make(chan fastaio.FastaRecord)
	cFRDone := make(chan bool)
	cTypeVariantsDone := make(chan bool)
	cNS := make(chan NodeStates)
	cNSResults := make(chan [][]byte)

	var wgTypeVariants sync.WaitGroup
	wgTypeVariants.Add(runtime.NumCPU())

	go fastaio.ReadAlignment(alignmentFile, cFR, cErr, cFRDone)

	for n := 0; n < runtime.NumCPU(); n++ {
		go func() {
			typeVariants(idx, length, characterStates, cFR, cNS, cErr)
			wgTypeVariants.Done()
		}()
	}

	go func() {
		wgTypeVariants.Wait()
		cTypeVariantsDone <- true
	}()

	go assignNodeStatesToStatesArray(t, length, cNS, cErr, cNSResults)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		case <-cTypeVariantsDone:
			close(cNS)
			n--
		}
	}

	var states [][]byte
	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return make([]CharacterStruct, 0), make([]StartStop, 0), make([][]byte, 0), err
		case states = <-cNSResults:
			n--
		}
	}

	return characterStates, idx, states, nil
}
