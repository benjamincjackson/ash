package annotation

import (
	"errors"
	"strconv"
	"strings"

	"github.com/cov-ert/gofasta/pkg/genbank"
)

// Region is a struct containing a part of the genome which might
// be a CDS or intergenic, for example
type Region struct {
	Whichtype   string // int(ergenic) or CDS
	Name        string // name of CDS, if it is one
	Start       int    // 1-based first position of region, inclusive
	Stop        int    // 1-based last position of region, inclusive
	Codonstarts []int  // a slice of the 1-based start positions of all its codons, if this region is a CDS
}

// from a string representation of all the nucleotides a state is (sorted in increasing order),
// return the appropriate IUPAC code
// TO DO - write a test for this
func GetIUPACMap() map[string]string {
	m := make(map[string]string)
	m["A"] = "A"
	m["T"] = "T"
	m["C"] = "C"
	m["G"] = "G"
	m["AC"] = "M"
	m["AG"] = "R"
	m["AT"] = "W"
	m["CG"] = "S"
	m["CT"] = "Y"
	m["GT"] = "K"
	m["ACG"] = "V"
	m["ACT"] = "H"
	m["AGT"] = "D"
	m["CGT"] = "B"
	m["ACGT"] = "N"
	return m
}

// get the positions of the CDS and the not-CDS from the genbank.
// return a slice of Region structs
func GetRegions(genbankFileIn string, nuc bool) ([]Region, error) {
	gb, err := genbank.ReadGenBank(genbankFileIn)
	if err != nil {
		return make([]Region, 0), err
	}

	if nuc {
		REGION := Region{Whichtype: "int", Start: 1, Stop: len(gb.ORIGIN)}
		return []Region{REGION}, nil
	}

	CDSFEATS := make([]genbank.GenbankFeature, 0)
	for _, F := range gb.FEATURES {
		if F.Feature == "CDS" {
			CDSFEATS = append(CDSFEATS, F)
		}
	}

	cdsregions := make([]Region, 0)

	// we get all the CDSes
	for _, feat := range CDSFEATS {
		REGION := Region{Whichtype: "CDS", Name: feat.Info["gene"], Codonstarts: make([]int, 0)}

		// these are genbank positions, so they are 1-based, inclusive
		positions, err := parsePositions(feat.Pos)
		if err != nil {
			return make([]Region, 0), err
		}
		REGION.Start = positions[0]
		REGION.Stop = positions[len(positions)-1]

		// how many codons we have already defined if this CDS is two ranges Join()ed together:
		previouscodons := 0
		for i := 0; i < len(positions); i = i + 2 {
			start := positions[i]
			stop := positions[i+1]

			// start-1 to get the start position as 0-based
			if (stop-(start-1))%3 != 0 {
				return make([]Region, 0), errors.New("CDS position range is not a multiple of 3")
			}
			length := (stop - (start - 1))

			for j := 0; j < length; j = j + 3 {
				// we add another integer codon start to the slice
				REGION.Codonstarts = append(REGION.Codonstarts, start+j)
			}

			previouscodons = previouscodons + (length / 3)
		}
		cdsregions = append(cdsregions, REGION)
	}

	// then we add the intergenic regions (between the CDSes)
	regions := make([]Region, 0)
	newstart := 1
	for i, cdsregion := range cdsregions {
		start := newstart
		stop := cdsregion.Start - 1

		// hopefully this deals with any cases where there isn't an intergenic region:
		if !((stop - start) > 0) {
			continue
		}
		REGION := Region{Whichtype: "int", Start: start, Stop: stop}
		regions = append(regions, REGION)
		regions = append(regions, cdsregion)
		newstart = cdsregion.Stop + 1
		if i == len(cdsregions) {
			start := newstart
			stop := len(gb.ORIGIN)
			if !((stop - start) > 0) {
				continue
			}
			REGION := Region{Whichtype: "int", Start: start, Stop: stop}
			regions = append(regions, REGION)
		}
	}

	return regions, nil
}

func parsePositions(position string) ([]int, error) {
	var A []int
	if position[0:4] == "join" {
		A = make([]int, 0)
		position = strings.TrimLeft(position, "join(")
		position = strings.TrimRight(position, ")")
		ranges := strings.Split(position, ",")
		for _, x := range ranges {
			y := strings.Split(x, "..")
			for _, z := range y {
				temp, err := strconv.Atoi(z)
				if err != nil {
					return []int{}, err
				}
				A = append(A, temp)
			}
		}
	} else {
		A = make([]int, 0)
		y := strings.Split(position, "..")
		for _, z := range y {
			temp, err := strconv.Atoi(z)
			if err != nil {
				return []int{}, err
			}
			A = append(A, temp)
		}
	}

	return A, nil
}
