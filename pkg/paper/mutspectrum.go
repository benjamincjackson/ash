package paper

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/benjamincjackson/gotree/tree"
)

func GetPrintMutSpec(t *tree.Tree) {
	m := getMutSpectrum(t)
	for k, v := range m {
		fmt.Println(k + "\t" + strconv.Itoa(v))
	}
}

func getMutSpectrum(t *tree.Tree) map[string]int {
	m := make(map[string]int)

	for i, a := range []string{"A", "C", "G", "T"} {
		for j, b := range []string{"A", "C", "G", "T"} {
			if i == j {
				continue
			}
			s := a + "->" + b
			m[s] = 0
		}
	}

	for _, e := range t.Edges() {
		for _, c := range e.GetComments() {
			temp := strings.Split(c, "=")[1]
			temp = strings.Split(temp, ",")[0]
			if _, ok := m[temp]; ok {
				m[temp]++
			}
		}
	}

	return m
}

func GetPrintSynNonsynMutSpec(t *tree.Tree) {
	m := getSynNonsynMutSpectrum(t)
	fmt.Println("change\tsyn\tnonSyn")

	for k, v := range m {
		fmt.Println(k + "\t" + strconv.Itoa(v["syn"]) + "\t" + strconv.Itoa(v["nonSyn"]))
	}
}

func getSynNonsynMutSpectrum(t *tree.Tree) map[string]map[string]int {
	m := make(map[string]map[string]int)

	for i, a := range []string{"A", "C", "G", "T"} {
		for j, b := range []string{"A", "C", "G", "T"} {
			if i == j {
				continue
			}
			s := a + "->" + b
			m[s] = make(map[string]int)
			m[s]["syn"] = 0
			m[s]["nonSyn"] = 0
		}
	}

	for _, e := range t.Edges() {
		for _, c := range e.GetComments() {
			SorN := strings.Split(c, "=")[0]
			trans := strings.Split(c, "=")[1]
			if _, ok := m[trans][SorN]; ok {
				m[trans][SorN]++
			}
		}
	}

	return m
}
