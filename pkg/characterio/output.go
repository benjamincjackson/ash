package characterio

import (
	"strconv"
	"strings"

	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/gotree/tree"
)

// Transition is a struct for one inferred change on an edge/branch
type Transition struct {
	Upnode     *tree.Node
	Downnode   *tree.Node
	Edge       *tree.Edge
	Upstate    string
	Downstate  string
	Number     int
	Transition string
	Label      string
}

// print some info about each transition's children
func SummarizeTransitions(threshold int, transitions []Transition, states [][]byte, idx StartStop, c CharacterStruct) []string {

	// the first lines are the name of the character
	output := make([]string, 0)
	output = append(output, c.Name)
	bs := make([]byte, len(c.Name))
	for i, _ := range bs {
		bs[i] = '-'
	}
	output = append(output, string(bs))

	// need to aggregate the changes at this site by transition type
	aggregated := make([][]Transition, 0)
	counter := 0
	previousTrans := ""
	for i, t := range transitions {
		if i == 0 {
			aggregated = append(aggregated, make([]Transition, 0))
			aggregated[counter] = append(aggregated[counter], t)
			previousTrans = t.Transition
		} else if t.Transition == previousTrans {
			aggregated[counter] = append(aggregated[counter], t)
		} else {
			aggregated = append(aggregated, make([]Transition, 0))
			counter++
			aggregated[counter] = append(aggregated[counter], t)
			previousTrans = t.Transition
		}
	}

	// how many difference sorts of changes are there for this state
	output = append(output, strconv.Itoa(len(aggregated))+" type(s) of transition:")

	var m map[string]int

	for i := range aggregated {
		// what is this change
		output = append(output, "\t"+aggregated[i][0].Transition)
		// how many times does it occur:
		output = append(output, "\t"+strconv.Itoa(len(aggregated[i]))+" occurence(s) on the tree:")
		// then for each one:
		for j := range aggregated[i] {
			m = summarizeTransitions(aggregated[i][j], states, idx, c)
			total := 0
			for _, v := range m {
				total += v
			}
			if !(total > threshold) {
				continue
			}
			inner := "\t\t" + strings.Split(aggregated[i][j].Label, ",")[1] + " has " + strconv.Itoa(total) + " child tip(s), with character counts at " + c.Name + " of: ("
			temp := make([]string, 0)
			for k, v := range m {
				temp = append(temp, k+": "+strconv.Itoa(v))
			}
			inner = inner + strings.Join(temp, ", ") + ")"
			output = append(output, inner)
		}
	}

	return output
}

func summarizeTransitions(trans Transition, states [][]byte, idx StartStop, cS CharacterStruct) map[string]int {
	m := make(map[string]int)
	updateMap(m, trans.Downnode, trans.Upnode, states, idx, cS)
	return m
}

func updateMap(m map[string]int, cur, prev *tree.Node, states [][]byte, idx StartStop, cS CharacterStruct) {

	// Base State if downnode is a tip, then we update the map and return
	if cur.Tip() {
		node_id := cur.Id()
		setbits := bitsets.GetSetBits(states[node_id][idx.Start:idx.Stop])
		characters := make([]string, 0)
		for _, b := range setbits {
			characters = append(characters, cS.StateKey[b-1])
		}
		var state string
		if len(characters) >= 1 {
			state = strings.Join(characters, "|")
		} else {
			state = "missing"
		}
		if _, ok := m[state]; ok {
			m[state]++
		} else {
			m[state] = 1
		}
		return
	}
	for _, n := range cur.Neigh() {
		if n == prev {
			continue
		}
		updateMap(m, n, cur, states, idx, cS)
	}
	return
}

// // print some info about each transition's children
// func SummarizeTransitions(cS []CharacterStruct, m map[string]map[string][]tree.TipBag) {
// 	fmt.Println()
// 	for _, c := range cS {
// 		fmt.Println(c.Name)
// 		bs := make([]byte, len(c.Name))
// 		for i, _ := range bs {
// 			bs[i] = '-'
// 		}
// 		fmt.Println(string(bs))
// 		fmt.Println(strconv.Itoa(len(m[c.Name])) + " type(s) of transition:")
// 		for transition := range m[c.Name] {
// 			fmt.Println("\t" + transition)
// 			n := len(m[c.Name][transition])
// 			fmt.Println("\t" + strconv.Itoa(n) + " occurence(s) on the tree:")
// 			for i, tb := range m[c.Name][transition] {
// 				tipNames := tb.TipNames()
// 				genoCount := aggregateCharacters(c, tipNames)

// 				nTips := tb.Size()
// 				fmt.Print("\t" + "\t" + transition + "#" + strconv.Itoa(i) + " has " + strconv.Itoa(nTips) + " child tip(s), ")
// 				fmt.Print("with character counts at " + c.Name + " of: (")
// 				temp := make([]string, 0)
// 				for k, v := range genoCount {
// 					if k == "" {
// 						k = "missing"
// 					}
// 					temp = append(temp, k+": "+strconv.Itoa(v))
// 				}
// 				fmt.Print(strings.Join(temp, ", "))
// 				fmt.Println(")")
// 			}
// 		}
// 	}
// 	fmt.Println()
// }

// func aggregateCharacters(c characterio.CharacterStruct, n []string) map[string]int {

// 	m := make(map[string]int)
// 	var allele string

// 	for _, name := range n {
// 		allele = c.Tipmap[name]
// 		if _, ok := m[allele]; ok {
// 			m[allele]++
// 		} else {
// 			m[allele] = 1
// 		}
// 	}

// 	return m
// }

// func WriteChildren(filename string, cS []characterio.CharacterStruct, m map[string]map[string][]tree.TipBag) error {

// 	f, err := os.Create(filename)
// 	if err != nil {
// 		return err
// 	}

// 	_, err = f.WriteString("character,transition,child,state\n")
// 	if err != nil {
// 		return err
// 	}

// 	var genotype string

// 	for _, c := range cS {
// 		char := c.Name
// 		for transition := range m[c.Name] {
// 			for i, tb := range m[c.Name][transition] {
// 				tipNames := tb.TipNames()
// 				for _, tip := range tipNames {
// 					genotype = c.Tipmap[tip]
// 					if genotype == "" {
// 						genotype = "missing"
// 					}
// 					_, err = f.WriteString(char + "," + transition + "#" + strconv.Itoa(i) + "," + tip + "," + genotype + "\n")
// 					if err != nil {
// 						return err
// 					}
// 				}
// 			}
// 		}
// 	}

// 	return nil
// }
