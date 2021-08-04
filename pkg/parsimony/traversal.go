package parsimony

import (
	"sort"
	"strconv"
	"strings"

	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/ash/pkg/characterio"
	"github.com/benjamincjackson/ash/pkg/tree"

	"github.com/cov-ert/gofasta/pkg/alphabet"
)

// LabelChanges traverses over the tree and labels inferred changes onto the branches
func LabelChanges(t *tree.Tree, characters []characterio.CharacterStruct, states [][]byte, idx []characterio.StartStop, transitions [][]characterio.Transition) {
	labelChanges(t.Root(), nil, characters, states, idx, transitions)
}

// recursive function to label branches with state changes
func labelChanges(cur, prev *tree.Node, characters []characterio.CharacterStruct, states [][]byte, idx []characterio.StartStop, transitions [][]characterio.Transition) {

	// Base state: if cur is a tip, there's nothing to do here
	if len(cur.Neigh()) == 1 {
		return
	}

	for i, n := range cur.Neigh() {
		if n != prev {
			// we see if we need to label this edge with a change
			labelEdge(cur, n, cur.Edges()[i], characters, states, idx, transitions)
			// and then we carry on recurring down the tree
			labelChanges(n, cur, characters, states, idx, transitions)
		}
	}
}

// type CharacterStruct struct {
// 	Name     string   // name of the variant
// 	V        variant  // detailed info about this variant
// 	StateKey []string // the slice might be ["A", "C", "G", "T"] (if all four nucs are present at this site in the alignment)
// }

// label the branch between two nodes with any inferred changes, and record the transition
func labelEdge(upnode, downnode *tree.Node, edge *tree.Edge, characters []characterio.CharacterStruct, states [][]byte, idx []characterio.StartStop, transitions [][]characterio.Transition) {

	var start, stop int

	// for every character
	up_id := upnode.Id()
	down_id := downnode.Id()
	for i := range idx {
		start = idx[i].Start
		stop = idx[i].Stop

		// if the states are different
		if bitsets.Different(states[up_id][start:stop], states[down_id][start:stop]) {
			// we don't care about transitios to missing data:
			if !bitsets.IsAnyBitSet(states[down_id][start:stop]) {
				continue
			}

			// we shouldn't care about transitions to ambiguous tips?
			if downnode.Tip() && bitsets.IsSubset(states[up_id][start:stop], states[down_id][start:stop]) {
				continue
			}

			// then we annotate the edge

			// the set bits for this character
			upstatebits := bitsets.GetSetBits(states[up_id][start:stop])
			downstatebits := bitsets.GetSetBits(states[down_id][start:stop])

			// then what these mean as states
			upstate := make([]string, 0)
			for _, b := range upstatebits {
				upstate = append(upstate, characters[i].StateKey[b-1])
			}

			downstate := make([]string, 0)
			for _, b := range downstatebits {
				downstate = append(downstate, characters[i].StateKey[b-1])
			}

			// then we can build the label and add it to the edge
			anc := strings.Join(upstate, "|")
			der := strings.Join(downstate, "|")
			trans := anc + "->" + der
			number := getTransitionNumber(transitions[i], trans)
			label := characters[i].Name + "=" + anc + "->" + der + "," + anc + "->" + der + "#" + strconv.Itoa(number)
			edge.AddComment(label)

			// and we can add the location of this transition to the array of transitions
			transition := characterio.Transition{Upnode: upnode, Downnode: downnode, Edge: edge, Upstate: anc, Downstate: der, Number: number, Transition: trans, Label: label}
			transitions[i] = append(transitions[i], transition)
		}
	}
}

func LabelChangesAnno(t *tree.Tree, regions []annotation.Region, characters []characterio.CharacterStruct, states [][]byte) {
	labelChangesAnno(t.Root(), nil, regions, characters, states)
}

// recursive function to label branches with state changes
func labelChangesAnno(cur, prev *tree.Node, regions []annotation.Region, characters []characterio.CharacterStruct, states [][]byte) {

	// Base state: if cur is a tip, there's nothing to do here
	if len(cur.Neigh()) == 1 {
		return
	}

	for i, n := range cur.Neigh() {
		if n != prev {
			// we see if we need to label this edge with a change
			labelEdgeAnno(cur, n, cur.Edges()[i], regions, characters, states)
			// and then we carry on recurring down the tree
			labelChangesAnno(n, cur, regions, characters, states)
		}
	}
}

// label an edge with annotated changes - amino acid changing versus neutral nucleotide change
func labelEdgeAnno(upnode, downnode *tree.Node, edge *tree.Edge, regions []annotation.Region, characters []characterio.CharacterStruct, states [][]byte) {
	up_id := upnode.Id()
	down_id := downnode.Id()

	IUPACMap := annotation.GetIUPACMap()
	codonDict := alphabet.MakeCodonDict()

	// the characters are all nucleotides so we don't need the idx of states
	// instead we use the regions slice to annotate things
	for _, region := range regions {
		switch region.Whichtype {
		case "int":
			for pos := region.Start - 1; pos < region.Stop; pos++ {
				if bitsets.Different(states[up_id][pos:pos+1], states[down_id][pos:pos+1]) {

					// we don't care about transitions to missing data:
					if !bitsets.IsAnyBitSet(states[down_id][pos : pos+1]) {
						continue
					}

					// we shouldn't care about transitions to ambiguous tips?
					if downnode.Tip() && bitsets.IsSubset(states[up_id][pos:pos+1], states[down_id][pos:pos+1]) {
						continue
					}

					// then we annotate the edge

					// the set bits for this character
					upstatebits := bitsets.GetSetBits(states[up_id][pos : pos+1])
					downstatebits := bitsets.GetSetBits(states[down_id][pos : pos+1])

					// then what these mean as states
					upstate := make([]string, 0)
					for _, b := range upstatebits {
						upstate = append(upstate, characters[pos].StateKey[b-1])
					}

					downstate := make([]string, 0)
					for _, b := range downstatebits {
						downstate = append(downstate, characters[pos].StateKey[b-1])
					}

					// then we can build the label and add it to the edge
					anc := strings.Join(upstate, "|")
					der := strings.Join(downstate, "|")
					// trans := anc + "->" + der
					// number := getTransitionNumber(transitions[i], trans)
					label := "nuc=" + anc + strconv.Itoa(pos+1) + der
					edge.AddComment(label)

					// Skipping this for now (means we can't summarise things). To do: re-implement this
					// // and we can add the location of this transition to the array of transitions
					// transition := characterio.Transition{Upnode: upnode, Downnode: downnode, Edge: edge, Upstate: anc, Downstate: der, Number: number, Transition: trans, Label: label}
					// transitions[i] = append(transitions[i], transition)
				}
			}
		case "CDS":
			// 1-based position of amino acid in this CDS:
			AACounter := 1

			// for every codon in this CDS:
			for _, codonstart := range region.Codonstarts {
				// if the bitsets for each codon's (three nucleotides') states are different:
				if bitsets.Different(states[up_id][codonstart-1:codonstart+2], states[down_id][codonstart-1:codonstart+2]) {
					// then we need to get the nucleotides and attempt to translate them
					nuclabels := make([]string, 0)

					upcodon := ""
					downcodon := ""

					// nucleotide position:
					pos := codonstart - 1

					for i := 0; i < 3; i++ {

						// the set bits for this character
						upstatebits := bitsets.GetSetBits(states[up_id][pos : pos+1])
						downstatebits := bitsets.GetSetBits(states[down_id][pos : pos+1])

						// then what these mean as states
						upstate := make([]string, 0)
						for _, b := range upstatebits {
							upstate = append(upstate, characters[pos].StateKey[b-1])
						}
						// sort them, for translating to the correct ambiguity code
						sort.Strings(upstate)

						downstate := make([]string, 0)
						for _, b := range downstatebits {
							downstate = append(downstate, characters[pos].StateKey[b-1])
						}
						// sort them, for translating to the correct ambiguity code
						sort.Strings(downstate)

						// build the ancestral codon
						switch {
						case len(upstate) == 1:
							upcodon = upcodon + upstate[0]
						case len(upstate) > 1:
							upcodon = upcodon + IUPACMap[strings.Join(upstate, "")]
						default:
							upcodon = upcodon + "N"
						}

						// build the derived codon
						switch {
						case len(downstate) == 1:
							downcodon = downcodon + downstate[0]
						case len(downstate) > 1:
							downcodon = downcodon + IUPACMap[strings.Join(downstate, "")]
						default:
							downcodon = downcodon + "N"
						}

						// we also build the nuc label and add it to the temporary slice of nuc changes,
						// if the nucs are different
						anc := strings.Join(upstate, "|")
						der := strings.Join(downstate, "|")
						if anc != der {
							// we don't care about transitions to missing data:
							if !bitsets.IsAnyBitSet(states[down_id][pos : pos+1]) {
								continue
							}

							// we shouldn't care about transitions to ambiguous tips?
							if downnode.Tip() && bitsets.IsSubset(states[up_id][pos:pos+1], states[down_id][pos:pos+1]) {
								continue
							}

							label := "nuc=" + anc + strconv.Itoa(pos+1) + der
							nuclabels = append(nuclabels, label)
						}
						// increment the nucleotide position:
						pos++
					}

					var upAA string
					var downAA string

					if _, ok := codonDict[upcodon]; ok {
						upAA = codonDict[upcodon]
					} else {
						upAA = "X"
					}

					if _, ok := codonDict[downcodon]; ok {
						downAA = codonDict[downcodon]
					} else {
						downAA = "X"
					}

					// we switch on whether we know both the ancestral and derived amino acids unambiguously
					ResolvedAAs := upAA != "X" && downAA != "X"
					switch ResolvedAAs {

					// if we can, we want to further switch on whether they're different
					case true:
						DifferentAAs := upAA != downAA
						switch DifferentAAs {

						// They are different. We label the edge with the AA change, we don't label any SNPs
						case true:
							label := "AA=" + region.Name + ":" + upAA + strconv.Itoa(AACounter) + downAA
							edge.AddComment(label)

						// They are the same. We label the edge with any nucleotide changes that there are
						case false:
							for _, label := range nuclabels {
								edge.AddComment(label)
							}
						}

					// if we can't, then we want to record the SNPs (SOMETHING FOR LATER- do we want to call them synonymous?)
					case false:
						for _, label := range nuclabels {
							edge.AddComment(label)
						}
					}
				}

				AACounter++
			}
		}
	}
}

func getTransitionNumber(ta []characterio.Transition, label string) int {
	count := 0
	for i := range ta {
		if label == ta[i].Transition {
			count++
		}
	}
	return count
}

func LabelNodes(t *tree.Tree, characters []characterio.CharacterStruct, states [][]byte, idx []characterio.StartStop) {
	for _, n := range t.Nodes() {
		id := n.Id()
		n.AddComment("nodenumber=" + strconv.Itoa(id))
		for i := range idx {
			start := idx[i].Start
			stop := idx[i].Stop

			statebits := bitsets.GetSetBits(states[id][start:stop])
			state := make([]string, 0)
			for _, b := range statebits {
				state = append(state, characters[i].StateKey[b-1])
			}
			statestring := strings.Join(state, "|")
			n.AddComment(characters[i].Name + "node=" + statestring)
		}
	}
}
