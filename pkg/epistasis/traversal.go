package epistasis

import (
	"sort"
	"strconv"
	"strings"

	"github.com/cov-ert/gofasta/pkg/alphabet"

	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/ash/pkg/characterio"
	"github.com/benjamincjackson/ash/pkg/tree"
)

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

					// then we can build the label and add it to the edge,
					// and increment the synonymous branch length
					anc := strings.Join(upstate, "|")
					der := strings.Join(downstate, "|")
					// trans := anc + "->" + der
					// number := getTransitionNumber(transitions[i], trans)
					label := "syn=" + anc + strconv.Itoa(pos+1) + der
					edge.AddComment(label)
					edge.SynLen++

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

							label := "syn=" + anc + strconv.Itoa(pos+1) + der
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
							label := "AA=" + region.Name + ":" + strconv.Itoa(AACounter) + ":" + upAA + downAA
							edge.AddComment(label)

						// They are the same. We label the edge with any nucleotide changes that there are,
						// and increment the synonymous branch length
						case false:
							for _, label := range nuclabels {
								edge.AddComment(label)
								edge.SynLen++
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
