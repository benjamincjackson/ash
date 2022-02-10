package main

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/spf13/cobra"

	"github.com/benjamincjackson/ash/pkg/ancestry"
	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/ash/pkg/characterio"
	"github.com/benjamincjackson/ash/pkg/epistasis"
	"github.com/benjamincjackson/ash/pkg/paper"
	"github.com/benjamincjackson/ash/pkg/parsimony"

	"github.com/benjamincjackson/gotree/newick"
	"github.com/benjamincjackson/gotree/tree"
)

func morethanone(ba ...bool) bool {
	counter := 0
	for _, b := range ba {
		if b {
			counter++
		}
		if counter > 1 {
			return true
		}
	}
	return false
}

// to do possibly - sanity check arguments if --civet is given
func checkArgs(treeFile string, alignmentFile string, variantsConfig string, genbankFile string, tipFile string,
	algorithmUp string, algorithmDown string, treeOut string,
	civet bool, nuc bool, p bool, epi bool, common_anc bool) (int, int, string, string, error) {

	algoUp := -1
	switch algorithmUp {
	case "hard":
		algoUp = 0
	case "soft":
		algoUp = 1
	default:
		return -1, -1, "", "", errors.New("unknown up-pass algorithm: choose one of soft or hard")
	}

	algoDown := -1
	switch algorithmDown {
	case "acctrans":
		algoDown = 0
	case "deltrans":
		algoDown = 1
	case "downpass":
		algoDown = 2
	default:
		return -1, -1, "", "", errors.New("unknown down-pass algorithm: choose one of acctrans, deltrans or downpass")
	}

	preset := "none"

	if morethanone(civet, nuc, p, epi, common_anc) {
		return algoUp, algoDown, "", "", errors.New("use one preset, not a combination")
	} else if civet {
		preset = "civet"
	} else if nuc {
		preset = "nuc"
	} else if p {
		preset = "paper"
	} else if common_anc {
		preset = "common_anc"
	} else if epi {
		preset = "epistasis"
	}

	// this is all wrong now that we have presets (civet and nuc)
	if len(alignmentFile) > 0 && len(tipFile) > 0 {
		return algoUp, algoDown, "", "", errors.New("either use a --tipfile OR an --alignment and a --variants-config file, not a mixture")
	}
	if len(variantsConfig) > 0 && len(tipFile) > 0 {
		return algoUp, algoDown, "", "", errors.New("either use a --tipfile OR an --alignment and a --variants-config file, not a mixture")
	}
	// if len(variantsConfig) > 0 && len(alignmentFile) == 0 || len(variantsConfig) == 0 && len(alignmentFile) > 0 {
	// 	return algoUp, algoDown, "", "", errors.New("if you provide an --alignment file you must provide a --variants-config file, and vice versa")
	// }

	if len(variantsConfig) > 0 {
		aa := false
		f, err := os.Open(variantsConfig)
		if err != nil {
			return algoUp, algoDown, "", "", err
		}
		defer f.Close()
		s := bufio.NewScanner(f)
		for s.Scan() {
			line := s.Text()
			fields := strings.Split(line, ":")
			if fields[0] == "aa" {
				aa = true
				break
			}
		}
		err = s.Err()
		if err != nil {
			return algoUp, algoDown, "", "", err
		}
		if aa && len(genbankFile) == 0 {
			return algoUp, algoDown, "", "", errors.New("you must provide a --genbank file if there are amino acids in your --variants-config file")
		}
	}

	var s string
	if len(alignmentFile) > 0 {
		s = "alignment"
	} else {
		s = "csv"
	}

	return algoUp, algoDown, s, preset, nil
}

func readTree(treeFile string) (*tree.Tree, error) {
	var t *tree.Tree
	var f *os.File
	var err error

	f, err = os.Open(treeFile)
	defer f.Close()
	if err != nil {
		return new(tree.Tree), err
	}

	t, err = newick.NewParser(f).Parse()
	if err != nil {
		return new(tree.Tree), err
	}

	// we find the max depth for each node, 'cos we want to sort on it
	t.MaxDepthRooted(t.Root(), nil)

	// then we sort by it
	t.SortNeighborsByDepth(t.Root(), nil)

	// Must update(/initiate?) the tip index so we can map the character states for the tips straight to the tree
	t.UpdateTipIndex()

	return t, nil
}

// func getRealSizeOf(v interface{}) (int, error) {
// 	b := new(bytes.Buffer)
// 	if err := gob.NewEncoder(b).Encode(v); err != nil {
// 		return 0, err
// 	}
// 	return b.Len(), nil
// }

func ash(treeIn string, alignmentFile string, variantsConfig string, genbankFile string, tipFile string,
	algorithmUp string, algorithmDown string, annotateNodes bool, annotateTips bool, threshold int,
	treeOut string, childrenOut string,
	summarize bool, civet bool, nuc bool, p bool, epi bool, common_anc bool, outgroup string, rescale bool,
	threads int) error {

	// algoUp, algoDown, input, err := checkArgs(treeIn, alignmentFile, variantsConfig, genbankFile, tipFile, algorithmUp, algorithmDown, treeOut)
	algoUp, algoDown, input, preset, err := checkArgs(treeIn, alignmentFile, variantsConfig, genbankFile, tipFile, algorithmUp, algorithmDown, treeOut, civet, nuc, p, epi, common_anc)
	if err != nil {
		return err
	}

	/*
		read in the tree
	*/
	t, err := readTree(treeIn)
	if err != nil {
		return err
	}

	if !t.Rooted() {
		return errors.New("the input tree is not rooted")
	}

	/*
		read in the tip states to the array of all nodes' states, and keep the characters around for looking up later
	*/
	var characterStates []characterio.CharacterStruct
	var idx []characterio.StartStop
	var states [][]byte

	// to do - incorporate the civet/nuc presets into the logic here?
	switch input {
	case "alignment":
		switch preset {
		case "none":
			characterStates, idx, states, err = characterio.TypeAlignment(t, alignmentFile, variantsConfig, genbankFile)
			if err != nil {
				return err
			}
		default:
			characterStates, idx, states, err = characterio.TypeAlignmentNuc(t, alignmentFile)
			if err != nil {
				return err
			}
		}
	case "csv":
		// TO DO- in tipfile columns that contain nucleotide data, IUPAC codes are treated as non-overlapping states, e.g. W != (A & T), which is different from the same data in an alignment input
		characterStates, idx, states, err = characterio.TypeTipfile(t, tipFile)
		if err != nil {
			return err
		}
	default:
		return errors.New("couldn't choose where the states are coming from")
	}

	// TO DO- check all the tree's tips are in the character state input (we do the converse already when we read the character states in)

	// TO DO- maybe just use the hard polytomies interpretation?
	switch algoUp {
	case 0: // hard polytomies
		parsimony.UpPass(t, 0, states, idx)
	case 1: // soft polytomies (resolve them [separately for each character!])
		parsimony.UpPass(t, 1, states, idx)
	}

	// for _, n := range t.Nodes() {
	// 	if n.Tip() {
	// 		continue
	// 	}
	// 	fmt.Println(states[n.Id()])
	// }

	// fmt.Println(algoDown)
	switch algoDown {
	case 0: // Acctrans
		parsimony.Acctrans(t, states, idx)
	case 1: // Deltrans
		parsimony.DownPass(t, states, idx)
		parsimony.Deltrans(t, states, idx)
	case 2: // Downpass only
		parsimony.DownPass(t, states, idx)
	}

	// for _, n := range t.Nodes() {
	// 	if n.Tip() {
	// 		continue
	// 	}
	// 	fmt.Println(states[n.Id()])
	// }

	switch preset {
	case "civet":
		// genbank annotation parsing:
		features, err := annotation.GetRegions(genbankFile, nuc)
		if err != nil {
			return err
		}

		parsimony.LabelChangesAnno(t, features, characterStates, states)

		if len(treeOut) > 0 {
			fout, err := os.Create(treeOut)
			if err != nil {
				return err
			}
			defer fout.Close()

			fout.WriteString(t.NexusOptionalComments(annotateNodes, annotateTips))
		}

	case "nuc":
		// genbank annotation parsing:
		features, err := annotation.GetRegions(genbankFile, nuc)
		if err != nil {
			return err
		}

		parsimony.LabelChangesAnno(t, features, characterStates, states)

		f, err := os.Create("branchlengths100k.tsv")
		if err != nil {
			return err
		}
		defer f.Close()
		var tip string
		f.WriteString("branch\tlength\tterminal\tnummuts\ttransitions\n")
		for i, e := range t.Edges() {
			if e.Right().Tip() {
				tip = "true"
			} else {
				tip = "false"
			}
			f.WriteString(strconv.Itoa(i) + "\t" + strconv.FormatFloat(e.Length(), 'f', 8, 64) + "\t" + tip + "\t" + strconv.Itoa(len(e.GetComments())) + "\t" + strings.Join(e.GetComments(), " ") + "\n")
		}

		// rescale the tree for JT
		if rescale {
			for _, e := range t.Edges() {
				e.SetLength(float64(len(e.GetComments())))
			}
		}

		// write the treefile...
		if len(treeOut) > 0 {
			fout, err := os.Create(treeOut)
			if err != nil {
				return err
			}
			defer fout.Close()

			// fout.WriteString(t.NewickOptionalComments(annotateNodes, annotateTips) + "\n")
			fout.WriteString(t.NexusOptionalComments(annotateNodes, annotateTips))
		}

	case "common_anc":
		// get the sequence at the node immediately ancestral to a set of samples
		// first step is as for "nuc"
		features, err := annotation.GetRegions(genbankFile, nuc)
		if err != nil {
			return err
		}

		parsimony.LabelChangesAnno(t, features, characterStates, states)

		// then get the ancestral node and print its sequence
		commonAncNodeID, err := ancestry.MRCA(t, outgroup)
		// _, err = ancestry.MRCA(t, outgroup)
		if err != nil {
			return err
		}
		// for _, tip := range t.Tips() {
		// 	fmt.Print(tip.Name() + " ")
		// 	fmt.Println(states[tip.Id()])
		// }
		// for i := range characterStates {
		// 	fmt.Println(characterStates[i])
		// }
		// for i := range characterStates {
		// 	fmt.Println(bitsets.GetSetBits(states[commonAncNodeID][i : i+1]))
		// }
		fmt.Println(">root")
		IUPACMap := annotation.GetIUPACMap()
		for i := range states[commonAncNodeID] {
			setBits := bitsets.GetSetBits(states[commonAncNodeID][i : i+1])
			nucstates := make([]string, 0)
			for _, b := range setBits {
				// fmt.Println(characterStates[i])
				nucstates = append(nucstates, characterStates[i].StateKey[b-1])
			}
			sort.Strings(nucstates)
			if nuc, ok := IUPACMap[strings.Join(nucstates, "")]; ok {
				fmt.Print(nuc)
			} else {
				fmt.Print("N")
			}
		}
		fmt.Println()

	case "paper":
		features, err := annotation.GetRegions(genbankFile, nuc)
		if err != nil {
			return err
		}
		transitions := make([][]characterio.Transition, len(characterStates), len(characterStates))
		for i := range transitions {
			transitions[i] = make([]characterio.Transition, 0)
		}
		paper.LabelChangesSynNonsyn(t, features, characterStates, states)
		paper.GetPrintSynNonsynMutSpec(t)

	case "epistasis":
		features, err := annotation.GetRegions(genbankFile, nuc)
		if err != nil {
			return err
		}
		// transitions := make([][]characterio.Transition, len(characterStates), len(characterStates))
		// for i := range transitions {
		// 	transitions[i] = make([]characterio.Transition, 0)
		// }
		// label the changes
		epistasis.LabelChangesAnno(t, features, characterStates, states)
		// calculate the epistasis statistic
		epistasis.Epistasis(t, features, threads)
	default:
		// we can keep all the transitions (structs containing pointers to the relevant nodes and edge) in an array
		// that matches the dimensions of the characters
		transitions := make([][]characterio.Transition, len(characterStates), len(characterStates))
		for i := range transitions {
			transitions[i] = make([]characterio.Transition, 0)
		}

		// to do- switch on whether we need to label the transitions or not
		parsimony.LabelChanges(t, characterStates, states, idx, transitions)

		// if we do, we should sort them
		for k := range transitions {
			sort.SliceStable(transitions[k], func(i, j int) bool {
				return transitions[k][i].Label < transitions[k][j].Label
			})
		}

		if annotateNodes {
			parsimony.LabelNodes(t, characterStates, states, idx)
		}

		if summarize {
			// TO DO: swap between stdout + a hard file
			fTrans := os.Stdout
			var lines []string
			for i := range transitions {
				lines = characterio.SummarizeTransitions(threshold, transitions[i], states, idx[i], characterStates[i])
				for _, l := range lines {
					fTrans.WriteString(l + "\n")
				}
			}
		}
	}

	// // write the treefile...
	// if len(treeOut) > 0 {
	// 	fout, err := os.Create(treeOut)
	// 	if err != nil {
	// 		return err
	// 	}
	// 	defer fout.Close()

	// 	// fout.WriteString(t.NewickOptionalComments(annotateNodes, annotateTips) + "\n")
	// 	fout.WriteString(t.NexusOptionalComments(annotateNodes, annotateTips))
	// }

	// f, err := os.Create("branchlengths.tsv")
	// if err != nil {
	// 	return err
	// }
	// defer f.Close()
	// var tip string
	// f.WriteString("branch\tlength\tterminal\tnummuts\ttransitions\n")
	// for i, e := range t.Edges() {
	// 	if e.Right().Tip() {
	// 		tip = "true"
	// 	} else {
	// 		tip = "false"
	// 	}
	// 	f.WriteString(strconv.Itoa(i) + "\t" + strconv.FormatFloat(e.Length()*29903, 'f', 4, 64) + "\t" + tip + "\t" + strconv.Itoa(len(e.GetComments())) + "\t" + strings.Join(e.GetComments(), " ") + "\n")
	// }

	// t.SortNeighborsByTips(t.Root(), nil)

	// // TO DO write the genotypes to file? (if so, do this in align.go)

	// // if we are going to write a CSV of the children of transitions:
	// if len(childrenOut) > 0 {
	// 	err = characterio.WriteChildren(childrenOut, characters, atlas)
	// 	if err != nil {
	// 		return err
	// 	}
	// }

	return nil
}

var treeFile string
var alignmentFile string
var variantsConfig string
var genbankFile string
var tipFile string
var algorithmUp string   // which algorithm to use for the uppass when there are polytomies (Madison 1989)
var algorithmDown string // which algorithm to use for resolving ties (Acctrans/Deltrans etc.)
var annotateNodes bool
var annotateTips bool
var treeOut string
var threshold int
var childrenOut string
var summarize bool
var civet bool
var nuc bool
var p bool
var common_anc bool
var epi bool
var outgroup string
var rescale bool
var threads int

var mainCmd = &cobra.Command{
	Use:   "ash",
	Short: "ancestral state helper",
	Long: `ancestral state helper

Example usage:

./ash --treefile tree.newick --alignment sequences.fasta --variants-config config --genbank MN908947.gb --algo-down deltrans --threshold 1
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = ash(treeFile, alignmentFile, variantsConfig, genbankFile, tipFile,
			algorithmUp, algorithmDown, annotateNodes, annotateTips, threshold,
			treeOut, childrenOut, summarize, civet, nuc, p, epi, common_anc, outgroup, rescale,
			threads)

		return
	},
}

func init() {

	mainCmd.Flags().StringVarP(&treeFile, "treefile", "", "", "Tree file to read - must be in newick format, must be rooted")
	mainCmd.Flags().StringVarP(&alignmentFile, "alignment", "", "", "Fasta format alignment to read")
	mainCmd.Flags().StringVarP(&variantsConfig, "config", "", "", "Variants to type in the alignment")
	mainCmd.Flags().StringVarP(&genbankFile, "genbank", "", "", "Genbank format annotation of a sequence in the same coordinates as the alignment")
	mainCmd.Flags().StringVarP(&tipFile, "tipfile", "", "", "CSV format table of tip to character relationships (instead of --alignment, --variants-config and --genbank)")
	mainCmd.Flags().StringVarP(&algorithmUp, "algo-up", "", "hard", "Algorithm to use for dealing with polytomies (choose one of soft/hard)")
	mainCmd.Flags().StringVarP(&algorithmDown, "algo-down", "", "", "Algorithm to use for breaking ties (choose one of acctrans/deltrans/downpass)")
	mainCmd.Flags().IntVarP(&threshold, "threshold", "", 0, "Threshold number of children, above which a transition will be included in the output (default: 0)")
	mainCmd.Flags().StringVarP(&treeOut, "tree-out", "", "", "Tree file to write (optionally) - will be in nexus format")
	mainCmd.Flags().BoolVarP(&annotateNodes, "annotate-nodes", "", false, "Annotate internal nodes of output tree with inferred states (default: false)")
	mainCmd.Flags().BoolVarP(&annotateTips, "annotate-tips", "", false, "Annotate tips of output tree with known states (default: false)")
	mainCmd.Flags().StringVarP(&childrenOut, "children-out", "", "", "CSV format file of the children of transitions to write (optionally)")
	mainCmd.Flags().BoolVarP(&summarize, "summarize-children", "", false, "Optionally summarize the counts of children with different states under each transition to stdout")
	mainCmd.Flags().BoolVarP(&civet, "civet", "", false, "annotate all amino acid changes + neutral nucleotide changes")
	mainCmd.Flags().BoolVarP(&nuc, "nuc", "", false, "annotate all nucleotide changes")
	mainCmd.Flags().BoolVarP(&p, "paper", "", false, "do papery things")
	mainCmd.Flags().BoolVarP(&epi, "epistasis", "", false, "do epistasis things")
	mainCmd.Flags().BoolVarP(&common_anc, "common_anc", "", false, "do common_anc things")
	mainCmd.Flags().StringVarP(&outgroup, "outgroup", "", "", "the outgroup")
	mainCmd.Flags().BoolVarP(&rescale, "rescale", "", false, "rescale --tree-out so branch lengths are inferred # nuc substitutions")
	mainCmd.Flags().IntVarP(&threads, "threads", "t", 1, "number of threads to use for epistasis")

	mainCmd.Flags().Lookup("annotate-nodes").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("annotate-tips").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("summarize-children").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("civet").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("nuc").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("paper").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("epistasis").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("common_anc").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("rescale").NoOptDefVal = "true"

	mainCmd.Flags().SortFlags = false
}

func main() {
	// f, err := os.Create("CPU.prof")
	// if err != nil {
	// 	log.Fatal("could not create CPU profile: ", err)
	// }
	// defer f.Close() // error handling omitted for example
	// if err := pprof.StartCPUProfile(f); err != nil {
	// 	log.Fatal("could not start CPU profile: ", err)
	// }
	// defer pprof.StopCPUProfile()

	mainCmd.Execute()

	// f2, err := os.Create("mem.prof")
	// if err != nil {
	// 	log.Fatal("could not create memory profile: ", err)
	// }
	// defer f2.Close() // error handling omitted for example
	// runtime.GC()     // get up-to-date statistics
	// if err := pprof.WriteHeapProfile(f2); err != nil {
	// 	log.Fatal("could not write memory profile: ", err)
	// }
}
