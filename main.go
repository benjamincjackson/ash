package main

import (
	"bufio"
	"bytes"
	"encoding/gob"
	"errors"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"

	"github.com/spf13/cobra"

	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/characterio"
	"github.com/benjamincjackson/ash/pkg/newick"
	"github.com/benjamincjackson/ash/pkg/parsimony"
	"github.com/benjamincjackson/ash/pkg/tree"
)

// to do possibly - sanity check arguments if --civet is given
func checkArgs(treeFile string, alignmentFile string, variantsConfig string, genbankFile string, tipFile string,
	algorithmUp string, algorithmDown string, treeOut string, civet bool) (int, int, string, error) {

	algoUp := -1
	switch algorithmUp {
	case "hard":
		algoUp = 0
	case "soft":
		algoUp = 1
	default:
		return -1, -1, "", errors.New("unknown up-pass algorithm: choose one of soft or hard")
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
		return -1, -1, "", errors.New("unknown down-pass algorithm: choose one of acctrans, deltrans or downpass")
	}

	if len(alignmentFile) > 0 && len(tipFile) > 0 {
		return algoUp, algoDown, "", errors.New("either use a --tipfile OR an --alignment and a --variants-config file, not a mixture")
	}
	if len(variantsConfig) > 0 && len(tipFile) > 0 {
		return algoUp, algoDown, "", errors.New("either use a --tipfile OR an --alignment and a --variants-config file, not a mixture")
	}
	if len(variantsConfig) > 0 && len(alignmentFile) == 0 || len(variantsConfig) == 0 && len(alignmentFile) > 0 {
		return algoUp, algoDown, "", errors.New("if you provide an --alignment file you must provide a --variants-config file, and vice versa")
	}

	if len(variantsConfig) > 0 {
		aa := false
		f, err := os.Open(variantsConfig)
		if err != nil {
			return algoUp, algoDown, "", err
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
			return algoUp, algoDown, "", err
		}
		if aa && len(genbankFile) == 0 {
			return algoUp, algoDown, "", errors.New("you must provide a --genbank file if there are amino acids in your --variants-config file")
		}
	}

	var s string
	if len(alignmentFile) > 0 {
		s = "alignment"
	} else {
		s = "csv"
	}

	return algoUp, algoDown, s, nil
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

func getRealSizeOf(v interface{}) (int, error) {
	b := new(bytes.Buffer)
	if err := gob.NewEncoder(b).Encode(v); err != nil {
		return 0, err
	}
	return b.Len(), nil
}

func ash(treeIn string, alignmentFile string, variantsConfig string, genbankFile string, tipFile string,
	algorithmUp string, algorithmDown string, annotateNodes bool, annotateTips bool, threshold int,
	treeOut string, childrenOut string, summarize bool, civet bool) error {

	// algoUp, algoDown, input, err := checkArgs(treeIn, alignmentFile, variantsConfig, genbankFile, tipFile, algorithmUp, algorithmDown, treeOut)
	algoUp, algoDown, input, err := checkArgs(treeIn, alignmentFile, variantsConfig, genbankFile, tipFile, algorithmUp, algorithmDown, treeOut, civet)
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

	switch input {
	case "alignment":
		characterStates, idx, states, err = characterio.TypeAlignment(t, alignmentFile, variantsConfig, genbankFile)
		if err != nil {
			return err
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
	// time.Sleep(10 * time.Second)

	// fmt.Println(characterStates)

	/*
		TO DO- check all the tree's tips are in the character state input (we do the converse already when we read the character states in)
	*/

	// fmt.Println("starting uppass")
	// TO DO- maybe just use the hard polytomies interpretation?
	switch algoUp {
	case 0: // hard polytomies
		parsimony.UpPass(t, 0, states, idx)
	case 1: // soft polytomies (resolve them [separately for each character!])
		parsimony.UpPass(t, 1, states, idx)
	}
	// fmt.Println("finished uppass")
	// time.Sleep(10 * time.Second)

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
	// fmt.Println("finished downpass")
	// time.Sleep(10 * time.Second)

	// fmt.Println(t.Root().Id())
	// fmt.Println(t.Root().Nneigh())
	// fmt.Println(t.Rooted())

	// t.PreOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) bool {
	// 	fmt.Print(cur.Name() + " ")
	// 	for i := range idx {
	// 		start := idx[i].Start
	// 		stop := idx[i].Stop
	// 		ba := downstates[cur.Id()][start:stop]
	// 		for _, k := range bitsets.GetSetBits(ba) {
	// 			fmt.Print(characterStates[i].StateKey[k-1])
	// 		}
	// 		fmt.Println()
	// 	}
	// 	return true
	// })

	switch civet {
	case true:
		// genbank annotation parsing:
		features, err := annotation.GetRegions(genbankFile)
		if err != nil {
			return err
		}
		// for _, f := range features {
		// 	fmt.Println(f)
		// }

		transitions := make([][]characterio.Transition, len(characterStates), len(characterStates))
		for i := range transitions {
			transitions[i] = make([]characterio.Transition, 0)
		}

		parsimony.LabelChangesAnno(t, features, characterStates, states, transitions)

	default:
		// we can keep all the transitions (structs containing pointers to the relevant nodes and edge) in an array
		// that matches the dimensions of the characters
		transitions := make([][]characterio.Transition, len(characterStates), len(characterStates))
		for i := range transitions {
			transitions[i] = make([]characterio.Transition, 0)
		}

		// for i := range transitions {
		// 	for j := range transitions[i] {
		// 		fmt.Println(transitions[i][j])
		// 	}
		// }

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

	f, err := os.Create("branchlengths.tsv")
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
		f.WriteString(strconv.Itoa(i) + "\t" + strconv.FormatFloat(e.Length()*29903, 'f', 4, 64) + "\t" + tip + "\t" + strconv.Itoa(len(e.GetComments())) + "\t" + strings.Join(e.GetComments(), " ") + "\n")
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

	// t.SortNeighborsByTips(t.Root(), nil)

	// fmt.Println()

	// // summarise the transitions (to stdout)
	// if summarize {
	// 	// TO DO: do this as we traverse the tree - don't keep everything in a map, if possible
	// 	characterio.SummarizeTransitions(characters, atlas)
	// }

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
			treeOut, childrenOut, summarize, civet)

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

	mainCmd.Flags().Lookup("annotate-nodes").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("annotate-tips").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("summarize-children").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("civet").NoOptDefVal = "true"

	mainCmd.Flags().SortFlags = false
}

func main() {
	f, err := os.Create("CPU.prof")
	if err != nil {
		log.Fatal("could not create CPU profile: ", err)
	}
	defer f.Close() // error handling omitted for example
	if err := pprof.StartCPUProfile(f); err != nil {
		log.Fatal("could not start CPU profile: ", err)
	}
	defer pprof.StopCPUProfile()

	mainCmd.Execute()

	f2, err := os.Create("mem.prof")
	if err != nil {
		log.Fatal("could not create memory profile: ", err)
	}
	defer f2.Close() // error handling omitted for example
	runtime.GC()     // get up-to-date statistics
	if err := pprof.WriteHeapProfile(f2); err != nil {
		log.Fatal("could not write memory profile: ", err)
	}
}
