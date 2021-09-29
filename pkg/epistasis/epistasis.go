package epistasis

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/tree"
)

// type Pair struct {
// }

// a struct containing info about one pair's consecutive mutations over the tree
type Pair struct {
	i     string    // the name of the first site
	j     string    // the name of the second site
	m_ij  int       // number of branches with changes at BOTH sites i and j on
	t_pi  []float64 // list of synonymous distances between pairs of nonsynonymous changes at sites i and j
	E_tau float64
}

func (p *Pair) add_m_ij() {
	p.m_ij++
}

func (p *Pair) append_t_pi(i float64) {
	p.t_pi = append(p.t_pi, i)
}

func (p *Pair) get_i() string {
	return p.i
}

func (p *Pair) get_j() string {
	return p.j
}

// func (p *Pair) get_t_pi() []float64 {
// 	return p.t_pi
// }

// func get_aa_pairs(features []annotation.Region) []Pair {
// 	justCDS := make([]annotation.Region, 0)
// 	for i := range features {
// 		if features[i].Whichtype == "CDS" {
// 			justCDS = append(justCDS, features[i])
// 		}
// 	}

// 	pairs := make([]Pair, 0)
// 	for _, cds1 := range justCDS {
// 		for _, cds2 := range justCDS {
// 			for i := range cds1.Codonstarts {
// 				for j := range cds2.Codonstarts {
// 					// not with itself:
// 					if i == j {
// 						continue
// 					}
// 					// otherwise:
// 					pairs = append(pairs, Pair{i: cds1.Name + ":" + strconv.Itoa(i+1), j: cds2.Name + ":" + strconv.Itoa(i+1)})
// 				}
// 			}
// 		}
// 	}

// 	return pairs
// }

// following Kryazhimskiy, Sergey, et al. "Prevalence of epistasis in the evolution of influenza A surface proteins." PLoS genetics 7.2 (2011): e1001301,
// get the mean synonymous distance between randomly chosen pairs of nonsynonymous substitutions from the tree
func get_tau(t *tree.Tree) float64 {
	tns := Tns{}

	tau_recur(nil, t.Root(), &tns)

	sum_distances := 0.0
	for i := range tns.distances {
		sum_distances += tns.distances[i]
	}

	sum_weights := 0
	for i := range tns.weights {
		sum_weights += tns.weights[i]
	}

	tau := sum_distances / float64(sum_weights)

	return tau
}

// for recording the synonymous distances between nonsynon substitutions
type Tns struct {
	distances []float64
	weights   []int
}

func (tns *Tns) add_distance(d float64) {
	tns.distances = append(tns.distances, d)
}

func (tns *Tns) add_weight(w int) {
	tns.weights = append(tns.weights, w)
}

func tau_recur(prevEdge *tree.Edge, curNode *tree.Node, tns *Tns) {

	// base state: if we have reached a tip, we don't want to do anything (per the original paper -
	// mutations on external edges are discounted)
	if curNode.Tip() {
		return
	}

	for i, e := range curNode.Edges() {

		// skip the rootward edge
		if e == prevEdge {
			continue
		}

		// first, if there are no non-synonymous mutations on this branch, we can skip the next two steps
		AAs := e.Get_AA_comments()
		switch {
		// if there is 1, we set the synonymous distance and collect the downstream muts starting from here:
		case len(AAs) == 1:
			synDist := e.SynLen / 2
			n := 1
			collect_distances(prevEdge, curNode, tns, n, synDist)

		// if there is more than one, we first add all possible pairs for this branch, the collect the downstream muts
		case len(AAs) > 1:

			for j := range AAs {
				for k := j + 1; k < len(AAs); k++ {
					_ = k
					tns.add_distance(float64(len(AAs)) * float64(len(AAs)) * (e.SynLen / 3))
					tns.add_weight(len(AAs) * len(AAs))
				}
			}

			synDist := e.SynLen / 2
			// n is the number of non-synon mutations on this branch - need this to get all the pairs
			n := len(AAs)
			collect_distances(prevEdge, curNode, tns, n, synDist)
		}

		// then we carry on recurring down the tree
		tau_recur(e, curNode.Neigh()[i], tns)
	}
}

func collect_distances(prevEdge *tree.Edge, curNode *tree.Node, tns *Tns, n int, synDist float64) {

	// base state: if we have reached a tip, we don't want to do anything (per the original paper -
	// mutations on external edges are discounted)
	if curNode.Tip() {
		return
	}

	for i, e := range curNode.Edges() {

		// skip the rootward edge
		if e == prevEdge {
			continue
		}

		AAs := e.Get_AA_comments()
		if len(AAs) > 0 {
			tns.add_distance(float64(n) * float64(len(AAs)) * (synDist + (e.SynLen / 2)))
			tns.add_weight(n * len(AAs))
		}

		d := synDist + e.SynLen

		// we carry on recurring down the tree
		collect_distances(e, curNode.Neigh()[i], tns, n, d)
	}
}

// Traverse the tree, storing the information about the changes present on branches.
// Filter the changes based on some criteria (not external branches, > 1 occurence).
// Make the pairs
func get_pairs_from_tree(t *tree.Tree) []Pair {

	m := make(map[string]int)
	for _, e := range t.Edges() {
		// skip external branches
		if e.Right().Tip() {
			continue
		}
		AAs := e.Get_AA_comments()
		for _, AA := range AAs {
			if _, ok := m[AA]; ok {
				m[AA]++
			} else {
				m[AA] = 1
			}
		}
	}

	tokeep := make([]string, 0)
	for key, value := range m {
		if value > 1 {
			tokeep = append(tokeep, key)
		}
	}

	sort.Strings(tokeep)

	pairs := make([]Pair, 0)
	for i := range tokeep {
		for j := i + 1; j < len(tokeep); j++ {
			aa1 := tokeep[i]
			aa2 := tokeep[j]
			pairs = append(pairs, Pair{i: aa1, j: aa2})
		}
	}

	return pairs
}

// the tree's branches are labelled with "syn=" for non-protein-changing nucleotide changes,
// and "AA=" for amino acid changes
// edge.SynLen is also set to the inferred synonymous branch length
func Epistasis(t *tree.Tree, features []annotation.Region, threads int) {
	// the mean synonymous distance between non-synonymous pairs:
	tau := get_tau(t)
	fmt.Println("tau is: " + strconv.FormatFloat(tau, 'f', 8, 64))
	fmt.Println()

	pairs := get_pairs_from_tree(t)
	// fmt.Println(len(pairs))
	// os.Exit(0)

	if threads > len(pairs) {
		threads = len(pairs)
	}

	runtime.GOMAXPROCS(threads)

	chunkSize := int(math.Floor(float64(len(pairs)) / float64(threads)))

	var wgPairs sync.WaitGroup
	wgPairs.Add(threads)

	for i := 0; i < threads; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if i == threads-1 {
			end = len(pairs)
		}
		go func() {
			processChunk(t, pairs[start:end])
			wgPairs.Done()
		}()
	}

	wgPairs.Wait()

	for i := range pairs {
		pairs[i].E_tau = calc_E_tau(pairs[i], tau)
		// fmt.Println(calc_E_tau(pair, tau))
		// fmt.Println(pair.E_tau)
		// fmt.Println()
	}

	sort.SliceStable(pairs, func(i, j int) bool {
		return pairs[i].E_tau > pairs[j].E_tau
	})

	for _, pair := range pairs {
		fmt.Println(pair)
	}

	// fmt.Println(pairs[72])
}

func processChunk(t *tree.Tree, pairs []Pair) {
	for i := range pairs {
		ij(t, &(pairs[i]))
	}

	return
}

// {AA=orf1ab:4387 AA=orf1ab:6485 0 []}

func ij(t *tree.Tree, p *Pair) {
	ij_recur(nil, t.Root(), p, [2]bool{false, false}, 0.0)
	// fmt.Println(*p)
}

func ij_recur(prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open [2]bool, synDist float64) {

	// base state: if we have reached a tip, we don't want to do anything (per the original paper -
	// mutations on external edges are discounted)
	if curNode.Tip() {
		return
	}

	for i, e := range curNode.Edges() {

		// skip the rootward edge
		if e == prevEdge {
			continue
		}

		// here is the predefined synonymous length of this branch:
		edgeSynDist := e.SynLen

		// check to see if there are any relevant changes on this branch
		// TO DO - move this to its own function?
		i_present := false
		j_present := false
		var site string
		for _, comment := range e.GetComments() {
			if strings.HasPrefix(comment, "syn=") {
				continue
			}
			site = strings.Join(strings.Split(comment, ":")[0:2], ":")
			if site == pair.get_i() {
				i_present = true
			}
			if site == pair.get_j() {
				j_present = true
			}
		}

		var d float64
		o := open
		// if this branch has a relevant mutation:
		if i_present || j_present {

			// then do the appropriate thing...
			switch whatToDo(open, i_present, j_present) {
			case -2:
				panic("got a -2 from whatToDo()")
			case -1: // i is not open, j is not open, but one or both are present here (for the first time ever in this tree)
				// set synonymous distance to half this branch's length,
				// update open appropriately
				d = float64(edgeSynDist) / 2
				o = [2]bool{i_present, j_present}
			case 0: // i is open, j is not open, i is present, j is not | 1010
				// reset synonymous distance to half this branch's length,
				// update open
				d = float64(edgeSynDist) / 2
				o = [2]bool{true, false}
			case 1: // i is not open, j is open, i is not present, j is present | 0101
				// reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				o = [2]bool{false, true}
			case 2: // i is open, j is not open, i is not present, j is present | 1001
				// record an i_j pair, reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				o = [2]bool{false, true}
			case 3: // i is not open, j is open, i is present, j is not present | 0110
				// record a j_i pair, reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				o = [2]bool{true, false}
			case 4: // i is open, j is not open, i is present, j is present | 1011
				// record two i_j pairs one with a long and one with a short distance, record a j_i pair with short distance
				// add a count to m_ij
				// reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				pair.append_t_pi(float64(edgeSynDist) / 3)
				pair.append_t_pi(float64(edgeSynDist) / 3)
				pair.add_m_ij()
				o = [2]bool{true, true}
			case 5: // i is not open, j is open, i is present, j is present | 0111
				// record two j_i pairs one with a long and one with a short distance, record an i_j pair with short distance
				// add a count to m_ij
				// reset synonymous distance to half this branch's length
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				pair.append_t_pi(float64(edgeSynDist) / 3)
				pair.append_t_pi(float64(edgeSynDist) / 3)
				pair.add_m_ij()
				o = [2]bool{true, true}
			case 6: // i is open, j is open, i is present, j is not present | 1110
				// record a j_i pair, reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				o = [2]bool{true, false}
			case 7: // i is open, j is open, i is not present, j is present | 1101
				// record an i_j pair, reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				o = [2]bool{false, true}
			case 8: // i is open, j is open, i is present, j is present | 1111
				// record two (short and long) i_j pairs, record two (short and long) j_i pairs
				// add a count to m_ij
				// reset synonymous distance to half this branch's length
				// update open
				d = float64(edgeSynDist) / 2
				pair.append_t_pi(synDist + d)
				pair.append_t_pi(synDist + d)
				pair.append_t_pi(float64(edgeSynDist) / 3)
				pair.append_t_pi(float64(edgeSynDist) / 3)
				pair.add_m_ij()
				o = [2]bool{true, true}
			}

		} else { // otherwise, we just increment the synonymous distance and retain what's open:
			d = synDist + float64(edgeSynDist)
			o = open
		}

		// func ij_recur(prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open [2]bool, synDist float64) {
		// then we carry on recurring down the tree
		ij_recur(e, curNode.Neigh()[i], pair, o, d)
	}
}

func whatToDo(open [2]bool, i_present, j_present bool) int {
	// nothing has occurred yet (e.g. at the root):
	if !open[0] && !open[1] {
		return -1
	}
	// i is open, j is not, i is present, j is not
	if open[0] && !open[1] && i_present && !j_present {
		return 0
	}
	// i is not open, j is open, i is not present, j is present
	if !open[0] && open[1] && !i_present && j_present {
		return 1
	}
	// i is open, j is not open, i is not present, j is present
	if open[0] && !open[1] && !i_present && j_present {
		return 2
	}
	// i is not open, j is open, i is present, j is not present
	if !open[0] && open[1] && i_present && !j_present {
		return 3
	}
	// i is open, j is not open, i is present, j is present
	if open[0] && !open[1] && i_present && j_present {
		return 4
	}
	// i is not open, j is open, i is present, j is present
	if !open[0] && open[1] && i_present && j_present {
		return 5
	}
	// i is open, j is open, i is present, j is not present
	if open[0] && open[1] && i_present && !j_present {
		return 6
	}
	// i is open, j is open, i is not present, j is present
	if open[0] && open[1] && !i_present && j_present {
		return 7
	}
	// i is open, j is open, i is present, j is present
	if open[0] && open[1] && i_present && j_present {
		return 8
	}

	return -2
}

// NB: don't do it this way - just take the mean of p.t_pi by dividing its sum by its length
// see equation (1) in Kryazhimskiy, Sergey, et al. "Prevalence of epistasis in the evolution of influenza A surface proteins." PLoS genetics 7.2 (2011): e1001301.
// func calc_E_tau(p Pair, tau float64) float64 {
// 	d := 1 / math.Pow(2, float64(p.m_ij))
// 	sum := 0.0
// 	for i := range p.t_pi {
// 		sum += math.Exp((p.t_pi[i] * -1) / tau)
// 	}
// 	E_tau := sum / d
// 	return E_tau
// }

func calc_E_tau(p Pair, tau float64) float64 {
	d := len(p.t_pi)
	if d == 0 {
		return 0.0
	}
	sum := 0.0
	for i := range p.t_pi {
		sum += math.Exp((p.t_pi[i] * -1) / tau)
	}
	E_tau := sum / float64(d)
	return E_tau
}
