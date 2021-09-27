package epistasis

import (
	"math"
	"strconv"
	"strings"

	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/tree"
)

// type Pair struct {
// }

// a struct containing info about one pair's consecutive mutations over the tree
type Pair struct {
	i    string    // the name of the first site
	j    string    // the name of the second site
	m_ij int       // number of branches with changes at BOTH sites i and j on
	t_pi []float64 // list of synonymous distances between pairs of nonsynonymous changes at sites i and j
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

func get_aa_pairs(features []annotation.Region) []Pair {
	justCDS := make([]annotation.Region, 0)
	for i := range justCDS {
		if justCDS[i].Whichtype == "CDS" {
			justCDS = append(justCDS, justCDS[i])
		}
	}

	pairs := make([]Pair, 0)
	for _, cds1 := range justCDS {
		for _, cds2 := range justCDS {
			for i := range cds1.Codonstarts {
				for j := range cds2.Codonstarts {
					// not with itself:
					if i == j {
						continue
					}
					// otherwise:
					pairs = append(pairs, Pair{i: cds1.Name + ":" + strconv.Itoa(i+1), j: cds2.Name + ":" + strconv.Itoa(i+1)})
				}
			}
		}

	}

	return pairs
}

// the tree's branches are labelled with "syn=" for non-protein-changing nucleotide changes,
// and "AA=" for amino acid changes
// edge.SynLen is also set to the inferred synonymous branch length
func Epistasis(t *tree.Tree, features []annotation.Region) {
	pairs := get_aa_pairs(features)
	// parallelise this eventually:
	for _, pair := range pairs {
		ij(t, &pair)
	}
}

func ij(t *tree.Tree, p *Pair) {
	ij_recur(nil, nil, t.Root(), p, [2]bool{false, false}, 0.0)
}

func ij_recur(curEdge, prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open [2]bool, synDist float64) {

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
				// reset synonymous distance to half this branch's length,
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

		// then we carry on recurring down the tree
		ij_recur(e, curEdge, curNode.Neigh()[i], pair, o, d)
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
	sum := 0.0
	for i := range p.t_pi {
		sum += math.Exp((p.t_pi[i] * -1) / tau)
	}
	E_tau := sum / float64(d)
	return E_tau
}
