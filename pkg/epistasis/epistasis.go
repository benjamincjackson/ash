package epistasis

import (
	"github.com/benjamincjackson/ash/pkg/annotation"
	"github.com/benjamincjackson/ash/pkg/tree"
)

// type Pair struct {
// }

// a struct containing
type Pair struct {
	m_ij int // number of

}

// the tree's branches are labelled with "syn=" for non-protein-changing nucleotide changes,
// and "AA=" for amino acid changes
// edge.SynLen is also set to the inferred synonymous branch length
func Epistasis(t *tree.Tree, features []annotation.Region) {

}

func ij(t *tree.Tree, i, j string) {

}

func ij_recur(curEdge, prevEdge *tree.Edge, curNode *tree.Node, synDist int, foundOne bool, pairs *Pair) {

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

		// update synDist & foundOne...

		// we carry on recurring down the tree
		ij_recur(e, curEdge, curNode.Neigh()[i], synDist, foundOne, pairs)
	}
}
