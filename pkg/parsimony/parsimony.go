package parsimony

import (
	"github.com/benjamincjackson/ash/pkg/bitsets"
	"github.com/benjamincjackson/ash/pkg/characterio"
	"github.com/benjamincjackson/ash/pkg/tree"
)

// First pass of the parsimony reconstruction
//
// Tree must be sorted by node depth, and then we can (reverse) post-order traverse
// over it to push the states up the tree from the tips to the root, in-place
func UpPass(t *tree.Tree, algoUp int, states [][]byte, idx []characterio.StartStop) {
	uppass(t.Root(), nil, algoUp, states, idx)
}

func uppass(cur, prev *tree.Node, algoUp int, states [][]byte, idx []characterio.StartStop) {

	// Base state: if cur is a tip, there's nothing to do here
	if len(cur.Neigh()) == 1 {
		return
	}

	cur_id := cur.Id()

	// outermost dimension is node, next dimension is bitset character array
	downnodestates := make([][]byte, 0)

	// iterating in reverse order gives us the reverse post-order traversal,
	// ([only] needed because of the way we sorted the tree previous to this)
	for i := len(cur.Neigh()) - 1; i > -1; i-- {
		n := cur.Neigh()[i]
		n_id := n.Id()
		if n != prev {
			uppass(n, cur, algoUp, states, idx)
			downnodestates = append(downnodestates, states[n_id])
		}
	}

	uppassMove(states[cur_id], downnodestates, algoUp, idx)
}

func uppassMove(upnodestates []byte, downnodestates [][]byte, algoUp int, idx []characterio.StartStop) {
	// each i represents one character
	for i := range idx {
		// start and stop define the elements for this character's states in the slice of all characters
		start := idx[i].Start
		stop := idx[i].Stop

		// Here we switch on a bifurcating/multifurcating node
		switch {
		// bifurcating:
		case len(downnodestates) == 2:
			// for each character, switch on intersection being empty/not empty
			switch {
			case bitsets.IsAnyBitSet(bitsets.Intersection(downnodestates[0][start:stop], downnodestates[1][start:stop])): // intersection not empty (take the intersection)
				bitsets.InPlaceIntersection(upnodestates[start:stop], downnodestates[0][start:stop], downnodestates[1][start:stop])
			default: // intersection empty (take the union)
				bitsets.InPlaceUnion(upnodestates[start:stop], downnodestates[0][start:stop], downnodestates[1][start:stop])
			}
		// multifurcating:
		case len(downnodestates) > 2:
			// subset the upstate slices to get only this individual character for now
			subset := make([][]byte, len(downnodestates), len(downnodestates))
			for j := range downnodestates {
				subset[j] = downnodestates[j][start:stop]
			}
			// switches on treating polytomies as hard/soft
			switch algoUp {
			case 0: // hard
				bitsets.InPlaceVarMax(upnodestates[start:stop], subset)
			case 1: // soft
				_ = bitsets.InPlaceVarCover(upnodestates[start:stop], subset)
			}
		default: // (length == 0) this is a tip, we don't need to do anything (actually we should never get here)
			break
		}
	}
}

// NOTE- unclear to me if you need to do a downpass first.
// Swofford & Maddison + Gotree say no, but Felsenstein (2007, ch6, pp70) seems to say yes
// Conclusion is that you don't have to. (But you could)
func Acctrans(t *tree.Tree, states [][]byte, idx []characterio.StartStop) {
	acctrans(t.Root(), nil, states, idx)
}

func acctrans(cur, prev *tree.Node, states [][]byte, idx []characterio.StartStop) {
	cur_id := cur.Id()
	for _, n := range cur.Neigh() {
		if n != prev && !n.Tip() {
			n_id := n.Id()
			acctransMove(states, cur_id, n_id, idx)
			acctrans(n, cur, states, idx)
		}
	}
}

func acctransMove(states [][]byte, upnode, downnode int, idx []characterio.StartStop) {
	// we want to maximize the difference between the two nodes' sets for the downnode's states
	// in our case this is the set difference Downside-Upside, otherwise we don't have to do anything
	for i := range idx {
		start := idx[i].Start
		stop := idx[i].Stop
		if bitsets.IsAnyBitSet(bitsets.SetDiff(states[downnode][start:stop], states[upnode][start:stop])) {
			bitsets.InPlaceSetDiff(states[downnode][start:stop], states[downnode][start:stop], states[upnode][start:stop])
		} // else we dont have to do anything
	}
}

// The root -> tips pass, to get the MPRs
func DownPass(t *tree.Tree, states [][]byte, idx []characterio.StartStop) {
	downpass(t.Root(), nil, states, idx)
}

// recur down the tree, calling downpassMove at each interior node, before moving on
func downpass(cur, prev *tree.Node, states [][]byte, idx []characterio.StartStop) {

	// if this is not the root
	if prev != nil {
		// base state- if this is a tip, there is nothing to do here:
		if cur.Tip() {
			return
		}

		// the index of this node in the slice of states:
		cur_id := cur.Id()

		// make a slice for the indices of its neighbours (both rootwards and tipwards)
		n_ids := make([]int, 0)
		// and populate it
		for _, n := range cur.Neigh() {
			n_ids = append(n_ids, n.Id())
		}

		// we calculate this node's MPR set:
		downpassMove(states, cur_id, states[cur_id], n_ids, idx)
	}

	// then we may move tipwards down the tree
	for _, n := range cur.Neigh() {
		if n != prev {
			downpass(n, cur, states, idx)
		}
	}
}

// we act as if we have rerooted the tree at each interior node by applying the parsimony method
// to this node's ancestor's MPR set (that has previously been defined by calling this function) +
// all of its tipward neighbours' first-pass state sets.
func downpassMove(states [][]byte, node_id int, nodestates []byte, neighbour_ids []int, idx []characterio.StartStop) {
	// then we calculate the MPR sets for each character:
	for i := range idx {
		start := idx[i].Start
		stop := idx[i].Stop
		// we get all the neighbours in one place (for this character):
		neighbour_states := make([][]byte, len(neighbour_ids), len(neighbour_ids))
		for i, n_id := range neighbour_ids {
			neighbour_states[i] = states[n_id][start:stop]
		}
		switch len(neighbour_states) {
		case 3:
			bitsets.InPlaceThreeSetMPR(nodestates[start:stop], neighbour_states[0], neighbour_states[1], neighbour_states[2])
		default:
			bitsets.InPlaceVarMax(nodestates[start:stop], neighbour_states)
		}
	}
}

//
func Deltrans(t *tree.Tree, states [][]byte, idx []characterio.StartStop) {
	deltrans(t.Root(), nil, states, idx)
}

//
func deltrans(cur, prev *tree.Node, states [][]byte, idx []characterio.StartStop) {
	cur_id := cur.Id()
	for _, n := range cur.Neigh() {
		n_id := n.Id()
		if n != prev && !n.Tip() {
			deltransMove(states, cur_id, n_id, idx)
			deltrans(n, cur, states, idx)
		}
	}
}

//
func deltransMove(states [][]byte, upnode, downnode int, idx []characterio.StartStop) {
	// we want to minimize the difference between the two nodes' sets for the downnode's states
	// in our case this is the intersection between up and down if this set is not empty, otherwise we don't have to do anything
	for i := range idx {
		start := idx[i].Start
		stop := idx[i].Stop
		if bitsets.IsAnyBitSet(bitsets.Intersection(states[upnode][start:stop], states[downnode][start:stop])) {
			bitsets.InPlaceIntersection(states[downnode][start:stop], states[upnode][start:stop], states[downnode][start:stop])
		} // else we do nothing
	}
}
