package tree

import (
	"bytes"
	"errors"
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// Tree structure having a root and a tip index, that maps tip names to their index
type Tree struct {
	root     *Node            // root node: If the tree is unrooted the root node should have 3 children
	tipIndex map[string]*Node // Map between tip name and Node
}

// Initialize a new empty Tree
func NewTree() *Tree {
	return &Tree{
		root:     nil,
		tipIndex: make(map[string]*Node, 0),
	}
}

// TEMPORARY:
func (t *Tree) PrintMap() {
	for k, v := range t.tipIndex {
		fmt.Print(k + ": ")
		fmt.Println(v.id)
	}
}

// Return the tip index if the tip with given name exists in the tree
// May return an error if tip index has not been initialized
// With UpdateTipIndex or if the tip does not exist
func (t *Tree) TipIndex(name string) (int, error) {
	if len(t.tipIndex) == 0 {
		return 0, errors.New("No tips in the index, tip name index is not initialized")
	}
	n, ok := t.tipIndex[name]
	if !ok {
		return 0, errors.New("No tip named " + name + " in the index")
	}
	return n.tipid, nil
}

func (t *Tree) TipId(name string) (int, error) {
	if len(t.tipIndex) == 0 {
		return 0, errors.New("No tips in the index, tip name index is not initialized")
	}
	n, ok := t.tipIndex[name]
	if !ok {
		return 0, errors.New("No tip named " + name + " in the index")
	}
	return n.id, nil
}

// Updates the tipindex which maps tip names to
// their index in the bitsets.
//
// Bitset indexes correspond to the position
// of the tip in the alphabetically ordered tip
// name list
func (t *Tree) UpdateTipIndex() (err error) {
	tips := t.SortedTips()
	for k := range t.tipIndex {
		delete(t.tipIndex, k)
	}
	for i, tip := range tips {
		if _, ok := t.tipIndex[tip.Name()]; ok {
			err = errors.New("Cannot create a tip index when several tips have the same name")
			return
		}
		t.tipIndex[tip.Name()] = tip
		tip.tipid = i
	}
	return
}

/* Tips, sorted by their order in the bitsets*/
func (t *Tree) SortedTips() []*Node {
	tips := t.Tips()
	sort.Slice(tips, func(i, j int) bool {
		return strings.Compare(tips[i].Name(), tips[j].Name()) < 0
	})
	return tips
}

// set the character states for a tip
func (t *Tree) SetTipStates(tip string, states [][]byte) error {
	if _, ok := t.tipIndex[tip]; ok {
		t.tipIndex[tip].SetUpstates(states)
		return nil
	} else {
		return errors.New("tip not found in tree: " + tip)
	}
}

// Initialize a new empty Node
func (t *Tree) NewNode() *Node {
	return &Node{
		name:       "",
		comment:    make([]string, 0),
		neigh:      make([]*Node, 0, 3),
		br:         make([]*Edge, 0, 3),
		depth:      NIL_DEPTH,
		id:         NIL_ID,
		tipid:      NIL_TIPID,
		Upstates:   make([][]byte, 0),
		Downstates: make([][]byte, 0),
	}
}

// Set a root for the tree. This does not check that the
// node is part of the tree. It may be useful to call
//	t.ReinitIndexes()
// After setting a new root, to update branch bitsets.
func (t *Tree) SetRoot(r *Node) {
	t.root = r
}

// Returns the current root of the tree
func (t *Tree) Root() *Node {
	return t.root
}

// Returns true if the tree is rooted (i.e. root node
// has 2 neighbors), and false otherwise.
func (t *Tree) Rooted() bool {
	return t.root.Nneigh() == 2
}

// Initialize a new empty Edge
func (t *Tree) NewEdge() *Edge {
	return &Edge{
		length:  NIL_LENGTH,
		support: NIL_SUPPORT,
		id:      NIL_ID,
		pvalue:  NIL_PVALUE,
	}
}

// Returns all the edges of the tree (do it recursively)
func (t *Tree) Edges() []*Edge {
	edges := make([]*Edge, 0, 2000)
	for _, e := range t.Root().br {
		edges = append(edges, e)
		t.edgesRecur(e, &edges)
	}
	return edges
}

// Recursive function to list all edges of the tree
func (t *Tree) edgesRecur(edge *Edge, edges *[]*Edge) {
	if len(edge.right.neigh) > 1 {
		for _, child := range edge.right.br {
			if child.left == edge.right {
				*edges = append((*edges), child)
				t.edgesRecur(child, edges)
			}
		}
	}
}

// Returns all the nodes of the tree (do it recursively)
func (t *Tree) Nodes() []*Node {
	nodes := make([]*Node, 0, 2000)
	t.nodesRecur(&nodes, nil, nil)
	return nodes
}

// recursive function that lists all nodes of the tree
func (t *Tree) nodesRecur(nodes *[]*Node, cur *Node, prev *Node) {
	if cur == nil {
		cur = t.Root()
	}
	*nodes = append((*nodes), cur)
	for _, n := range cur.neigh {
		if n != prev {
			t.nodesRecur(nodes, n, cur)
		}
	}
}

// Returns all the tips of the tree (do it recursively)
func (t *Tree) Tips() []*Node {
	tips := make([]*Node, 0)
	t.tipsRecur(&tips, nil, nil)
	return tips
}

// recursive function that lists all tips of the tree
func (t *Tree) tipsRecur(tips *[]*Node, cur *Node, prev *Node) {
	if cur == nil {
		cur = t.Root()
	}
	if cur.Tip() {
		*tips = append((*tips), cur)
	}
	for _, n := range cur.neigh {
		if n != prev {
			t.tipsRecur(tips, n, cur)
		}
	}
}

// Returns all the tip name in the tree
// Starts with n==nil (root)
func (t *Tree) AllTipNames() []string {
	names := make([]string, 0, 1000)
	t.allTipNamesRecur(&names, nil, nil)
	return names
}

// Returns all the tip name in the tree
// Starts with n==nil (root)
// It is an internal recursive function
func (t *Tree) allTipNamesRecur(names *[]string, n *Node, parent *Node) {
	if n == nil {
		n = t.Root()
	}
	// is a tip
	if len(n.neigh) == 1 {
		*names = append(*names, n.name)
	} else {
		for _, child := range n.neigh {
			if child != parent {
				t.allTipNamesRecur(names, child, n)
			}
		}
	}
}

// Connects the two nodes in argument by an edge that is returned.
func (t *Tree) ConnectNodes(parent *Node, child *Node) *Edge {
	newedge := t.NewEdge()
	newedge.setLeft(parent)
	newedge.setRight(child)
	parent.addChild(child, newedge)
	child.addChild(parent, newedge)
	return newedge
}

func (t *Tree) MaxDepthRooted(n *Node, prev *Node) int {
	maxdepth := 0
	// this is the base condition (if the current node is a tip, its depth is 0):
	if n.Tip() {
		// set its depth
		n.depth = 0
		// and return the amount to add to the previous node (for the recurrence)
		return 0
		// otherwise, we cycle through its neighbours
	} else {
		for _, neighbour := range n.neigh {
			// if this neighbour isn't its parent...
			if neighbour != prev {
				// calculate the depth (+1 because it's depth for n, not neighbour)
				depth := t.MaxDepthRooted(neighbour, n) + 1
				// if this depth is greater than the nodes current max depth, then set maxdepth to it
				if depth > maxdepth {
					maxdepth = depth
				}
			}
		}
	}
	// finally, we finish by setting this node's depth to maxdepth, and returning it (for the recurrence)
	n.depth = maxdepth
	return n.depth
}

// Sort neighbours of all nodes by their maximum depth
// so that when we post-order recur to come back up through the tree,
// we start with the tip(s) furthest from the root and each internal node
// is visited immediately after all of its descendants
func (t *Tree) SortNeighborsByDepth(cur, prev *Node) {
	// max depth at the end of each neighbor of cur
	neighbours := make([]struct {
		depth int
		neigh *Node
		br    *Edge
	}, len(cur.Neigh()))

	for i, n := range cur.Neigh() {
		neighbours[i].neigh = n
		neighbours[i].br = cur.Edges()[i]
		if n != prev {
			neighbours[i].depth = n.depth
			t.SortNeighborsByDepth(n, cur)
		} else {
			neighbours[i].depth = -1
		}
	}
	// we sort neighbor slice
	sort.SliceStable(neighbours, func(i, j int) bool { return neighbours[i].depth < neighbours[j].depth })
	// replace the original slices with the sorted ones
	for i, _ := range cur.Neigh() {
		cur.neigh[i] = neighbours[i].neigh
		cur.br[i] = neighbours[i].br
	}
}

// Sorts the neighbors of current node after having recursively
// sorted the neigbors of its children
func (t *Tree) SortNeighborsByTips(cur, prev *Node) int {
	// Number of tips at the end of each neighbor of cur
	neighbors := make([]struct {
		ntips int
		neigh *Node
		br    *Edge
	}, len(cur.Neigh()))
	total := 0
	for i, c := range cur.Neigh() {
		neighbors[i].neigh = c
		neighbors[i].br = cur.Edges()[i]
		if c != prev {
			neighbors[i].ntips = t.SortNeighborsByTips(c, cur)
			total += neighbors[i].ntips
		} else {
			neighbors[i].ntips = 0
		}
	}
	// we sort neighbor slice
	sort.SliceStable(neighbors, func(i, j int) bool { return neighbors[i].ntips < neighbors[j].ntips })
	//
	for i, _ := range cur.Neigh() {
		cur.neigh[i] = neighbors[i].neigh
		cur.br[i] = neighbors[i].br
	}

	if cur.Tip() {
		return 1
	} else {
		return total
	}
}

// Ben edit:
func (t *Tree) NexusOptionalComments(annotate_nodes bool, annotate_tips bool) string {
	newick := t.NewickOptionalComments(annotate_nodes, annotate_tips)
	var buffer bytes.Buffer
	buffer.WriteString("#NEXUS\n")
	buffer.WriteString("BEGIN TAXA;\n")
	tips := t.Tips()
	buffer.WriteString(" DIMENSIONS NTAX=")
	buffer.WriteString(strconv.Itoa(len(tips)))
	buffer.WriteString(";\n")
	buffer.WriteString(" TAXLABELS")
	for _, tip := range tips {
		buffer.WriteString(" " + tip.Name())
	}
	buffer.WriteString(";\n")
	buffer.WriteString("END;\n")
	buffer.WriteString("BEGIN TREES;\n")
	buffer.WriteString("  TREE tree1 = ")
	buffer.WriteString(newick)
	buffer.WriteString("\n")
	buffer.WriteString("END;\n")
	return buffer.String()
}
