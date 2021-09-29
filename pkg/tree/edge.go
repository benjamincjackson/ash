package tree

import (
	"strings"

	fbitset "github.com/fredericlemoine/bitset"
)

// Structure of an edge
// TO DO - add fields to contain synonymous, nonsynonymous (deletion?) changes, etc
type Edge struct {
	left, right *Node    // Left and right nodes
	length      float64  // length of branch
	comment     []string // Comment if any in the newick file
	support     float64  // -1 if no support
	pvalue      float64  // -1 if no pvalue
	// a Bit at index i in the bitset corresponds to the position of the tip i
	//left:0/right:1 .
	// i is the index of the tip in the sorted tip name array
	bitset        *fbitset.BitSet // Bitset of length Number of taxa each
	hashcoderight uint64          // HashCode related to Taxa on the left
	hashcodeleft  uint64          // HashCode related to Taxa on the right
	ntaxright     int             // Number of taxa below : Initialized with hashes / ReinitIndexes
	ntaxleft      int             // Number of taxa above : Initialized with hashes / ReinitIndexes
	id            int             // this field is used at discretion of the user to store information
	SynLen        float64         // the length of this branch in units of synonymous snps
}

// Constant for uninitialized values
const (
	NIL_SUPPORT = -1.0
	NIL_LENGTH  = -1.0
	NIL_PVALUE  = -1.0
	NIL_ID      = -1
	NIL_TIPID   = 0
)

// Sets left node (parent)
func (e *Edge) setLeft(left *Node) {
	e.left = left
}

// Sets right node (child)
func (e *Edge) setRight(right *Node) {
	e.right = right
}

// Sets the id of the branch
func (e *Edge) SetId(id int) {
	e.id = id
}

// Adds a comment to the edge. It will be coded by a list of []
// In the Newick format.
func (e *Edge) AddComment(comment string) {
	e.comment = append(e.comment, comment)
}

func (e *Edge) GetComments() []string {
	return e.comment
}

// get the comments that start "AA="
func (e *Edge) Get_AA_comments() []string {
	comments := e.GetComments()
	AAs := make([]string, 0)
	for _, c := range comments {
		info := strings.Split(c, ":")
		// if it's not an amino acid, skip it
		if strings.Split(info[0], "=")[0] != "AA" {
			continue
		}
		name := strings.Join(info[0:2], ":")
		AAs = append(AAs, name)
	}

	return AAs
}

// Sets the length of the branch
func (e *Edge) SetLength(length float64) {
	e.length = length
}

// returns the length of the branch
func (e *Edge) Length() float64 {
	return e.length
}

// Sets the branch support
func (e *Edge) SetSupport(support float64) {
	e.support = support
}

// Returns the support of that branch
func (e *Edge) Support() float64 {
	return e.support
}

// Sets the pvalue of this edge (if not null, pvalue
// is stored/parsed as "/pvalue" in the bootstrap value
// field.
func (e *Edge) SetPValue(pval float64) {
	e.pvalue = pval
}

// Returns the node at the right side of the edge (child)
func (e *Edge) Right() *Node {
	return e.right
}

// Returns the node at the left side of the edge (parent)
func (e *Edge) Left() *Node {
	return e.left
}
