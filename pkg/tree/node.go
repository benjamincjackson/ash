package tree

import (
	"bytes"
	"errors"
	"fmt"
	"strconv"
	"strings"
)

// Node structure
type Node struct {
	name       string   // Name of the node
	comment    []string // Comment if any in the newick file
	neigh      []*Node  // neighbors array
	br         []*Edge  // Branches array (same order than neigh)
	depth      int      // Depth of the node
	id         int      // This field is used at discretion of the user to store information
	tipid      int      // This is used by TipIndex to match tip names of different trees
	Upstates   [][]byte // character states will be encoded here
	Downstates [][]byte // character states will be encoded here
}

// Uninitialized depth is coded as -1
const (
	NIL_DEPTH = -1
)

func (n *Node) SetUpstates(states [][]byte) {
	n.Upstates = states
}

func (n *Node) SetDownstates(states [][]byte) {
	n.Downstates = states
}

// Adds a child n to the node p, connected with edge e
func (p *Node) addChild(n *Node, e *Edge) {
	p.neigh = append(p.neigh, n)
	p.br = append(p.br, e)
}

// Returns the name of the node
func (n *Node) Name() string {
	return n.name
}

// Sets the name of the node. No verification if another node
// has the same name
func (n *Node) SetName(name string) {
	n.name = name
}

// Returns the Id of the node. Id==NIL_ID means that
// it has not been set yet.
func (n *Node) Id() int {
	return n.id
}

// Sets the id of the node
func (n *Node) SetId(id int) {
	n.id = id
}

// List of neighbors of this node
func (n *Node) Neigh() []*Node {
	return n.neigh
}

// Number of neighbors of this node
func (n *Node) Nneigh() int {
	return len(n.neigh)
}

// Is a tip or not?
func (n *Node) Tip() bool {
	return len(n.neigh) == 1
}

// Adds a comment to the node. It will be coded by a list of []
// In the Newick format.
func (n *Node) AddComment(comment string) {
	n.comment = append(n.comment, comment)
}

// Returns the string of comma separated comments
// surounded by [].
func (n *Node) CommentsString() string {
	var buf bytes.Buffer
	buf.WriteRune('[')
	for i, c := range n.comment {
		if i > 0 {
			buf.WriteRune(',')
		}
		buf.WriteString(c)
	}
	buf.WriteRune(']')
	return buf.String()
}

// List of edges going from this node
func (n *Node) Edges() []*Edge {
	return n.br
}

// Returns the depth of the node. Returns an error if the depth has
// not been set yet.
func (n *Node) Depth() (int, error) {
	if n.depth == NIL_DEPTH {
		return n.depth, errors.New("Node depth has not been computed")
	}
	return n.depth, nil
}

// aggregate k=v pairs into the correct format for metadata comments:
// AA={"S:L1234I", "...", ...}
func aggregateComments(comments []string) string {
	m := make(map[string][]string)
	for _, comment := range comments {
		sa := strings.Split(comment, "=")
		k := sa[0]
		v := sa[1]
		if _, ok := m[k]; ok {
			m[k] = append(m[k], "\""+v+"\"")
		} else {
			m[k] = []string{"\"" + v + "\""}
		}
	}
	newstring := "&"
	list := make([]string, 0)
	for k, vs := range m {
		temp := k + "={"
		temp = temp + strings.Join(vs, ",")
		temp = temp + "}"
		list = append(list, temp)
	}
	newstring = newstring + strings.Join(list, ",")

	return newstring
}

// Recursive function that outputs newick representation
// from the current node
func (n *Node) Newick(parent *Node, newick *bytes.Buffer) {
	if len(n.neigh) > 0 {
		if len(n.neigh) > 1 {
			newick.WriteString("(")
		}
		nbchild := 0
		for i, child := range n.neigh {
			if child != parent {
				if nbchild > 0 {
					newick.WriteString(",")
				}
				child.Newick(n, newick)
				if n.br[i].support != NIL_SUPPORT && child.Name() == "" {
					newick.WriteString(strconv.FormatFloat(n.br[i].support, 'f', -1, 64))
					if n.br[i].pvalue != NIL_PVALUE {
						newick.WriteString(fmt.Sprintf("/%s", strconv.FormatFloat(n.br[i].pvalue, 'f', -1, 64)))
					}
				}
				if len(child.comment) != 0 {
					for _, c := range child.comment {
						newick.WriteString("[")
						newick.WriteString(c)
						newick.WriteString("]")
					}
				}
				if n.br[i].length != NIL_LENGTH {
					newick.WriteString(":")
					newick.WriteString(strconv.FormatFloat(n.br[i].length, 'f', -1, 64))
				}
				if len(n.br[i].comment) != 0 {
					newick.WriteString("[")
					newick.WriteString(aggregateComments(n.br[i].comment))
					newick.WriteString("]")
					// for _, c := range n.br[i].comment {
					// 	newick.WriteString("[")
					// 	newick.WriteString(c)
					// 	newick.WriteString("]")
					// }
				}
				nbchild++
			}
		}
		if len(n.neigh) > 1 {
			newick.WriteString(")")
		}
	}
	newick.WriteString(n.name)
}

func (n *Node) NewickOptionalComments(parent *Node, newick *bytes.Buffer, annotate_nodes bool, annotate_tips bool) {
	if len(n.neigh) > 0 {
		if len(n.neigh) > 1 {
			newick.WriteString("(")
		}
		nbchild := 0
		for i, child := range n.neigh {
			if child != parent {
				if nbchild > 0 {
					newick.WriteString(",")
				}
				child.NewickOptionalComments(n, newick, annotate_nodes, annotate_tips)
				if n.br[i].support != NIL_SUPPORT && child.Name() == "" {
					newick.WriteString(strconv.FormatFloat(n.br[i].support, 'f', -1, 64))
					if n.br[i].pvalue != NIL_PVALUE {
						newick.WriteString(fmt.Sprintf("/%s", strconv.FormatFloat(n.br[i].pvalue, 'f', -1, 64)))
					}
				}
				if len(child.comment) != 0 {
					if annotate_nodes && !child.Tip() {
						newick.WriteString("[&")
						newick.WriteString(strings.Join(child.comment, ","))
						newick.WriteString("]")
					}
					if annotate_tips && child.Tip() {
						newick.WriteString("[&")
						newick.WriteString(strings.Join(child.comment, ","))
						newick.WriteString("]")
					}
				}
				if n.br[i].length != NIL_LENGTH {
					newick.WriteString(":")
					newick.WriteString(strconv.FormatFloat(n.br[i].length, 'f', -1, 64))
				}
				if len(n.br[i].comment) != 0 {
					newick.WriteString("[")
					newick.WriteString(aggregateComments(n.br[i].comment))
					newick.WriteString("]")
				}
				nbchild++
			}
		}
		if len(n.neigh) > 1 {
			newick.WriteString(")")
		}
	}
	newick.WriteString(n.name)
}
