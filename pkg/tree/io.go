package tree

import (
	"bytes"
	"strings"
)

// Modified from Newick()
func (t *Tree) NewickOptionalComments(annotate_nodes bool, annotate_tips bool) string {
	var buffer bytes.Buffer
	t.root.NewickOptionalComments(nil, &buffer, annotate_nodes, annotate_tips)
	if len(t.root.comment) != 0 {
		if annotate_nodes {
			buffer.WriteString("[&")
			buffer.WriteString(strings.Join(t.root.comment, ","))
			buffer.WriteString("]")
		}
	}
	buffer.WriteString(";")
	return buffer.String()
}

// Returns a newick string representation of this tree
func (t *Tree) Newick() string {
	var buffer bytes.Buffer
	t.root.Newick(nil, &buffer)
	if len(t.root.comment) != 0 {
		for _, c := range t.root.comment {
			buffer.WriteString("[")
			buffer.WriteString(c)
			buffer.WriteString("]")
		}
	}
	buffer.WriteString(";")
	return buffer.String()
}
