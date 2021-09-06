package tree

func (t *Tree) PostOrder(f func(cur *Node, prev *Node, e *Edge) (keep bool)) {
	t.postOrderRecur(t.Root(), nil, nil, f)
}

func (t *Tree) postOrderRecur(cur *Node, prev *Node, e *Edge, f func(cur *Node, prev *Node, e *Edge) (keep bool)) (keep bool) {
	keep = true
	for i, n := range cur.neigh {
		if n != prev {
			if keep = t.postOrderRecur(n, cur, cur.Edges()[i], f); !keep {
				return
			}
		}
	}
	keep = keep && f(cur, prev, e)
	return
}

// reverse postorder traversal
func (t *Tree) PostOrderRev(f func(cur *Node, prev *Node, e *Edge) (keep bool)) {
	t.postOrderRecurRev(t.Root(), nil, nil, f)
}

func (t *Tree) postOrderRecurRev(cur *Node, prev *Node, e *Edge, f func(cur *Node, prev *Node, e *Edge) (keep bool)) (keep bool) {
	keep = true
	for i := len(cur.neigh) - 1; i > -1; i-- {
		n := cur.neigh[i]
		if n != prev {
			if keep = t.postOrderRecurRev(n, cur, cur.Edges()[i], f); !keep {
				return
			}
		}
	}
	keep = keep && f(cur, prev, e)
	return
}

func (t *Tree) PreOrder(f func(cur *Node, prev *Node, e *Edge) (keep bool)) {
	t.preOrderRecur(t.Root(), nil, nil, f)
}

func (t *Tree) preOrderRecur(cur *Node, prev *Node, e *Edge, f func(cur *Node, prev *Node, e *Edge) bool) bool {
	keepgoing := true
	if keepgoing = f(cur, prev, e); !keepgoing {
		return keepgoing
	}
	for i, n := range cur.neigh {
		if n != prev {
			if keepgoing = t.preOrderRecur(n, cur, cur.Edges()[i], f); !keepgoing {
				return keepgoing
			}
		}
	}
	return keepgoing
}

// func (t *Tree) PreOrderRev(f func(cur *Node, prev *Node, e *Edge) (keep bool)) {
// 	t.preOrderRecurRev(t.Root(), nil, nil, f)
// }

// func (t *Tree) preOrderRecurRev(cur *Node, prev *Node, e *Edge, f func(cur *Node, prev *Node, e *Edge) bool) bool {
// 	keepgoing := true
// 	if keepgoing = f(cur, prev, e); !keepgoing {
// 		return keepgoing
// 	}
// 	for i := len(cur.neigh) - 1; i > -1; i-- {
// 		n := cur.neigh[i]
// 		if n != prev {
// 			if keepgoing = t.preOrderRecurRev(n, cur, cur.Edges()[i], f); !keepgoing {
// 				return keepgoing
// 			}
// 		}
// 	}
// 	return keepgoing
// }
