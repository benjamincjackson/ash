package ancestry

import (
	"errors"

	"github.com/benjamincjackson/ash/pkg/tree"
)

// get the MRCA (node # of) of all the sequences that AREN'T the outgroup
func MRCA(t *tree.Tree, outgroup string) (int, error) {
	for _, node := range t.Root().Neigh() {
		if node.Tip() {
			if !(node.Name() == outgroup) {
				return -1, errors.New("didn't find the outgroup where I expected to")
			}
			continue
		} else {
			return node.Id(), nil
		}
	}

	return -1, errors.New("didn't find the outgroup where I expected to")
}
