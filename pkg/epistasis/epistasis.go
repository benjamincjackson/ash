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
	"github.com/benjamincjackson/gotree/tree"
)

// a reduced-size struct which just the result for one pair in it
type Result struct {
	i     string
	j     string
	E_tau float64
}

// a struct containing info about one pair's consecutive mutations over the tree
type Pair struct {
	i             string                   // the name of the first site
	j             string                   // the name of the second site
	m_ij          int                      // number of branches with changes at BOTH sites i and j on
	t_pi          []float64                // list of synonymous distances between pairs of nonsynonymous changes at sites i and j
	order         map[*tree.Edge][2]string // what is the order of the temporally unresolved mutations in this permutation
	E_tau_tempsum float64
	E_tau         float64
}

// func (p *Pair) add_m_ij() {
// 	p.m_ij++
// }

func (p *Pair) set_m_ij(m_ij int) {
	p.m_ij = m_ij
}

// see equation (1) in Kryazhimskiy, Sergey, et al. "Prevalence of epistasis in the evolution of influenza A surface proteins." PLoS genetics 7.2 (2011): e1001301.
func (p *Pair) add_E_tau(tau float64) {
	E_tau_temp := 0.0
	for i := range p.t_pi {
		E_tau_temp += math.Exp((p.t_pi[i] * -1) / tau)
	}
	p.E_tau_tempsum += E_tau_temp
}

func (p *Pair) calc_mean_E_tau() {
	d := 1 << p.m_ij
	p.E_tau = p.E_tau_tempsum / float64(d)
}

// n is the nth possible order, i is a value of t_pi
func (p *Pair) append_t_pi(i float64) {
	p.t_pi = append(p.t_pi, i)
}

func (p *Pair) reset_t_pi() {
	p.t_pi = make([]float64, 0)
}

func (p *Pair) get_i() string {
	return p.i
}

func (p *Pair) get_j() string {
	return p.j
}

// n is the nth possible order
func (p *Pair) get_order() map[*tree.Edge][2]string {
	return p.order
}

// // 1 << m_ij is the same as 2^^m_ij
// func (p *Pair) set_order_length(m_ij int) {
// 	if m_ij > 0 {
// 		p.order = make([]map[*tree.Edge][2]string, 1<<m_ij, 1<<m_ij)
// 		for i := 0; i < 1<<m_ij; i++ {
// 			p.order[i] = make(map[*tree.Edge][2]string)
// 		}
// 	}
// }

func (p *Pair) set_order(m map[*tree.Edge][2]string) {
	p.order = m
}

// func (p *Pair) initialise_t_ij_slice(m_ij int) {
// 	p.t_pi = make([][]float64, 1<<m_ij)
// 	for i := 0; i < 1<<m_ij; i++ {
// 		p.t_pi[i] = make([]float64, 0)
// 	}
// }

// the tree's branches are labelled with "syn=" for non-protein-changing nucleotide changes,
// and "AA=" for amino acid changes
// edge.SynLen is also set to the inferred synonymous branch length
func Epistasis(t *tree.Tree, features []annotation.Region, threads int) {
	// get_m_total(t)
	// os.Exit(0)
	// the mean synonymous distance between non-synonymous pairs:
	tau := get_tau(t)
	fmt.Println("tau is: " + strconv.FormatFloat(tau, 'f', 8, 64))

	aas_to_keep := get_aas_to_keep(t)
	fmt.Println("number of amino acids to analyse is: " + strconv.Itoa(len(aas_to_keep)))
	// fmt.Println()

	runtime.GOMAXPROCS(threads)

	cPair := make(chan Pair, threads)
	cResults := make(chan Result, threads)
	cResultsAgg := make(chan []Result)

	// cWriteDone := make(chan bool)

	go func() {
		for i := range aas_to_keep {
			for j := i + 1; j < len(aas_to_keep); j++ {
				aa1 := aas_to_keep[i]
				aa2 := aas_to_keep[j]
				// the original paper defined the pairs as ordered, so i -> j is a different proposition to j -> i.
				// hence, we have a separate struct for each order
				cPair <- Pair{i: aa1, j: aa2}
				cPair <- Pair{i: aa2, j: aa1}
			}
		}

		close(cPair)
	}()

	var wgPairs sync.WaitGroup
	wgPairs.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			process_pairs(cPair, cResults, t, tau)
			wgPairs.Done()
		}()
	}

	go aggregateResults(cResults, cResultsAgg)

	// go printResults(cResults, cWriteDone)

	wgPairs.Wait()
	close(cResults)

	// _ = <-cWriteDone

	results := <-cResultsAgg

	sort.SliceStable(results, func(i, j int) bool {
		return results[i].E_tau > results[j].E_tau
	})

	for i := range results {
		fmt.Println(results[i])
	}
}

// possible temporary function to summarise all the temporarily-unresolved mutation-branch combinations
func get_m_total(t *tree.Tree) {
	m := make(map[string]int)
	for _, e := range t.Edges() {
		// skip external branches
		if e.Right().Tip() {
			continue
		}
		AAs := e.Get_AA_residues()
		for i := range AAs {
			for j := i + 1; j < len(AAs); j++ {
				p := []string{AAs[i], AAs[j]}
				sort.Strings(p)
				s := strings.Join(p, "_")
				m[s]++
			}
		}
	}

	type simple struct {
		name  string
		value int
	}

	key_val_pair := make([]simple, 0)

	for k, v := range m {
		key_val_pair = append(key_val_pair, simple{k, v})
	}

	sort.Slice(key_val_pair, func(i, j int) bool {
		return key_val_pair[i].value > key_val_pair[j].value
	})

	for _, kvp := range key_val_pair {
		fmt.Println(kvp.name + " " + strconv.Itoa(kvp.value))
	}
}

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

// traverse the tree, and return a list of residues where changes occur at least once
func get_aas_to_keep(t *tree.Tree) []string {
	// m is just a map from residue to counts of it
	m := make(map[string]int)
	for _, e := range t.Edges() {
		// skip external branches
		if e.Right().Tip() {
			continue
		}
		// AAs
		AAs := e.Get_AA_residues()
		for _, AA := range AAs {
			if _, ok := m[AA]; ok {
				m[AA]++
			} else {
				m[AA] = 1
			}
		}
	}

	// here we only keep residues where mutations occur at least once on the tree (follows original paper)
	tokeep := make([]string, 0)
	for residue, count := range m {
		if count >= 1 {
			tokeep = append(tokeep, residue)
		}
	}

	return tokeep
}

func process_pairs(cPairIn chan Pair, cPairOut chan Result, t *tree.Tree, tau float64) {

	for p := range cPairIn {
		cPairOut <- process_one_pair(p, t, tau)
	}

}

func process_one_pair(p Pair, t *tree.Tree, tau float64) Result {

	i := p.get_i()
	j := p.get_j()
	m_ij, unresolved_edges := get_temporally_unresolved(t, &p)

	p.set_m_ij(m_ij)

	// if there are no unresolved edges for this pair, we can
	// recur down the tree without too many cares in the world
	if m_ij == 0 {
		ij_recur(nil, t.Root(), &p, false, 0.0)
		p.add_E_tau(tau)
		p.reset_t_pi()
	} else {
		// otherwise, let's do the order 2^^n stuff in a loop to
		// prevent an explosion of memory use
		for k := 0; k < 1<<m_ij; k++ {
			order := order(m_ij, k)
			om := make(map[*tree.Edge][2]string)
			for l := range unresolved_edges {
				switch order[l] {
				case 0:
					om[unresolved_edges[l]] = [2]string{i, j}
				case 1:
					om[unresolved_edges[l]] = [2]string{j, i}
				}
			}
			p.set_order(om)
			ij_recur_unresolved(nil, t.Root(), &p, false, 0.0)
			p.add_E_tau(tau)
			p.reset_t_pi()
		}
	}

	// get_set_temporally_unresolved(t, &p)
	// ij(t, &p)

	p.calc_mean_E_tau()

	return Result{i: p.i, j: p.j, E_tau: p.E_tau}
}

func aggregateResults(cResultsIn chan Result, cResultsOut chan []Result) {

	results := make([]Result, 0)

	for r := range cResultsIn {
		results = append(results, r)
	}

	cResultsOut <- results
}

func printResults(cResultsIn chan Result, cWriteDone chan bool) {

	fmt.Println("i,j,e_tau")

	for r := range cResultsIn {
		fmt.Println(r.i + "," + r.j + "," + strconv.FormatFloat(r.E_tau, 'f', 6, 64))
	}

	cWriteDone <- true
}

// // Traverse the tree, storing the information about the changes present on branches.
// // Filter the changes based on some criteria (not external branches, >= 1 occurence).
// // Make the pairs
// func get_pairs_from_tree(t *tree.Tree) []Pair {

// 	// m is just a map from residue to counts of it
// 	m := make(map[string]int)
// 	for _, e := range t.Edges() {
// 		// skip external branches
// 		if e.Right().Tip() {
// 			continue
// 		}
// 		// AAs
// 		AAs := e.Get_AA_residues()
// 		for _, AA := range AAs {
// 			if _, ok := m[AA]; ok {
// 				m[AA]++
// 			} else {
// 				m[AA] = 1
// 			}
// 		}
// 	}

// 	// here we only keep residues where mutations occur at least once on the tree (follows original paper)
// 	tokeep := make([]string, 0)
// 	for key, value := range m {
// 		if value >= 1 {
// 			tokeep = append(tokeep, key)
// 		}
// 	}

// 	sort.Strings(tokeep)

// 	pairs := make([]Pair, 0)
// 	for i := range tokeep {
// 		for j := i + 1; j < len(tokeep); j++ {
// 			aa1 := tokeep[i]
// 			aa2 := tokeep[j]
// 			// the original paper defined the pairs as ordered, so i -> j is a different proposition to j -> i.
// 			// hence, we have a separate struct for each order
// 			pairs = append(pairs, Pair{i: aa1, j: aa2})
// 			pairs = append(pairs, Pair{i: aa2, j: aa1})
// 		}
// 	}

// 	return pairs
// }

// /*
// modified from: https://github.com/mxschmitt/golang-combinations/blob/master/combinations.go

// MIT License

// Copyright (c) 2018 Max Schmitt

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// */
// func orders(n int) (subsets [][]int) {
// 	if n == 0 {
// 		return make([][]int, 0)
// 	}
// 	length := uint(n)

// 	for subsetBits := 0; subsetBits < (1 << length); subsetBits++ {
// 		subset := make([]int, n)

// 		for object := uint(0); object < length; object++ {
// 			// checks if object is contained in subset
// 			// by checking if bit 'object' is set in subsetBits
// 			if (subsetBits>>object)&1 == 1 {
// 				// add object to subset
// 				subset[object] = 1
// 			}
// 		}
// 		// add subset to subsets
// 		subsets = append(subsets, subset)
// 	}
// 	return subsets
// }

func order(n, m int) []int {
	if n == 0 {
		return make([]int, 0)
	}
	length := uint(n)
	subsetBits := uint(m)
	order := make([]int, n)

	for object := uint(0); object < length; object++ {
		if (subsetBits>>object)&1 == 1 {
			order[object] = 1
		}
	}
	return order
}

func get_temporally_unresolved(t *tree.Tree, p *Pair) (int, []*tree.Edge) {
	var i_present bool
	var j_present bool
	var site string

	i := p.get_i()
	j := p.get_j()

	m_ij := 0

	unresolved := make([]*tree.Edge, 0)

	for _, e := range t.Edges() {
		i_present = false
		j_present = false
		for _, comment := range e.GetComments() {
			if strings.HasPrefix(comment, "syn=") {
				continue
			}
			site = strings.Join(strings.Split(comment, ":")[0:2], ":")
			if site == i {
				i_present = true
			}
			if site == j {
				j_present = true
			}
			if i_present && j_present {
				m_ij++
				unresolved = append(unresolved, e)
				break
			}
		}
	}

	return m_ij, unresolved
}

// func get_set_temporally_unresolved(t *tree.Tree, p *Pair) {
// 	var i_present bool
// 	var j_present bool
// 	var site string

// 	i := p.get_i()
// 	j := p.get_j()

// 	m_ij := 0

// 	unresolved := make([]*tree.Edge, 0)

// 	for _, e := range t.Edges() {
// 		i_present = false
// 		j_present = false
// 		for _, comment := range e.GetComments() {
// 			if strings.HasPrefix(comment, "syn=") {
// 				continue
// 			}
// 			site = strings.Join(strings.Split(comment, ":")[0:2], ":")
// 			if site == i {
// 				i_present = true
// 			}
// 			if site == j {
// 				j_present = true
// 			}
// 			if i_present && j_present {
// 				m_ij++
// 				unresolved = append(unresolved, e)
// 				break
// 			}
// 		}
// 	}

// 	p.set_m_ij(m_ij)
// 	p.set_order_length(m_ij)
// 	p.initialise_t_ij_slice(m_ij)

// 	os := orders(m_ij)
// 	for n := range os {
// 		om := make(map[*tree.Edge][2]string)
// 		for m := range unresolved {
// 			switch os[n][m] {
// 			case 0:
// 				om[unresolved[m]] = [2]string{i, j}
// 			case 1:
// 				om[unresolved[m]] = [2]string{j, i}
// 			}

// 		}
// 		p.set_order(n, om)
// 	}
// }

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
		AAs := e.Get_AA_residues()
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

		AAs := e.Get_AA_residues()
		if len(AAs) > 0 {
			tns.add_distance(float64(n) * float64(len(AAs)) * (synDist + (e.SynLen / 2)))
			tns.add_weight(n * len(AAs))
		}

		d := synDist + e.SynLen

		// we carry on recurring down the tree
		collect_distances(e, curNode.Neigh()[i], tns, n, d)
	}
}

// func processChunktemporal(t *tree.Tree, pairs []Pair) {
// 	for i := range pairs {
// 		get_set_temporally_unresolved(t, &(pairs[i]))
// 	}
// }

// func processChunkij(t *tree.Tree, pairs []Pair) {
// 	for i := range pairs {
// 		ij(t, &(pairs[i]))
// 	}
// }

// {AA=orf1ab:4387 AA=orf1ab:6485 0 []}

// func ij(t *tree.Tree, p *Pair) {
// 	switch p.m_ij {
// 	case 0:
// 		ij_recur(nil, t.Root(), p, false, 0.0)
// 	default:
// 		for i := range p.order {
// 			ij_recur_unresolved(nil, t.Root(), p, i, false, 0.0)
// 		}
// 	}
// 	// fmt.Println(*p)
// }

func ij_recur(prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open bool, synDist float64) {

	// base state: if we have reached a tip, we don't want to do anything (per the original paper -
	// mutations on external edges are discounted)
	if curNode.Tip() {
		return
	}

	i := pair.get_i()
	j := pair.get_j()

	for k, e := range curNode.Edges() {

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
			if site == i {
				i_present = true
			}
			if site == j {
				j_present = true
			}
		}

		if i_present && j_present {
			panic("i and j both present but this is the function for temporally resolved pairs")
		}

		// d is the synonymous distance that will be passed to the next call of the function.
		// it will be reset if there are any relevant nonsynonymous mutations on this branch,
		// or it will be the current accrued synonymous distance + the edge synonymous distance
		// if there aren't.
		var d float64

		// o is a boolean array of length 2, [i is present (last) on the parent branch, j is present (last) on the parent branch]
		var o bool

		// we switch on whether there is an ancestral mutation at site i, and whether there are mutations on this branch at sites i or j
		switch {
		case !open && !i_present && !j_present: // 000
			o = false
		case !open && !i_present && j_present: // 001
			o = false
		case !open && i_present && !j_present: // 010
			o = true
			d = float64(edgeSynDist) / 2
		case open && !i_present && !j_present: // 100
			o = true
			d = synDist + float64(edgeSynDist)
		case open && !i_present && j_present: // 101
			pair.append_t_pi(synDist + (float64(edgeSynDist) / 2))
			o = false
			d = 0.0
		case open && i_present && !j_present: // 110
			o = true
			d = float64(edgeSynDist) / 2
		}

		// func ij_recur(prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open [2]bool, synDist float64) {
		// then we carry on recurring down the tree
		ij_recur(e, curNode.Neigh()[k], pair, o, d)
	}
}

func ij_recur_unresolved(prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open bool, synDist float64) {

	// base state: if we have reached a tip, we don't want to do anything (per the original paper -
	// mutations on external edges are discounted)
	if curNode.Tip() {
		return
	}

	i := pair.get_i()
	j := pair.get_j()

	for k, e := range curNode.Edges() {

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
			if site == i {
				i_present = true
			}
			if site == j {
				j_present = true
			}
		}

		// d is the synonymous distance that will be passed to the next call of the function.
		// it will be reset if there are any relevant nonsynonymous mutations on this branch,
		// or it will be the current accrued synonymous distance + the edge synonymous distance
		// if there aren't.
		var d float64

		// o is a boolean array of length 2, [i is present (last) on the parent branch, j is present (last) on the parent branch]
		var o bool

		// we switch on whether there is an ancestral mutation at site i, and whether there are mutations on this branch at sites i or j for this order
		switch {
		case !open && !i_present && !j_present: // 000
			o = false
		case !open && !i_present && j_present: // 001
			o = false
		case !open && i_present && !j_present: // 010
			o = true
			d = float64(edgeSynDist) / 2
		case open && !i_present && !j_present: // 100
			o = true
			d = synDist + float64(edgeSynDist)
		case open && !i_present && j_present: // 101
			pair.append_t_pi(synDist + (float64(edgeSynDist) / 2))
			o = false
			d = 0.0
		case open && i_present && !j_present: // 110
			o = true
			d = float64(edgeSynDist) / 2

		case !open && i_present && j_present: // 011

			first := pair.get_order()[e][0]
			// second := pair.get_order(order)[e][0]

			if first == i { // xij
				pair.append_t_pi((float64(edgeSynDist) / 3))
				o = false
				d = 0.0
			} else { // xji
				o = true
				d = float64(edgeSynDist) / 3
			}
		case open && i_present && j_present: // 111

			first := pair.get_order()[e][0]
			// second := pair.get_order(order)[e][0]

			if first == i { // iij
				pair.append_t_pi((float64(edgeSynDist) / 3))
				o = false
				d = 0.0
			} else { // iji
				pair.append_t_pi(synDist + (float64(edgeSynDist) / 3))
				o = true
				d = float64(edgeSynDist) / 3
			}
		}

		// func ij_recur(prevEdge *tree.Edge, curNode *tree.Node, pair *Pair, open [2]bool, synDist float64) {
		// then we carry on recurring down the tree
		ij_recur_unresolved(e, curNode.Neigh()[k], pair, o, d)
	}
}
