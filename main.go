package main

import (
	"fmt"
	"math/rand"
	"time"
)

type Populations [2]Population

type Population []Cell

type Cell struct {
	strategy string
	status   string
	partner  *Cell
}

//FreqAnalysis takes a list of host and parasite populations, and returns the
//population characteristics at each step: number of individuals in host and
//parasite, and frequency of killers and diplomats
func FreqAnalysis(poplist []Populations) string {
	for i := range poplist { //each step of simulation
		if len(poplist[i][0]) == 0 {
			return "Oops! Everyone is dead."
		}
		if len(poplist[i][1]) == 0 {
			return "Yay! We eradicated the virus."
		}
		khcount := 0
		kpcount := 0
		fmt.Println("At step", i, "host count is ", len(poplist[i][0]))
		for j := range poplist[i][0] { //host killer percentage
			//fmt.Println(poplist[i][0][j].strategy)
			if poplist[i][0][j].strategy == "k" {
				khcount++
			}
		}
		fmt.Println("At step", i, "host killer percentage is", float64(khcount)/float64(len(poplist[i][0])))
		fmt.Println("At step", i, "parasite count is ", len(poplist[i][1]))
		for k := range poplist[i][1] { //parasite killer percentage
			if poplist[i][1][k].strategy == "k" {
				kpcount++
			}
		}
		fmt.Println("At step", i, "parasite killer percentage is", float64(kpcount)/float64(len(poplist[i][1])))
		fmt.Println()
	}
	return "Done! :)"
}

//Popnlist takes initial host and parasite population and updates it numGens times
//by simulating host-parasite interactions between among all members of the population.
func Popnlist(inithpop, initppop Population, numGens int) []Populations {
	poplist := make([]Populations, numGens)
	poplist[0][0], poplist[0][1] = inithpop, initppop
	for j := 1; j < numGens; j++ {
		//fmt.Println("Previous host population: ", poplist[j-1][0])
		newhpop, newppop := MakePairs(poplist[j-1][0], poplist[j-1][1])
		//fmt.Println("Made pairs: ", newhpop)
		poplist[j][0], poplist[j][1] = UpdatePopns(newhpop, newppop)
		//fmt.Println("Updated host population: ", poplist[j][0])
		//fmt.Println()
	}
	return poplist
}

//UpdatePopns takes hostpop, parpop and updates the individuals in the populations
//as per the results of a pairwise host-parasite competition.
func UpdatePopns(hostpop, parpop Population) (Population, Population) {
	var newhostpop Population
	var newparpop Population
	for i := range hostpop {
		hostpop[i], *hostpop[i].partner = UpdateStatus(hostpop[i], *hostpop[i].partner)
	}
	for i := range hostpop {
		if hostpop[i].status == "alive" {
			newhostpop = append(newhostpop, hostpop[i])
		}
	}
	for i := range parpop {
		if parpop[i].status == "alive" {
			newparpop = append(newparpop, parpop[i])
			//fmt.Println(parpop[i])
		}
	}
	return newhostpop, newparpop
}

//UpdateStatus takes a pair of cells, one host and one parasite, and returns the
//status of each cell after a simulated interaction.
func UpdateStatus(host, parasite Cell) (Cell, Cell) {
	if host.strategy == "d" {
		if parasite.strategy == "k" {
			host.status = "dead"
			parasite.status = "alive"
		} else if parasite.strategy == "d" {
			host.status = "alive"
			parasite.status = "alive"
		}
	} else if host.strategy == "k" {
		if parasite.strategy == "k" {
			probhostwin := rand.Float64()
			//fmt.Println(probhostwin)
			if probhostwin <= 0.6 {
				host.status = "alive"
				parasite.status = "dead"
			} else {
				host.status = "dead"
				parasite.status = "alive"
			}
		} else if parasite.strategy == "d" {
			host.status = "alive"
			parasite.status = "dead"
		}
	}
	return host, parasite
}

//MakePairs takes two Populations (host and parasite) and matches each parasite
//Cell to a host Cell
func MakePairs(hpop, ppop Population) (Population, Population) {
	hostindices := make([]int, 0) //array for storing host indices
	for i := range hpop {
		hostindices = append(hostindices, i)
	}
	for j := range ppop {
		if len(hostindices) == 0 {
			return hpop, ppop
		}
		index := rand.Intn(len(hostindices))
		match := hostindices[index]
		//fmt.Println(hpop[match])
		ppop[j].partner = &hpop[match]
		hpop[match].partner = &ppop[j]
		hostindices = Delete(hostindices, index)
		//fmt.Println(hostindices)
	}
	return hpop, ppop
}

func Delete(a []int, item int) []int {
	a = append(a[:item], a[item+1:]...)
	return a
}

//Initialize takes numbers of hosts and parasites and the percentage of killers
//in each, and creates host and parasite populations accordingly.
func Initialize(hostnum, parnum int, hkillers, pkillers float64) (Population, Population) {
	hpop := make(Population, hostnum)
	ppop := make(Population, parnum)
	for i := range hpop { //initializing host cells
		hpop[i].status = "alive"
		strategy := rand.Float64()
		if strategy <= hkillers {
			hpop[i].strategy = "k"
			//fmt.Println(hcell.strategy)
		} else {
			hpop[i].strategy = "d"
		}
		//fmt.Println(hpop)
	}
	for j := range ppop {
		ppop[j].status = "alive"
		strategy := rand.Float64()
		if strategy <= pkillers {
			ppop[j].strategy = "k"
		} else {
			ppop[j].strategy = "d"
		}
		//fmt.Println(pcell)
	}
	return hpop, ppop
}

func main() {
	rand.Seed(time.Now().UnixNano())
	hostnum := 500
	parnum := 1000
	hkillers := 0.7
	pkillers := 0.4
	hpop, ppop := Initialize(hostnum, parnum, hkillers, pkillers)
	mypop := Popnlist(hpop, ppop, 100)
	fmt.Println(FreqAnalysis(mypop))
	//fmt.Println(hpop)
}
