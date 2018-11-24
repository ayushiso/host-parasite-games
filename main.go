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

//NewPopnlist takes initial host and parasite population and updates it numGens times
//by simulating host-parasite interactions between among all members of the population.
func NewPopnlist(inithpop, initppop Population, numGens int, hostrep, parrep float64) Cell {
	poplist := make([]Populations, numGens+1)

	poplist[0][0], poplist[0][1] = inithpop, initppop

	//starting parameters from args
	fmt.Println("Starting with host count: ", len(inithpop))
	fmt.Println("Starting with host killer frequency: ", CountKillFreq(inithpop))
	fmt.Println("Starting with parasite count: ", len(initppop))
	fmt.Println("Starting with parasite killer frequency: ", CountKillFreq(initppop))

	//Update all populations, and add to list
	for j := 1; j <= numGens; j++ {
		fmt.Println("entering loop...")

		//population size limit due to RAM constraints
		if len(poplist[j-1][0]) >= 10000000 || len(poplist[j-1][1]) >= 10000000 {
			fmt.Println("Populations out of control!")
			break
		}

		var newhpop Population
		var newppop Population

		//Making pairs, depending on which population is more
		if len(poplist[j-1][0]) > len(poplist[j-1][1]) { //hostpopn more
			newhpop, newppop = HMakePairs(poplist[j-1][0], poplist[j-1][1])
		} else { //both equal, or parpopn more
			newhpop, newppop = PMakePairs(poplist[j-1][0], poplist[j-1][1])
		}

		//update populations by simulating one round of killer-diplomat game
		uphpop, uppop := UpdatePopns(newhpop, newppop)

		//replicate populations as per respective replication rate
		poplist[j][0], poplist[j][1] = ReplPopns(uphpop, uppop, hostrep, parrep)

		fmt.Println("At step", j, "host count is", len(poplist[j][0]))
		fmt.Println("At step", j, "host killer percentage is", CountKillFreq(poplist[j][0]))
		fmt.Println("At step", j, "parasite count is", len(poplist[j][1]))
		fmt.Println("At step", j, "parasite killer percentage is", CountKillFreq(poplist[j][1]))

		if len(poplist[j][0]) == 0 { //all hosts are dead
			fmt.Println("Oops! Everyone is dead.")
			break
		}

		if len(poplist[j][1]) == 0 { //all parasites are dead
			fmt.Println("Nice! We eradicated the virus.")
			break
		}
	}
	return poplist[0][0][0]
}

//CountKillFreq takes a population and returns the proportion of Cells in the
//population with "killer" strategy
func CountKillFreq(popn Population) float64 {
	kcount := float64(CountKill(popn))
	return kcount / float64(len(popn))
}

//CountKill takes a population and returns the number of killer cells in that population
func CountKill(popn Population) int {
	killers := 0
	for i := range popn {
		if popn[i].strategy == "k" {
			killers++
		}
	}
	return killers
}

//ReplPopns takes hostpop, parpop and returns the new population according to the
//reproduction rate of the host and parasite.
func ReplPopns(hostpop, parpop Population, hrep, prep float64) (Population, Population) {

	newhostpop := make(Population, int(float64(len(hostpop))*hrep))
	newparpop := make(Population, int(float64(len(parpop))*prep))

	hkillfreq := CountKillFreq(hostpop)
	pkillfreq := CountKillFreq(parpop)

	//assign strategies in new population as per parent population

	for k := range newhostpop { //for all host cells
		newhostpop[k].status = "alive"
		strategy := rand.Float64()
		if strategy <= hkillfreq {
			newhostpop[k].strategy = "k"

		} else {
			newhostpop[k].strategy = "d"
		}

	}
	for l := range newparpop { //for all parasite cells
		newparpop[l].status = "alive"
		strategy := rand.Float64()
		if strategy <= pkillfreq {
			newparpop[l].strategy = "k"
		} else {
			newparpop[l].strategy = "d"
		}
	}
	return newhostpop, newparpop
}

//UpdatePopns takes hostpop, parpop and updates the individuals in the populations
//as per the results of a pairwise host-parasite competition.
func UpdatePopns(hostpop, parpop Population) (Population, Population) {
	var newhostpop Population
	var newparpop Population

	//match all cells in one population to cells in the other population
	if len(hostpop) > len(parpop) {
		for i := range parpop {
			parpop[i], *parpop[i].partner = UpdateStatus(parpop[i], *parpop[i].partner)
		}
	} else {
		for i := range hostpop {
			hostpop[i], *hostpop[i].partner = UpdateStatus(hostpop[i], *hostpop[i].partner)
		}
	}

	//add to populations if status is "alive" after matching
	for i := range hostpop { //for all host cells
		if hostpop[i].status == "alive" {
			newhostpop = append(newhostpop, hostpop[i])
		}
	}
	for i := range parpop { //for all parasite cells
		if parpop[i].status == "alive" {
			newparpop = append(newparpop, parpop[i])

		}
	}
	return newhostpop, newparpop
}

//UpdateStatus takes a pair of cells, one host and one parasite, and returns the
//status of each cell after a simulated interaction.
func UpdateStatus(host, parasite Cell) (Cell, Cell) {

	//bottom row of payoff matrix
	if host.strategy == "d" {
		if parasite.strategy == "k" {
			host.status = "dead"
			parasite.status = "alive"
		} else if parasite.strategy == "d" {
			host.status = "alive"
			parasite.status = "alive"
		}

		//top row of payoff matrix
	} else if host.strategy == "k" {
		if parasite.strategy == "k" {
			probhostwin := rand.Float64()
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

//MakePairs takes two Populations (host and parasite) and matches each host
//Cell to a parasite Cell (hostnum > parnum)
func HMakePairs(hpop, ppop Population) (Population, Population) {
	if len(ppop) > len(hpop) {
		panic("Parasite popn more: wrong subroutine")
	}

	parindices := make([]int, 0) //array for storing parasite indices
	for i := range ppop {
		parindices = append(parindices, i)
	}

	for j := range hpop {
		if len(parindices) == 0 {
			return hpop, ppop
		}
		index := rand.Intn(len(parindices))
		match := parindices[index]

		hpop[j].partner = &ppop[match]
		ppop[match].partner = &hpop[j]
		parindices = Delete(parindices, index)
	}
	return hpop, ppop
}

//MakePairs takes two Populations (host and parasite) and matches each parasite
//Cell to a host Cell (hostnum < parnum)
func PMakePairs(hpop, ppop Population) (Population, Population) {
	if len(hpop) > len(ppop) {
		panic("Host popn more: Wrong subroutine")
	}
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

		ppop[j].partner = &hpop[match]
		hpop[match].partner = &ppop[j]
		hostindices = Delete(hostindices, index)
	}
	return hpop, ppop
}

func Delete(a []int, item int) []int {
	a = append(a[:item], a[item+1:]...)
	return a
}

//*********************************************************************
//remove randomness in frequency assignment:
//*********************************************************************
//Initialize takes numbers of hosts and parasites and the percentage of killers
//in each, and creates host and parasite populations accordingly.
func Initialize(hostnum, parnum int, hkillers, pkillers float64) (Population, Population) {
	hpop := make(Population, hostnum)
	ppop := make(Population, parnum)

	hknum := int(float64(hostnum) * hkillers)
	pknum := int(float64(parnum) * pkillers)

	for i := range hpop {
		hpop[i].status = "alive"
	}
	for j := 0; j < hknum; j++ {
		hpop[j].strategy = "k"
	}
	for j := hknum; j < len(hpop); j++ {
		hpop[j].strategy = "d"
	}

	for i := range ppop {
		ppop[i].status = "alive"
	}
	for j := 0; j < pknum; j++ {
		ppop[j].strategy = "k"
	}
	for j := pknum; j < len(ppop); j++ {
		ppop[j].strategy = "d"
	}

	return hpop, ppop
}

func main() {

	rand.Seed(time.Now().UnixNano())

	hostnum := 100
	parnum := 100
	hkillers := 0.3
	pkillers := 0.4
	//hostrate := 1.3
	//parrate := 1.4
	hpop, ppop := Initialize(hostnum, parnum, hkillers, pkillers)
	reppop := DynPopnList(hpop, ppop, 100)
	fmt.Println(reppop)
	/*
		nhpop, nppop := DynUpdatePopn(hpop, ppop)
		fmt.Println(hpop)
		fmt.Println(nhpop)
		fmt.Println()
		fmt.Println(ppop)
		fmt.Println(nppop)
	*/
	/* TESTING FOR PART1
	mypop := NewPopnlist(hpop, ppop, 100, hostrate, parrate)
	fmt.Println(mypop)
	*/
}

/*
********************************************************************************
REPLICATOR DYNAMIC PROGRAM BEGINS HERE
********************************************************************************
*/

//DynPopnList takes host and parasite populations and updates the populations
//numGens times, and returns a list of updates populations.
func DynPopnList(inithpop, initppop Population, numGens int) Cell {
	poplist := make([]Populations, numGens+1)

	poplist[0][0], poplist[0][1] = inithpop, initppop

	//starting parameters from args
	fmt.Println("Starting with host count: ", len(inithpop))
	fmt.Println("Starting with host killer frequency: ", CountKillFreq(inithpop))
	fmt.Println("Starting with parasite count: ", len(initppop))
	fmt.Println("Starting with parasite killer frequency: ", CountKillFreq(initppop))

	for j := 1; j <= numGens; j++ {
		fmt.Println("entering loop...")

		//population size limit due to RAM constraints
		if len(poplist[j-1][0]) >= 10000000 || len(poplist[j-1][1]) >= 10000000 {
			fmt.Println("Populations out of control!")
			break
		}
		//update populations using replicator dynamic
		poplist[j][0], poplist[j][1] = DynUpdatePopn(poplist[j-1][0], poplist[j-1][1])

		fmt.Println("At step", j, "host count is", len(poplist[j][0]))
		fmt.Println("At step", j, "host killer percentage is", CountKillFreq(poplist[j][0]))
		fmt.Println("At step", j, "parasite count is", len(poplist[j][1]))
		fmt.Println("At step", j, "parasite killer percentage is", CountKillFreq(poplist[j][1]))

		if len(poplist[j][0]) == 0 { //all hosts are dead
			fmt.Println("Oops! Everyone is dead.")
			break
		}

		if len(poplist[j][1]) == 0 { //all parasites are dead
			fmt.Println("Nice! We eradicated the virus.")
			break
		}
	}
	return poplist[0][0][0]
}

//DynUpdatePopn takes host and parasite populations
//and returns an updated population in accordance with the replicator dynamic.
func DynUpdatePopn(hostpop, parpop Population) (Population, Population) {
	//finding frequencies in initial host population
	killhfreq := CountKillFreq(hostpop)
	diphfreq := 1.0 - killhfreq

	//finding frequencies in initial parasite population
	killpfreq := CountKillFreq(parpop)
	dippfreq := 1.0 - killpfreq

	//computing fitness values for killer and diplomat hosts
	fitness_hk := killpfreq*0.6 - killpfreq*0.4 + dippfreq*1.0
	fitness_hd := killpfreq*(-1.0) + dippfreq*1.0

	//computing fitness values for killer and diplomat parasites
	fitness_pk := killhfreq*0.4 - killhfreq*0.6 + diphfreq*1.0
	fitness_pd := killhfreq*(-1.0) + diphfreq*1.0

	//computing population fitness for host and parasite popns
	host_fitness := killhfreq*fitness_hk + diphfreq*fitness_hd
	par_fitness := killpfreq*fitness_pk + dippfreq*fitness_pd

	//computing rate of increase for killer and diplomat hosts
	hkrate := killhfreq * (fitness_hk - host_fitness)
	hdrate := diphfreq * (fitness_hd - host_fitness)

	//computing rate of increase for killer and diplomat parasites
	pkrate := killpfreq * (fitness_pk - par_fitness)
	pdrate := dippfreq * (fitness_pd - par_fitness)

	//making and adding two new host populations corresponding to killers and diplomats
	khosts := CountKill(hostpop)
	dhosts := len(hostpop) - khosts
	khnum := Replicate(khosts, hkrate)
	dhnum := Replicate(dhosts, hdrate)

	var newkhosts Population
	var newdhosts Population

	if khnum > 0 {
		newkhosts = make(Population, khnum)
	} else {
		newkhosts = make(Population, 0)
	}
	if dhnum > 0 {
		newdhosts = make(Population, dhnum)
	} else {
		newdhosts = make(Population, 0)
	}

	for h := range newkhosts {
		newkhosts[h].strategy = "k"
	}
	for h := range newdhosts {
		newdhosts[h].strategy = "d"
	}
	newdhosts = append(newdhosts, newkhosts...)

	//making and adding two new parasite populations corresponding to killers and diplomats
	kpars := CountKill(parpop)
	dpars := len(parpop) - khosts
	kpnum := Replicate(kpars, pkrate)
	dpnum := Replicate(dpars, pdrate)

	var newkpars Population
	var newdpars Population
	if kpnum > 0 {
		newkpars = make(Population, kpnum)
	} else {
		newkpars = make(Population, 0)
	}
	if dpnum > 0 {
		newdpars = make(Population, dpnum)
	} else {
		newdpars = make(Population, 0)
	}

	for p := range newkpars {
		newkpars[p].strategy = "k"
	}
	for p := range newdpars {
		newdpars[p].strategy = "d"
	}
	newdpars = append(newdpars, newkpars...)

	return newdhosts, newdpars
}

//Replicate takes a population size and returns the new population size based on
//replication rate.
func Replicate(popsize int, rate float64) int {
	if rate < 0 { //population size decreasing
		dec := (-1.0) * rate
		if int(float64(popsize)*dec) > popsize {
			return 0
		}
		return popsize - int(float64(popsize)*dec)
	} else {
		return popsize + int(float64(popsize)*rate)
	}
}
