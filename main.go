package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"strconv"
	"time"
)

type Populations [2]Population

type Population []Cell

type Cell struct {
	strategy string
	status   string
	partner  *Cell
}

func main() {

	//variables common to both algorithms
	var numGens int
	var hostnum, parnum, hkillers, pkillers, hostrate, parrate, mortality float64

	//variables for replicator algorithm
	var vircost, rescost, gifthost, giftpar float64

	rand.Seed(time.Now().UnixNano())

	//taking inputs from the user
	reader := bufio.NewReader(os.Stdin)

	fmt.Println("\tWelcome to the host-parasite game theory simulator!")
	fmt.Println("\tChoose a mode, or [e]xit: \n\t[c]ompare | [a]nalyze | [s]andbox")

	text, _, err := reader.ReadRune()
	if err != nil {
		fmt.Println("There was some kind of error during input, please try again.")
	}

	if text == 'c' { //comparison stage

		fmt.Println("\tComparing Matching and Replicator algorithms on the same input.")

		//getting user input for all parameters
		//since no variable is excluded, -1 is a placeholder for "variable"
		GetParameters(-1, &hostnum, &parnum, &hkillers, &pkillers, &hostrate, &parrate, &mortality, &vircost, &rescost, &gifthost, &giftpar)

		//getting number of generations
		fmt.Print("Number of generations:")
		_, err1 := fmt.Scan(&numGens)
		if err1 != nil {
			fmt.Println(err1)
		}

		//map of all parameters used in Replicator algorithm
		params := map[string]float64{
			"hostrate":  hostrate,
			"parrate":   parrate,
			"mortality": mortality,
			"vircost":   vircost,
			"rescost":   rescost,
			"gifthost":  gifthost,
			"giftpar":   giftpar,
		}

		hpop, ppop := Initialize(int(hostnum), int(parnum), hkillers, pkillers)
		hpop1, ppop1 := Initialize(int(hostnum), int(parnum), hkillers, pkillers)

		//evolving population using both algorithms
		matchpop := MatchPopnlist(hpop, ppop, numGens, hostrate, parrate, mortality)
		reppop := DynPopnList(hpop1, ppop1, numGens, params)

		//writing to separate files
		WritePopns("matching.txt", "Evolution by matching algorithm", matchpop)
		WritePopns("replicator.txt", "Evolution by replicator algorithm", reppop)

	} else if text == 'a' { //analysis stage

		//getting number ofgenerations first
		fmt.Print("Number of generations:")
		_, err1 := fmt.Scan(&numGens)
		if err1 != nil {
			fmt.Println(err1)
		}

		var choice int
		var n int

		//getting choice of parameter to be analyzed
		fmt.Println("\tWhich parameter are you studying?")
		fmt.Println("\tChoose a parameter:")
		fmt.Println("\t[0]host population | [1]parasite population | [2]host killer frequency | [3]parasite killer frequency | [4]host replication rate | [5]parasite replication rate | [6]parasite-induced mortality |\n\t[7]cost of virulence | [8]cost of resistance | [9]host energy loss | [10]parasite energy loss")

		_, err2 := fmt.Scan(&choice)
		if err2 != nil {
			fmt.Println(err2)
		}

		fmt.Print("\tNumber of parameters(maximum 10, minimum 2):")
		_, err3 := fmt.Scan(&n)
		if err3 != nil {
			fmt.Println(err3)
		}
		if n > 10 || n < 2 {
			panic("Wrong number of parameters! Try again")
		}

		//making a list of the parameter values being analyzed
		pars := GetValues(n)

		//get user input for all other parameters
		GetParameters(choice, &hostnum, &parnum, &hkillers, &pkillers, &hostrate, &parrate, &mortality, &vircost, &rescost, &gifthost, &giftpar)

		//map of all parameters used in Replicator algorithm
		params := map[string]float64{
			"hostrate":  hostrate,
			"parrate":   parrate,
			"mortality": mortality,
			"vircost":   vircost,
			"rescost":   rescost,
			"gifthost":  gifthost,
			"giftpar":   giftpar,
		}

		//generates different populations depending on the variable being studied
		Analyze(choice, numGens, pars, params, hostnum, parnum, hkillers, pkillers)

	} else if text == 's' { //sandbox (simple simulation mode)
		fmt.Println("\tTinkering with the replicator dynamic! Enter your own values:")

		//getting parameters from user
		//Since all parameters are being analyzed, -1 is a placeholder variable
		GetParameters(-1, &hostnum, &parnum, &hkillers, &pkillers, &hostrate, &parrate, &mortality, &vircost, &rescost, &gifthost, &giftpar)

		//getting number of generations
		fmt.Print("Number of generations:")
		_, err1 := fmt.Scan(&numGens)
		if err1 != nil {
			fmt.Println(err1)
		}

		//map of all parameters used in Replicator algorithm
		params := map[string]float64{
			"hostrate":  hostrate,
			"parrate":   parrate,
			"mortality": mortality,
			"vircost":   vircost,
			"rescost":   rescost,
			"gifthost":  gifthost,
			"giftpar":   giftpar,
		}

		//generate population with Replicator algorithm, and write to file
		hpop, ppop := Initialize(int(hostnum), int(parnum), hkillers, pkillers)
		reppop := DynPopnList(hpop, ppop, numGens, params)
		WritePopns("sandbox.txt", "Evolution of host and parasite populations", reppop)
	}
}

//Analyze takes the index of the variable under study, along with the values given
//by the user, and generates the corresponding populations and writes them to file.
func Analyze(variable, numGens int, pars []float64, params map[string]float64, hostnum, parnum, hkillers, pkillers float64) {

	//list of variables which can be analyzed, used for writing to file
	vals := []string{"Host population", "Parasite population", "Host killfreq", "Parasite killfreq", "Host reprate", "Parasite reprate", "Host mortality", "Virulence cost", "Resistance cost", "Host energy lost", "Parasite energy lost"}

	//default values of populations from args
	hpop, ppop := Initialize(int(hostnum), int(parnum), hkillers, pkillers)

	//initialize populations or parameters differently depending on variable under study
	for i := range pars {
		switch variable {
		case 0:
			hpop, ppop = Initialize(int(pars[i]), int(parnum), hkillers, pkillers)
		case 1:
			hpop, ppop = Initialize(int(hostnum), int(pars[i]), hkillers, pkillers)
		case 2:
			hpop, ppop = Initialize(int(hostnum), int(parnum), pars[i], pkillers)
		case 3:
			hpop, ppop = Initialize(int(hostnum), int(parnum), hkillers, pars[i])
		case 4:
			params["hostrate"] = pars[i]
		case 5:
			params["parrate"] = pars[i]
		case 6:
			params["mortality"] = pars[i]
		case 7:
			params["vircost"] = pars[i]
		case 8:
			params["rescost"] = pars[i]
		case 9:
			params["gifthost"] = pars[i]
		case 10:
			params["giftpar"] = pars[i]
		default:
			panic("Wrong variable!")
		}

		//generate the population corresponding to this parameter value, and write to file
		popnlist := DynPopnList(hpop, ppop, numGens, params)

		//filename format for parsing by graphing.py
		name := "parameter" + strconv.Itoa(i+1) + ".txt"
		//header for plot in graphing.py
		header := vals[variable] + " = " + fmt.Sprint(pars[i])
		WritePopns(name, header, popnlist)
	}
}

//GetParameters takes all the user input for different variables (except the value
//chosen for analysis, if any)
func GetParameters(variable int, hostnum, parnum, hkillers, pkillers, hostrate, parrate, mortality, vircost, rescost, gifthost, giftpar *float64) {
	//inputs stores the user prompts
	inputs := []string{"Host number:", "Parasite number:", "Killer host percentage:", "Killer parasite percentage:", "Host replication rate:", "Parasite replication rate:", "Parasite-induced host mortality:", "Cost of virulence:", "Cost of resistance:", "Host energy loss:", "Parasite energy loss:"}
	//vars stores pointers to the variables which will be assigned values as per user input
	vars := []*float64{hostnum, parnum, hkillers, pkillers, hostrate, parrate, mortality, vircost, rescost, gifthost, giftpar}

	for char := range inputs {
		if char != variable { //gets only those parameters which are not the variable being analyzed
			fmt.Print(inputs[char])        //ask for value
			_, err := fmt.Scan(vars[char]) //assign to variable
			if err != nil {
				fmt.Println(err)
			}
		}
	}
}

//GetValues gets user values for the parameter being investigated
func GetValues(n int) []float64 {

	//getting a list of parameters
	pars := make([]float64, 0)
	for i := 1; i <= n; i++ {
		var value float64
		fmt.Print("Parameter", i, ":")
		fmt.Scan(&value)
		pars = append(pars, value)
	}
	return pars
}

/*
********************************************************************************
MATCHING ALGORITHM
********************************************************************************
*/

//MatchPopnlist takes initial host and parasite Population and updates it numGens times
//using the matching algorithm, and returns a list of evolving Populations
func MatchPopnlist(inithpop, initppop Population, numGens int, hostrep, parrep, mortality float64) []Populations {
	poplist := make([]Populations, 1)

	//first element in population is the tuple of initial host and parasite populations
	poplist[0][0], poplist[0][1] = inithpop, initppop

	//Update all populations, and add to list
	for j := 1; j <= numGens; j++ {

		//population size limit due to RAM constraints
		if len(poplist[j-1][0])+len(poplist[j-1][1]) >= 10000000 {
			fmt.Println("Total population size reached 10^7")
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
		uphpop, uppop := UpdatePopns(newhpop, newppop, mortality)

		//replicate populations as per respective replication rate
		var newpops Populations
		newpops[0], newpops[1] = ReplPopns(uphpop, uppop, hostrep, parrep)
		poplist = append(poplist, newpops)

		if len(poplist[j][0]) == 0 { //all hosts are dead
			fmt.Println("Oops! Every host is dead.")
			break
		}

		if len(poplist[j][1]) == 0 { //all parasites are dead
			fmt.Println("Nice! We eradicated the parasite.")
			break
		}
	}
	return poplist
}

//ReplPopns takes hostpop, parpop and returns the new population according to the
//reproduction rate of the host and parasite.
func ReplPopns(hostpop, parpop Population, hrep, prep float64) (Population, Population) {

	newhostpop := make(Population, int(float64(len(hostpop))*hrep))
	newparpop := make(Population, int(float64(len(parpop))*prep))

	hkillfreq := CountKillFreq(hostpop)
	pkillfreq := CountKillFreq(parpop)

	//assign strategies in new population as per parent population (randomized)

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
func UpdatePopns(hostpop, parpop Population, mortality float64) (Population, Population) {
	var newhostpop Population
	var newparpop Population

	//match all cells in one population to cells in the other population
	if len(hostpop) > len(parpop) {
		for i := range parpop {
			parpop[i], *parpop[i].partner = UpdateStatus(parpop[i], *parpop[i].partner, mortality)
		}
	} else {
		for i := range hostpop {
			hostpop[i], *hostpop[i].partner = UpdateStatus(hostpop[i], *hostpop[i].partner, mortality)
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
func UpdateStatus(host, parasite Cell, mortality float64) (Cell, Cell) {

	//in a K-K competition, probability of host winning is 1- mortality rate for the
	//prticular disease
	hostwin := 1.0 - mortality

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
			//introducing randomness in choosing whether host or parasite wins
			probhostwin := rand.Float64()

			if probhostwin <= hostwin {
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

/*
*******************************************************************************
HELPER FUNCTIONS
*******************************************************************************
*/

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

//Deletes an element from a list of ints (basic delete function covered in class)
func Delete(a []int, item int) []int {
	a = append(a[:item], a[item+1:]...)
	return a
}

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

	for i := range ppop {
		ppop[i].status = "alive"
	}

	//assign host killer strategies
	if hknum != 0 {
		for j := 0; j < hknum; j++ {
			hpop[j].strategy = "k"
		}
	}
	//assign host diplomat strategies
	if hknum != hostnum {
		for j := hknum; j < len(hpop); j++ {
			hpop[j].strategy = "d"
		}
	}

	//assign parasite killer strategies
	if pknum != 0 {
		for j := 0; j < pknum; j++ {
			ppop[j].strategy = "k"
		}
	}
	//assign host diplomat strategies
	if pknum != parnum {
		for j := pknum; j < len(ppop); j++ {
			ppop[j].strategy = "d"
		}
	}

	return hpop, ppop
}

/*
********************************************************************************
REPLICATOR ALGORITHM
********************************************************************************
*/

//DynPopnList takes host and parasite populations and updates the populations
//numGens times using Replicator algorithm, and returns a list of updated populations.
func DynPopnList(inithpop, initppop Population, numGens int, params map[string]float64) []Populations {
	poplist := make([]Populations, 1)

	//first element in population is the tuple of initial host and parasite populations
	poplist[0][0], poplist[0][1] = inithpop, initppop

	//computing Wh and Wp values for each cell of payoff matrix
	payoffs := ComputePayoffs(params)

	for j := 1; j <= numGens; j++ {

		//population size limit due to RAM constraints
		if len(poplist[j-1][0])+len(poplist[j-1][1]) >= 1000000 {
			fmt.Println("Populations out of control!")
			break
		}

		var newpops Populations
		//update once, form a new Population tuple, add to growing list
		newpops[0], newpops[1] = DynUpdatePopn(poplist[j-1][0], poplist[j-1][1], payoffs)
		poplist = append(poplist, newpops)

		if len(poplist[j][0]) == 0 { //all hosts are dead
			fmt.Println("Oops! Every host is dead.")
			break
		}

		if len(poplist[j][1]) == 0 { //all parasites are dead
			fmt.Println("Nice! We eradicated the parasite.")
			break
		}
	}
	return poplist
}

//DynUpdatePopn takes host and parasite populations
//and returns an updated population in accordance with the replicator dynamic.
func DynUpdatePopn(hostpop, parpop Population, payoffs map[string]float64) (Population, Population) {

	//compute the replication rate for each strategy subpopulation using fitness equations
	rates := ComputeFitness(hostpop, parpop, payoffs)

	//making two new host populations corresponding to killers and diplomats
	khosts := CountKill(hostpop) //getting existing size
	dhosts := len(hostpop) - khosts
	khnum := Replicate(khosts, rates["hkrate"]) //getting new size
	dhnum := Replicate(dhosts, rates["hdrate"])

	//calculating killer frequency in new host population
	totalh := khnum + dhnum
	hkfreq := float64(khnum) / float64(totalh)

	//making two new parasite populations corresponding to killers and diplomats
	kpars := CountKill(parpop) //getting existing size
	dpars := len(parpop) - kpars
	kpnum := Replicate(kpars, rates["pkrate"]) //getting new size
	dpnum := Replicate(dpars, rates["pdrate"])

	//calculating killer frequency in new parasite population
	totalp := kpnum + dpnum
	pkfreq := float64(kpnum) / float64(totalp)

	//making the actual Populations
	newhpop, newppop := Initialize(totalh, totalp, hkfreq, pkfreq)

	return newhpop, newppop
}

//Replicate takes a population size and returns the new population size based on
//replication rate.
func Replicate(popsize int, rate float64) int {
	if popsize == 0 {
		return 0
	}
	if rate < 0 { //population size decreasing
		dec := (-1.0) * rate
		if int(float64(popsize)*dec) > popsize { //population number cannot be negative
			return 0
		}
		return popsize - int(float64(popsize)*dec)
	} else { //population size increasing
		return popsize + int(float64(popsize)*rate)
	}
}

//ComputeParameters takes population parameters as a map of strings to floats
//and returns calculated payoffs for each possible interaction
func ComputePayoffs(params map[string]float64) map[string]float64 {
	payoffs := make(map[string]float64)

	//calculate payoffs as per model (Renaud & Meeus, 1991)
	payoffs["Whkk"] = (1.0 - params["mortality"]) * (params["hostrate"] - params["rescost"])
	payoffs["Whkd"] = params["hostrate"] - params["rescost"]
	payoffs["Whdd"] = params["hostrate"] - params["gifthost"]
	payoffs["Wpkk"] = (params["mortality"]) * (params["parrate"] - params["vircost"])
	payoffs["Wpkd"] = params["parrate"] - params["vircost"]
	payoffs["Wpdd"] = params["parrate"] - params["giftpar"]

	return payoffs
}

//ComputeFitness takes a payoff matrix and current population frequency, and
//returns the host and parasite rate of increase according to replicator dynamic
//equations (Hoeffman & Yoeli, 2013 (MIT OCW)) as a map of strings to floats
func ComputeFitness(hostpop, parpop Population, payoffs map[string]float64) map[string]float64 {
	rates := make(map[string]float64)

	//finding frequencies in initial host population
	killhfreq := CountKillFreq(hostpop)
	diphfreq := 1.0 - killhfreq

	//finding frequencies in initial parasite population
	killpfreq := CountKillFreq(parpop)
	dippfreq := 1.0 - killpfreq

	//computing fitness values for killer and diplomat hosts
	fitness_hk := killpfreq*payoffs["Whkk"] + dippfreq*payoffs["Whkd"]
	fitness_hd := killpfreq*0 + dippfreq*payoffs["Whdd"]

	//computing fitness values for killer and diplomat parasites
	fitness_pk := killhfreq*payoffs["Wpkk"] + diphfreq*payoffs["Wpkd"]
	fitness_pd := killhfreq*0 + diphfreq*payoffs["Wpdd"]

	//computing population fitness for host and parasite popns
	host_fitness := killhfreq*fitness_hk + diphfreq*fitness_hd
	par_fitness := killpfreq*fitness_pk + dippfreq*fitness_pd

	//computing rate of increase for killer and diplomat hosts
	rates["hkrate"] = killhfreq * (fitness_hk - host_fitness)
	rates["hdrate"] = diphfreq * (fitness_hd - host_fitness)

	//computing rate of increase for killer and diplomat parasites
	rates["pkrate"] = killpfreq * (fitness_pk - par_fitness)
	rates["pdrate"] = dippfreq * (fitness_pd - par_fitness)

	return rates
}

/*
********************************************************************************
WRITING FUNCTION
********************************************************************************
*/

//WritePopns writes a poplist to file, with host size, population size, host frequency, and
//population frequency in separate lines in a .txt file
func WritePopns(filename, header string, poplist []Populations) {
	outfile, err := os.Create(filename)
	if err != nil {
		fmt.Println("Error! couldn’t create", filename)
	}
	defer outfile.Close()

	//print plot heading first
	fmt.Fprintf(outfile, "%v\n", header)

	//printing host sizes
	for i := range poplist {
		fmt.Fprint(outfile, len(poplist[i][0]), " ")
	}
	fmt.Fprintf(outfile, "\n")

	//printing host kill frequencies
	for i := range poplist {
		fmt.Fprint(outfile, CountKillFreq(poplist[i][0]), " ")
	}
	fmt.Fprintf(outfile, "\n")

	//printing parasite sizes
	for i := range poplist {
		fmt.Fprint(outfile, len(poplist[i][1]), " ")
	}
	fmt.Fprintf(outfile, "\n")

	//printing parasite kill frequencies
	for i := range poplist {
		fmt.Fprint(outfile, CountKillFreq(poplist[i][1]), " ")
	}
	fmt.Fprintf(outfile, "\n")
	fmt.Println("Finished printing to file!")
}
