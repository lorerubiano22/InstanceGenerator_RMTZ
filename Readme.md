
$\LARGE{Instance}$ $\LARGE{generator}$ $\LARGE{for}$ $\LARGE{Revenue}$ $\LARGE{Maximizing}$ $\LARGE{Tariff}$ $\LARGE{Zone}$ $\LARGE{Planning}$ $\LARGE{for}$ $\LARGE{Public}$ $\LARGE{Transportation}$

The size of an instance is determined by the number of nodes in $\left|\mathcal{I}\right|$ and the connectivity of $\mathcal{G}^{\text{PT}}$, measured here as $\left|\mathcal{A}\right|/\left|\mathcal{I}\right|$. Each instance is a square with uneven $\left|\mathcal{I}\right|$ such that the district in the middle of the square represents the central business district or the city center attracting most of the trips.
 
Problem Statement

The problem consists in generating artificial instances that integrate realistic characteristics of the transportation system. The objective is to generate artificial instances that allow to test and validate models dedicated to design a zone counting tariff system. The notation used for the instance generator is presented below:

Notation

$\mathcal{I}$ as the set of nodes (stop) indexed by $i,j,k$ and $l$

$\mathcal{A}$ as the set of edges $\mathcal{A} \subseteq \mathcal{I} \times \mathcal{I}$ representing available direct connections without intermediate nodes.

$\mathcal{OD}_{ij}$ as the sequence of edges $(k,l) \in \mathcal{A}$ on the shortest-path from $i \in \mathcal{I}$ to $j \in \mathcal{I}$.

We designed the instance generator with these characteristics:

$Demand$: number of public transport trips determined by mode choice behaviour of customers.

Endogeneity: demand depends on price and tariff zones.

$Revenue$: demand and the price of a tariff zone determine the revenue.

Passengers always select the shortest route.

Artificial cities with one or more than one central business district: $\left|\mathcal{I}\right| / 2 \notin \mathbb{Z}$

Service areas with rectangular shapes: $\sqrt{\left|\mathcal{I}\right|} \in \mathbb{Z}$

Network connectivity (NC) defines the size of the set $\mathcal{A} $ as $\mid \mathcal{A} \mid$ = $\mid \mathcal{I} \mid \cdot (\mid \mathcal{I} \mid -\text{1}) \cdot$ NC , where NC >0  and NC $\leq$ 1


Number of the center node: $bcd$ = $\mid \mathcal{I} \mid$ - $\frac{ \mid \mathcal{I} \mid - \text{1}}{\text{2}}$

Distance ($\tau_{ij}$) from stop $i \in \mathcal{I}$ and stop $j \in \mathcal{I}$ is computed by: $\tau_{ij} = \frac{\text{10 + }\mathcal{U}\text{(-1, 1)}}{\sqrt{\left| \mathcal{I} \right|}}$

Distance from center to boundary of the city is defined by: $\varphi = \sqrt{\left|\mathcal{I}\right|} - \text{1}$

The X coordinates of the nodes of node $i$ are defined by: $a_i = \text{mod}\left(i + \varphi, \sqrt{\left|\mathcal{I}\right|}\right) + \text{1}$

The Y coordinates of the nodes of node $i$ are defined by: $b_i = \left\lceil i / \sqrt{\left|\mathcal{I}\right|}\right\rceil$

Rectilinear distance from $i$ to $bcd$: $\delta_i$

The graph $\mathcal{G}^{\text{PT}}$ of each instance can be described as a square matrix with $\mid \mathcal{I} \mid = 2g+1$ with $g \in \mathbb {Z}$

The node $i \in \mathcal{I}$ in the center of the matrix indicates the business center of each instance



The population density of each district determines the total travel demand. We assume that the population density of each district $i$ of the artificial city follows a uniform distribution $pop_i$= $\mathcal{U}$(10000, 20000), $pop_i$ is skewed to integer values and represents the population density in-district $i$. We perform a simple tuning to determine the parametrization of the population density in each district of the artificial city. From the above, the total population density of the artificial city is given by: $Pop= \sum_{i} pop_i $

Using the gravity model as a reference, we calculate the total number of trips demanded to the business center of the artificial city. The number of trips is determined from the trip attractivity coefficient from each node $i \in \mathcal{I}$ to the central business district. We define the trip attractivity $trip_{ij}$ as: 

![eq1](https://user-images.githubusercontent.com/39961021/202230459-7e7f3aa0-5764-4540-bf45-ac9189a59b3c.PNG)

$trip_{ij} =  pop_i \cdot A_j$
	
Then, we determine the proportion of total demand ($trips_{ij}$) that will be made using public transport. The demand for public transport depends on the cost of travel and travel times.

$\LARGE{Multinomial}$ $\LARGE{logit}$

![eq3](https://user-images.githubusercontent.com/39961021/202230896-01beaba1-c697-4e8f-a12c-c45f31b1498e.PNG)
														
Then the number of trips (total demand) from $i$ to $j$ are given by:

![eq4](https://user-images.githubusercontent.com/39961021/202231548-d2e41708-50ec-4fcd-b532-668d1c3776b8.PNG)

Expected Public transport revenue:

![eq5](https://user-images.githubusercontent.com/39961021/202231655-085cb890-f3b6-4e77-8e24-573b42063d7c.PNG)

$\LARGE{Inputs}$

Sets

$\mathcal{I}$: Public transport stops (nodes), indexed $i,j,k,l$

$\mathcal{N}$: District border nodes, indexed $n,m$

$\mathcal{P}$: Price systems with $p \in \mathcal{P}$

Parameters

$T_{ij}$: Maximum number of tariff zones along the shortest path from $i \in \mathcal{I}$ to $j \in \mathcal{I}$

$\pi_{pt}$: Price per tariff zone, if $t$ tariff zones are visited given price system $p$

$xdim$: Number of stops in a level of the $\mathcal{G}^{\text{PT}}$

NC: Network connectivity

Input example

![exampleInputs](https://user-images.githubusercontent.com/39961021/202231982-5d4e142d-ef2e-439c-ab8f-2912bbc6af34.PNG)

$\LARGE{Output}$

Sets

$\mathcal{I}$: Public transport stops (nodes), indexed $i,j,k,l$

$\mathcal{N}$: District border nodes, indexed $n,m$

$\mathcal{P}$: Price systems with $p \in \mathcal{P}$

$\mathcal{A}$: Public transport arcs 

$\mathcal{B}$: District border arcs 

$D_{ij}$: District border arcs $(n,m)$ corresponding to $\left(j,i\right)\in \mathcal{A}$, with $i \in \mathcal{I}$ and $j \in \mathcal{I}$ corresponding to the stops

$\mathcal{OD}_{ij}$: District border arcs $(n,m)$ along shortest path from stop $i \in \mathcal{I}$ to stop $j \in \mathcal{I}$

Parameters

$T_{ij}$: Maximum number of tariff zones along the shortest path from $i \in \mathcal{I}$ to $j \in \mathcal{I}$

$\pi_{pt}$: Price per tariff zone, if $t$ tariff zones are visited given price system $p$


$r_{ijt}\left( \pi_{pt} \right)$: Expected revenue corresponding to $\left(j,i\right)\in \mathcal{A}$, if $t= 1, \ldots, T_{ij}$ tariff zones are visited on the shortest path from $i \in \mathcal{I}$ to $j \in \mathcal{I}$

Instance visualization

![exampleOutput](https://user-images.githubusercontent.com/39961021/202222178-404ed0ec-ecfa-4751-9a40-5dde41fff07e.png)


												
