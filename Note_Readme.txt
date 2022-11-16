Instance generator for Revenue Maximizing Tariff Zone Planning for Public Transportation

The size of an instance is determined by the number of nodes in $\left|\mathcal{I}\right|$ and the connectivity of $\mathcal{G}^{\text{PT}}$, measured here as $\left|\mathcal{A}\right|/\left|\mathcal{I}\right|$. Each instance is a square with uneven $\left|\mathcal{I}\right|$ such that the district in the middle of the square represents the central business district or the city center attracting most of the trips.
 
Problem Statement

The problem consists in generating artificial instances that integrate realistic characteristics of the transportation system. The objective is to generate artificial instances that allow to test and validate models dedicated to design a zone counting tariff system. The notation used for the instance generator is presented below:

Notation

$\mathcal{I}$ as the set of nodes (stop) indexed by $i,j,k$ and $l$

$\mathcal{A}$ as the set of edges $\mathcal{A} \subseteq \mathcal{I} \times \mathcal{I}$ representing available direct connections without intermediate nodes.

$\mathcal{OD}_{ij}$ as the sequence of edges $(k,l) \in \mathcal{A}$ on the shortest-path from $i \in \mathcal{I}$ to $j \in \mathcal{I}$.

We designed the instance generator with these characteristics:


Demand: number of public transport trips determined by mode choice behaviour of customers.

Endogeneity: demand depends on price and tariff zones.

Revenue: demand and the price of a tariff zone determine the revenue.

Passengers always select the shortest route.

Artificial cities with one or more than one central business district: $\left|\mathcal{I}\right| / 2 \notin \mathbb{Z}$

Service areas with rectangular shapes: $\sqrt{\left|\mathcal{I}\right|} \in \mathbb{Z}$

Network connectivity (NC) defines the size of the set $\mathcal{A} $ as $\mid \mathcal{A} \mid$ = $\mid \mathcal{I} \mid \cdot (\mid \mathcal{I} \mid -\text{1}) \cdot$ NC , where NC is given by $\{0.25, 0.5, 0.75, 1\}$

Number of the center node: $bcd$ = $\mid \mathcal{I} \mid$ - $\frac{ \mid \mathcal{I} \mid - \text{1}}{\text{2}}$

Distance ($\tau_{ij}$) from stop $i \in \mathcal{I}$ and stop $j \in \mathcal{I}$ is computed by: $\tau_{ij} = \frac{\text{10 + }\mathcal{U}\text{(-1, 1)}}{\sqrt{\left| \mathcal{I} \right|}}$

Distance from center to boundary of the city is defined by: $\varphi = \sqrt{\left|\mathcal{I}\right|} - \text{1}$

The X coordinates of the nodes of node $i$ are defined by: $a_i = \text{mod}\left(i + \varphi, \sqrt{\left|\mathcal{I}\right|}\right) + \text{1}$

The Y coordinates of the nodes of node $i$ are defined by: $b_i = \left\lceil i / \sqrt{\left|\mathcal{I}\right|}\right\rceil$

Rectilinear distance from $i$ to $bcd$: $\delta_i$



The population density of each district determines the total travel demand. We assume that the population density of each district $i$ of the artificial city follows a uniform distribution $pop_i$= $\mathcal{U}$(10000, 20000), $pop_i$ is skewed to integer values and represents the population density in-district $i$. We perform a simple tuning to determine the parametrization of the population density in each district of the artificial city. From the above, the total population density of the artificial city is given by: $Pop= \sum_{i} pop_i $

Using the gravity model as a reference, we calculate the total number of trips demanded to the business center of the artificial city. The number of trips is determined from the trip attractivity coefficient from each node $i \in \mathcal{I}$ to the central business district. We define the trip attractivity $trip_{ij}$ as: 

$A_j$ := $j$ =bcd,     if  $Pop= \frac{\sum_{i} pop_i }{\delta_i}$; 0, otherwise

$trip_{ij}$ =  $pop_i \cdot A_j$


Then, we determine the proportion of total demand ($trips_{ij}$) that will be made using public transport. The demand for public transport depends on the cost of travel and travel times.

Multinomial logit



$PT_{ij}^{\text{\tiny public transport}} = \frac{\exp \left(\text{DetUtil}_{ij}^{\text{\tiny public transport}}\right)}{\sum_{\text{\tiny mode}}\text{DetUtil}_{ij}^{\text{\tiny mode}}}$																

Then the number of trips (total demand) from $i$ to $j$ are given by:

$Demand_{ij}^{\text{\tiny public transport}} = \text{Trips}_{ij} \cdot PT_{ij}^{\text{\tiny public transport}}\left({\beta^*}\right)$

Expected Public transport revenue:

$Revenue_{ij}^{\text{\tiny public transport}} = Price_{ij}^{\text{\tiny public transport}} \cdot Demand_{ij}^{\text{\tiny public transport}} =Price_{ij}^{\text{\tiny public transport}} \cdot \text{Trips}_{ij} \cdot P_{ij}^{\text{\tiny public transport}}\left({\beta}\right)$

${\beta} \ $ utility coefficients: weighting attributes																