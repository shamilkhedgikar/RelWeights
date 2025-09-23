# Setting up RelWeights

This repo provides a demonstration the concept of creating RelWeights.

## Setup

Download/Clone this repo using [https://github.com/shamilkhedgikar/RelWeights.git]

## Theory

RelWeights *(Relational Weights)* are a new class of Spatial Weights that allow us to incoporate multi-layer interactions into spatial econometric models.This class of Spatial Weights were identified were first proposed in an attempt to understand the effects of administrative boundaries on groundwater extraction.

RelWeights (Relational Weights) generalize the idea of contiguity by constructing neighborhood linkages across two or more overlayed layers. The approach treats inter-layer intersections as “inherited contiguities,” such that adjacency in one layer can induce connections in another. For example, two districts (Layer α) may not be contiguous in the Queen sense, but if they share sub-basins or infrastructure corridors (Layer β), RelWeights encode a relationship between them by creating a relationship between α and β.
Mathematically, if W is an n×k incidence matrix representing the intersection of polygons in layer α (n) and layer β (k). The RelWeights matrix is defined as:

R = WW′− diag(WW′)

Row-standardization produces a spatial weights matrix of dimension *n x n* that captures both direct and inherited adjacency. Compared with single-layer Queen weights, RelWeights matrices are typically denser, yielding higher variance and skewness in contiguity histograms. RelWeights allow us to:

- Model irregular but non-arbitrary spatial lags induced by overlay structures (e.g., flood risk propagating across watersheds into administrative districts).

- Capture fragmentation effects, where small overlay units (sub-basins, land parcels) create multiple inherited linkages.

- Generate richer variance structures in neighborhood density, with implications for detecting clustering, spillovers, and heterogeneous risk exposure.

For more theoretical details see:

[Khedgikar, S. S. (2023). Fragmented Futures: Understanding the Role of Spatial Boundaries on Groundwater in India.](https://knowledge.uchicago.edu/record/7136?ln=en&v=pdf)
