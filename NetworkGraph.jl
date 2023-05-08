using CSV
using LightGraphs
using Plots

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Convert the SMILES strings to molecular graphs
molgraphs = [smiles_to_molgraph(smiles) for smiles in tox21_data[:smiles]]

# Extract the molecular descriptors
descriptors = [extract_descriptors(molgraph) for molgraph in molgraphs]

# Compute the pairwise similarity between the compounds based on their descriptors
similarity_matrix = pairwise_similarity(descriptors)

# Convert the similarity matrix to a graph
g = Graph(similarity_matrix)

# Create the network graph
graphplot(g, layout=GraphPlot.circular_layout, nodefillc=:lightblue, nodeborderc=:white, nodelabel=tox21_data[:id], nodelabelfontsize=6, edgelinewidth=
