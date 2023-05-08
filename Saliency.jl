using MoleculeGraphs
using Plots
using Flux

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Select a compound to visualize
compound_index = 1
smiles = tox21_data[compound_index, :smiles]
molgraph = smiles_to_molgraph(smiles)

# Train the neural network model on the Tox21 dataset
model = train_model(inputs, outputs, n_epochs=10)

# Compute the saliency map for the compound
saliency_map = compute_saliency_map(model, molgraph)

# Create the saliency map plot
plot_molecule_saliency_map(molgraph, saliency_map, size=(500, 500))
