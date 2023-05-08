using RDKit
using DataFrames
using CSV
using Flux
using Statistics
using ROCAnalysis
using Plots
using HTTP

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Convert the SMILES strings to molecular graphs
molgraphs = [smiles_to_molgraph(smiles) for smiles in tox21_data[:smiles]]

# Extract the existing molecular descriptors
descriptors = [extract_descriptors(molgraph) for molgraph in molgraphs]

# Extract the Morgan fingerprint descriptor
morgan = [morgansimilarity(molgraph, radius=2) for molgraph in molgraphs]

# Add the Morgan fingerprint descriptor to the existing set of descriptors
for i in 1:length(descriptors)
    descriptors[i]["Morgan"] = morgan[i]
end

# Convert the descriptors to a matrix
inputs = Matrix([values(desc) for desc in descriptors])

# Train the neural network model on the Tox21 dataset with the augmented descriptors
model = train_model(inputs, outputs, n_epochs=10)

# Evaluate the performance of the model on the Tox21 dataset with the augmented descriptors
evaluate_model(model, inputs, outputs)
