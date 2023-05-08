using CSV
using Flux
using Statistics
using ROCAnalysis
using Plots

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Extract the molecular weight and boiling point for each compound
molweights = [430.5, 356.9, 405.5, 505.6, 315.8, 401.4, 494.5, 412.6, 383.5, 305.4]
boilingpoints = [203.5, 259.2, 237.5, 402.9, 63.5, 364.6, 601.4, 387.4, 319.8, 210.9]

# Convert the SMILES strings to molecular graphs
molgraphs = [smiles_to_molgraph(smiles) for smiles in tox21_data[:smiles]]

# Extract the molecular descriptors
descriptors = [extract_descriptors(molgraph) for molgraph in molgraphs]

# Add the molecular weight and boiling point as additional features
for i in 1:length(descriptors)
    descriptors[i]["MolWeight"] = molweights[i]
    descriptors[i]["BoilingPoint"] = boilingpoints[i]
end

# Convert the descriptors to a matrix
inputs = Matrix([values(desc) for desc in descriptors])

# Train the neural network model on the Tox21 dataset with the augmented descriptors
model = train_model(inputs, outputs, n_epochs=10)

# Evaluate the performance of the model on the Tox21 dataset with the augmented descriptors
evaluate_model(model
