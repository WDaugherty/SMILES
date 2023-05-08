using RDKit
using DataFrames
using CSV
using Flux
using Statistics
using ROCAnalysis
using Plots
using HTTP

# Define a function to convert SMILES strings to molecular graphs
function smiles_to_molgraph(smiles::AbstractString)
    mol = Mol(smiles)
    atom_indices = [get_idx(mol, atom) for atom in atoms(mol)]
    bond_indices = [get_idx(mol, bond) for bond in bonds(mol)]
    adj_matrix = adjacency_matrix(mol)
    node_features = Dict("atom_type" => [get_atomic_num(atom) for atom in atoms(mol)])
    edge_features = Dict("bond_type" => [get_bond_type(bond) for bond in bonds(mol)])
    return (adj_matrix, atom_indices, bond_indices, node_features, edge_features)
end

# Define a function to extract molecular descriptors from a molecular graph
function extract_descriptors(molgraph::Tuple)
    adj_matrix, atom_indices, bond_indices, node_features, edge_features = molgraph
    mol = GraphMol(adj_matrix)
    rdfeats = RDKit.featFactory.GetFeaturesForMol(mol)
    descriptors = OrderedDict()
    for feat in rdfeats
        if RDKit.featTypeToString(feat.GetType()) == "2DPharm2D"
            continue
        end
        name = RDKit.featTypeToString(feat.GetType())
        val = feat.GetFamily() == RDKit.ATT_NumHs ? feat.GetNumHs() : feat.GetVal()
        descriptors[name] = val
    end
    return descriptors
end

# Define the neural network architecture
function create_model(n_inputs, n_hidden=50, dropout=0.5)
    Chain(
        Dense(n_inputs, n_hidden, relu),
        Dropout(dropout),
        Dense(n_hidden, n_hidden, relu),
        Dropout(dropout),
        Dense(n_hidden, 12, sigmoid)
    )
end

# Define the loss function
function binary_cross_entropy(y_pred, y_true)
    -mean(y_true .* log.(y_pred) .+ (1 .- y_true) .* log.(1 .- y_pred))
end

# Define the optimizer and the training function
function train_model(model, inputs, outputs, n_epochs=100, batch_size=128)
    optimizer = ADAM()
    for epoch in 1:n_epochs
        # Shuffle the data
        indices = shuffle(1:size(inputs, 1))
        # Train the model in mini-batches
        for batch_start in 1:batch_size:size(inputs, 1)-batch_size+1
            batch_end = batch_start+batch_size-1
            batch_indices = indices[batch_start:batch_end]
            x_batch = inputs[batch_indices, :]
            y_batch = outputs[batch_indices, :]
            loss = binary_cross_entropy(model(x_batch), y_batch)
            Flux.back!(loss)
            Flux.update!(optimizer, model)
        end
        # Print the loss after each epoch
        println("Epoch $epoch: Loss $(binary_cross_entropy(model(inputs), outputs))")
    end
end

# Define a function to evaluate the performance of the model
function evaluate_model(model, inputs, outputs)
    y_pred = model(inputs)
    auc_scores = [rocauc(outputs[:, i], y_pred[:, i]) for i in 1:size(outputs, 2)]
    # Print the AUC scores for each toxicity endpoint
    for (i, endpoint) in enumerate(["NR
