using RDKit
using DataFrames
using CSV
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

# Download and load the Tox21 dataset
url = "https://tripod.nih.gov/tox21/challenge/download?id=tox21_10k_data_csv.zip"
filename = "tox21_10k_data_csv.zip"
data_dir = "tox21_data"
mkdir(data_dir, exist_ok=true)
download(url, joinpath(data_dir, filename))
CSV.File(joinpath(data_dir, "tox21_10k_data.csv"), header=true) |> DataFrame

# Convert the SMILES strings to molecular graphs
tox21_data = CSV.File(joinpath(data_dir, "tox21_10k_data.csv"), header=true) |> DataFrame
molgraphs = [smiles_to_molgraph(smiles) for smiles in tox21_data.smiles]

# Extract molecular descriptors from the molecular graphs and add them to the DataFrame
descriptors = [extract_descriptors(molgraph) for molgraph in molgraphs]
for (name, values) in Iterators.flatten(descriptors)
    tox21_data[Symbol(name)] = values
end
