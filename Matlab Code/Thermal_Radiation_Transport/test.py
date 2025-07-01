import pickle as pkl
import numpy as np
# import tensor_networks as tn
import pytens as tn
import scipy.io


timesteps = [0, 2, 3, 4, 5, 1e4, 5e4, 99992]

# Load data from pickle files
with open("pre_rounding1.pkl", "rb") as f:
    data1 = pkl.load(f)

with open("pre_rounding2.pkl", "rb") as f:
    data2 = pkl.load(f)

def create_dict(data):
    """Creates a dictionary dict where dict[timestep] represents a Tensor-Train
    by a list of 3 numpy arrays (TT-cores)"""
    data_dict = {}
    for i, network in enumerate(data):
        data_dict[timesteps[i]] = [network.value(j) for j in range(len(network.shape()))]
    return data_dict

def print_info(data_dict):
    """Prints the shape of the TT-cores at each timestep"""
    for key, value in data_dict.items():
        print("Timestep : ", key, ", Shape : ", [val.shape for val in value])


data1_dict = create_dict(data1)
data2_dict = create_dict(data2)
# for k in timesteps:
#     name = 'data1_'+str(int(k))+'.mat'
#     scipy.io.savemat(name,{'core': data1_dict[k],'N':3})
#     name = 'data2_'+str(int(k))+'.mat'
#     scipy.io.savemat(name,{'core': data2_dict[k],'N':3})


# Uncomment to print the shape of TT-cores at each timestep
print("TT before rounding after transport solver update:")
print_info(data1_dict)

print()

print("TT before rounding after source term solver update:")
print_info(data2_dict)

# # Uncomment to save the dictionaries to pickle files
# with open("data1_dict.pkl", "wb") as f:
#     pkl.dump(data1_dict, f)
# with open("data2_dict.pkl", "wb") as f:
#     pkl.dump(data2_dict, f)
