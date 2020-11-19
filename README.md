This is a folder containing codes used in the paper "Exact Distribution of Linkage Disequilibrium in the Presence of Mutation, Selection or Minor Allele Frequency Filtering".

In the folder of src, the julia script named "functions.jl" contains all the functions necessary for generating the exact distribution of LD over generations. The julia script named "Recursive_pvec.jl" is the script used to generate the p vectors every 50 generations (where output interval can be redefined in the script). Values of Ne, c, and u are required from the users. Selections are not considered in this script. The script named "Recursive_pvec_selection.jl"is the script used to generate the p vectors over generations when selection occurs. 

