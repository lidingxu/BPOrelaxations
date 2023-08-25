#=========================================================
 benchmark test
=========================================================#

#include("../src/BPOrelaxations.jl")
using BPOrelaxations

function main(args)
    instance_dir = args[1]
    instance_name = args[2]
    algorithm = args[3]
    output_dir = args[4]

    problem = readmc(instance_dir, instance_name)
    option = Option(parseOption(algorithm)) 
    stat = relax(problem, option)
    abs_output = string(output_dir , "/" , instance_name , "." , algorithm) 
    stat_info = string("instance: ", instance_name, "\n", "obj: ", string(-stat.val), "\n", "time: ", string(stat.time), "\n", 
        "algo: ", algorithm)
    io = open(abs_output, "w")
    print(io, stat_info)
    close(io)
end

main(ARGS)
