#!/bin/bash
algorithms=("SheraliAdams1" "Lasserre1" "Submodular1" "Submodular2" "Submodular3" "Submodular4")
project_dir=$(pwd)
datapath="${project_dir}/benchmark"
resultpath="${project_dir}/result"
jualia_bin="julia"

runInstance() {
    benchmark_dir=$1
    instance=$2
    algo=$3
    result_dir=$4

    ${jualia_bin}  "${project_dir}/experiment/runbenchmark.jl" "$benchmark_dir" "$instance" "$algo" "$result_dir"
                
}
export -f runInstance


benchmarks=$(ls ${datapath})

echo $benchmarks


for benchmark in $benchmarks
do
    if [ $benchmark == "test" ]
    then
        continue
    fi
    
    instances=$(ls $datapath/$benchmark)

    for instance in  $instances
    do
        for algo in ${algorithms[@]}
        do
            runInstance  "$datapath/$benchmark"  "$instance"  "$algo"  "$resultpath" 
        done
    done

    #find   $datapath/$benchmark -name *.cbp
done

