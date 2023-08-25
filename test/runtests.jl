using BPOrelaxations
using Test
using Printf


@testset "test" begin    
   problem = readmc("../benchmark/test", "bqp50.mc")
   #option = Option(parseOption("SheraliAdams1")) 
   #option = Option(parseOption("Lasserre1")) 
   option = Option(parseOption("Submodular")) 
   #option = Option(parseOption("Solve")) 
   relax(problem, option)
end 