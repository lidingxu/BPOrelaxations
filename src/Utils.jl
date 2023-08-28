

@enum Algorithm::Int begin
   SheraliAdams   =   1   # the first level of Sherali-Adams
   Lasserre       =   2   # the first level of Lasserre
   Submodular     =   3   # the first level of submodular
   Solve          =   4
end


struct Option
   algorithm::Algorithm    # relaxation algorithm
   level::Int              # level
   time_limit::Float64     # time limit
   log_level::Int          # log level
   thread::Int             # thread number
   silent::Bool            # silent model

   function Option(algorithm::Algorithm, level::Int=1, time_limit::Float64=3600.0, log_level::Int=1, thread::Int=1, silent::Bool=true)
      new(algorithm, level, time_limit, log_level, thread, silent)
   end
end

function parseOption(string::String)
   level = 0
   if string == "SheraliAdams1"
      algorithm = SheraliAdams
      level = 1
   elseif  string == "Lasserre1"
      algorithm = Lasserre
      level = 1
   elseif  string == "Submodular1"
      algorithm = Submodular
      level = 1
   elseif  string == "Submodular2"
      algorithm = Submodular
      level = 2
   elseif  string == "Submodular3"
      algorithm = Submodular
      level = 3
   elseif  string == "Submodular4"
      algorithm = Submodular
      level = 4
   elseif string == "Solve"
      algorithm = Solve
   else
      "algorithm is not supported!"
      exit(1)
   end
   return algorithm, level
end

mutable struct Stat
   termintaion_status
   val::Float64
   gap::Float64
   time::Float64
   algo::Algorithm
   instance::String
   
   function Stat()
      new()
   end
end


function getValues(m::JuMP.Model)
   av = JuMP.all_variables(m)
   d = [(string(v), value(v)) for v in av]
   return d
end


function set_optimal_start_values(model::Model)
   # Store a mapping of the variable primal solution
   variable_primal = Dict(x => value(x) for x in all_variables(model))
   # In the following, we loop through every constraint and store a mapping
   # from the constraint index to a tuple containing the primal and dual
   # solutions.
   constraint_solution = Dict()
   for (F, S) in list_of_constraint_types(model)
       # We add a try-catch here because some constraint types might not
       # support getting the primal or dual solution.
       try
           for ci in all_constraints(model, F, S)
               constraint_solution[ci] = (value(ci), dual(ci))
           end
       catch
           @info("Something went wrong getting $F-in-$S. Skipping")
       end
   end
   # Now we can loop through our cached solutions and set the starting values.
   for (x, primal_start) in variable_primal
       set_start_value(x, primal_start)
   end
   for (ci, (primal_start, dual_start)) in constraint_solution
       set_start_value(ci, primal_start)
       set_dual_start_value(ci, dual_start)
   end
   return
end
