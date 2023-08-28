function relax(problem::Problem, option::Option)
   model = createModel(option, option.algorithm)
   if option.algorithm == SheraliAdams && option.level == 1
      stat = relaxSheraliAdams1(model, problem, option)
   elseif option.algorithm == Lasserre && option.level == 1
      stat = relaxLasserre1(model, problem, option)
   elseif option.algorithm == Submodular
      stat = relaxSubmodular(model, problem, option)
   elseif option.algorithm == Solve
      stat = solveSheraliAdams1(model, problem, option)
   else 
      print("unsuppoted algorithm")
      exit(1)
   end
   print(stat)
   return stat
end

# the following relaxations only work for degee 2 BPO
function relaxSheraliAdams1(model, problem::Problem, option::Option)::Stat
   # variables
   @variables(model, begin 
      x[j in 1:problem.num_vars] 
      y[i in 1:problem.num_terms_nl]
   end) 

   # constraints
   @constraints(model, begin 
      [i in 1:problem.num_terms_nl, j in problem.terms_nl[i][1]], y[i] <= x[j]
      [i in 1:problem.num_terms_nl, j in problem.terms_nl[i][1]], y[i] >= sum(x[j] for j in problem.terms_nl[i][1]) - 1
   end)

   # objective
   @objective(model, Min, sum(y[i] * problem.terms_nl[i][2] for i in 1:problem.num_terms_nl) + sum(x[j] * problem.terms_l[j][2] for j in 1:problem.num_vars) )
   optimize!(model)

   stat = Stat();
   stat.termintaion_status = termination_status(model)
   stat.val = objective_value(model)
   stat.time = solve_time(model)
   stat.algo = option.algorithm
   stat.instance = problem.instance

   return stat
end


# 
function solveSheraliAdams1(model, problem::Problem, option::Option)::Stat
   # variables
   @variables(model, begin 
      x[j in 1:problem.num_vars], Bin
      y[i in 1:problem.num_terms_nl], Bin
   end) 

   # constraints
   @constraints(model, begin 
      [i in 1:problem.num_terms_nl, j in problem.terms_nl[i][1]], y[i] <= x[j]
      [i in 1:problem.num_terms_nl, j in problem.terms_nl[i][1]], y[i] >= sum(x[j] for j in problem.terms_nl[i][1]) - 1
   end)

   # objective
   @objective(model, Min, sum(y[i] * problem.terms_nl[i][2] for i in 1:problem.num_terms_nl) + sum(x[j] * problem.terms_l[j][2] for j in 1:problem.num_vars) )
   optimize!(model)

   stat = Stat();
   stat.termintaion_status = termination_status(model)
   stat.val = objective_value(model)
   stat.time = solve_time(model)
   stat.algo = option.algorithm
   stat.instance = problem.instance

   return stat
end


function relaxLasserre1(model, problem::Problem, option::Option)::Stat
   matsize = 1+problem.num_vars
   # variables
   @variables(model, begin 
       x[j in 1:problem.num_vars] 
       y[i in 1:problem.num_terms_nl] 
      X[1:matsize, 1:matsize]
   end) 

   # constraints
   @constraints(model, begin 
      X[1,1] == 1
      [j in 1:problem.num_vars], X[1 + j, 1 + j] == x[j] 
      [j in 1:problem.num_vars], X[1, 1+j] == x[j] 
      [j in 1:problem.num_vars], X[1+j, 1] == x[j] 
      [i in 1:problem.num_terms_nl], X[collect(problem.terms_nl[i][1])[1] + 1, collect(problem.terms_nl[i][1])[2] + 1] == y[i]
      X >= 0, PSDCone()
   end)

   # objective
   @objective(model, Min, sum(y[i] * problem.terms_nl[i][2] for i in 1:problem.num_terms_nl) + sum(x[j] * problem.terms_l[j][2] for j in 1:problem.num_vars) )
   optimize!(model)

   stat = Stat();
   stat.termintaion_status = termination_status(model)
   stat.val = objective_value(model)
   stat.time = solve_time(model)
   stat.algo = option.algorithm
   stat.instance = problem.instance

   return stat   
end


function createModel(option::Option, algorithm::Algorithm)
   model =  JuMP.Model(() -> Mosek.Optimizer())
   set_optimizer_attribute(model, "MSK_IPAR_NUM_THREADS", option.thread)
   set_optimizer_attribute(model, "MSK_DPAR_OPTIMIZER_MAX_TIME", option.time_limit) 
   return model
end
