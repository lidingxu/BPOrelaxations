function relaxSubmodular(model, problem::Problem, option::Option)::Stat
   # get decomposition polynomials, the shared negative nonlinear terms
   terms_neg_nl, inds_neg_nl, polynomials = preprocessSubmodular(problem, option)

   print("\n num of polynomials:", length(polynomials), "\n") 

   stat = extendedFormulation!(model, problem, terms_neg_nl, inds_neg_nl, polynomials, option)
   
   return stat
end

function extendedFormulation!(model, problem::Problem,  terms_neg_nl::Vector{Set{Int}}, inds_neg_nl::Vector{Int}, polynomials::Vector{Polynomial}, option::Option)::Stat
   # createsubmodular model
   poly_vars, lambda = createSubmodularModel!(model, problem, terms_neg_nl, inds_neg_nl, polynomials)

   # max flow based extended formulation
   constructExtended!(model, problem, terms_neg_nl, inds_neg_nl, poly_vars, polynomials)

   stat = Stat()
   optimize!(model)
   stat.termintaion_status = termination_status(model)
   stat.val = objective_value(model)
   stat.time = solve_time(model)
   stat.algo = option.algorithm
   stat.instance = problem.instance

   return stat
end

function constructExtended!(model, problem, terms_neg_nl, inds_neg_nl, poly_vars, polynomials)
   separator = Separator(problem.num_vars, terms_neg_nl, inds_neg_nl)
   num_polys = length(poly_vars)
   # construct for each polynomial
   for i in 1:num_polys
      poly_var = poly_vars[i]
      polynomial = polynomials[i]
      # construct for each seletor
      for selector in polynomial.selectors
         addMaxFlowFormulation!(model, problem, terms_neg_nl, inds_neg_nl, poly_var, polynomial, selector, separator)
      end
   end
end

function projectedFormulation!(model, problem::Problem,  terms_neg_nl::Vector{Set{Int}}, inds_neg_nl::Vector{Int}, polynomials::Vector{Polynomial}, option::Option)::Stat
   # createsubmodular model
   poly_vars, lambda = createSubmodularModel!(model, problem, terms_neg_nl, inds_neg_nl, polynomials)

   # add random some cuts to prevent dual infeasibility (primal unbounded)
   x_global = Vector{Bool}(undef, problem.num_vars)
   xs = Vector{Vector{Bool}}()

   is_initi = true

   if is_initi
      push!(xs, falses(problem.num_vars))
      print("\n", problem.num_vars, "\n")
      for j in 1:problem.num_vars
         x = falses(problem.num_vars)
         x[j] = true
         push!(xs, x)
         x = trues(problem.num_vars)
         x[j] = false
         push!(xs, x)
      end
      push!(xs, trues(problem.num_vars))
   end

   # the vector that holds the Boolean values of negative nonlinear term 
   terms_neg_nl_values = Vector{Bool}(undef, length(terms_neg_nl))
   shrinkSubmodularModel!(model, problem, terms_neg_nl, terms_neg_nl_values, inds_neg_nl, poly_vars, polynomials, xs)

   #print(model)
   # initilize the separator' data
   separator = Separator(problem.num_vars, terms_neg_nl, inds_neg_nl)

   is_print_model = false
   stat = Stat()
   counter = 1
   is_zero_shrink = false
   while true
      optimize!(model)
      print("\n obj:", objective_value(model), " status:", termination_status(model), "\n")
      #print(getValues(model))
      if termination_status(model) != OPTIMAL || (abs(objective_value(model)) < 1e-4 && is_zero_shrink)
         if termination_status(model) == DUAL_INFEASIBLE || (abs(objective_value(model)) < 1e-4 && is_zero_shrink)
            print("\n shrink \n")
            xs = Vector{Vector{Bool}}()
            for i in 1:problem.num_vars
               x = rand(Bool, problem.num_vars)
               push!(xs, x)
            end
            shrinkSubmodularModel!(model, problem, terms_neg_nl, terms_neg_nl_values, inds_neg_nl, poly_vars, polynomials, xs)
            continue
         else
            stat.termintaion_status = termination_status(model)
            stat.val = objective_value(model)
            stat.time = solve_time(model)
            stat.algo = option.algorithm
            stat.instance = problem.instance
            break
         end
      end
      is_separated = separate!(model, problem, terms_neg_nl, terms_neg_nl_values, inds_neg_nl, poly_vars, polynomials, separator, x_global)
      if ! is_separated
         stat.termintaion_status = termination_status(model)
         stat.val = objective_value(model)
         stat.time = solve_time(model)
         stat.algo = option.algorithm
         stat.instance = problem.instance
         break
      end
   end
   return stat
end

# create a submodular model
function createSubmodularModel!(model, problem::Problem, terms_neg_nl::Vector{Set{Int}}, inds_neg_nl::Vector{Int}, polynomials::Vector{Polynomial})
   poly_vars = Vector{PolyVar}() 
   for (poly_id, polynomial) in enumerate(polynomials)
      poly_var =  PolyVar()
      poly_var.poly_id = poly_id
      poly_var.term_constant = @variable(model, base_name = string("const", poly_id))
      poly_var.terms_l = Vector{VariableRef}()
      poly_var.terms_nl  = Dict{Int, VariableRef}()
      for i in inds_neg_nl
         poly_var.terms_nl[i] =  @variable(model, base_name = string("w", poly_id, "-_", i))
         @constraint(model, poly_var.terms_nl[i] <= 0)
      end
      for i in polynomial.inds_pos_nl
         poly_var.terms_nl[i] =  @variable(model, base_name = string("w", poly_id, "+_", i))
         @constraint(model, poly_var.terms_nl[i] >= 0)
      end
      for j in 1:problem.num_vars
         push!(poly_var.terms_l, @variable(model, base_name = string("x", poly_id, "_", j)))
      end
      poly_var.terms_l_vals = Vector{Float64}(undef, problem.num_vars)
      poly_var.terms_nl_vals = Dict{Int, Float64}()
      push!(poly_vars, poly_var)
   end

   lambda = @variable(model, base_name = "lambda")

   #decompose
   @constraints(model, begin
      -lambda >= sum(poly_var.term_constant for poly_var in poly_vars) # constant
      [j in 1:problem.num_vars], problem.terms_l[j][2] >= sum(poly_var.terms_l[j] for poly_var in poly_vars) # linear terms
      [i in 1:problem.num_terms_nl], problem.terms_nl[i][2] >= sum( haskey(poly_var.terms_nl, i) ? poly_var.terms_nl[i] : 0 for poly_var in poly_vars) # nonlinear terms
   end)

   @objective(model, Max, lambda)

   return poly_vars, lambda
end

# shrk a submodular model by adding some random cuts
function shrinkSubmodularModel!(model, problem::Problem, terms_neg_nl::Vector{Set{Int}}, terms_neg_nl_values::Vector{Bool}, inds_neg_nl::Vector{Int}, poly_vars::Vector{PolyVar}, polynomials::Vector{Polynomial}, xs::Vector{Vector{Bool}})
   num_polys = length(poly_vars)
   for x in xs
      for (ind, term_neg_nl) in enumerate(terms_neg_nl)
         terms_neg_nl_values[ind] =  true 
         for j in term_neg_nl
            if x[j] == false
               terms_neg_nl_values[ind] =  false
               break
            end
         end
      end
      for i in 1:num_polys
         poly_var = poly_vars[i]
         polynomial = polynomials[i]
         for selector in polynomial.selectors
            addSubmodularConstraint!(model, problem, terms_neg_nl_values, inds_neg_nl, poly_var, polynomial, x, selector)
         end
      end   
   end
end

function separate!(model, problem::Problem, terms_neg_nl::Vector{Set{Int}}, terms_neg_nl_values::Vector{Bool}, inds_neg_nl::Vector{Int}, poly_vars::Vector{PolyVar}, polynomials::Vector{Polynomial}, separator::Separator, x_global::Vector{Bool})
   num_polys = length(poly_vars)
   is_separated = false
   # query the values of variables
   for i in 1:num_polys
      poly_var = poly_vars[i]
      poly_var.term_constant_val = value(poly_var.term_constant)
      for j in 1:problem.num_vars
         poly_var.terms_l_vals[j] = value(poly_var.terms_l[j])
      end
      for ind in keys(poly_var.terms_nl)
         poly_var.terms_nl_vals[ind] = value(poly_var.terms_nl[ind])
      end
      #print(poly_var.term_constant_val, poly_var.terms_l_vals, poly_var.y_vals, poly_var.terms_nl_vals, "\n")
   end   
   # separate
   for i in 1:num_polys
      poly_var = poly_vars[i]
      polynomial = polynomials[i]
      InitFlowNetwork!(separator)
      partialSetFlowNetwork!(separator, inds_neg_nl, poly_var, polynomial)
      for selector in polynomial.selectors
         fullSetFlowNetwork!(separator, inds_neg_nl, poly_var, polynomial, selector)
         minCut!(separator, x_global)
         
         poly_val = polyVal!(problem, terms_neg_nl, terms_neg_nl_values, inds_neg_nl, poly_var, polynomial, x_global, selector)
         # test the f(x) < 0, i.e., the constraint f(x) >= 0 is violated
         if poly_val < -epsilon
            #print("\n sep ", poly_val, "\n")
            is_separated = true
            addSubmodularConstraint!(model, problem, terms_neg_nl_values, inds_neg_nl, poly_var, polynomial, x_global, selector)
         end
      end
   end
   return is_separated
end


function polyVal!(problem::Problem, terms_neg_nl::Vector{Set{Int}}, terms_neg_nl_values::Vector{Bool},  inds_neg_nl::Vector{Int}, poly_var::PolyVar, polynomial::Polynomial, x::Union{Vector{Bool}, BitVector}, selector::Dict{Int, Int})
   for (ind, term_neg_nl) in enumerate(terms_neg_nl)
      terms_neg_nl_values[ind] =  true 
      for j in term_neg_nl
         if x[j] == false
            terms_neg_nl_values[ind] =  false
            break
         end
      end
   end
   f_val = poly_var.term_constant_val
   for j in 1:problem.num_vars
      f_val += x[j] ? poly_var.terms_l_vals[j] : 0
   end
   for (ind, i) in enumerate(inds_neg_nl)
      f_val += terms_neg_nl_values[ind] ? poly_var.terms_nl_vals[i] : 0
   end
   for i in polynomial.inds_pos_nl
      f_val += x[selector[i]] ? poly_var.terms_nl_vals[i] : 0
   end
   return f_val
end

# add a constraint  f(x) >= 0
function addSubmodularConstraint!(model, problem::Problem, terms_neg_nl_values::Vector{Bool}, inds_neg_nl::Vector{Int}, poly_var::PolyVar, polynomial::Polynomial, x::Vector{Bool}, selector::Dict{Int, Int})
   # values of negative nonlinear terms
   ex = @expression(model, sum( x[j] ?  poly_var.terms_l[j] : 0 for j in 1:problem.num_vars) + poly_var.term_constant)
   ex = @expression(model, ex + sum( terms_neg_nl_values[i] ? poly_var.terms_nl[ind_neg_nl] : 0 for (i, ind_neg_nl) in enumerate(inds_neg_nl) ) )
   @constraint(model, ex + sum( x[selector[inds_pos_nl]] ? poly_var.terms_nl[inds_pos_nl] : 0 for inds_pos_nl in polynomial.inds_pos_nl) >= 0)
end

# get decomposition polynomials, the shared negative nonlinear terms
function preprocessSubmodular(problem::Problem, option::Option)
   # partition the positive nonlinear terms
   collectors = partition(problem, option.level)
   # use the multilinear technique to extend the positive nonlinear terms
   transforms = extendMultilinear!(problem, collectors)

   # record negative nonlinear polynomials
   terms_neg_nl = Vector{Set{Int}}()
   inds_neg_nl = problem.inds_nl_negcoef
   for ind in problem.inds_nl_negcoef
      push!(terms_neg_nl, problem.terms_nl[ind][1])
   end

   # build the polynomials, whose positive nonlinear terms are recorded
   polynomials = Vector{Polynomial}()
   # the first deal with none of positive nonlinear terms
   if length(transforms) == 0
      push!(polynomials, Polynomial())
   end
   for (i, selectors) in enumerate(transforms)
      terms_pos_nl = Vector{Set{Int}}()
      inds_pos_nl = collect(collectors[i])
      for ind in inds_pos_nl
         push!(terms_pos_nl, problem.terms_nl[ind][1])
      end
      push!(polynomials, Polynomial(terms_pos_nl, inds_pos_nl, selectors))
      if is_debug
         for ind_pos_nl in polynomials[end].inds_pos_nl
            for selector in selectors
               if !haskey(selector, ind_pos_nl)
                  print("\n", ind_pos_nl, " ", inds_pos_nl, " ", polynomials[i].inds_pos_nl, " ", collectors[i], " ", selector, "\n")
               end
               @assert(haskey(selector, ind_pos_nl))
            end
         end
      end
   end

   return terms_neg_nl, inds_neg_nl, polynomials
end


# the parition algorithm prefers to partition the nonlinear terms with the same variables
function partition(problem::Problem, card::Int)
   inds = Set(problem.inds_nl_poscoef)
   collectors = Vector{Set{Int}}()
   while !isempty(inds)
      collector = Set{Int}()
      # add an ind of a positive nonlinear term
      ind = first(inds)
      push!(collector, ind)
      delete!(inds, ind)
      # get vars of the positive nonlinear term
      vars = problem.terms_nl[ind][1] 
      # the cardinality of this collector
      card_ = min(card - 1, length(inds))
      # a heuristic search to find terms with less variables
      for _ in 1:card_
         ind_add = first(inds)
         max_inter_size = -1
         for ind_search in inds
            inter_size = length(intersect(problem.terms_nl[ind_search][1], vars))
            if inter_size > max_inter_size
               ind_add = ind_search
               max_inter_size = inter_size
            end
         end
         push!(collector, ind_add)
         delete!(inds, ind_add)
         vars = union(vars, problem.terms_nl[ind_add][1])
      end
      # add the collector
      push!(collectors, collector)
   end

   if is_debug
      print("debug in partition...\n")
      inds_pos = Set(problem.inds_nl_poscoef)
      inds_union = Set{Int}()
      for collector in collectors
         inds_union = union(inds_union, collector)
         @assert( !isempty(intersect(inds_union, collector)) )
      end
      @assert( isempty(setdiff(inds_pos, inds_union)) )
   end

   return collectors
end

# build all selectors
function recursiveExtend!(terms_vars::Vector{Vector{Int}}, terms_inds::Vector{Int}, selector::Dict{Int, Int}, selectors::Vector{Dict{Int, Int}}, cur_ind::Int)
   for term_var in terms_vars[cur_ind]
      tmp_selector = copy(selector)
      term_ind = terms_inds[cur_ind]
      tmp_selector[term_ind] = term_var
      if cur_ind == length(terms_vars)
         push!(selectors, tmp_selector)
      else
         recursiveExtend!(terms_vars, terms_inds, tmp_selector, selectors, cur_ind + 1)
      end
   end
end

# use the multilinear technique to extend the positive nonlinear terms of each polynomial
function extendMultilinear!(problem::Problem, collectors::Vector{Set{Int}})
   num_collector = length(collectors)
   transforms = Vector{Vector{Dict{Int, Int}}}()
   for i in 1:num_collector
      # get the positive nonlinear terms
      collector = collectors[i]
      terms_vars = Vector{Vector{Int}}()
      terms_inds = Vector{Int}()
      # get all the variables
      for ind in collector
         push!(terms_inds, ind)
         push!(terms_vars, collect(problem.terms_nl[ind][1]))
      end
      # construct the selectors that map a nonlinear term to a variable
      selectors = Vector{Dict{Int, Int}}()
      recursiveExtend!(terms_vars, terms_inds, Dict{Int, Int}(), selectors, 1)

      if is_debug
         num_selectors = 1
         for ind in collector
            num_selectors *= length(problem.terms_nl[ind][1])
         end
         if i == 1
            print(selectors)
         end
         @assert(length(selectors) == num_selectors)
      end
      push!(transforms, selectors)
   end
   return transforms
end