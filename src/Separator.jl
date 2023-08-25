

mutable struct Separator
   # original binary polynomial data
   flow_net::DiGraph # node id 1 (source) and 2 (sink), 2 + k (k-th terms) , 2 + num_terms + k (k-th variables)
   in_edges_terms_l::Vector{Vector{Int}}
   capacity_matrix::Matrix{Float64}
   inf_entries::Vector{Tuple{Int, Int}}
   base_var_values::Vector{Float64}
   neg_indicators::Vector{Bool}
   original_caps::Vector{Float64}
   num_vars::Int
   num_terms::Int
   terms_nl::Vector{Set{Int}} # nonlinear terms
   inds_nl::Vector{Int} # nonlinear terms

   # constructor
   function Separator(num_vars::Int, terms_nl::Vector{Set{Int}}, inds_nl::Vector{Int}) 
      num_terms = length(inds_nl)
      flow_net = DiGraph(2 + num_terms + num_vars)
      in_edges_terms_l = Vector{Vector{Int}}()
      inf_entries = Vector{Tuple{Int, Int}}()
      for i in 1:num_terms
         add_edge!(flow_net, 1, 2 + i)
      end
      for j in 1:num_vars
         add_edge!(flow_net, 2 + num_terms + j, 2)
         push!(in_edges_terms_l, Vector{Int}())
      end
      for (i, term_nl) in enumerate(terms_nl)
         for j in term_nl
            add_edge!(flow_net, 2 + i, 2 + num_terms + j)
            push!(inf_entries, (2 + i, 2 + num_terms + j))
            push!(in_edges_terms_l[j], i)
         end
      end
      capacity_matrix = zeros(Float64, nv(flow_net), nv(flow_net))
      neg_indicators = Vector{Bool}(undef, num_vars)
      original_caps =  Vector{Float64}(undef, num_vars)
      base_var_values = zeros(Float64, num_vars)
      for (e1, e2) in inf_entries
         capacity_matrix[e1,e2] = typemax(Float64)
      end
      new(flow_net, in_edges_terms_l, capacity_matrix, inf_entries, base_var_values, neg_indicators, original_caps, num_vars, num_terms, terms_nl, inds_nl)
   end
end

# run every time to initilize a flow network of poly var
function InitFlowNetwork!(separator::Separator)
   fill!(separator.capacity_matrix, 0.0)
   for (e1, e2) in separator.inf_entries
      separator.capacity_matrix[e1, e2] = typemax(Float64)
   end
end

# run every time to partially set an a-priori flow network of poly var 
function partialSetFlowNetwork!(separator::Separator, inds_neg_nl::Vector{Int}, poly_var::PolyVar, polynomial::Polynomial)
   InitFlowNetwork!(separator)
   for j in 1:separator.num_vars
      separator.base_var_values[j] = poly_var.terms_l_vals[j]
   end
   for i in 1:separator.num_terms
      separator.capacity_matrix[1, 2 + i] = -poly_var.terms_nl_vals[inds_neg_nl[i]]
   end
end

#  run every time to fully set  flow network of poly var and a selector
function fullSetFlowNetwork!(separator::Separator, inds_neg_nl::Vector{Int}, poly_var::PolyVar, polynomial::Polynomial, selector::Dict{Int, Int})
   for i in polynomial.inds_pos_nl
      j = selector[i]
      separator.capacity_matrix[2 + separator.num_terms + j, 2] = poly_var.terms_nl_vals[i]
   end
   for j in 1:separator.num_vars
      separator.original_caps[j] = separator.base_var_values[j] + separator.capacity_matrix[2 + separator.num_terms + j, 2]
      separator.neg_indicators[j] = separator.original_caps[j] <= 0
      separator.capacity_matrix[2 + separator.num_terms + j, 2] = max(separator.original_caps[j], 0)
   end
end


# run min cut and return the value of x for which x_j = 1 indicates j is in the part1
function minCut!(separator::Separator, x::Vector{Bool})
   part1, part2, cut_val = GraphsFlows.mincut(separator.flow_net, 1, 2, separator.capacity_matrix, BoykovKolmogorovAlgorithm())   
   fill!(x, false)
   for j in 1:separator.num_vars
      if (2 + separator.num_terms + j) in part1
         x[j] = true
      elseif separator.neg_indicators[j]
         x[j] = true
      end
   end

   if is_debug
      verifyCut!(separator, x, part1, part2, cut_val)
   end
end

# verify the cut
function verifyCut!(separator::Separator, x::Vector{Bool}, part1, part2, cut_val)
   @assert(1 in part1)
   @assert(2 in part2)
   for j in 1:separator.num_vars
      if separator.capacity_matrix[2 + separator.num_terms + j, 2] <= 0
         @assert(x[j] == true)
      end 
   end
   for i in 1:separator.num_terms
      term_nl = separator.terms_nl[i]
      val = true
      for j in term_nl
         val = x[j] && val
      end 
      if val == true
         if !((2 + i) in part1)
            @assert(abs(separator.capacity_matrix[1, 2 + i]) <  epsilon)
         end
      else 
         @assert( (2 + i) in part2)
      end
   end
end


function addMaxFlowFormulation!(model, problem::Problem, terms_neg_nl::Vector{Set{Int}}, inds_neg_nl::Vector{Int}, poly_var::PolyVar, polynomial::Polynomial, selector::Dict{Int, Int}, separator::Separator)
   # node id 1 (source) and 2 (sink), 2 + k (k-th terms) , 2 + num_terms + k (k-th variables) 
   flow_vars = Dict{Tuple{Int, Int}, VariableRef}()

   num_terms = separator.num_terms
   num_vars = separator.num_vars
   terms_nl = separator.terms_nl
   in_edges_terms_l = separator.in_edges_terms_l
   # flow capacity constraints
   for i in 1:num_terms
      e = (1, 2 + i)
      flow_vars[e] = @variable(model)
      @constraint(model, 0 <= flow_vars[e])
      @constraint(model, flow_vars[e] <= -poly_var.terms_nl[inds_neg_nl[i]] )
   end
   for j in 1:num_vars
      e = (2 + num_terms + j, 2)
      flow_vars[e] = @variable(model)
      @constraint(model, flow_vars[e] <= poly_var.terms_l[j] +  sum( selector[i] == j ? poly_var.terms_nl[i] : 0 for i in polynomial.inds_pos_nl) )
      e_ = (1, 2 + num_terms + j)
      flow_vars[e_] = @variable(model)
      @constraint(model, flow_vars[e_] <= 0 )
   end
   for (i, term_nl) in enumerate(terms_nl)
      for j in term_nl
         e = (2 + i, 2 + num_terms + j)
         flow_vars[e] = @variable(model)
         @constraint(model, 0 <= flow_vars[e])
      end
   end

   # flow balance constraints
   for (i, term_nl) in enumerate(terms_nl)
      e_in = (1, 2 + i)
      es_out = [(2 + i, 2 + num_terms + j) for j in term_nl]
      @constraint(model, flow_vars[e_in] == sum(flow_vars[e_out] for e_out in es_out))
   end
   for j in 1:num_vars
      e_out = (2 + num_terms + j, 2)
      @constraint(model, sum(flow_vars[(2 + i, 2 + num_terms + j)] for i in in_edges_terms_l[j]) + flow_vars[(1, 2 + num_terms + j)] == flow_vars[e_out])
   end

   # nonnegativity constraint
   @constraint(model, poly_var.term_constant + sum( poly_var.terms_nl[i] for i in inds_neg_nl)  + sum(flow_vars[(2 + num_terms + j, 2)] for j in 1:num_vars) >= 0) 

end


