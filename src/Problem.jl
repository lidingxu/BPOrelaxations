# problem data
struct Problem
   # parameters
   instance::String # name of instance
   num_vars::Int # number of variables
   num_terms_nl::Int # number of nonlinear terms

   # binary polynomial data
   terms_l::Vector{Tuple{Int, Float64}} # vector of var index, coefficient
   terms_nl::Vector{Tuple{Set{Int}, Float64}} # vector of vars indices, coefficient
   inds_nl_poscoef::Vector{Int}
   inds_nl_negcoef::Vector{Int}
   map_nl::Dict{Set{Int}, Int}

   # construtor
   function Problem(instance::String, num_vars::Int, num_terms_nl::Int, terms_l::Vector{Tuple{Int, Float64}}, terms_nl::Vector{Tuple{Set{Int}, Float64}},
      inds_nl_poscoef::Vector{Int}, inds_nl_negcoef::Vector{Int})
      map_nl = Dict{Set{Int}, Int}()
      for (i, term_nl) in enumerate(terms_nl)
         map_nl[term_nl[1]] = i
      end
      problem = new(instance::String, num_vars::Int, num_terms_nl::Int, terms_l::Vector{Tuple{Int, Float64}}, terms_nl::Vector{Tuple{Set{Int}, Float64}},
      inds_nl_poscoef::Vector{Int}, inds_nl_negcoef::Vector{Int}, map_nl::Dict{Set{Int}, Int})
      return problem
   end
end

# a polynomial f(x)
mutable struct Polynomial
   # original binary polynomial data
   terms_pos_nl::Vector{Set{Int}} # nonlinear terms
   inds_pos_nl::Vector{Int} # nonlinear terms
   selectors::Vector{Dict{Int, Int}} #

   # constructor
   Polynomial(terms_pos_nl::Vector{Set{Int}}, inds_pos_nl::Vector{Int}, selectors::Vector{Dict{Int, Int}}) =
      new(terms_pos_nl, inds_pos_nl, selectors)
   
   function Polynomial()
      terms_pos_nl=Vector{Set{Int}}()
      inds_pos_nl=Vector{Int}()
      selectors=Vector{Dict{Int, Int}}()
      push!(selectors, Dict{Int, Int}())
      new(terms_pos_nl, inds_pos_nl, selectors)
   end
end

# variables of g(x):= yx - f(x)
mutable struct PolyVar
   poly_id::Int
   term_constant::VariableRef # constant var of f(x)
   terms_l::Vector{VariableRef} # vars of linear term of f(x)
   terms_nl::Dict{Int, VariableRef} # vars of nonlinear terms of f(x)

   term_constant_val::Float64
   terms_l_vals::Vector{Float64}
   terms_nl_vals::Dict{Int, Float64}

   PolyVar() = new()
end