

function readmc(datapath, filename)
   file_lines = readlines(datapath * "/" * filename)
   meta_data = file_lines[1]
   meta_data = split(meta_data)
   num_nodes = parse(Int, meta_data[1])
   num_edges = parse(Int, meta_data[2])
   
   terms_l = Vector{Float64}()
   for _ in  1:num_nodes
      push!(terms_l, 0.0)
   end

   terms_nl = Vector{Tuple{Set{Int}, Float64}}()
   ind = 1
   inds_nl_poscoef = Vector{Int}()
   inds_nl_negcoef = Vector{Int}()
   for line in file_lines[2:end]
      edge = split(line)
      v1 = parse(Int, edge[1])
      v2 = parse(Int, edge[2])
      @assert(v1 != v2)
      w = parse(Float64, edge[3])
      # w ( (1 - x_v1) x_v2 + (1 - x_v2) x_v1 ) = w ( x_v1 + x_v2  - 2 x_v1 x_v2)
      terms_l[v1] -= w
      terms_l[v2] -= w
      push!(terms_nl, (Set{Int}((v1, v2)), 2*w))
      if w > 0
         push!(inds_nl_poscoef, ind)
      else
         push!(inds_nl_negcoef, ind)
      end
      ind += 1
   end

   terms_l_ = Vector{Tuple{Int,Float64}}()
   for (j, w) in enumerate(terms_l)
      push!(terms_l_, (j,w))
   end
   
   num_vars = num_nodes
   num_terms_nl = num_edges

   return Problem(filename, num_vars, num_terms_nl, terms_l_, terms_nl, inds_nl_poscoef, inds_nl_negcoef)
end



function readBQP(datapath, filename)
   file_lines = readlines(datapath * "/" * filename)
   meta_data = file_lines[1]
   meta_data = split(meta_data)
   num_vars = parse(Int, meta_data[1])
   num_terms_nl = parse(Int, meta_data[2])
   
   terms_l = Vector{Float64}()
   for j in  1:num_vars
      push!(terms_l, j, 0.0)
   end

   terms_nl = Vector{Tuple{Set{Int}, Float64}}()
   ind = 1
   inds_nl_poscoef = Vector{Int}()
   inds_nl_negcoef = Vector{Int}()
   for line in file_lines[2:end]
      line_data = split(line)
      if length(line_data) == 2
         j = parse(Int, line_data[1])
         w = parse(Float64, line_data[2])
         terms_l[j] = w
      elseif length(line_data) == 3
         j1 = parse(Int, line_data[1])
         j2 = parse(Int, line_data[2])
         w = parse(Float64, line_data[3])
         push!(terms_nl, (Set{Int}((j1, j2)), w))
         if w > 0
            push!(inds_nl_poscoef, ind)
         else
            push!(inds_nl_negcoef, ind)
         end
         ind += 1
      end
   end

   terms_l_ = Vector{Tuple{Int,Float64}}()
   for (j, w) in enumerate(terms_l)
      push!(terms_l_, (j,w))
   end
   return Problem(filename, num_vars, num_terms_nl, terms_l_, terms_nl, inds_nl_poscoef, inds_nl_negcoef)
end