

function random_f() 
   return rand() * (1-0.1) + 0.1
end



ns = [40, 50]

nposs = [20, 30]

cuid = 1
for n in ns
   for npos in nposs
      m = n * (n-1) / 2
      pos_inds = rand(1:m, npos)
      ind = 0
      wvs = [(j, random_f()) for j in 1:n]
      wes = []
      for j1 in 1:n
         for j2 in (j1 + 1):n
            val = random_f()
            wval = ind in pos_inds ? val : ( val > 0.2 ? -val : 0)
            push!(wes, (j1, j2, wval))
            ind += 1
         end
      end
      m = length(wes)
      abs_output =  string("randBQPM." , cuid) 
      data = string(n, " ", m, "\n")
      for (j, wval) in wvs
         data *= string(j, " ", wval, "\n")
      end
      for (j1, j2, wval) in wes
         data *= string(j1, " ", j2, " ", wval, "\n")
      end
      io = open(abs_output, "w")
      print(io, data)
      close(io)
      global cuid += 1
   end
end