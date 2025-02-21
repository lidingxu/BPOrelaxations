using JuMP
using MosekTools

# f = - c1 x2 x3 - c2 x1 x3 x4 - c3 x3 x5 +  \sum_{j} fj xj + f0

f0 = 7
c = [1,2,5]
f = [0,1,1,-1,0]
terms = [([2,3], -c[1]), ([1,3,4], -c[2]), ([3,5], -c[3])]

# print f
f_string = "f:"

# Add the terms from `terms`
for (i, term) in enumerate(terms)
    global f_string
    vars = ["x_" * string(j) for j in term[1]]
    coef = term[2]
    term_str = string( coef == 1 ? "" : (coef == -1 ? "-" : coef) ) * " " * join(vars, " ")
    f_string *= (i > 1 ? " " : "") * term_str
end

for (i, term) in enumerate(f)
    global f_string
    if term == 0
        continue
    end
    if term > 0
        f_string *= " + "
    end
    f_string *= string( term == 1 ? "" : (term == -1 ? "-" : term) ) * "x_" * string(i)
end

f_string *= " + " * string(f0)
fa = f0 - sum(c)


# print fb
fb_string = "fb: "

# Add the terms from `terms`
for (i, term) in enumerate(terms)
    global fb_string
    vars = ["x_" * string(j) for j in term[1]]
    coef = -term[2]
    term_str = (coef == 1 ? "+" : "+" * string(coef)  ) * " (1-" * join(vars, " ") * ")"
    fb_string *= (i > 1 ? " " : "") * term_str
end

for (i, term) in enumerate(f)
    global fb_string
    if term == 0
        continue
    end
    if term > 0
        fb_string *= " + "
    end
    fb_string *= string( term == 1 ? "" : (term == -1 ? "-" : term) ) * "x_" * string(i)
end

# print fc
fc_string = "fc:"

# Add the terms from `terms`
for (i, term) in enumerate(terms)
    global fc_string
    vars = ["x_" * string(j) for j in term[1]]
    coef = -term[2]
    term_str = (coef == 1 ? "+" : "+" * string(coef)  ) * " (1-" * join(vars, " ") * ")"
    fc_string *= (i > 1 ? " " : "") * term_str
end

for (i, term) in enumerate(f)
    global fc_string
    term = max(term, 0)
    if term == 0
        continue
    end
    if term > 0
        fc_string *= " + "
    end
    fc_string *= string( term == 1 ? "" : (term == -1 ? "-" : term) ) * "x_" * string(i)
end

print(f_string, "\n")
print("fa:", fa, "\n")
print(fb_string, "\n")
print(fc_string, "\n")

using IterTools

# Generate all vectors in {0,1}^5
binary_vectors = collect(product(0:1, 0:1, 0:1, 0:1, 0:1))

# Print the binary vectors
minval = Inf
minvec = -1
for (i, vec) in enumerate(binary_vectors)
    global minval
    global minvec
    val = sum(vec[j] * f[j] for j in 1:5) + f0
    val += sum(-c[i] * prod(vec[j] for j in terms[i][1]) for i in 1:3)
    if val < minval
        minval = val
        minvec = i
    end
end
print("minimum and value of f by enumeration:", (minval, binary_vectors[minvec]), "\n")

modelcut = Model()

# Binary variables for the cut: 1 if on source side, 0 if on sink side
@variable(modelcut, edgeleft[1:3] >= 0)
@variable(modelcut, edgemiddle[i in 1:3, j in terms[i][1]] >= 0)
@variable(modelcut, edgeright[j in 1:5] >= 0)

@variable(modelcut, labelleft[1:3] >= 0)  # For left nodes
@variable(modelcut, labelright[1:5] >= 0)  # For right nodes

# Objective: minimize the capacity of the cut
@objective(modelcut, Min, sum(c[i] * edgeleft[i] for i in 1:3) + sum(edgemiddle[i, j] * 1000 for i in 1:3 for j in terms[i][1])
                      +  sum(max(f[j], 0) * edgeright[j] for j in 1:5))

# Cut constraints: if source node is in S, destination must also be in S
@constraints(modelcut, begin
    [j in 2:2], edgeright[j] == 0
    [i in 1:3], edgeleft[i] - 1  +  labelleft[i] >= 0
    [j in 1:5], edgeright[j] - labelright[j] >= 0
    [i in 1:3, j in terms[i][1]], edgemiddle[i, j] - labelleft[i] + labelright[j] >= 0
end)



set_silent(modelcut)
set_optimizer(modelcut, Mosek.Optimizer)
optimize!(modelcut)

print("network cut model for min(fc), min cut value:", objective_bound(modelcut), "\n")
print("label left:", value.(labelleft), "edge left:", value.(value.(edgeleft)),
 "label right:", value.(labelright),  "edge right:", value.(value.(edgeright)),"\n")

model = Model()


@variable(model, flowleft[1:3] >= 0)
@variable(model, flowmiddle[i in 1:3, j in terms[i][1]] >= 0)
@variable(model, flowright[j in 1:5])
@variable(model, shortleft[j in 1:5])

@constraints(model, begin
    [i in 1:3], flowleft[i] <= c[i]
    [j in 1:5], flowright[j] <= f[j]
    [i in 1:3], flowleft[i] == sum(flowmiddle[i, j] for j in  terms[i][1])
    [j in 1:5], shortleft[j] + sum(flowmiddle[i, j] for i in  1:3 if j in terms[i][1]) == flowright[j]
end)


@objective(model, Max,  sum(flowright[j] for j in  1:5) + fa)
set_silent(model)
set_optimizer(model, Mosek.Optimizer)
optimize!(model)
print("network flow model for min(f), max flow value:",objective_bound(model),"\n")
print("network flow model: ", model)
print("network cut model: ", modelcut)

#print(termination_status(model))
