using JuMP
using MosekTools

# f = - c1 x2 x3 - c2 x1 x3 x4 - c3 x3 x5 +  \sum_{j} fj xj + f0

f0 = 7
c = [1,2,5]
f = [0,1,1,-1,0]
terms = [([2,3], -c[1]), ([1,3,4], -c[2]), ([3,5], -c[3])]

# print f
f_string = ""

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
fb_string = ""

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
fc_string = ""

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
print(fa, "\n")
print(fb_string, "\n")
print(fc_string, "\n")

model = Model()


@variable(model, rholeft[1:3] >= 0)
@variable(model, rhomiddle[i in 1:3, j in terms[i][1]] >= 0)
@variable(model, rhoright[j in 1:5])
@variable(model, shortleft[j in 1:5])

@constraints(model, begin
    [i in 1:3], rholeft[i] <= c[i]
    [j in 1:5], rhoright[j] <= f[j]
    [i in 1:3], rholeft[i] == sum(rhomiddle[i, j] for j in  terms[i][1])
    [j in 1:5], shortleft[j] + sum(rhomiddle[i, j] for i in  1:3 if j in terms[i][1]) == rhoright[j]
end)


@objective(model, Max,  sum(rhoright[j] for j in  1:5) + fa)

set_optimizer(model, Mosek.Optimizer)
optimize!(model)
print(objective_bound(model))
print(value.(rhoright), value.(shortleft))
print(model)
print(termination_status(model))


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
print((minval, binary_vectors[minvec]))
