function _swap_degree!(Y::BinaryNetwork)
   iy = interactions(Y)
   i1, i2 = StatsBase.sample(iy, 2, replace=false)
   n1, n2 = (from=i1.from, to=i2.to), (from=i2.from, to=i1.to)

   while (n2 ∈ iy)|(n1 ∈ iy)
      i1, i2 = StatsBase.sample(iy, 2, replace=false)
      n1, n2 = (from=i1.from, to=i2.to), (from=i2.from, to=i1.to)
   end

   Y[i1.from, i1.to] = false
   Y[i2.from, i2.to] = false
   Y[n1.from, n1.to] = true
   Y[n2.from, n2.to] = true
end

function _swap_fill!(Y::BinaryNetwork)
   iy = interactions(Y)
   i1 = StatsBase.sample(iy)
   n1 = (from=StatsBase.sample(species(Y; dims=1)), to=StatsBase.sample(species(Y; dims=2)))

   while (n1 ∈ iy)
      n1 = (from=StatsBase.sample(species(Y; dims=1)), to=StatsBase.sample(species(Y; dims=2)))
   end

   Y[i1.from, i1.to] = false
   Y[n1.from, n1.to] = true
end

function _swap_vulnerability!(Y::BinaryNetwork)
   iy = interactions(Y)
   i1 = StatsBase.sample(iy)
   n1 = (from=StatsBase.sample(species(Y; dims=1)), to=i1.to)

   while (n1 ∈ iy)
      n1 = (from=StatsBase.sample(species(Y; dims=1)), to=i1.to)
   end

   Y[i1.from, i1.to] = false
   Y[n1.from, n1.to] = true
end

function _swap_generality!(Y::BinaryNetwork)
   iy = interactions(Y)
   i1 = StatsBase.sample(iy)
   n1 = (from=i1.from, to=StatsBase.sample(species(Y; dims=2)))

   while (n1 ∈ iy)
      n1 = (from=i1.from, to=StatsBase.sample(species(Y; dims=2)))
   end

   Y[i1.from, i1.to] = false
   Y[n1.from, n1.to] = true
end


"""
    shuffle(N::BinaryNetwork; constraint::Symbol=:degree)

Return a shuffled copy of the network (the original network is not modified).
See `shuffle!` for a documentation of the keyword arguments.

#### References

- Fortuna, M.A., Stouffer, D.B., Olesen, J.M., Jordano, P., Mouillot, D.,
  Krasnov, B.R., Poulin, R., Bascompte, J., 2010. Nestedness versus modularity
  in ecological networks: two sides of the same coin? Journal of Animal Ecology
  78, 811–817. https://doi.org/10.1111/j.1365-2656.2010.01688.x
"""
function Random.shuffle(N::BinaryNetwork; constraint::Symbol=:degree)
   Y = copy(N)
   shuffle!(Y; constraint=constraint)
   return Y
end

"""
    shuffle!(N::BinaryNetwork; constraint::Symbol=:degree)

Shuffles interactions inside a network (the network is *modified*), under the
following `constraint`:

- :degree, which keeps the degree distribution intact
- :generality, which keeps the out-degree distribution intact
- :vulnerability, which keeps the in-degree distribution intact
- :fill, which moves interactions around freely

The function will take two interactions, and swap the species establishing them.
By repeating the process a large enough number of times, the resulting network
should be relatively random. Note that this function will conserve the degree
(when appropriate under the selected constraint) of *every* species. Calling the
function will perform **a single** shuffle. If you want to repeat the shuffling
a large enough number of times, you can use something like:

    [shuffle!(n) for i in 1:10_000]

If the keyword arguments are invalid, the function will throw an
`ArgumentError`.

#### References

- Fortuna, M.A., Stouffer, D.B., Olesen, J.M., Jordano, P., Mouillot, D.,
  Krasnov, B.R., Poulin, R., Bascompte, J., 2010. Nestedness versus modularity
  in ecological networks: two sides of the same coin? Journal of Animal Ecology
  78, 811–817. https://doi.org/10.1111/j.1365-2656.2010.01688.x
"""
function Random.shuffle!(N::BinaryNetwork; constraint::Symbol=:degree)
   constraint ∈ [:degree, :generality, :vulnerability, :fill] || throw(ArgumentError("The constraint argument you specificied ($(constraint)) is invalid -- see ?shuffle! for a list."))

   f = EcologicalNetworks._swap_degree!
   constraint == :generality && (f = EcologicalNetworks._swap_generality!)
   constraint == :vulnerability && (f = EcologicalNetworks._swap_vulnerability!)
   constraint == :fill && (f = EcologicalNetworks._swap_fill!)

   f(N)
end
