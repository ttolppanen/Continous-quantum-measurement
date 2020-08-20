using LinearAlgebra
function te(v, k)
    k*sum(v)
end

v = [[[1,2,3], [3,4], [3,4]],[[1,9,3], [3,4], [3,4]]]
a = [1,2]
b = [1,2]
show(te.(v[2], 2))
