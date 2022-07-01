using LinearAlgebra, Printf

function get_particular_svd_soln(E,y)
    threshold = 1.0e-10 
    svd_E = svd(E,full=true)

    p = findall(x -> abs.(x) > threshold,svd_E.S)
    MM = size(svd_E.V,1)
    NN = size(p,1)
    x = zeros(MM)       # Generates a (column) *vector* of MMx1 zeros (unlike MATLAB!)

    # Assemble solution
    for ii = 1:NN
       x = x .+ (svd_E.U[:,ii]' * y / svd_E.S[ii]) .* svd_E.V[:,ii]
    end # ii

    # Test solution:
    Base.@assert y ≈ E*x "Fail y = Ex"

    x2 = x
    for ii = NN+1:MM
        x2 = x2 .+ rand(1,1) .* svd_E.V[:,ii]
     end # ii
    Base.@assert y ≈ E*x2 "Fail null space test!"

    # All done
    return x
end

N = 8
M = N + 4
y = rand(N,1)
E = rand(N,M) + im*rand(N,M)

x = get_particular_svd_soln(E,y)

E2 = (hcat(real(E),imag(E)))
x2 = get_particular_svd_soln(E2,y)
x2 = x2[1:M] + im.*x2[M+1:2*M]

@printf("Solution to y = E*x:\n ")
print(x)

@printf("\n\nSolution to y = Re[ E*x ]:\n ")
print(x2)
