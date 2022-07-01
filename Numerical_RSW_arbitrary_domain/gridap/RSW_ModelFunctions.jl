# Function definitions

function DefineProblem(MeshName)

    @printf("Reading Gmsh discrete model mesh [%s]\n", MeshName)
    model = GmshDiscreteModel(MeshName)
    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2 * order)
    Γ = BoundaryTriangulation(model)
    dΓ = Measure(Γ, 2 * order)
    n_Γ = get_normal_vector(Γ)
    r((n1, n2)) = VectorValue(-n2, n1)       # rotate normal vector
    t_Γ = r ∘ n_Γ       # tangent vector using \circ = ∘ function composition

    # For modes:
    V = TestFESpace(model, reffe, vector_type = Vector{ComplexF64})
    U = TrialFESpace(V)

    # For steady solution:
    V∞ = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "edges")
    ηbc(x) = 0.0       # \eta boundary condition
    U∞ = TrialFESpace(V∞, ηbc)

    # Number of nodes
    n_nodes = size(model.grid.node_coordinates, 1)

    return U, V, U∞, V∞, Ω, dΩ, dΓ, t_Γ, n_nodes
end

function solveModesProblem(F)

    # Sort the eigenvalue/vector results
    function sortEigen()
        p = sortperm(abs.(ω))
        ωsort = ω[p]
        ηmatsort = ηmat[:, p]
        return ωsort, ηmatsort
    end

    # Check eigenvalues and eigenvectors satisfy the equation
    function checkEigen()
        nev = size(ω, 1)
        @printf("Matrix problem size                = [%6d x %6d]\n", size(A0Mat, 1), size(A0Mat, 2))
        @printf("Number of eigenvalues              = [%6d]\n", nev)
        @printf("Eigenvalue range                   = [%6.2e -- %6.2e]\n", minimum(abs.(ω)), maximum(abs.(ω)))
        @printf("Size of eigenvector matrix         = [%6d x %6d]\n", size(η, 1), size(η, 2))
        res = zeros(nev, 1)
        for nn = 1:nev
            this_ω = ω[nn]
            this_η = ηmat[:, nn]
            tmp = (A0Mat + this_ω * A1Mat + (this_ω^2) * A2Mat + (this_ω^3) * A3Mat) * this_η
            res[nn] = maximum(abs.(tmp))
        end # nn
        tmp, ind = findmax(res)
        @printf("Maximum residual for eigenvalue nn = [%6.0d], ω = [%6.3f + %6.3fim], max residual = [%6.3e]\n", ind[1], real(ω[ind]), imag(ω[ind]), tmp)
        return nothing
    end

    @printf("\nSolving modal problem for (ω,η) eigenpairs with F = [%8.4f]:\n", F)
    t = @elapsed begin
        # Weak form of PDE
        A3(ψ, η) = ∫(ψ * η) * dΩ                         # cubic omega term
        A1(ψ, η) = ∫(-∇(ψ) ⋅ ∇(η) - F * ψ * η) * dΩ          # linear omega term
        A0(ψ, η) = ∫(im * sqrt(F) * ψ * (∇(η) ⋅ t_Γ)) * dΓ     # constant term. IS THIS SIGN RIGHT?

        # Assemble matrices explicitly
        A0Mat = assemble_matrix(A0, U, V)  # ω^0 term
        A1Mat = assemble_matrix(A1, U, V)    # ω^1 term
        A2Mat = 0.0 * A1Mat                  # ω^2 term is zero
        A3Mat = assemble_matrix(A3, U, V)    # ω^3 term

        # Solve problem
        nep = PEP([A0Mat, A1Mat, A2Mat, A3Mat])
        #nep = PEP([A0Mat,A1Mat])
        #A3Mat = 0.0*A1Mat                  # ω^3 term
        ω, ηmat = polyeig(nep)         # Similar to MATLAB function except output variables are flipped
        ω, ηmat = sortEigen()
        η, ηnorm = computeNorm_assembleOutput(ω, ηmat)
    end
    @printf("Done in [%6.2f]s.\n\n", t)

    # Diagnostics
    if (checkEigenFlag)
        checkEigen()
    end #if    

    # All done    
    return ω, ηmat, η, ηnorm
end

function computeNorm_assembleOutput(ω, ηmat)
    el2(e) = real(sqrt(sum(∫(e * conj(e)) * dΩ)))       # L2 norm of the FE field e
    nev = size(ω, 1)
    η = Vector{CellField}(undef, nev)
    ηnorm = Vector{Float64}(undef, nev)
    for m = 1:nev
        η[m] = FEFunction(U, ηmat[:, m])     # This is the gridap FE form of the eigenvector
        ηnorm[m] = el2(η[m])
    end
    return η, ηnorm
end

# Trim the eigenvalue/vector results.
function trimEigen(ω, ηmat, η, ηnorm, trimFn)
    p = findall(trimFn, ω)
    ωtrim = ω[p]
    ηmattrim = ηmat[:, p]
    ηtrim = η[p]
    ηnormtrim = ηnorm[p]
    @printf("\nTrimmed [%d] eigenvalues to leave [%d = 2 x %d] significant eigenvalues with [%d] cells.\n", size(ω, 1) - size(ωtrim, 1), size(ωtrim, 1), size(ωtrim, 1) / 2, size(η, 1) / 3)
    @printf("Trimmed eigenvalue range           = [%6.2e -- %6.2e]\n", minimum(abs.(ωtrim)), maximum(abs.(ωtrim)))
    return ωtrim, ηmattrim, ηtrim, ηnormtrim
end

function solveSteadyProblem()
    @printf("\nSolving steady problem for η∞.\n")
    t = @elapsed begin

        # Weak form of PDE
        a(ψ, η) = ∫(∇(ψ) ⋅ ∇(η) + F * ψ * η) * dΩ
        b(ψ) = ∫(F * ψ * ηᵢ) * dΩ

        # Solve problem
        op = AffineFEOperator(a, b, U∞, V∞)
        ls = LUSolver()
        solver = LinearFESolver(ls)
        η∞ = solve(solver, op)
    end
    @printf("Done in [%6.2f]s.\n\n", t)

    # All done    
    return η∞
end

function solveInertialProblem()
    @printf("\nSolving inertial problem for η (this mode has frequency ω = √F, but a mode structure independent of F).\n")
    t = @elapsed begin

        # Weak form of PDE
        A(ψ, η) = ∫(∇(ψ) ⋅ ∇(η)) * dΩ + ∫(im * ψ * (∇(η) ⋅ t_Γ)) * dΓ

        # Assemble matrices explicitly
        AMat = assemble_matrix(A, U, V)

        # Solve problem
        λ, ϕmat = eigs(AMat * AMat', nev = 8, which = :SM)
        ϕ, ϕnorm = computeNorm_assembleOutput(λ, ϕmat)

    end
    @printf("Done in [%6.2f]s.\n\n", t)

    # All done    
    return λ, ϕmat, ϕ, ϕnorm
end

function findExpansionCoefficients(Ein)

    # Setup
    @printf("Finding expansion coefficients with [%d] modes on a mesh with [%d] nodes:\n", size(Ein, 2), size(Ein, 1))
    @printf("This ASSUMES the initial divergence and vorticity vanish.\n")
    SE = size(Ein)
    y = zeros(ComplexF64,SE[1]*3)
    E = zeros(ComplexF64,SE[1]*3,SE[2])
    
    # \eta at t=0
    inds = 1:SE[1]
    fᵢ = interpolate_everywhere(ηᵢ - η∞, U)
    yᵢ = get_free_dof_values(fᵢ)
    y[inds] = yᵢ
    E[inds,:] = Ein

    # d/dt \eta = 0 at t=0 (ASSUMING initial divergence vanishes)
    inds = SE[1]+1:SE[1]*2
    E[inds,:] = Ein * Diagonal(-im * ω)

    # d2/dt2 \eta - \nabla^2= 0 at t=0 (ASSUMING initial vorticity vanishes)
    inds = 2*SE[1]+1:SE[1]*3
    E[inds,:] = Ein * Diagonal(-ω.^2)
    # fᵢ = interpolate_everywhere(∇²(ηᵢ), U)
    fᵢ = interpolate_everywhere(∇⋅∇(ηᵢ), U)
    yᵢ = get_free_dof_values(fᵢ)
    y[inds] = yᵢ
    
    # Find solution
    αo = solve_yeqReEx(E, y)

    return αo
end

function solve_yeqEx(E, y)
    # Compute particular SVD solution to y = E*x
    @printf("Computing eigenvector expansion of initial condition using SVD. Solves y = E*x:\n")
    @printf("size(E) = [%d x %d], size(y) = [%d] => size(x) = [%d]\n", size(E, 1), size(E, 2), size(y, 1), size(E, 2))
    t = @elapsed begin
        # Solve y = E*x
        x, svd_E, p = get_particular_svd_soln(E, y)
        res = (y - E * x)' * (y - E * x)
        @printf("Residual of (y - E*x)'*(y - E*x) = [%7.3e].\n", abs(res))
        #@assert real(res) < 1e-12 "solve_yeqEx: Fail y = Ex"
        plotExpansionCoefficients(x, svd_E, p)
    end
    @printf("Done in [%6.2f]s.\n\n", t)
    return x
end

function solve_yeqReEx(E, y)
    # Compute particular SVD solution to y = Re[ E*x ]
    @printf("Computing eigenvector expansion of initial condition using SVD. Solves y = Re[ E*x ]:\n")
    @printf("size(E) = [%d x %d], size(y) = [%d] => size(x) = [%d]\n", size(E, 1), size(E, 2), size(y, 1), size(E, 2))
    t = @elapsed begin
        # Solve y = Re[ E*x ]
        M = size(E, 2)
        E2 = (hcat(real(E), imag(E)))
        x, svd_E, p = get_particular_svd_soln(E2, y)
        x = x[1:M] - im .* x[M+1:2*M]
        res = (y - real(E * x))' * (y - real(E * x))
        @printf("Residual of (y - real(E*x))'*(y - real(E*x)) = [%7.3e].\n", res)
        #@assert real(res) < 1e-12 "solve_yeqReEx: Fail y = Re[ Ex ]"
        plotExpansionCoefficients(x, svd_E, p)
    end
    @printf("Done in [%6.2f]s.\n\n", t)
    return x
end

function get_particular_svd_soln(E, y)
    threshold = 1.0e-10
    svd_E = svd(E, full = true)

    p = findall(x -> abs.(x) > threshold, svd_E.S)
    MM = size(svd_E.V, 1)
    NN = size(p, 1)
    x = zeros(MM)
    @printf("get_particular_svd_soln: Using [%d] vectors to construct solution with [%d] nullspace vectors.\n", NN, MM - NN)

    # Assemble solution
    for ii = 1:NN
        x = x .+ (svd_E.U[:, ii]' * y / svd_E.S[ii]) .* svd_E.V[:, ii]
    end # ii

    # Test solution:
    # 1. Does this x satisfy y = E*x ?
    #@assert y ≈ E * x "get_particular_svd_soln: 1. Fail y = Ex"

    # 2. Does this x plus a random nullspace contribution satisfy y = E*x ?
    x2 = x
    for ii = NN+1:MM
        x2 = x2 .+ rand(1, 1) .* svd_E.V[:, ii]
    end # ii
    #@assert y ≈ E * x2 "get_particular_svd_soln: 2. Fail null space test!"

    # All done
    return x, svd_E, p
end

function plotExpansionCoefficients(αₙ, svd_E, p)
    fig = Figure(resolution = (600, 1200), font = :sans)
    ax1 = Axis(fig, title = "Modal coefficients", xlabel = "mode", ylabel = "coefficient", yscale = log10)
    scatter!(ax1, abs.(αₙ), color = :black, markersize = 5, label = L"\alpha_n")
    fig[1, 1] = ax1
    ax2 = Axis(fig, xlabel = "mode", ylabel = "eigenvalue ω", title = "Raw eigenvalues", yscale = log10)
    scatter!(ax2, abs.(ωraw), color = :black, markersize = 5, label = L"\omega")
    lines!(ax2, [0; size(ωraw, 1)], [trimThreshold; trimThreshold], label = "Threshold")
    lines!(ax2, [0; size(ωraw, 1)], [sqrt(F); sqrt(F)], label = L"$\sqrt{F}$")
    axislegend(ax2, position = :rb)
    fig[2, 1] = ax2
    ax3 = Axis(fig, xlabel = "eigenvalue ω", ylabel = "coefficient", title = "Coefficients vs. trimmed eigenvalues", xscale = log10, yscale = log10)
    scatter!(ax3, abs.(ω), abs.(αₙ), color = :black, markersize = 5, label = L"\alpha_n")
    lines!(ax3, [sqrt(F); sqrt(F)], [minimum(abs.(αₙ)); maximum(abs.(αₙ))], label = L"$\sqrt{F}$")
    axislegend(ax3, position = :rt)
    fig[3, 1] = ax3
    ax4 = Axis(fig, xlabel = "Index", ylabel = "singular value", title = "Singular values", yscale = log10)
    scatter!(ax4, svd_E.S, color = :grey, markersize = 5)
    scatter!(ax4, svd_E.S[p], color = :black, markersize = 5)
    fig[4, 1] = ax4
    display(fig)
end

function computeη(ω, ηmat, αₙ, thistime)

    thisη = 0.0 * ηmat[:, 1]
    for mm = 1:size(ω, 1)
        timefac = exp(-im * ω[mm] * thistime)
        thisη = thisη + αₙ[mm] * timefac * ηmat[:, mm]
    end # mm

    return thisη
end

function loadη(filename, timesteps)
    @printf("Reading Oceananigans datafile [%s]\n", filename)
    # Grid coordinates at cell centers:
    grid_x, Lx = load_grid(filename, "x", "ᶜᵃᵃ")
    grid_y, Ly = load_grid(filename, "y", "ᵃᶜᵃ")

    # Loop over timesteps
    t_dim = []
    interp_linear = []
    for timestep ∈ timesteps
        fld_key = "timeseries/η/" * string(timestep)
        time_key = "timeseries/t/" * string(timestep)
        η = load(filename, fld_key)
        η = η[:, :, 1]
        this_t_dim = load(filename, time_key)
        push!(t_dim, this_t_dim)
        @printf(" timestep [%d] <=> time = [%8.2f]s\n", timestep, this_t_dim)

        # Build interpolator
        push!(interp_linear, LinearInterpolation((grid_x / Lx, grid_y / Lx), η, extrapolation_bc = Flat()))

        # Plot
        kwargs = (
            fillrange = true,
            levels = range(-2.2, 2.2, length = 64),
            colormap = :seismic
        )
        fig, ax, plt = contourf(grid_x / Lx, grid_y / Lx, η, figure = (resolution = (700, 600),), axis = (title = "Oceananigans η at step " * string(timestep), xlabel = L"x/L_x", ylabel = L"y/L_x",); kwargs...)
        ax.aspect = DataAspect()
        Colorbar(fig[1, 2], plt)
        display(fig)
    end

    return interp_linear, t_dim / timescale
end

function load_grid(filename, type, supscript)
    grid = load(filename, "grid/" * type * supscript)
    grid_L = load(filename, "grid/L" * type)
    grid_N = load(filename, "grid/N" * type)
    grid_N_halo = Int((size(grid, 1) - grid_N) / 2)
    @printf(" load_grid: for [%s]: size(grid) = [%d], grid_N = [%d], halo = [%d]\n", type * supscript, size(grid, 1), grid_N, grid_N_halo)
    grid = grid[grid_N_halo+1:size(grid, 1)-grid_N_halo]
    return grid, grid_L
end

function plotModes(ω, η, ηnorm, fileprefix)
    # Write modes    
    if (writeparaView)
        nmodes = size(η, 1)
        @printf("\nWriting [%d] modes to paraView files [%s].\n", nmodes, fileprefix)
        for ind = 1:nmodes
            writevtk(Ω, "output/$(fileprefix)_$(ind)", cellfields = ["Real" => real(η[ind]), "Imag" => imag(η[ind])])
        end
    end # if

    # Plot frequencies
    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel = "eigenvalue ω", ylabel = "eigenvector norm", title = "Eigenvalues and norms", xscale = log10, yscale = log10)
    scatter!(fig[1, 1], abs.(ω), ηnorm, color = :black, markersize = 5)

    ylimmax = maximum([maximum(ηnorm); maximum(ηnorm)])
    ylimmin = minimum([minimum(ηnorm); minimum(ηnorm)])
    ylimits = [ylimmin; ylimmax]

    if (F > 1e-6)
        lines!(fig[1, 1], sqrt.([F; F]), ylimits, label = L"$\omega = \sqrt{F}$")
    end # if

    for n = 1:4
        lab_text = latexstring("\$ \\omega = $(n) \\pi \$")
        lines!(fig[1, 1], n * pi * [1; 1], ylimits, label = lab_text)
    end # n 

    axislegend(ax1, position = :lt)
    display(fig)

    # Plot modes
    fig = Figure(resolution = (600, 1200), font = :sans)
    for mm = 1:4
        thisω = round(real(ω[mm]), digits = 3)

        tit_txt1 = L"$\Re ( \eta_{%$(mm)} )$, $\omega = %$(thisω)$"
        ax1 = Axis(fig[mm, 1], xlabel = L"x", ylabel = L"y", title = tit_txt1)
        ax1.aspect = DataAspect()
        plot!(Ω, real(η[mm]))
        hidedecorations!(ax1)

        tit_txt2 = L"$\Im ( \eta_{%$(mm)} )$, $\omega = %$(thisω)$"
        ax2 = Axis(fig[mm, 2], xlabel = L"x", ylabel = L"y", title = tit_txt2)
        ax2.aspect = DataAspect()
        plot!(Ω, imag(η[mm]))
        hidedecorations!(ax2)

        tit_txt3 = L"$\left| \eta_{%$(mm)} \right|$, $\omega = %$(thisω)$"
        ax3 = Axis(fig[mm, 3], xlabel = L"x", ylabel = L"y", title = tit_txt3)
        ax3.aspect = DataAspect()
        plot!(Ω, abs(η[mm]))
        hidedecorations!(ax3)
    end # mm

    display(fig)
    return nothing
end

function plotSteadySoln(fileprefix)
    if (writeparaView)
        @printf("Writing steady solution to paraView file [%s].\n\n", fileprefix)
        writevtk(Ω, "output/$(fileprefix)", cellfields = ["eta_infty_re" => real(η∞), "eta_infty_im" => imag(η∞)])
    end # if

    # Plot:
    kwargs = (
        aspect = DataAspect(),
        fill = true,
        levels = 20,
        linewidth = 0,
        colormap = :seismic,
        colorrange = (-1.1, 1.1) .* 2,
    )
    fig = Figure(fontsize = 18)
    Axis(fig[1, 1], xlabel = L"x", ylabel = L"y", title = L"Steady state $\eta_\infty$ solution"; kwargs...)
    plt = plot!(Ω, η∞)
    Colorbar(fig[1, 2], plt)
    display(fig)

    return nothing
end

function plotInertialSoln(fileprefix)
    nmodes = 4
    if (writeparaView)
        @printf("\nWriting [%d] modes to paraView files [%s].\n", nmodes, fileprefix)
        for ind = 1:nmodes
            writevtk(Ω, "output/$(fileprefix)_$(ind)", cellfields = ["Real" => real(ϕ[ind]), "Imag" => imag(ϕ[ind])])
        end
    end # if

    # Plot:
    kwargs = (
        aspect = DataAspect(),
        fill = true,
        levels = 20,
        linewidth = 0,
        colormap = :seismic,
        colorrange = (-1.1, 1.1) .* 2,
    )
    fig = Figure(resolution = (600, 1200), font = :sans)
    for mm = 1:nmodes
        @printf("Singular value (should be zero): [%8.4e + %8.4e im].\n", real(λ[mm]), imag(λ[mm]))
        thisλ = round(real(λ[mm]), digits = 3)

        tit_txt1 = L"$\Re ( \eta_{%$(mm)} )$, $\lambda = %$(thisλ)$"
        ax1 = Axis(fig[mm, 1], xlabel = L"x", ylabel = L"y", title = tit_txt1)
        ax1.aspect = DataAspect()
        plot!(Ω, real(ϕ[mm]))
        hidedecorations!(ax1)

        tit_txt2 = L"$\Im ( \eta_{%$(mm)} )$, $\lambda = %$(thisλ)$"
        ax2 = Axis(fig[mm, 2], xlabel = L"x", ylabel = L"y", title = tit_txt2)
        ax2.aspect = DataAspect()
        plot!(Ω, imag(ϕ[mm]))
        hidedecorations!(ax2)

        tit_txt3 = L"$\left| \eta_{%$(mm)} \right|$, $\lambda = %$(thisλ)$"
        ax3 = Axis(fig[mm, 3], xlabel = L"x", ylabel = L"y", title = tit_txt3)
        ax3.aspect = DataAspect()
        plot!(Ω, abs(ϕ[mm]))
        hidedecorations!(ax3)
    end # mm

    display(fig)

    return nothing
end

function animateSolution(nsteps, Δt, Tf, ω, ηmat, αₙ, fileprefix)
    t = @elapsed begin
        videofile = "output/$(fileprefix)_results.mp4"
        @printf("Making animation of time dependent solution:\nWriting output to [%s] with [%d] steps of length [%8.6f] and final time = [%6.4f].\n", videofile, nsteps, Δt, Tf)

        kwargs = (
            fill = true,
            levels = 20,
            linewidth = 1,
            colormap = :seismic,
            colorrange = (-1.1, 1.1) .* 2,
        )

        # Define observables:
        thistime = Observable(0.0)
        u = lift(thistime) do thistime
            thisη = computeη(ω, ηmat, αₙ, thistime)
            ηₒ = FEFunction(U, thisη)     # This is the gridap FE form of the eigenvector
            u = real(ηₒ) + real(η∞)
            #u = ηₒ + η∞
        end

        tit_txt = lift(thistime) do thistime
            tit_txt = @sprintf("η at t = %4.2f = %3.1f hours", thistime, thistime * timescale / 3600)
        end

        # Make plot:
        fig, ax, plt = plot(Ω, u, figure = (resolution = (700, 600),), axis = (title = tit_txt, xlabel = L"x/L_x", ylabel = L"y/L_x",); kwargs...)
        ax.aspect = DataAspect()
        Colorbar(fig[1, 2], plt)

        # Control animation:
        framerate = 30
        timestamps = Tf * collect(1:nsteps) / nsteps
        record(fig, videofile, timestamps; framerate = framerate) do this_t
            thistime[] = this_t
        end
    end
    @printf("Done in [%6.2f]s.\n\n", t)
    return nothing

end