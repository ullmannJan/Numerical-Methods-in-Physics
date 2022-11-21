### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ ab68a0f9-32bf-400a-85dd-45092730658e
using LinearAlgebra

# ╔═╡ f4bd1716-61d4-446a-8a33-40b90983af48
using Statistics

# ╔═╡ d5fbd469-1230-4de5-8a9c-e4d58409deca
md"""
# Assignment 3 - Jan Ullmann

## Problem I.1
"""

# ╔═╡ 80b9b02e-67a9-4824-858a-de5759f6351c
function matrix_mult(matrix, vector)
	@assert size(matrix,2) == size(vector,1) "no valid size of vector and matrix"
	output = Vector(undef,size(matrix,1))
	for r ∈ 1:size(matrix,1)
		sums = 0
    	for c ∈ 1:size(vector,1)
			sums += matrix[r,c]*vector[c]
		end
		output[r] = sums
	end
	return output
end

# ╔═╡ 75600936-2ae0-462a-be14-1414f7fce220
m = [1 2 3; 4 5 6; 7 8 9]

# ╔═╡ b9a734b1-99d7-4cb6-9bc7-e81b11c736ed
b = [10:12;]

# ╔═╡ 9385c80a-567b-4267-8394-3f6d262e5097
m*b

# ╔═╡ ee1e4521-fb2a-4222-bb27-e87983445837
matrix_mult(m,b)

# ╔═╡ ec45c4c0-c223-4c67-8dcf-36f2ae8be76a
md"""
## Problem I.2
"""

# ╔═╡ 639dd38a-e987-475f-91e3-ef8763c0ea7f
function solve_triang(matrix::UpperTriangular, vector::Vector; eps::Number=1e-14)
	sol = zeros(size(vector))
	for i ∈ size(vector,1):-1:1
		sums = 0
		for k ∈ (i+1):size(vector,1)
			sums += matrix[i,k]*sol[k]
		end
		if matrix[i,i] != 0
			sol[i] = (vector[i]-sums)/matrix[i,i]
		else
			if vector[i] != sums
				return "This system of equations has no solution"
			else
				return "This system of equations has more than one solution"
			end
		end
	end
	if norm(matrix*sol - vector) < eps
		return sol
	else
		return "No solution found"
	end
end

# ╔═╡ e64ec780-97f3-4681-9d27-11a22296f6f5
m2 = UpperTriangular(m)

# ╔═╡ 4c7c99fb-2cef-4c46-9211-6d0aa76f795b
solve_triang(m2, b)

# ╔═╡ 17f8f021-46c4-436e-8ec7-63a642396634
solve_triang(UpperTriangular([1 2 3; 0 0 2; 0 0 2]),[1,2,1])

# ╔═╡ 6bdf7099-cac0-4862-822d-225a08fdc577
md"""
## Problem II.1
"""

# ╔═╡ 498334b5-f237-4517-88cc-3ed1966beffc
function solve_gauss_wo_pivoting(matrix::Matrix, vector::Vector)
	sys = hcat(matrix, vector)
	for i ∈ 1:size(vector,1)
		for j in (i+1):size(vector,1)
			sys[j,:] = sys[j,:] - sys[i,:] / sys[i,i] * sys[j,i] 
		end
	end
	tri = UpperTriangular(sys[:,1:end-1])
	vector_new = sys[:,end]
	display(sys)
	return solve_triang(tri, vector_new)
end

# ╔═╡ 728071c0-563c-41c9-a0e5-130f6cc30dea
md"""
## Problem II.2
"""

# ╔═╡ 3e9d433b-3783-40b8-aa76-5785b184b3f2
m3 = [2.0 0.1 −0.2
0.05 4.2 0.032
0.12 −0.07 5.0]

# ╔═╡ b6a1b4f2-dd2e-4951-b1fa-758117c1cf50
s = solve_gauss_wo_pivoting(m3, b)

# ╔═╡ 2c7c854b-82e2-4f41-aca0-7c6e233a454e
m3*s

# ╔═╡ ca77e3ca-64a5-407b-a461-284455325687
md"""
## Problem II.3
"""

# ╔═╡ ea9d5e85-39f7-404e-836e-8e98beae0297
m4 = [1 1 0
2 2 −2
0 3 15]

# ╔═╡ a7f819e4-2aee-4db9-ba73-e5da67279bde
b2 = [1, -2, 33]

# ╔═╡ 961bda15-b22d-47d1-b4ae-635595ac851f
solve_gauss_wo_pivoting(m4, b2)

# ╔═╡ 92326547-87c4-4743-b818-65166d7dfdb5
md"""
This error makes sense as we are dividing by 0. Without pivoting we do not get a triangular matrix, because the second row turns to 0 before we can subtract it from the third. Thus, we have to reorder the rows first before putting it into the triangular shape.
"""

# ╔═╡ d38bf061-f0b2-4235-97d0-3c02b080af26
function solve_gauss(matrix::Matrix, vector::Vector)
	sys = hcat(matrix, vector)
	for i ∈ 1:size(vector,1)
		if sys[i,i] == 0
			#find correct shift value for swapping
			shift=1
			while sys[i+shift,i] == 0
				@assert	shift < size(vector,1)-i "system of equations not solvable"
				shift +=1
			end
			#swapping rows
			temp = sys[i,:]
			sys[i,:] = sys[i+shift,:]
			sys[i+shift,:] = temp
		end
		for j in (i+1):size(vector,1)
			sys[j,:] = sys[j,:] - sys[i,:] / sys[i,i] * sys[j,i] 
		end
	
	end
	tri = UpperTriangular(sys[:,1:end-1])
	vector_new = sys[:,end]
	display(sys)
	return solve_triang(tri, vector_new)
end

# ╔═╡ 98d8d4ba-261f-4f04-adba-9d8e76f29be0
s2 = solve_gauss(m4, b2)

# ╔═╡ 9ee68f85-1663-4a74-888c-059e59b718c1
m4*s2

# ╔═╡ 2a19ce6c-e767-4f4b-b3fc-0246a99365e4
md"""
## Problem III.1
"""

# ╔═╡ dcba17dc-0f46-4aa3-b7c9-0162edc8e3c3
m5 = [6 5 −5
2 6 −2
2 5 −1]

# ╔═╡ dd7b92ef-8cb5-42e5-8587-ea53bd69cd5d
eigen(m5)

# ╔═╡ 3159490f-d2e5-473b-81ac-8d598cdd1826
function power_iteration(matrix::AbstractMatrix;initial_vector::AbstractVector=[1,0.5,0.5], eps::Number=1e-6, iter_max::Integer=100)
	iteration = 0
	v_old = matrix*initial_vector
	v_new = matrix*v_old
	e_val = 0
	e_val_new = 1
	while abs(e_val_new-e_val) > eps && iteration < iter_max
		iteration += 1
		v_old = v_new
		v_new = matrix*v_new
		e_val = e_val_new
		e_val_new = mean(v_new./v_old)
	end
	return e_val_new, normalize!(v_new), iteration
end

# ╔═╡ b2f4559d-151f-4f12-89f5-328708a28726
power_iteration(m5)

# ╔═╡ d0128d19-40f3-4815-ae93-38e2bcd8bc7e
function power_iteration_improved(matrix::AbstractMatrix;initial_vector::AbstractVector=[1,0.5,0.5], eps::Number=1e-6, iter_max::Integer=100)
	iteration = 0
	v_old = matrix*initial_vector
	v_new = matrix*v_old
	e_val_old = -1
	e_val = 0
	e_val_new = 1
	while abs(e_val_new-e_val) > eps && iteration < iter_max
		iteration += 1
		v_old = v_new
		v_new = matrix*v_new
		e_val_old = e_val
		e_val = e_val_new
		e_val_new = mean(v_new./v_old)
	end
	e_val = (e_val_old * e_val_new - e_val*e_val)/(e_val_new - 2*e_val + e_val_old)
	print(e_val)
	return e_val, normalize!(v_new), iteration
end

# ╔═╡ 7c8b9e93-1dd5-41e4-988a-625a70cfa4be
power_iteration_improved(m5)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "8c3efc2a8bec7327ac93e0344d9faedc429a153e"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╟─d5fbd469-1230-4de5-8a9c-e4d58409deca
# ╠═ab68a0f9-32bf-400a-85dd-45092730658e
# ╠═f4bd1716-61d4-446a-8a33-40b90983af48
# ╠═80b9b02e-67a9-4824-858a-de5759f6351c
# ╠═75600936-2ae0-462a-be14-1414f7fce220
# ╠═b9a734b1-99d7-4cb6-9bc7-e81b11c736ed
# ╠═9385c80a-567b-4267-8394-3f6d262e5097
# ╠═ee1e4521-fb2a-4222-bb27-e87983445837
# ╟─ec45c4c0-c223-4c67-8dcf-36f2ae8be76a
# ╠═639dd38a-e987-475f-91e3-ef8763c0ea7f
# ╠═e64ec780-97f3-4681-9d27-11a22296f6f5
# ╠═4c7c99fb-2cef-4c46-9211-6d0aa76f795b
# ╠═17f8f021-46c4-436e-8ec7-63a642396634
# ╟─6bdf7099-cac0-4862-822d-225a08fdc577
# ╠═498334b5-f237-4517-88cc-3ed1966beffc
# ╟─728071c0-563c-41c9-a0e5-130f6cc30dea
# ╠═3e9d433b-3783-40b8-aa76-5785b184b3f2
# ╠═b6a1b4f2-dd2e-4951-b1fa-758117c1cf50
# ╠═2c7c854b-82e2-4f41-aca0-7c6e233a454e
# ╟─ca77e3ca-64a5-407b-a461-284455325687
# ╠═ea9d5e85-39f7-404e-836e-8e98beae0297
# ╠═a7f819e4-2aee-4db9-ba73-e5da67279bde
# ╠═961bda15-b22d-47d1-b4ae-635595ac851f
# ╟─92326547-87c4-4743-b818-65166d7dfdb5
# ╠═d38bf061-f0b2-4235-97d0-3c02b080af26
# ╠═98d8d4ba-261f-4f04-adba-9d8e76f29be0
# ╠═9ee68f85-1663-4a74-888c-059e59b718c1
# ╟─2a19ce6c-e767-4f4b-b3fc-0246a99365e4
# ╠═dcba17dc-0f46-4aa3-b7c9-0162edc8e3c3
# ╠═dd7b92ef-8cb5-42e5-8587-ea53bd69cd5d
# ╠═3159490f-d2e5-473b-81ac-8d598cdd1826
# ╠═b2f4559d-151f-4f12-89f5-328708a28726
# ╠═d0128d19-40f3-4815-ae93-38e2bcd8bc7e
# ╠═7c8b9e93-1dd5-41e4-988a-625a70cfa4be
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
