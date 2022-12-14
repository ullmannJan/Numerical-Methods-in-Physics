### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 6221bca8-731c-11ed-29e1-63f63ac51d11
md"""
# Assignments 5

## Problem 4.I
To show: 

$$y''_0 = \frac{-y_2 + 16 y_1 - 30 y_2 + 16 y_3 - y_2}{12 h**2}$$

First we will look at the taylor series of y(x) at h, 2h, -h, -2h

```math
\begin{align}
y(x_0+h) &= y_0 + hy_0' + \frac{h^2}{2}y_0''+ \frac{h^3}{6}y_0'''+ \frac{h^4}{24}y_0^{(4)}\\
y(x_0-h) &= y_0 - hy_0' + \frac{h^2}{2}y_0''- \frac{h^3}{6}y_0'''+ \frac{h^4}{24}y_0^{(4)}\\
y(x_0+2h) &= y_0 + 2hy_0' + 2h^2y_0''+ \frac{4h^3}{3}y_0'''+ \frac{2h^4}{3}y_0^{(4)}\\
y(x_0-2h) &= y_0 - 2hy_0' + 2h^2y_0''- \frac{4h^3}{3}y_0'''+ \frac{2h^4}{3}y_0^{(4)}\\
\end{align}
```

Now we calculate

```math
\begin{align}
&-y(x_0+2h)-y(x_0-2h) + 16 \left(y(x_0+h)+y(x_0-h)\right)\\
&= -2y_0 - 4h^2 y_0'' - \frac{4}{3} h^4 y_0^{(4)} + 16 \left( 2y_0 + h^2y_0'' + \frac{1}{12} h^4 y_0^{(4)} \right)\\
&= 30 y_0 + 12 h^2 y_0'' + 0 h^4y_0^{(4)}\\
\\
\Rightarrow y_0'' &= \frac{-y_2 + 16y_1 - 30 y_2 + 16y_3 - y_2}{12 h^2} + \mathcal{O}(h^6) 
\end{align}
```
"""

# ╔═╡ 0924dac5-08bd-4459-840f-6de0692d95f8
md"""
## Problem 4.II a)
"""

# ╔═╡ a5666d0b-b9ec-48de-a4c5-ad8c6ce4f379
y(x) = atan(x)

# ╔═╡ 7ca77f51-b4d5-431a-a00b-9fc1c6e0ca9a
x_0 = 0

# ╔═╡ 953fbd88-a1c7-46ab-b081-5aad0d9b45ae
h = π/12

# ╔═╡ 9d4b6258-655b-4080-bfbd-9f698b9589b3
y(x_0)

# ╔═╡ 06d05906-d427-43d8-a06d-27fa4bc7274c
function F(h, func, x) 
	return 0.5/h*(func(x+h) - func(x-h))
end

# ╔═╡ 867d9a7d-e8c9-4ba1-ae33-0f6c943932b7
function F_tilde(h, func, x)
	return 4/3*F(h/2, func, x)-F(h, func, x)/3
end

# ╔═╡ c5deb1b5-7b15-4f2f-b38e-e9660941c199
F_tilde(h, y, x_0)

# ╔═╡ d36fc59e-d5fb-485f-8535-402007137bf4
function F_tilde_2(h, func, x)
	return 16/15*F_tilde(h/2, func, x)-F_tilde(h, func, x)/15
end

# ╔═╡ c8620d65-ef48-48e2-b1fd-7a7e13202a2a
F_tilde_2(h, y, x_0)

# ╔═╡ d3b6b2f4-e86a-468b-bace-a105f344f2af
md"""
## Problem 5.I

Simpson 1/3 rule:

$$\int_a^b f(x) dx = \sum_{i=0,i = i+2}^{n-2}\frac{h}{3} (f_i + 4 f_{i+1} + f_{i+2})$$

Romberg integration

```math
\begin{align}
A &= I_1 + c h^n\\
A &= I_2 + c(kh)^n\\
A &= I_2 + \frac{I_2 - I_1}{2^n-1}
\end{align}
```

This means in our case

```math
\begin{align}
I_1 &= \frac{2h}{3} (f_0+4f_2+f_4) + \mathcal{O}(h^4)\\
I_2 &= \frac{h}{3} (f_0+4f1+f_2) + \frac{h}{3} (f_2+4f_3+f_4) + \mathcal{O}(h^4)\\
\Rightarrow A &= I_2 + \frac{I_2 - I_1}{2^n-1}\\
&= I_2 + \frac{h}{45} (-f_0 +4f_1 - 6f_2+4f_3 - f_4)\\
&= \frac{h}{45} (15f_0 - f_0 + 60f_1 +4f_1 +30f_2 - 6f_2 + +60f_3 + 4f_3 +15f_4 - f_4)\\
&= \frac{h}{45} (14f_0 + 64f_1 +24f_2 + 64f_3 +14f_4 )\\
&= \frac{2h}{45} (7f_0 + 32f_1 +12f_2 + 32f_3 +7f_4 )\\
&\Rightarrow \text{4th order interpolating polynomial integration}
\end{align}
```

## Problem 5.II
"""

# ╔═╡ 64e92b6b-ffde-4794-91fb-2b3568a8ff03
function f(x) 
	if x == 0
		return 1
	end
	return 2^x*sin(x)/x
end

# ╔═╡ f456c1b0-db74-4135-bba2-df9786e3ef8f
function simpson(func, a, b, steps)
	x = range(a,b,length=steps)
	h = (b-a)/(steps-1)
	I = 0
	for i in 1:2:(steps-2)
		I += h/3*(func(x[i]) + 4*func(x[i+1]) + func(x[i+2]))
	end
	return I
end

# ╔═╡ 0a6e0e68-0b7a-4a35-b304-6ba620fc7257
simpson(f, 0,1, 4000)

# ╔═╡ 63297be8-d6be-4960-b000-f3bdc8387d89
function romberg_simpson(func, a, b, steps)
	x = range(a,b,length=steps)
	h = (b-a)/(steps-1)
	I = 0
	for i in 1:4:(steps-4)
		I += 2*h/45*(7*func(x[i]) + 32*func(x[i+1]) + 12*func(x[i+2])+ 32*func(x[i+3])+ 7*func(x[i+4]))
	end
	return I
end

# ╔═╡ 7f8ad0ec-39ba-4737-b9a2-5a4fcfb773a3
romberg_simpson(f, 0,1, 4000)

# ╔═╡ 26b46b81-2102-4f21-a2d1-699e6f6c9813
# not working yet

function gaussian_int(func, a::Number, b::Number, steps::Int)
	x = range(a,b,length=steps)
	h = (b-a)/(steps-1)
	
	root_0 = 0.3394810436 *h*0.5
	root_1 = 0.8611363116 *h*0.5
	A_1 = 0.6521451549 *h*0.5
	A_0 = 0.3478548451 *h*0.5

	I = 0
	for i in 1:(steps-1)
		mid = (x[i+1]+x[i])*0.5
		I += A_0*(func(mid+root_0)+func(mid-root_0)) + A_1*(func(mid+root_1)+func(mid-root_1))
	end
	return I
end

# ╔═╡ facf2feb-da9a-4014-bf21-3593a7f15deb
gaussian_int(f, 0,1, 4000)

# ╔═╡ 9c96f1ed-d645-4d4e-a6fe-f65db0b4101f
md"""
## Problem 5.III
"""

# ╔═╡ 572c421f-a682-4e45-a466-f49834d29f9e
f2(x, y) = x*y*y

# ╔═╡ e78077e1-7510-43d0-b812-e4b5fe617dcc
integral(func, a, b) = gaussian_int(func, a, b, Int(trunc(abs(b-a)*1000)+2))

# ╔═╡ ff397ef5-9237-4f33-b2a0-a6cdda055729
2/3

# ╔═╡ b3ce2c4b-c9a5-4503-b129-61ba9726072a
integral(y -> integral(x -> f2(x, y),0,2), 0,1)

# ╔═╡ bd7e4a24-826d-4ec0-b8d0-08cfcd5cbace
4/15

# ╔═╡ a1346fd3-901c-4ef8-b4cf-7a64760bb4af
integral(y -> integral(x -> f2(x, y),y*2,2), 0,1)

# ╔═╡ ded665ff-093a-4eac-9002-8c5a2947c913
integral(x -> integral(y -> f2(x, y),0,x*0.5), 0,2)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╟─6221bca8-731c-11ed-29e1-63f63ac51d11
# ╟─0924dac5-08bd-4459-840f-6de0692d95f8
# ╠═a5666d0b-b9ec-48de-a4c5-ad8c6ce4f379
# ╠═7ca77f51-b4d5-431a-a00b-9fc1c6e0ca9a
# ╠═953fbd88-a1c7-46ab-b081-5aad0d9b45ae
# ╠═9d4b6258-655b-4080-bfbd-9f698b9589b3
# ╠═06d05906-d427-43d8-a06d-27fa4bc7274c
# ╠═867d9a7d-e8c9-4ba1-ae33-0f6c943932b7
# ╠═c5deb1b5-7b15-4f2f-b38e-e9660941c199
# ╠═d36fc59e-d5fb-485f-8535-402007137bf4
# ╠═c8620d65-ef48-48e2-b1fd-7a7e13202a2a
# ╟─d3b6b2f4-e86a-468b-bace-a105f344f2af
# ╠═64e92b6b-ffde-4794-91fb-2b3568a8ff03
# ╠═f456c1b0-db74-4135-bba2-df9786e3ef8f
# ╠═0a6e0e68-0b7a-4a35-b304-6ba620fc7257
# ╠═63297be8-d6be-4960-b000-f3bdc8387d89
# ╠═7f8ad0ec-39ba-4737-b9a2-5a4fcfb773a3
# ╠═26b46b81-2102-4f21-a2d1-699e6f6c9813
# ╠═facf2feb-da9a-4014-bf21-3593a7f15deb
# ╟─9c96f1ed-d645-4d4e-a6fe-f65db0b4101f
# ╠═572c421f-a682-4e45-a466-f49834d29f9e
# ╠═e78077e1-7510-43d0-b812-e4b5fe617dcc
# ╠═ff397ef5-9237-4f33-b2a0-a6cdda055729
# ╠═b3ce2c4b-c9a5-4503-b129-61ba9726072a
# ╠═bd7e4a24-826d-4ec0-b8d0-08cfcd5cbace
# ╠═a1346fd3-901c-4ef8-b4cf-7a64760bb4af
# ╠═ded665ff-093a-4eac-9002-8c5a2947c913
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
