### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ df165d16-bf34-11ec-3a94-bff9ab323db7
"""
    orbit_plain(S, x)
Compute the orbit of `x` under the action of a group `G` generated by set `S`.

It is assumed that elements of `G` act on `x` via operation `^`.

### Input
 * `x` - a point 
 * `S` - a finite generating set for `G = ⟨S⟩`
### Output
 * `{xᵍ | g ∈ G}` - the orbit of `x` under the action of `G`
"""
function orbit_plain(x, S)
	@assert !isempty(S) # groups need generators
	Δ = [x]
	for δ in Δ
		for s in S
			γ = δ^s
			if γ ∉ Δ # the expensive check
				push!(Δ, γ)
			end
		end
	end
	return Δ
end

# ╔═╡ f763036d-18d0-4635-af2e-104007bff48f
md"
But before we can actually call it we need to implement some permutation groups.
Let's begin then!

# Permutations
Each permutation is defined by the images of points under the natural action of permutations on set `1..n`. We could define `Permutation`s as follows:
```julia
struct Permutation
    images::Vector{Int}
end
```
Here we make promise that images of permutations will be stored as vectors of `Int`s and `n = length(images)` determines group `Sym(n)` a permutation belongs to. 

One may ask if `UInt`s would be better, or if creating a parametrized type
```julia
struct Permutation{T}
    images::Vector{T}
end
```
would be more flexible, but we let those questions rest for now.

The important thing now is to set conventions for permutations, namely:
* if ``\sigma \in Sym(n)``, then `σ.images = [σ(1), …, σ(n)]`, and
* ``(1,2)\cdot(2,3) = (1,3,2)``, according to our convention of **acting from right**.

According to those conventions we represent ``(1,2) \in S_3`` as `Permutation` with images `[2,1,3]`.

----

Before implementing this we need to decide how to access `σ(i)`. There are two possibilities here:
* elements of `Vectors`, `Matrices` etc. are accessed through `v[n]` notation, i.e. if
`images = [2,1,3]`, then `v[1] == 2`, `v[2] == 1` and `v[3] == 3` (all other accesses will result in `BoundsError`).
* functions are called simply by `σ(n)` (remember -- permutations are just bijections!)

Here we will make an assumption that outside of its images we each bijection is the identity, that is if `σ` is defined by `[σ(1), …, σ(n)]`, then this permutation (as a map) is the identity elsewhere. To say it differently we maintain implicit embedings 

$$\operatorname{Sym}(n) → \operatorname{Sym}(n+1)$$

by extending each permutation ``\sigma \in \operatorname{Sym}(n)`` to ``\sigma' \in  \operatorname{Sym}(n+1)`` so that ``\sigma'(n+1) = n+1``.

----
"

# ╔═╡ b8d09e00-eece-45f2-bbf1-d244e12e38c5
md"
### Technical Note
On a more technical note we should probably add some checks that the vector of images indeed defines a permutation (protect from shooting ourselves in the foot :). Since such checks are usually quite expensive we will also hide them behind a \"flag\".

This way if we're absolutely sure that `v` is a vectors of permutation images we can simply call `Permutation(v, false)` to skip the expensive check.
"

# ╔═╡ b8d7bbe9-7066-47b5-b002-a6b6500bc2a6
begin
	# struct definition
	struct Permutation
		images::Vector{Int}
		
		function Permutation(v::AbstractVector{<:Integer}, check=true)
			if check
				@assert sort(v) == 1:length(v) "Image vector doesn't define a permutation"
			end
			return new(v) # calls convert(Vector{Int}, v)
		end
	end
	# accessing images of a permutation
	function (σ::Permutation)(n::Integer)
		if n > length(σ.images)
			return convert(Int, n)
		else
			return σ.images[n]
		end
	end
end

# ╔═╡ a1b8b8dc-2137-49de-a5a3-e40af09d152c
σ = Permutation(BigInt[2,1,3])

# ╔═╡ 3c5b4a0b-b883-44b1-8822-671e0ff36c70
@code_warntype σ(UInt(4))

# ╔═╡ 07ca94be-1f15-46b8-98fe-8f7c6ddcfa1f
Permutation([1,2,4])

# ╔═╡ a0991f77-c787-48bd-896f-371a4378d181
(σ(1), σ(2), σ(3), σ(6))

# ╔═╡ a8b15ba7-6d76-42c8-82f1-e4408c2ff718
md"
Finally we need a method to determine the (possible) support of a permutation, i.e. ``n`` such that ``σ`` actually belongs to ``\operatorname{Sym}(n)``. We will call this function (somewhat non-canonically) the `degree` of a permutation.
"

# ╔═╡ 5f323c62-576b-4707-997b-1d4b79cdb029
function degree(σ::Permutation)
	return length(σ.images)
end

# ╔═╡ 0bfdfb23-7132-4314-8b38-86c3b9c7a25f
md"
>**Exercise 1**: Enhance `degree` so that it will always return the minimal `n`.
"

# ╔═╡ c567c037-e60e-4f4e-9916-f7dd6b24528f
function Base.:(==)(σ::Permutation, τ::Permutation)
	deg = max(degree(σ), degree(τ))

	for i in 1:deg
		if σ(i) ≠ τ(i)
			return false
		end
	end
	return true
end

# ╔═╡ c073f1df-b6a6-440c-9385-19a4e1e19a20
md"
## Group Structure
Of course every group needs three things:
* the identity element -- we will implement `Base.one`
* inversion -- we will implement `Base.inv`
* binary operation -- we will use `Base.:*` for this
"

# ╔═╡ 35c9fa94-1cb0-4528-8dcb-6c1bf03b4a03
md"
>**Exercise 2**:
>Using syntax `σ(i)` and it `degree` **only** implement:
> * `Base.one(σ::Permutation)`
> * `Base.inv(σ::Permutation)` 
> * `Base.:*(σ::Permutation, τ::Permutation)`
> (Remember **not to** access `σ.images` directly!)
"

# ╔═╡ db9c7876-427a-4e45-8636-4098ec4d2f69
function Base.one(σ::Permutation)
   
end

# ╔═╡ 0a12688b-6684-4905-9c33-84948db6eb83
@assert degree(one(σ)) == 1

# ╔═╡ b49dd619-e944-4507-8c38-2a1abec02303
function Base.inv(σ::Permutation)

end

# ╔═╡ b7d1dcd3-d0d2-47d4-ad01-3a9efb40ced3
begin
	@assert inv(σ) == one(σ)
	@assert inv(Permutation([2,3,1])) == Permutation([3,1,2])
end

# ╔═╡ 419ab7e8-a2f1-4d4a-a4bc-6b48991a3760
function Base.:(*)(σ::Permutation, τ::Permutation)

end

# ╔═╡ dd280601-8ec7-46a4-9583-074027937c5a
τ = Permutation([1,3,2])

# ╔═╡ 48d0b6ef-a739-43f7-83e7-11f457c4d40a
σ * τ

# ╔═╡ a6bc98a0-4d3c-442e-b5a5-c5edfbffb10c
let σ = Permutation([2,1,3]), τ = Permutation([1,3,2])
	@assert inv(one(σ)) == one(σ)
	@assert inv(σ)*σ == one(σ)
	@assert τ*inv(τ) == one(τ)
	@assert inv(σ*τ) == inv(τ)*inv(σ)
	# (1,2)·(2,3) == (1,3,2)
	@assert σ*τ == Permutation([3,1,2])
end

# ╔═╡ 387806a7-3a6f-4961-9198-480137ca69de
md"
## Cycle decomposition

We would like to display permutations in more human-friendly manner, namely as a cycle decomposition, i.e. instead

	Permutation([3, 2, 1, 4])

we would like to see

	(1,3,2)(4)

(possibly even omitting cycles of length `1`). But how to achieve it? How do we compute the cycles? We start with a fixed number and observe where does it travel under repeated applications of our permutation... That's just computing decomposing `1:degree(σ)` into orbits under the action of cyclic group generated by `σ`!

Since we already have the code for orbit computation lying around let's use it! We need just to define the `^`-action of permutations on integers:
"

# ╔═╡ 0ca057c2-220a-4054-ac9d-7adb4db5ef22
σ

# ╔═╡ 2619d02f-cec6-4e94-910e-84dd85be9af3
# Base.:^(i::Integer, σ::Permutation) = σ(i)

# ╔═╡ 007d963d-c7cc-4360-a9a9-aa4a8bf5576a
orbit_plain(1, [σ])

# ╔═╡ 83799213-014f-4bb7-9d84-07a68af3e767
function cycle_decomposition(σ::Permutation)
	visited = falses(degree(σ))
	cycles = Vector{Vector{Int}}()
	# each cycle will be a Vector{Int} and we have a whole bunch of them
	for i in 1:degree(σ)
		if visited[i]
			# if we have already seen this point there is no point in computing
			# the same orbit twise
			continue # i.e. skip the rest of the body and continue with i+1
		end
		Δ = orbit_plain([σ], i)
		visited[Δ] .= true # modify the `visited` along the whole orbit
		push!(cycles, Δ) # add obtained orbit to cycles
	end
	return cycles
end

# ╔═╡ 10f67cdd-706c-472c-8e18-7638d977499f
cycle_decomposition(σ)

# ╔═╡ ea2a53f4-0738-41d9-9233-a8060c60def8
cycle_decomposition(τ)

# ╔═╡ bb405a6a-4117-4734-8076-b2cd76330364
cycle_decomposition(σ*τ)

# ╔═╡ 557de2a4-e366-48c2-9ade-485c76c888c6
cycle_decomposition(one(σ))

# ╔═╡ ae7cb768-6882-46a5-bfbd-af6ce1579f85
md"
Finally we can define nice visuals for `Permutations` by neatly printing its cycles. To do so we will juse `join` function. Have a look at its help to understand it better.
"

# ╔═╡ aa46864a-c656-448c-8896-f62b485bfa4c
function Base.show(io::IO, σ::Permutation)
	if isone(σ)
		print(io, "()")
	end
	
	for cycle in cycle_decomposition(σ)
		if length(cycle) == 1
			continue
		else
			print(io, "(")
			join(io, cycle, ",")
			print(io, ")")
		end
	end
end

# ╔═╡ 447a4f69-5c9a-499f-8f02-917fab664d00
σ

# ╔═╡ 0860e4c4-9075-4ce6-ac20-e549dda125ff
τ

# ╔═╡ 9dc57b08-ac07-4c1a-bee4-2d9be075ac50
σ*τ

# ╔═╡ ced3a00e-64db-4f13-86a9-8ef5c0856ce0
inv(σ)*σ

# ╔═╡ 68adc013-5ac1-4793-82ed-4e1cc03f9698


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╠═df165d16-bf34-11ec-3a94-bff9ab323db7
# ╟─f763036d-18d0-4635-af2e-104007bff48f
# ╟─b8d09e00-eece-45f2-bbf1-d244e12e38c5
# ╠═b8d7bbe9-7066-47b5-b002-a6b6500bc2a6
# ╠═a1b8b8dc-2137-49de-a5a3-e40af09d152c
# ╠═3c5b4a0b-b883-44b1-8822-671e0ff36c70
# ╠═07ca94be-1f15-46b8-98fe-8f7c6ddcfa1f
# ╠═a0991f77-c787-48bd-896f-371a4378d181
# ╟─a8b15ba7-6d76-42c8-82f1-e4408c2ff718
# ╠═5f323c62-576b-4707-997b-1d4b79cdb029
# ╟─0bfdfb23-7132-4314-8b38-86c3b9c7a25f
# ╠═c567c037-e60e-4f4e-9916-f7dd6b24528f
# ╟─c073f1df-b6a6-440c-9385-19a4e1e19a20
# ╟─35c9fa94-1cb0-4528-8dcb-6c1bf03b4a03
# ╠═db9c7876-427a-4e45-8636-4098ec4d2f69
# ╠═0a12688b-6684-4905-9c33-84948db6eb83
# ╠═b49dd619-e944-4507-8c38-2a1abec02303
# ╠═b7d1dcd3-d0d2-47d4-ad01-3a9efb40ced3
# ╠═419ab7e8-a2f1-4d4a-a4bc-6b48991a3760
# ╠═dd280601-8ec7-46a4-9583-074027937c5a
# ╠═48d0b6ef-a739-43f7-83e7-11f457c4d40a
# ╠═a6bc98a0-4d3c-442e-b5a5-c5edfbffb10c
# ╟─387806a7-3a6f-4961-9198-480137ca69de
# ╠═0ca057c2-220a-4054-ac9d-7adb4db5ef22
# ╠═2619d02f-cec6-4e94-910e-84dd85be9af3
# ╠═007d963d-c7cc-4360-a9a9-aa4a8bf5576a
# ╠═83799213-014f-4bb7-9d84-07a68af3e767
# ╠═10f67cdd-706c-472c-8e18-7638d977499f
# ╠═ea2a53f4-0738-41d9-9233-a8060c60def8
# ╠═bb405a6a-4117-4734-8076-b2cd76330364
# ╠═557de2a4-e366-48c2-9ade-485c76c888c6
# ╟─ae7cb768-6882-46a5-bfbd-af6ce1579f85
# ╠═aa46864a-c656-448c-8896-f62b485bfa4c
# ╠═447a4f69-5c9a-499f-8f02-917fab664d00
# ╠═0860e4c4-9075-4ce6-ac20-e549dda125ff
# ╠═9dc57b08-ac07-4c1a-bee4-2d9be075ac50
# ╠═ced3a00e-64db-4f13-86a9-8ef5c0856ce0
# ╠═68adc013-5ac1-4793-82ed-4e1cc03f9698
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
