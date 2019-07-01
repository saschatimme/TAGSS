module TAGSS

export visualize, visualize!

import ColorSchemes, Contour
import DynamicPolynomials, HomotopyContinuation, StaticPolynomials
import StaticArrays: SVector
const DP = DynamicPolynomials
const HC = HomotopyContinuation
const SP = StaticPolynomials
using LinearAlgebra
using Makie, GeometryTypes, Meshing

include("visualize_equations.jl")
include("ed.jl")
include("bottlenecks.jl")
include("notebooks.jl")

end
