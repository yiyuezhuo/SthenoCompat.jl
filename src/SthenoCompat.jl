module SthenoCompat

# Write your package code here.

using Stheno
using Zygote
import Stheno: AV, ew, pw

include("transformed_linear_kernel.jl")

end
