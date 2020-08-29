# SthenoCompat

Extend some function to [Stheno](https://github.com/willtebbutt/Stheno.jl), a project transfers many feature development work to something like `KernelFunctions` and `AbstractGPs` while itself doesn't support them (:<). So I backport some features I need to `Stheno`.

## Features

* `LinearWithTransform` and `Transform`: LinearKernel but with an extra transformation.
* Disable `In slow method` check, so `f.(ColVec(A))` can be used in AD.