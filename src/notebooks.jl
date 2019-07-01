export notebooks

import IJulia

function notebooks()
    IJulia.notebook(; dir=joinpath(@__DIR__, "..", "notebooks"))
end
