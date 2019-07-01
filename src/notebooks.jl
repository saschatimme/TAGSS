export ed_notebook

import IJulia

function ed_notebook()
    IJulia.notebook(; dir=joinpath(@__DIR__, "..", "notebooks", "ed-degree.ipynb"))
end
