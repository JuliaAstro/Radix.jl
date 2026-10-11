# `RADIX_QUICK=1` skips the tests that take minutes with the database (the thermal equilibrium of zones, the thick slab, three passes: about 10 of the 15 minutes
# of a laptop), for the continuous integration; the whole suite is what is run otherwise
macro slow(name, body)
    :(isempty(get(ENV, "RADIX_QUICK", "")) ? $(esc(:(@testset $name $body))) : @info("Skipping the slow test (RADIX_QUICK): " * $name))
end
