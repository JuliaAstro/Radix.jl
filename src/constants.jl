# Physical constants

const mc² = 5.11e5

# exp() with the argument clamped to ±60, as in XSTAR's expo()
expo(x) = exp(clamp(x, -60, 60))
