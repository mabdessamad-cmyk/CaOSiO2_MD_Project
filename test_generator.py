from generator import GlassGenerator

g = GlassGenerator(
    x=0.30,
    formula_units=2400,
    density=1.50,
    seed=12345,
)

g.summary()
