from generator import GlassGenerator
from writer import write_header

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
)

cfg = g.generate()

with open("output/header_test.data", "w") as f:
    write_header(f, cfg)

print("Header written successfully.")
