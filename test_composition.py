from composition import calculate_composition, total_charge

comp = calculate_composition(
    x=0.30,
    formula_units=2400,
    density=1.50,
)

print()

print("Composition")
print("-----------")
print("Formula units :", comp["formula_units"])
print("Si            :", comp["Si"])
print("Ca            :", comp["Ca"])
print("O             :", comp["O"])
print("Total atoms   :", comp["natoms"])

print()

print("Molecular weight :", comp["molecular_weight"])

print("Box length       :", comp["box_length"])

print()

print("Total charge :", total_charge(comp))
