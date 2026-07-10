from writer import write_masses

with open("output/masses_test.data", "w") as f:
    write_masses(f)

print("Masses written successfully.")
