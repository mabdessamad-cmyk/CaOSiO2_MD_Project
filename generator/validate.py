from constants import CHARGE
from geometry import distance

def validate_composition(cfg):
    """
    Check composition and total charge.
    """

    n_si = cfg.species.count("Si")
    n_o = cfg.species.count("O")
    n_ca = cfg.species.count("Ca")

    total_charge = sum(cfg.charges)

    print()
    print("Validation")
    print("--------------------------")
    print(f"Si : {n_si}")
    print(f"O  : {n_o}")
    print(f"Ca : {n_ca}")
    print(f"Total atoms  : {cfg.natoms}")
    print(f"Total charge : {total_charge:.6f}")

    ok = (
        n_si == 1680 and
        n_o == 4080 and
        n_ca == 720 and
        cfg.natoms == 6480 and
        abs(total_charge) < 1.0e-10
    )

    if ok:
        print("Composition validation PASSED")
    else:
        print("Composition validation FAILED")

    return ok

def validate_distances(cfg):
    """
    Compute the minimum distance for every pair type.
    """
    pair_min = {
        ("Ca", "Ca"): float("inf"),
        ("Ca", "O"):  float("inf"),
        ("Ca", "Si"): float("inf"),
        ("O", "O"):   float("inf"),
        ("O", "Si"):  float("inf"),
        ("Si", "Si"): float("inf"),
    }

    for i in range(cfg.natoms - 1):

        spi = cfg.get_species(i)
        posi = cfg.get_position(i)

        for j in range(i + 1, cfg.natoms):

            spj = cfg.get_species(j)
            posj = cfg.get_position(j)

            pair = tuple(sorted((spi, spj)))

            d = distance(posi, posj, cfg.box_length)

            if d < pair_min[pair]:
                pair_min[pair] = d

    print()
    print("Minimum distances (Å)")
    print("--------------------------")

    for pair in [
        ("Si", "Si"),
        ("O", "Si"),
        ("Ca", "Si"),
        ("O", "O"),
        ("Ca", "O"),
        ("Ca", "Ca"),
    ]:
        print(f"{pair[0]}-{pair[1]:<2} : {pair_min[pair]:7.3f}")
