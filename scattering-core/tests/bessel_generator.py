import numpy as np
from scipy.special import jv, hankel1, kv, iv
import platform
from datetime import datetime, UTC
import textwrap
from importlib.metadata import version
import sys

this_system = platform.uname()

orders_max = 20  # Number of bessel orders

shared_info = textwrap.dedent(f"""\
    # Hardware details:
    #   date: {datetime.now(tz=UTC)}
    #   system: {this_system.system}
    #   release: {this_system.release}
    #   version: {this_system.version}
    #   machine: {this_system.machine}
    #   processor: {this_system.processor}
    #
    # Software details:
    #   python version: {sys.version}
    #   scipy version: {version("scipy")}
    #
    """)


def create_bessel_j_file():

    bessel_j_values = []  # (order, z_r, z_theta, eval_real, eval_imag)
    for order in range(0, orders_max + 1):
        for r in np.logspace(-3, 2, 30):
            for theta in np.linspace(0, 2 * np.pi, 20):
                value = jv(order, r * np.exp(1j * theta))
                bessel_j_values.append(
                    (
                        order,
                        r.item(),
                        theta.item(),
                        value.real.item(),
                        value.imag.item(),
                    )
                )

    with open("truth_data/bessel_j_evaluations.txt", "w") as f:
        f.write(
            textwrap.dedent(
                """\
            # Complex evaluations of Bessel function of the first kind generated on a square grid.
            #
            # Function representation: J(ν, z)
            #   J --> [complex] Bessel function of the first kind
            #   z --> [complex] argument
            #   ν --> [int] bessel order
            #
            """
            )
        )
        f.write(shared_info)
        f.write("# Data shape: (order, z_r, z_theta, real(J), imag(J)) \n")

        for v in bessel_j_values:
            f.write(str(v))
            f.write("\n")


def create_hankel1_file():

    bessel_h_values = []  # (order, z_real, z_imag, eval_real, eval_imag)
    for order in range(0, orders_max + 1):
        for r in np.logspace(-3, 2, 30):
            for theta in np.linspace(0, 2 * np.pi, 20):
                value = hankel1(order, r * np.exp(1j * theta))
                bessel_h_values.append(
                    (
                        order,
                        r.item(),
                        theta.item(),
                        value.real.item(),
                        value.imag.item(),
                    )
                )

    with open("truth_data/hankel1_evaluations.txt", "w") as f:
        f.write(
            textwrap.dedent(
                """\
            # Complex evaluations of Hankel function of the first find generated on a square grid.
            #
            # Function representation: H1(ν, z)
            #   H --> [complex] Hankel function of the first kind
            #   z --> [real] argument   
            #   ν --> [int] bessel order
            #
            """
            )
        )
        f.write(shared_info)
        f.write("# Data shape: (order, real(z), real(J), imag(J)) \n")

        for v in bessel_h_values:
            f.write(str(v))
            f.write("\n")


def create_bessel_k_file():
    """K_n(z) for z in the RHP (Re(z) > 0). Uses same grid as J/H1 but
    restricted to angles (-pi/2, pi/2) so Re(z) > 0."""

    values = []
    for order in range(0, orders_max + 1):
        for r in np.logspace(-3, 2, 30):
            # RHP only: theta in (-pi/2, pi/2)
            for theta in np.linspace(-np.pi / 2 + 0.01, np.pi / 2 - 0.01, 20):
                z = r * np.exp(1j * theta)
                value = kv(order, z)
                values.append(
                    (
                        order,
                        r.item(),
                        theta.item(),
                        value.real.item(),
                        value.imag.item(),
                    )
                )

    with open("truth_data/bessel_k_evaluations.txt", "w") as f:
        f.write(
            textwrap.dedent(
                """\
            # Complex evaluations of modified Bessel function K_n(z) for z in the RHP.
            #
            # Function representation: K(ν, z)
            #   K --> [complex] modified Bessel function of the second kind
            #   z --> [complex] argument (Re(z) > 0)
            #   ν --> [int] bessel order
            #
            """
            )
        )
        f.write(shared_info)
        f.write("# Data shape: (order, z_r, z_theta, real(K), imag(K)) \n")

        for v in values:
            f.write(str(v))
            f.write("\n")


def create_bessel_i_file():
    """I_n(z) for z in the RHP (Re(z) > 0). Same grid as K."""

    values = []
    for order in range(0, orders_max + 1):
        for r in np.logspace(-3, 2, 30):
            # RHP only: theta in (-pi/2, pi/2)
            for theta in np.linspace(-np.pi / 2 + 0.01, np.pi / 2 - 0.01, 20):
                z = r * np.exp(1j * theta)
                value = iv(order, z)
                values.append(
                    (
                        order,
                        r.item(),
                        theta.item(),
                        value.real.item(),
                        value.imag.item(),
                    )
                )

    with open("truth_data/bessel_i_evaluations.txt", "w") as f:
        f.write(
            textwrap.dedent(
                """\
            # Complex evaluations of modified Bessel function I_n(z) for z in the RHP.
            #
            # Function representation: I(ν, z)
            #   I --> [complex] modified Bessel function of the first kind
            #   z --> [complex] argument (Re(z) > 0)
            #   ν --> [int] bessel order
            #
            """
            )
        )
        f.write(shared_info)
        f.write("# Data shape: (order, z_r, z_theta, real(I), imag(I)) \n")

        for v in values:
            f.write(str(v))
            f.write("\n")


def _edge_case_angles():
    """Angles targeting implementation branch points and blind spots.

    The main grid uses 20 evenly-spaced angles (~0.33 rad apart), which misses
    narrow regions near the axes where our implementations switch code paths.
    These edge-case angles provide dense coverage at each critical boundary.

    Returns deduplicated, sorted angles in [0, 2π).
    """
    eps_angles = [0.01, 0.05]  # offsets from each axis
    angles = set()

    # Near θ=0 / θ=2π (positive real axis):
    #   hankel1() dispatches J+iY for Im(z)==0 vs hankel1_via_k for Im(z)!=0.
    #   Small angles just above/below the real axis stress this branch.
    for eps in eps_angles:
        angles.add(eps)
        angles.add(2 * np.pi - eps)

    # Near θ=π/2 and θ=3π/2 (imaginary axis):
    #   hankel1_via_k computes w = -iz, so z near +imag maps w near +real,
    #   and z near -imag maps w near -real. This is the RHP/LHP (ZBKNU vs
    #   ZACON) branch boundary in bessel_k. Small offsets from the imaginary
    #   axis put w.re ≈ 0, stressing that decision.
    for eps in eps_angles:
        angles.add(np.pi / 2 - eps)
        angles.add(np.pi / 2 + eps)
        angles.add(3 * np.pi / 2 - eps)
        angles.add(3 * np.pi / 2 + eps)

    # Near θ=π (negative real axis):
    #   bessel_j uses reflection J_n(-z) = (-1)^n J_n(z) for Re(z)<0.
    #   bessel_k hits the ZACON analytic continuation branch cut here.
    for eps in eps_angles:
        angles.add(np.pi - eps)
        angles.add(np.pi + eps)

    return sorted(angles)


def _edge_case_radii():
    """Radii targeting the series/Miller handoff boundary.

    bessel_k switches from power series (|z|<=2) to Miller backward
    recurrence (|z|>2). Points straddling this boundary catch discontinuities
    at the handoff. The logspace points cover the full 0.001–100 domain.
    """
    boundary_radii = [1.9, 1.99, 2.0, 2.01, 2.1]
    logspace_radii = np.logspace(-3, 2, 30).tolist()
    combined = sorted(set(boundary_radii + logspace_radii))
    return combined


def _generate_edge_case_data(func, angles, radii, orders):
    """Evaluate a scipy Bessel function on the edge-case grid."""
    values = []
    for order in orders:
        for r in radii:
            for theta in angles:
                z = r * np.exp(1j * theta)
                value = func(order, z)
                values.append((order, r, theta, value.real.item(), value.imag.item()))
    return values


def _write_edge_case_file(filename, header_comment, values):
    """Write edge-case evaluation data to an artifact file."""
    with open(f"truth_data/{filename}", "w") as f:
        f.write(header_comment)
        f.write(shared_info)
        f.write("# Data shape: (order, z_r, z_theta, real, imag)\n")
        for v in values:
            f.write(str(v))
            f.write("\n")


def create_edge_case_files():
    """Generate edge-case reference data for all four Bessel functions.

    These complement the main grid tests by densely sampling near
    implementation branch points: the real/imaginary axes and the
    |z|=2 series/Miller boundary.
    """
    angles_full = _edge_case_angles()
    radii = _edge_case_radii()
    orders = list(range(0, orders_max + 1))

    # K and I are only defined in the RHP (Re(z) > 0), so restrict angles
    # to (-π/2, π/2), i.e. angles in (0, π/2) ∪ (3π/2, 2π).
    angles_rhp = [a for a in angles_full if a < np.pi / 2 or a > 3 * np.pi / 2]

    funcs = [
        (
            "bessel_j_edge_cases.txt",
            jv,
            angles_full,
            textwrap.dedent("""\
            # Edge-case evaluations of Bessel J near implementation branch points.
            #
            # Angles cluster near the real axis (θ≈0, 2π), imaginary axis (θ≈π/2, 3π/2),
            # and negative real axis (θ≈π). Radii include points straddling the |z|=2
            # series/Miller boundary.
            #
            """),
        ),
        (
            "hankel1_edge_cases.txt",
            hankel1,
            angles_full,
            textwrap.dedent("""\
            # Edge-case evaluations of Hankel H1 near implementation branch points.
            #
            # Angles cluster near the real axis (J+iY vs K-path branch in hankel1),
            # imaginary axis (RHP/LHP branch in bessel_k via w=-iz), and negative real
            # axis (ZACON continuation). Radii include the |z|=2 boundary.
            #
            """),
        ),
        (
            "bessel_k_edge_cases.txt",
            kv,
            angles_rhp,
            textwrap.dedent("""\
            # Edge-case evaluations of Bessel K in the RHP near branch points.
            #
            # Angles cluster near the real axis (small Im(z), stressing Miller recurrence
            # convergence) and near the imaginary axis (Re(z)≈0, the RHP/LHP boundary).
            # Radii include the |z|=2 series/Miller handoff.
            #
            """),
        ),
        (
            "bessel_i_edge_cases.txt",
            iv,
            angles_rhp,
            textwrap.dedent("""\
            # Edge-case evaluations of Bessel I in the RHP near branch points.
            #
            # Same grid as K edge cases. I_n is computed via J_n(iz), so angles near the
            # imaginary axis test the J_n reflection for Re(iz)≈0.
            #
            """),
        ),
    ]

    for filename, func, angles, header in funcs:
        values = _generate_edge_case_data(func, angles, radii, orders)
        _write_edge_case_file(filename, header, values)
        print(f"  {filename}: {len(values)} points")


if __name__ == "__main__":
    create_bessel_j_file()
    create_hankel1_file()
    create_bessel_k_file()
    create_bessel_i_file()
    print("Edge-case files:")
    create_edge_case_files()
