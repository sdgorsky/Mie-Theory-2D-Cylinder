import numpy as np
from scipy.special import jv, hankel1, kv, iv
import platform
from datetime import datetime, UTC
import textwrap
from importlib.metadata import version
import sys

this_system = platform.uname()

orders_max = 20 # Number of bessel orders

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

    bessel_j_values = [] # (order, z_r, z_theta, eval_real, eval_imag)
    for order in range(0,orders_max+1):
        for r in np.logspace(-3,2,30):
            for theta in np.linspace(0, 2*np.pi, 20):
                value = jv(order, r * np.exp(1j*theta) )
                bessel_j_values.append((order, r.item(), theta.item(), value.real.item(), value.imag.item()))

    with open("artifacts/bessel_j_evaluations.txt", "w") as f:

        f.write(textwrap.dedent(
            """\
            # Complex evaluations of Bessel function of the first kind generated on a square grid.
            #
            # Function representation: J(ν, z)
            #   J --> [complex] Bessel function of the first kind
            #   z --> [complex] argument
            #   ν --> [int] bessel order
            #
            """))
        f.write(shared_info)
        f.write("# Data shape: (order, z_r, z_theta, real(J), imag(J)) \n")

        for v in bessel_j_values:
            f.write(str(v))
            f.write("\n")

def create_hankel1_file():

    bessel_h_values = [] # (order, z_real, z_imag, eval_real, eval_imag)
    for order in range(0,orders_max+1):
        for r in np.logspace(-3,2,30):
            for theta in np.linspace(0, 2*np.pi, 20):
                value = hankel1(order, r * np.exp(1j*theta) )
                bessel_h_values.append((order, r.item(), theta.item(), value.real.item(), value.imag.item()))

    with open("artifacts/hankel1_evaluations.txt", "w") as f:

        f.write(textwrap.dedent(
            """\
            # Complex evaluations of Hankel function of the first find generated on a square grid.
            #
            # Function representation: H1(ν, z)
            #   H --> [complex] Hankel function of the first kind
            #   z --> [real] argument   
            #   ν --> [int] bessel order
            #
            """))
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
                values.append((order, r.item(), theta.item(), value.real.item(), value.imag.item()))

    with open("artifacts/bessel_k_evaluations.txt", "w") as f:
        f.write(textwrap.dedent(
            """\
            # Complex evaluations of modified Bessel function K_n(z) for z in the RHP.
            #
            # Function representation: K(ν, z)
            #   K --> [complex] modified Bessel function of the second kind
            #   z --> [complex] argument (Re(z) > 0)
            #   ν --> [int] bessel order
            #
            """))
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
                values.append((order, r.item(), theta.item(), value.real.item(), value.imag.item()))

    with open("artifacts/bessel_i_evaluations.txt", "w") as f:
        f.write(textwrap.dedent(
            """\
            # Complex evaluations of modified Bessel function I_n(z) for z in the RHP.
            #
            # Function representation: I(ν, z)
            #   I --> [complex] modified Bessel function of the first kind
            #   z --> [complex] argument (Re(z) > 0)
            #   ν --> [int] bessel order
            #
            """))
        f.write(shared_info)
        f.write("# Data shape: (order, z_r, z_theta, real(I), imag(I)) \n")

        for v in values:
            f.write(str(v))
            f.write("\n")

if __name__ == "__main__":
    create_bessel_j_file()
    create_hankel1_file()
    create_bessel_k_file()
    create_bessel_i_file()