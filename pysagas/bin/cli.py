import click
import pysagas


@click.group()
def cli():
    """PySAGAS command line interface."""
    pass


@click.command()
def version():
    """Shows the installed version."""
    print(f"\x1B[3mPySAGAS\x1B[0m v{pysagas.__version__} installed")


@click.command()
@click.option("-m", "--method", default="newtonian", help="The solver method to use.")
@click.option("-g", "--geometry", default="Components.i.tri", help="The geometry file.")
@click.option(
    "-v", "--velocity", default="1,0,0", help="The direction of the freestream."
)
def aero(method, geometry, velocity: str):
    """Calculate the aerodynamic coefficients."""
    from pysagas import Vector
    from pysagas.utilities import newtonian_impact_solver
    from pysagas.wrappers.cart3d.utilities import parse_tri_file

    # Print banner
    pysagas.banner()

    # # Define flow properties
    # gamma = 1.4
    # M_inf = 6
    # A_ref = 1  # m2
    # a_inf = 299.499  # m/s
    # rho_inf = 0.0308742  # kg/m3

    # V_inf = M_inf * a_inf  # m/s
    # p_inf = rho_inf * a_inf**2 / gamma
    # q_inf = 0.5*rho_inf*V_inf**2

    # cells = parse_tri_file()
    # v = Vector(1, 0, 0)

    # F, M = newtonian_impact_solver(cells=cells, v=v, p_inf=p_inf, q_inf=q_inf)

    # C_f = F / q_inf

    # print("Coming soon!")
    v = Vector.from_coordinates([float(d) for d in velocity.split(",")])
    print(v)


# Add commands to CLI
cli.add_command(version)
cli.add_command(aero)
