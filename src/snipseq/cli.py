import typer
import snipseq
from .bsc import bsc
from .adv import adv

app = typer.Typer(help="snipseq: demultiplex sequencing reads by barcode pairs.",  invoke_without_command=True)

# longhand, without decorators
app.command()(bsc)
app.command()(adv)

# callback without decorator
def main(version: bool = typer.Option(False, "--version", "-v", help="Show version")):
    if version:
        print(f"snipseq version {snipseq.__version__}")
        raise typer.Exit()

app.callback()(main)

if __name__ == "__main__":
    app()