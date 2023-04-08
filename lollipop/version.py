try:
    import importlib.metadata

    __version__ = importlib.metadata.version(__package__ or __name__)
except:
    __version__ = "0.0.0"

if __version__ == "0.0.0":
    import subprocess
    import os

    __version__ = (
        subprocess.check_output(
            ["git", "describe", "--tag"], cwd=os.path.dirname(os.path.abspath(__file__))
        )
        .decode("ascii")
        .strip()
    )
