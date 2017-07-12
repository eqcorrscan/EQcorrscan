# Over-ride the -n auto in travis as this goes of the wall
import sys
import os


def pytest_cmdline_preparse(args):
    is_travis = 'TRAVIS' in os.environ
    if 'xdist' in sys.modules and is_travis:
        if "-n" in args:
            args[args.index("-n") + 1] = "2"
        else:
            args[:] = args + ["-n", "2"]
    print(args)


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",
                     help="run slow tests")
