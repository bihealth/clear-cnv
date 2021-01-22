import os, pathlib, yaml


def f():
    x = pathlib.Path(__file__).absolute().parent
    print(str(x))
    y = x / "file.txt"
    print(str(y))


f()
