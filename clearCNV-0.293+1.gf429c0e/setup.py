import os.path
import setuptools

import versioneer


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


with open("README.md") as readme_file:
    readme = readme_file.read()

with open("HISTORY.md") as history_file:
    history = history_file.read()


install_requirements = parse_requirements("requirements/base.txt")


setuptools.setup(
    name="clearCNV",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Vinzenz May",
    author_email="Vinzenz.May@mdc-berlin.de",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="Clinical sequencing panel CNV caller and visualizer",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://github.com/bihealth/clear-CNV",
    packages=setuptools.find_packages(),
    entry_points={
        "console_scripts": [
            "clearCNV = clearCNV.__main__:main",
        ],
    },
    install_requires=install_requirements,
    python_requires=">=3.7",
    zip_safe=False,
)
