import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="clearCNV", # Replace with your own username
    version="0.0.1",
    author="Vinzenz May",
    author_email="Vinzenz.May@mdc-berlin.de",
    description="Clinical sequencing panel CNV caller and visualizer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bihealth/clear-CNV",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'clearCNV = clearCNV.__main__:main',
        ],
    },
    install_requires=["numpy>=1.18.1", "pandas>=1.0.1", "scipy>=1.4.1",
    "scikit-learn>=0.22.1", "seaborn>=0.10.0", "hmmlearn>=0.2.3",
    "matplotlib>=3.1.3", "plotly>=4.5.4"],
    python_requires='>=3.7'
)
