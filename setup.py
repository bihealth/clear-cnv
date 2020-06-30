import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="clear_CNV-vinzenzmay", # Replace with your own username
    version="0.0.1",
    author="Vinzenz May",
    author_email="Vinzenz.May@mdc-berlin.de",
    description="Clinical sequencing panel CNV caller and visualizer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vinzenzmay/clear_CNV",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'clear_CNV = clear_CNV:main',
        ],
    },
    python_requires='>=3.7.6',
)
