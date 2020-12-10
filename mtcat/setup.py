import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name='mtcat',
    version='0.0.1',
    author='Francesca-Zhoufan Li, Evan Mun, and Noah Epsteing',
    author_email='fzl@caltech.edu',
    description='Packages for Microtubule Catastrophe Data Polishing.',
    long_description=long_description,
    long_description_content_type='ext/markdown',
    packages=setuptools.find_packages(),
    install_requires=["numpy","pandas", "bokeh>=1.4.0"],
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ),
)