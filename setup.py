
from setuptools import setup, find_packages 

# Extract the long description from the README file
with open("README.md", "r") as fh:
    long_description = fh.read()

# Extract the required files from the requirements text file
with open("requirements.txt") as f:
    required = f.read().splitlines()

# Generate the setup function
setup(
    name = "bmdrc",
    version = "0.0.1",
    description = "Benchmark Dose Response Curve Calculations",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/PNNL-CompBio/bmdrc",
    author = ["Degnan, David", "Gosline, Sara", "Waters, Katrina"],
    author_email = ["david.degnan@pnnl.gov", "sara.gosline@pnnl.gov", "katrina.waters@pnnl.gov"],
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent"
    ],
    package_data = {'external': ['disclaimer.txt']},
    packages = find_packages(),
    include_package_data = True,
    install_requires = required,
    setup_requires = ['pytest-runner', 'wheel'],
    test_suite = 'pytest',
    tests_require = ['pytest']
)

