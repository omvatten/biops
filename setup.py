from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["pandas", "numpy", "matplotlib", "scipy"]

setup(
    name="biops",
    version="0.2.1",
    author="Oskar Modin",
    author_email="omvatten@gmail.com",
    description="Bioprocess simulator",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/omvatten/biops/",
    packages=find_packages(),
    install_requires=requirements,
	include_package_data=True,
	package_data={'': ['data/*']},
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)