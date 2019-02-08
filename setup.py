import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="MERCluster",
	version="0.0.1",
	author="Stephen W. Eichhorn",
	author_email="stephen_eichhorn@fas.harvard.edu",
	description="A package for MERFISH/SC-seq clustering and comparisons",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/seichhorn/MERCluster",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU GPLv3",
		"Operating System :: OS Independent",
	],
)