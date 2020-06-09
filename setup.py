import os
import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()


install_requires = [line.rstrip() for line in open(
    os.path.join(os.path.dirname(__file__), "requirements.txt"))]


setuptools.setup(
	name="mercluster",
	version="0.0.1",
	author="Stephen W. Eichhorn",
	author_email="stephen_eichhorn@fas.harvard.edu",
	description="MERFISH analysis software",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/seichhorn/MERCluster",
	packages=setuptools.find_packages(),
	install_requires=install_requires,
	entry_points={'console_scripts': ["mercluster=mercluster.mercluster:mercluster"]}
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
)