import setuptools
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


thelibFolder = os.path.dirname(os.path.realpath(__file__))

requirementPath = thelibFolder + '/lib/requirements.txt'
install_requires = [] # Examples: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
print(requirementPath)
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()
        
setuptools.setup(
    name="PyImFCS", # Replace with your own username
    version="0.0.1",
    author="Aurélien Barbotin",
    description="Useful stuff for imFCS data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    py_modules=[],
    python_requires=">=3.6",
    install_requires=install_requires,
)
