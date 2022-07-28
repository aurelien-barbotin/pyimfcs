import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

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
)
