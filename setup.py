import setuptools
 
setuptools.setup(
    name='archer',
    version='1.0',
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={"archer": ["etc/*"]},
)
