import mypy
from setuptools import setup, find_packages  # type: ignore

# Function to read the list of dependencies from requirements.txt
def read_requirements():
    with open('requirements.txt') as req:
        return req.read().strip().split('\n')

setup(
    name='my_package',
    version='0.1',
    packages=find_packages(),
    install_requires=read_requirements(),
    author='Gustave Ronteix',
    author_email='gustave.ronteix@orakl-oncology.com',
    description='Save us from RedCap',
    license='MIT',
)
