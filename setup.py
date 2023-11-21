import mypy
from setuptools import setup, find_packages  # type: ignore

# Function to read the list of dependencies from requirements.txt
def read_requirements():
    with open('requirements.txt') as req:
        return req.read().strip().split('\n')

setup(
    name='redcappreprocessing',
    version='0.1',
    packages=find_packages(),
    install_requires=read_requirements(),
    entry_points={
        'console_scripts': [
            'start-myapp=redcap_preprocessing.app:run',
        ],
    },
    author='Gustave Ronteix',
    author_email='gustave.ronteix@orakl-oncology.com',
    description='Save us from RedCap',
    license='MIT',
)
