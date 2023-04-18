from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name='example',
    version='0.1dev0',
    author='WU Yuzhuo',
    author_email='wyzbelinda@gmail.com',
    packages=find_packages(),
    long_description=open('README.md').read()
)