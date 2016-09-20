from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

setup(
    name='crnpy',

    version='0.0.1',

    description='A Python package for the analysis of Chemical Reaction Networks.',
    long_description='A Python package for the analysis of Chemical Reaction Networks.',

    url='https://github.com/etonello/crnpy',

    author='Elisa Tonello',
    author_email='elisa@tonello.me',

    license='BSD',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science',
        'Intended Audience :: Research',
        'License :: BSD',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='chemical reaction networks',

    packages=find_packages(exclude=['tests']),

    # run-time dependencies that will be installed by pip
    install_requires=['python-libsbml', 'numpy', 'scipy', 'sympy'],
)

