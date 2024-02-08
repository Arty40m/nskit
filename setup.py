from setuptools import setup, Extension, find_packages
from pathlib import Path



with open(Path(__file__).parent/'README.md') as f:
    long_description = f.read()
    

setup(
    name = 'nskit',
    version = '1.8.0',
    author = 'Artem Mukanov',
    url = 'https://github.com/Arty40m/nskit',
    description = 'Library for processing nucleic acid secondary structure',
    long_description = long_description,
    license = 'LGPLv3', 
    classifiers = [
        'Programming Language :: Python :: 3 :: Only', 
        'Topic :: Scientific/Engineering', 
        'Topic :: Scientific/Engineering :: Bio-Informatics', 
        'Topic :: Software Development :: Libraries', 
        'Topic :: Software Development :: Libraries :: Python Modules', 
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)'
                  ],

    zip_safe = False,
    include_package_data = True, 
    packages = find_packages(), 
    python_requires = '>=3.8', 
    install_requires = ['numpy>=1.18'], 
    extras_require = {'dev':['pytest']}, 
    
    ext_modules=[
        Extension('nskit.algo.levenshtein._levenshtein', ['nskit/algo/levenshtein/_levenshtein.c']),
                ]
)