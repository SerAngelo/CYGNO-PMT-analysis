# setup.py
from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='pmt_analysis_library',
    version='0.1',
    author='SerAngelo',
    author_email='angeloserrecchia.phys@gmail.com',
    description='Una libreria per l\'analisi dati dei PMT per l'esperimento CYGNO',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/il_tuo_username/pmt_analysis_library',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

