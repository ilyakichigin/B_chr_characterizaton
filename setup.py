
from setuptools import setup, find_packages

setup(
    name='dopseq',
	version='1.0.0-alpha',
	description='Bioinformatics pipeline for chromosome sequencing analysis',
	author='Alex Makunin',
	author_email='alex@mcb.nsc.ru',
	license='',
	packages=find_packages(exclude=['anolis']),
	install_requires=['pysam>=0.10.0', 'cutadapt>=1.8.3', 'pyyaml>=3.12'],
    include_package_data=True,
	entry_points={
        'console_scripts': [
            'dopseq_pipeline=dopseq.dopseq_pipeline:main',
            'variation_pipeline=dopseq.variation_pipeline:main',
        ],
	},
	)