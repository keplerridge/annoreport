from setuptools import setup

setup(
    name='annoreport',
    version='0.1.0',
    description='Summary and visualization tool for MAG gene annotation workflows',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='keplerridge',
    url='https://github.com/keplerridge/annoreport',
    license='MIT',
    py_modules=['annotation_report'],
    entry_points={
        'console_scripts': [
            'annoreport=annotation_report:main',
        ],
    },
    python_requires='>=3.9',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
    ],
)
