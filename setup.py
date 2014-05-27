#! /usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Shantanu H. Joshi \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

from distutils.core import setup
# from setuptools import setup
base_dir = 'ccshape'
setup(
    name='ccshape',
    version='0.1dev',
    packages=['ccshape'],
    package_dir={'ccshape': base_dir, },
    test_suite='nose.collector',
    scripts=['bin/corpus_callosum_thickness.py', 'bin/corpus_callosum_analyze.py'],
    license='MIT/TBD',
    exclude_package_data={'': ['.gitignore', '.idea']},
    author='Brandon Ayers, Shantanu H. Joshi',
    author_email='ayersb@ucla.edu, s.joshi@ucla.edu',
    description='corpus callosum shape analysis',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: MIT/TBD',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    keywords='corpus callosum brain shape statistics',
)
