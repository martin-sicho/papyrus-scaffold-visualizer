"""
setup.py

Created by: Martin Sicho
On: 10.10.22, 17:20
"""

from setuptools import setup
from distutils.util import convert_path

main_ns = {}
ver_path = convert_path('src/scaffviz/__init__.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)

setup(version=main_ns['VERSION'])

