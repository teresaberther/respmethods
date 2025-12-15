from setuptools import setup, Extension
from numpy import get_include

mod_phase_ext = Extension(
    'PhaseExtraction',
    sources=['./src/phase-extraction.c'],
    include_dirs=[get_include()],
)

setup(
    name='respymethods',
    version='0.1',
    ext_modules=[mod_phase_ext],
)
