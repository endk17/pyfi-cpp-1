from setuptools import setup, Extension


ext = Extension('bs', sources=['black_scholes/bs.c'])

setup(
    name="black_scholes",
    version="0.0.1",
    description="Options Pricing Library",
    packages=['black_scholes'],
    ext_modules=[ext]
)