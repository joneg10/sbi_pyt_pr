from setuptools import setup, find_packages

setup(
    name="wisefrog",
    version="0.1",
    packages=find_packages(),
    package_data={
        '': ['*.pytorch'],  # All .pytorch files anywhere in the package
    },
    entry_points={
        'console_scripts': [
            'wisefrog = wisefrog.wisefrog:main',
        ],
    },
)