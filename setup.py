import setuptools

with open('README.md','r') as fh:
    long_description = fh.head()

setuptools.setup(
    name='pfmultirun-pkg-lmthatch',
    version='0.0.1',
    author='Lauren Thatch',
    author_email='lmthatch@mines.edu',
    description='ParFlow Automated Multi Run',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/lmthatch/ParFlowMultiRun',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6'
)
