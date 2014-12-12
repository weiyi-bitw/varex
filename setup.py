from distutils.core import setup

setup(
    name='varex',
    version='0.2',
    author='Wei-Yi Cheng',
    author_email='wei-yi.cheng@mssm.edu',
    packages=['varex'],
    package_dir={'varex':'varex'},
    include_package_data = True,
    scripts=['bin/explainVCF', 'bin/createVarMatrix'],
    url='http://hidysabc.com/blog',
    license='LICENSE.txt',
    description='Explain the causal variant of phenotype',
    long_description=open('README.md').read(),
    install_requires=[
        "varnorm >= 0.2",
        "db2util >= 0.2",
        "MySQL-python >= 1.2.5",
        "plotly >= 1.4.9",
        "argparse >= 1.2.1"
    ]
)

