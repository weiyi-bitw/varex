from distutils.core import setup

setup(
    name='varex',
    version='0.2',
    author='Wei-Yi Cheng',
    author_email='wei-yi.cheng@mssm.edu',
    packages=['varex', 'varex.commons'],
    include_package_data = True,
    scripts=['bin/explainVCF', 'bin/createVarMatrix', 'bin/plotHeatmap', 'bin/elementize_vcf', 'bin/elementize_ped'],
    url='https://github.com/weiyi-bitw/varex',
    license='LICENSE.txt',
    description='Explain the causal variant of phenotype',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy==1.9.1",
        "scikit-learn >= 0.15.1",
        "scipy==0.14.0",
        "varnorm >= 0.2",
        "db2util >= 0.2",
        "MySQL-python >= 1.2.5",
        "plotly >= 1.4.9",
        "argparse >= 1.2.1"
    ]
)

