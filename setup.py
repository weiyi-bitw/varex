from distutils.core import setup

setup(
    name="varex",
    version="0.2",
    author="Wei-Yi Cheng",
    author_email="ninpy.weiyi@gmail.com",
    packages=["varex", 
        "varex.commons"
        ],
    include_package_data = True,
    scripts=["bin/explain_vcf", 
        "bin/create_var_matrix", 
        "bin/elementize_vcf", 
        "bin/elementize_ped",
        "bin/gwas",
        "bin/gwas_hypergeo",
        "bin/normalize_gwas_out",
        "bin/count_allele",
        "bin/query_genotype_matrix",
        "bin/reinstall"],
    url="https://github.com/weiyi-bitw/varex",
    license="LICENSE.txt",
    description="Explain the causal variant of phenotype",
    long_description=open("README.md").read(),
    install_requires=[
        "numpy==1.9.1",
        "scikit-learn >= 0.15.1",
        "scipy==0.14.0",
        "varnorm >= 0.2",
        "db2util >= 0.2",
        "MySQL-python >= 1.2.5",
        "argparse >= 1.2.1"
    ]
)

