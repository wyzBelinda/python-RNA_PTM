import pandas as pd
from pyteomics import mgf, fasta

from example.utils.Element import create_elements_from_file
from example.utils.Modification import create_modifications_from_table
from example.utils.Nucleotide import create_nucleotides_from_table


def something():
    '''
    运行
    :return:
    '''
    # 创建 element 对象列表
    elements = create_elements_from_file('data/element.ini')

    # 读入表格
    df = pd.read_csv('data/nts.ini')

    # 创建 nucleotide 对象列表
    nucleotides = create_nucleotides_from_table(df)
    # 计算每个基团的质量并赋值
    for nt in nucleotides:
        nt.gen_mass_by_elements(elements)

    # 载入所有已知的修饰数据
    # 读入表格
    df = pd.read_csv('data/mods.csv')

    # 创建 nucleotide 对象列表
    mods = create_modifications_from_table(df)



if __name__ == '__main__':

    # Open the MGF file and read the first spectrum
    with mgf.read('example.mgf') as spectra:
        spectrum = next(spectra)

    # Print the m/z and intensity values for each peak in the first spectrum
    for mz, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
        print('m/z:', mz, 'Intensity:', intensity)

    # Open the FASTA file and read the first protein
    with fasta.read('example.fasta') as f:
        header, sequence = next(f)

    # Print the header and sequence of the first protein
    print('Header:', header)
    print('Sequence:', sequence)
