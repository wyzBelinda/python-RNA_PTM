import pandas as pd
from typing import List

from example.utils.Element import cal_mass_by_formula
from example.utils.Oligo import Oligo


class Nucleotide:
    def __init__(self, name, base, sugar, phosphate, base_mass=0, sugar_mass=0, phosphate_mass=0, whole_mass=0):
        self.name = name
        self.base = base
        self.sugar = sugar
        self.phosphate = phosphate
        self.base_mass = base_mass
        self.sugar_mass = sugar_mass
        self.phosphate_mass = phosphate_mass
        self.whole_mass = whole_mass

    def gen_mass_by_elements(self, elements):
        self.base_mass = cal_mass_by_formula(elements, self.base)
        self.sugar_mass = cal_mass_by_formula(elements, self.sugar)
        self.phosphate_mass = cal_mass_by_formula(elements, self.phosphate)
        self.whole_mass = self.base_mass + self.sugar_mass + self.phosphate_mass


def create_nucleotides_from_table(table):
    nts = {}
    for _, row in table.iterrows():
        nucleotide = Nucleotide(
            name=row['name'],
            base=row['base'],
            sugar=row['sugar'],
            phosphate=row['phosphate'],
        )
        nts[row['name']] = nucleotide
    return nts


def cal_mass_by_sequence(nucleos, sequence:str):
    mass = 0
    for nt in sequence:
        mass += nucleos[nt].whole_mass
    mass -= nucleos[sequence[-1]].phosphate_mass
    return mass


if __name__ == '__main__':
    # 读入表格
    df = pd.read_csv('data/nts.ini')

    # 创建 nucleotide 对象列表
    nucleotides = create_nucleotides_from_table(df)

    # 打印 nucleotide 对象列表中第一个对象的名字和磷酸根
    print(nucleotides['A'].name, nucleotides['A'].phosphate)
