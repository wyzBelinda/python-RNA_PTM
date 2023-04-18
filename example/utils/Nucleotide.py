import pandas as pd

from example.utils.Element import cal_mass_by_formula


class Nucleotide:
    def __init__(self, name, base, sugar, phosphate, base_mass=0, sugar_mass=0, phosphate_mass=0):
        self.name = name
        self.base = base
        self.sugar = sugar
        self.phosphate = phosphate
        self.base_mass = base_mass
        self.sugar_mass = sugar_mass
        self.phosphate_mass = phosphate_mass

    def gen_mass_by_elements(self, elements):
        self.base_mass = cal_mass_by_formula(elements, self.base)
        self.sugar_mass = cal_mass_by_formula(elements, self.sugar)
        self.phosphate_mass = cal_mass_by_formula(elements, self.phosphate)


def create_nucleotides_from_table(table):
    nts = []
    for _, row in table.iterrows():
        nucleotide = Nucleotide(
            name=row['name'],
            base=row['base'],
            sugar=row['sugar'],
            phosphate=row['phosphate'],
        )
        nts.append(nucleotide)
    return nts


if __name__ == '__main__':
    # 读入表格
    df = pd.read_csv('data/nts.ini')

    # 创建 nucleotide 对象列表
    nucleotides = create_nucleotides_from_table(df)

    # 打印 nucleotide 对象列表中第一个对象的名字和磷酸根
    print(nucleotides[0].name, nucleotides[0].phosphate)
