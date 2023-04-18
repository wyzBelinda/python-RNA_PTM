import pandas as pd


class Modification:
    def __init__(self, name, symbol, origin, nucleotide, base, delta_mass_1, sugar, delta_mass_2, phosphate,
                 delta_mass_3):
        self.name = name
        self.symbol = symbol
        self.origin = origin
        self.nucleotide = nucleotide
        self.base = base
        self.delta_mass_1 = delta_mass_1
        self.sugar = sugar
        self.delta_mass_2 = delta_mass_2
        self.phosphate = phosphate
        self.delta_mass_3 = delta_mass_3


def create_modifications_from_table(table):
    mods = []
    for _, row in table.iterrows():
        mod = Modification(
            name=row['Name'],
            symbol=row['Symbol'],
            origin=row['origin'],
            nucleotide=row['nucleotide'],
            base=row['Base'],
            delta_mass_1=row['delta_mass'],
            sugar=row['sugar'],
            delta_mass_2=row['delta_mass.1'],
            phosphate=row['phosphate'],
            delta_mass_3=row['delta_mass.2']
        )
        mods.append(mod)
    return mods


if __name__ == '__main__':
    # 读入表格
    df = pd.read_csv('data/mods.csv')

    # 创建 nucleotide 对象列表
    mods = create_modifications_from_table(df)

    # 打印 nucleotide 对象列表中第一个对象的名字和符号
    print(mods[0].name, mods[0].symbol)
