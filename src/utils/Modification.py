import pandas as pd


class Modification:
    def __init__(self, name, symbol, origin, nucleotide, base, base_delta_mass, sugar, sugar_delta_mass, phosphate,
                 phosphate_delta_mass):
        self.name = name
        self.symbol = symbol
        self.origin = origin
        self.nucleotide = nucleotide
        self.base = base
        self.base_delta_mass = base_delta_mass
        self.sugar = sugar
        self.sugar_delta_mass = sugar_delta_mass
        self.phosphate = phosphate
        self.phosphate_delta_mass = phosphate_delta_mass


def create_modifications_from_table(table):
    mods = {}
    for _, row in table.iterrows():
        mod = Modification(
            name=row['Name'],
            symbol=row['Symbol'],
            origin=row['origin'],
            nucleotide=row['nucleotide'],
            base=row['Base'],
            base_delta_mass=row['delta_mass'],
            sugar=row['sugar'],
            sugar_delta_mass=row['delta_mass.1'],
            phosphate=row['phosphate'],
            phosphate_delta_mass=row['delta_mass.2']
        )
        mods[row['Symbol']] = mod
    return mods


def get_known_mod(known_mod_path):
    df = pd.read_csv(known_mod_path, delimiter=' ')
    known_mod_data = []
    for _, row in df.iterrows():
        known_mod_data.append({
            'Molecule': row['Molecule'],
            'Position': row['Position'],
            'ID': row['ID'],
            'ID_ext': row['ID_ext'],
            'Include': row['Include']
        })

    print(known_mod_data)
    return known_mod_data


if __name__ == '__main__':
    # 读入表格
    df = pd.read_csv('data/mods.csv')

    # 创建 nucleotide 对象列表
    mods = create_modifications_from_table(df)

    # 打印 nucleotide 对象列表中第一个对象的名字和符号
    print(mods['Ap'].name, mods['Ap'].symbol)
