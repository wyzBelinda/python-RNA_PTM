import pandas as pd
from Bio.Restriction import RestrictionType


class Enzyme:
    def __init__(self, species, seq_regex, cut_3end, cut_5end):
        self.species = species
        self.seq_regex = seq_regex
        self.cut_3end = cut_3end
        self.cut_5end = cut_5end

    def to_dna_restriction_type(self):
        site,cut=find_star_position(self.seq_regex)
        # 自定义限制酶
        my_enzyme = RestrictionType(
            self.species,
            site="GATC",
            cut=(0,)
        )


def find_star_position(seq):
    star_pos = seq.index('*')
    seq_before_star = ''.join(seq[:star_pos])
    seq_after_star = ''.join(seq[star_pos + 1:])
    return seq_before_star + seq_after_star, star_pos


def create_enzymes_from_table(table):
    enzymes = []
    for _, row in table.iterrows():
        enzyme = Enzyme(
            species=row['species'],
            seq_regex=row['seq_regex'],
            cut_3end=row['cut_3end'],
            cut_5end=row['cut_5end']
        )
        enzymes.append(enzyme)
    return enzymes


if __name__ == '__main__':
    # 读入表格
    df = pd.read_csv('data/enzyme.csv')

    # 创建 enzyme 对象列表
    enzymes = create_enzymes_from_table(df)

    # 打印 enzyme 对象列表中第一个对象的物种和序列正则表达式
    print(enzymes[0].species, enzymes[0].seq_regex)
