import pandas as pd


class Enzyme:
    def __init__(self, name, seq_regex, cut_3end, cut_5end):
        self.name = name
        self.seq_regex = seq_regex
        self.cut_3end = cut_3end
        self.cut_5end = cut_5end
        self.cut_list = []

    def analysis(self):
        regex_list = self.seq_regex.split("|")
        for regex in regex_list:
            cut_signal, cut_pos = find_star_position(regex)
            self.cut_list.append((cut_signal, cut_pos))



def find_star_position(seq):
    star_pos = seq.index('*')
    seq_before_star = ''.join(seq[:star_pos])
    seq_after_star = ''.join(seq[star_pos + 1:])
    return seq_before_star + seq_after_star, star_pos


def create_enzymes_from_table(table):
    enzymes = {}
    for _, row in table.iterrows():
        enzyme = Enzyme(
            name=row['name'],
            seq_regex=row['seq_regex'],
            cut_3end=row['cut_3end'],
            cut_5end=row['cut_5end']
        )
        enzyme.analysis()
        enzymes[row['name']] = enzyme
    return enzymes


def get_used_enzymes(enzyme_path, enzymes_name):
    if enzymes_name =='':
        return {}

    # 读入表格
    df = pd.read_csv(enzyme_path)

    # 创建 enzyme 对象列表
    all_enzymes = create_enzymes_from_table(df)
    used_enzymes = {}

    for ezy_name in enzymes_name.split(','):
        used_enzymes[ezy_name] = all_enzymes[ezy_name]
    return used_enzymes


if __name__ == '__main__':
    # 读入表格
    df = pd.read_csv('data/enzyme.csv')

    # 创建 enzyme 对象列表
    enzymes = create_enzymes_from_table(df)

    # 打印 enzyme 对象列表中第一个对象的物种和序列正则表达式
    print(enzymes['T1'].name, enzymes['T1'].seq_regex)
