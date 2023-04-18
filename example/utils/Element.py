class Element:
    def __init__(self, name, masses, probabilities):
        self.name = name
        self.masses = masses
        self.probabilities = probabilities


def try_parse_float(x):
    try:
        return float(x)
    except ValueError:
        return -1


def create_elements_from_file(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    elements = {}

    for line in lines:
        parts = line.strip().split('=')
        values = parts[1].split('|')
        name = values[0]
        masses = tuple(try_parse_float(x) for x in values[1].split(','))
        probabilities = tuple(try_parse_float(x) for x in values[2].split(','))
        elements[name] = (masses, probabilities)

    # elements = []
    #
    # for line in lines:
    #     parts = line.strip().split('=')
    #
    #     values = parts[1].split('|')
    #     name = values[0]
    #     masses = tuple(try_parse_float(x) for x in values[1].split(','))
    #     probabilities = tuple(try_parse_float(x) for x in values[2].split(','))
    #     element = Element(
    #         name=name,
    #         masses=masses,
    #         probabilities=probabilities
    #     )
    #     elements.append(element)

    return elements


def cal_mass_by_formula(elems, formula):
    mass = 0
    i = 0
    while i < len(formula):
        if formula[i].isalpha():
            j = i + 1
            while j < len(formula) and formula[j].isdigit():
                j += 1
            ele_name = formula[i:j]
            print(ele_name)
            k = j + 1
            while k < len(formula) and formula[k].isdigit():
                k += 1
            a = formula[j + 1:k]
            count = float(a)
            if j - 1 - i == 0:
                mass += elems[ele_name][0][0] * count
            else:
                mass += max(elems[ele_name[0]][0]) * count
            i = k
        else:
            i += 1

    return mass


if __name__ == '__main__':
    # 创建 element 对象列表
    elements = create_elements_from_file('data/element.ini')

    # 打印元素名称为 "B" 的同位素的不同质量和它们在自然界中的概率
    print(elements['B'])
    print(cal_mass_by_formula(elements, "C13(1)H(2)"))
    # # 打印 element 对象列表中第一个对象的名称和相对丰度
    # print(elements[0].name, elements[0].probabilities)
