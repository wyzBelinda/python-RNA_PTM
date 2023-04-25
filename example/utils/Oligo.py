from example.utils.Element import cal_mass_by_formula


class Oligo:
    def __init__(self, sequence, miss, chemistry_3, chemistry_5, sequence_location, mod="", mass=0,
                 masses: list = None):
        if masses is None:
            masses = []
        self.sequence = sequence
        self.mod = mod
        self.miss = miss
        self.chemistry_3 = chemistry_3
        self.chemistry_5 = chemistry_5
        self.sequence_location = sequence_location
        self.mass = mass
        self.masses = masses

    @staticmethod
    def from_string(line: str) -> 'Oligo':
        data = Oligo()
        fields = line.strip().split()
        data.sequence = fields[0]
        data.mod = fields[1]
        data.miss = int(fields[2])
        data.chemistry_3 = fields[3]
        data.chemistry_5 = fields[4]
        data.num_copy = int(fields[5])
        data.sequence_location = fields[6]
        data.molecule = fields[7]
        data.residue_start = int(fields[8])
        data.residue_end = int(fields[9])
        return data

    def to_string(self) -> str:
        return f"{self.sequence} {self.mod} {self.miss} {self.chemistry_3} {self.chemistry_5} {self.sequence_location} {self.mass} {self.masses}"
        # return f"{self.sequence} {self.mod} {self.miss} {self.chemistry_3} {self.chemistry_5} {self.num_copy} {self.sequence_location} {self.molecule} {self.residue_start} {self.residue_end}"

    def cal_oligo_mass(self, elements, nucls):
        from example.utils.Nucleotide import cal_mass_by_sequence
        mass = cal_mass_by_sequence(nucls, self.sequence)
        mass5 = cal_mass_by_formula(elements, self.chemistry_5)
        mass3 = cal_mass_by_formula(elements, self.chemistry_3)
        # self.mass = mass5 + mass + mass3
        return mass5, mass, mass3

    def cal_fragments_mass(self, elements: dict, nucleotides, mods, max_charge=5):
        # (elements["P"] + elements["O"] * 4 + elements["H"][0][0] * 2)
        # TODO: 中性碎片片段转多电荷版未完成，考虑mods修饰版本未完成（只将c++原始版本与肽段版本分离）
        mass5, mass, mass3 = self.cal_oligo_mass(elements, nucleotides)
        print(mass5, mass, mass3)

        # cur_mass = mass5
        # # for info in self.sequence_location['mod_infos']:
        #
        # mods_index = [""] * len(self.sequence)
        #
        # for info in self.sequence_location[0]['mod_infos']:  # 多个母离子同一个self的时候，mods位置一样
        #     mods_index[info['pos']] = info['mod_symbol']
        # print(mods_index)
        # # mods_poses = set(zip(self.sequence_location[0]['mod_infos']['pos']))
        # mods_poses = set(d['pos'] for d in self.sequence_location[0]['mod_infos'])
        #
        # ion_masses = []
        # for i in range(len(self.sequence)):
        #     cur_ion_masses = []
        #     icon_mass = [nucleotides[self.sequence[i]].base_mass, nucleotides[self.sequence[i]].sugar_mass,
        #                  nucleotides[self.sequence[i]].phosphate_mass]
        #     if i in mods_poses:
        #         icon_mass += [mods[mods_index[i]].base_delta_mass, mods[mods_index[i]].sugar_delta_mass,
        #                       mods[mods_index[i]].phosphate_delta_mass]
        #
        #     cur_mass += icon_mass[0]
        #     cur_ion_masses.append(cur_mass)  # a-B
        #     cur_mass += icon_mass[1]
        #     cur_ion_masses.append(cur_mass)  # a
        #     cur_mass += elements['O']
        #     cur_ion_masses.append(cur_mass)  # b
        #     cur_mass += cal_mass_by_formula(elements, "P(1)O(2)H(1)")
        #     cur_ion_masses.append(cur_mass)  # c
        #     cur_mass += elements['O']
        #     cur_ion_masses.append(cur_mass)  # d

        # if self.mod == "":
        #     # 无需考虑修饰

        #     -----------------

        masses = []
        cur_mass_5 = mass5 + elements["e"][0][0]
        label_id = 1
        seq_len = len(self.sequence)
        whole_mass = mass + mass5 + mass3
        self.mass = whole_mass

        for B in self.sequence:
            if label_id == seq_len:
                break

            cur_mass_5 += nucleotides[B].sugar_mass - elements["H"][0][0] * 3
            flag = ord('a')

            for i in range(1, max_charge + 1):
                label1 = "a_{}-{}^{}".format(label_id, B, i)
                tmp_mass = (cur_mass_5 - elements["A"][0][0] * (i - 1)) / i
                masses.append((label1, tmp_mass))

            cur_mass_5 += nucleotides[B].base_mass + elements["H"][0][0]

            for i in range(1, max_charge + 1):
                label1 = "{}_{}^{}".format(chr(flag), label_id, i)
                masses.append((label1, (cur_mass_5 - elements["A"][0][0] * (i - 1)) / i))

            cur_mass_3 = whole_mass - cur_mass_5 - elements["A"][0][0] * 2 - elements["e"][0][0]

            for i in range(1, max_charge + 1):
                label2 = "{}_{}^{}".format(chr(flag + 22), seq_len - label_id, i)
                masses.append((label2, (cur_mass_3 - elements["A"][0][0] * (i - 1)) / i))

            cur_mass_5 += cal_mass_by_formula(elements, "H(2)O(1)")

            for i in range(1, max_charge + 1):
                label1 = "{}_{}^{}".format(chr(flag + 1), label_id, i)
                masses.append((label1, (cur_mass_5 - elements["A"][0][0] * (i - 1)) / i))

            cur_mass_3 = whole_mass - cur_mass_5 - elements["A"][0][0] * 2 - elements["e"][0][0]

            for i in range(1, max_charge + 1):
                label2 = "{}_{}^{}".format(chr(flag + 22 + 1), seq_len - label_id, i)
                masses.append((label2, (cur_mass_3 - elements["A"][0][0] * (i - 1)) / i))

            cur_mass_5 += cal_mass_by_formula(elements, "H(-1)O(2)P(1)")

            for i in range(1, max_charge + 1):
                label1 = "{}_{}^{}".format(chr(flag + 2), label_id, i)
                masses.append((label1, (cur_mass_5 - elements["A"][0][0] * (i - 1)) / i))

            cur_mass_3 = whole_mass - cur_mass_5 - elements["A"][0][0] * 2 - elements["e"][0][0]

            for i in range(1, max_charge + 1):
                label2 = "{}_{}^{}".format(chr(flag + 22 + 2), seq_len - label_id, i)
                masses.append((label2, (cur_mass_3 - elements["A"][0][0] * (i - 1)) / i))

        masses.sort(key=lambda x: x[1])

        self.masses = masses
        print(self.chemistry_5, self.sequence, self.chemistry_3)
        print(masses)

        return masses
