import re

nuclei_dict = {
    'H': '1H',
    'N': '15N',
    'C': '13C',
}


class Peak(object):
    def __init__(self, assignment):
        self.assignment = assignment
        self.resonances = [
            Resonance(name) for name in split_assignment(self.assignment)
            ]


class Resonance(object):
    def __init__(self, name):
        group, atom = split_resonance_name(name)

        self.name = ''.join([group, atom])
        self.group = Group(group)
        self.atom = Atom(atom)


class Group(object):
    def __init__(self, name):
        symbol, number, suffix = split_group_name(name)

        self.name = name
        self.symbol = symbol
        self.number = int(number)
        self.suffix = suffix


class Atom(object):
    def __init__(self, name):
        self.name = name
        self.nucleus = nuclei_dict[name[0].upper()]


def split_group_name(name):
    symbol, number, suffix = re.split('([0-9]+)', name, 1)
    return symbol, number, suffix


def split_resonance_name(name):
    group_atom = None

    if name == '?':
        group_atom = '', ''
    else:
        split_index = 0
        match_digit = re.search('[0-9]', name)
        if match_digit:
            split_index += match_digit.start()
            match_atom = re.search('[hHcCnNqQmM]', name[split_index:])
            if match_atom:
                split_index += match_atom.start()
                group_atom = name[:split_index], name[split_index:].upper()

    return group_atom


def split_assignment(name):
    resonances = []
    resonance_names = name.split('-')
    last_group = None
    for resonance_name in resonance_names:
        group_atom = split_resonance_name(resonance_name)
        if group_atom is None:
            if last_group is not None:
                group_atom = last_group, resonance_name
        resonances.append(''.join(group_atom))
        last_group = group_atom[0]
    return resonances


# If an assignment along an axis has the same group as the previous axis
# then the group name is not shown.

def format_assignment(resonances):
    assignment = []
    last_group = None
    for resonance in resonances:
        if resonance:
            if resonance.group == last_group:
                assignment.append(resonance.atom.name)
            else:
                assignment.append(resonance.name)
            last_group = resonance.group
        else:
            assignment.append('?')
    return '-'.join(assignment)
