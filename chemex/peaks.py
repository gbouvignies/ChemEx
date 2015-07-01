import functools
import re

re_peak_name = re.compile(
    '''
        (^\s*|\-)
        (                                # group name
            (?P<symbol>\D?)              # one letter amino acid (optional)
            0*(?P<number>\d+)            # residue number
            (?P<suffix>[abd-gi-mopr-z]*) # suffix (optional)
        )?
        (?P<nucleus>                     # nucleus name (for example: CA, HG, ...)
            (?P<atom>[hncq])             # nucleus type
            [a-z0-9]*                    # nucleus name - nucleus type
        )
        (\s*$)?
    ''',
    re.IGNORECASE | re.VERBOSE
)


@functools.total_ordering
class Peak(object):
    def __init__(self, assignment):
        self.resonances = get_resonances(assignment.upper())
        self.assignment = get_assignment(self.resonances)

    def __repr__(self):
        return self.assignment

    def __eq__(self, other):
        return self.assignment == other.assignment

    def __lt__(self, other):
        if isinstance(other, basestring):
            other = Peak(other)

        self_tuple = tuple((resonance['nucleus'], int(resonance['number'])) for resonance in self.resonances)
        other_tuple = tuple((resonance['nucleus'], int(resonance['number'])) for resonance in other.resonances)

        return self_tuple < other_tuple


def get_resonances(assignment):
    resonances = []
    last_resonance = re_peak_name.match('H').groupdict()

    for match in re.finditer(re_peak_name, assignment.lower()):

        resonance = match.groupdict()

        for key, value in resonance.items():
            if value is None:
                resonance[key] = last_resonance[key]

        resonance['group'] = ''.join(resonance[_] for _ in ('symbol', 'number', 'suffix') if resonance[_] is not None)
        resonance['name'] = ''.join(resonance[_] for _ in ('group', 'nucleus'))

        resonances.append(resonance)
        last_resonance = resonance

    return resonances


def get_assignment(resonances):
    parts = []
    last_group = None

    for resonance in resonances:
        part = []
        if resonance['group'] is not None and resonance['group'] != last_group:
            part.append(resonance['group'])
        part.append(resonance['nucleus'])
        parts.append(''.join(part))
        last_group = resonance['group']

    return '-'.join(parts)
