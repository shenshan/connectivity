import datajoint as dj
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
schema = dj.schema('shan_literature', locals())


@schema
class ExcelFile(dj.Lookup):
    definition = """
    excel_file_id       : tinyint
    ---
    filename             : varchar(128)
    sheetname            : varchar(128)
    """
    contents = [
        (0, "Connectivity studies in the last 3 decades.xlsx", "Sheet5")
    ]


@schema
class Study(dj.Imported):
    definition = """
    # studies considered
    -> ExcelFile
    citation_key    : char(64)
    ---
    row_idx         : int           # row index in the loaded pandas array
    type            : enum("review", "original")
    year            : decimal(4,0)
    """

    def _make_tuples(self, key):
        filename, sheetname = (ExcelFile() & key).fetch1['filename', 'sheetname']

        df = pd.read_excel(filename, sheetname=sheetname)

        regexp = re.compile("""^(?P<authors>.+)(?P<year>[0-9]{4})(\s(?P<review>\(review\)))?""")
        for i, row in df.iterrows():
            k = row['paper id']

            pat = re.match(regexp, k)

            if pat is None:
                print(k, ' could not be parsed')
            else:
                gr = pat.groupdict()
                key.update(dict(citation_key="{authors}{year}".format(**gr), year=int(gr['year']),
                           type='original' if gr['review'] is None else 'review',
                           row_idx=i))
                self.insert1(key)



@schema
class AgeRange(dj.Computed):
    definition = """
    -> Study
    range_idx     : tinyint
    ---
    low=null          : float
    high=null         : float
    unit_low=null     : enum("P","E","g","yrs",'kg')
    unit_high=null    : enum("P","E","g","yrs",'kg')
    unit_type         : enum("mass","time")
    age_category=null : enum("juvenile","adult","unspecified","both","prenatal")
    range_type              : enum("range","singular")
    """

    def _make_tuples(self, key):
        filename, sheetname = (ExcelFile() & key).fetch1['filename', 'sheetname']
        df = pd.read_excel(filename, sheetname=sheetname)
        row = df.iloc[(Study() & key).fetch1['row_idx']]
        if not pd.isnull(row.age):

            # try postnatal range
            for range_idx, r in enumerate(row.age.split(',')):
                key['range_idx'] = range_idx
                range_pat = re.compile("""\s*((?P<unit_low>[P,E]?)(?P<low>([0-9]+\.?[0-9]*))\s*(?P<si_low>(?:kg|yrs|g))?)?\s*(?P<range>-)?\s*((?P<unit_high>[P,E]?)(?P<high>([0-9]+\.?[0-9]*))\s*(?P<si_high>(?:kg|yrs|g))?)?\s*""")
                match = re.match(range_pat, r).groupdict()
                # print(range_idx, r, match)
                match = {k:v if k not in ('high','low') else float(v) for k,v in match.items() if v is not None}
                key['age_category'] = row['juvenile or adult?'].lower() if not pd.isnull(row['juvenile or adult?']) else 'unspecified'


                if not match:
                    print('Could not parse', row.age)
                else:
                    if 'range' in match:
                        key['range_type'] = 'range'
                        match.pop('range')
                    else:
                        key['range_type'] = 'singular'

                    if 'si_high' not in match and 'si_low' not in match:
                        key.update(match)

                        key['unit_type'] = 'time'
                    else:
                        if not 'si_low' in match:
                            match['si_low'] = match['si_high']
                        if not 'si_high' in match:
                            match['si_high'] = match['si_low']
                        if match['si_high'] in ('g','kg'):
                            key['unit_type'] = 'mass'
                        else:
                            key['unit_type'] = 'time'

                        if 'high' in match:
                            key['high'] = match['high']
                            key['unit_high'] = match['si_high']
                        if 'low' in match:
                            key['low'] = float(match['low'])
                            key['unit_low'] = match['si_low']
                    self.insert1(key)

@schema
class Species(dj.Imported):
    definition = """
    -> Study
    species     : enum("monkey","mouse","rat","cat")
    """

    def _make_tuples(self, key):
        filename, sheetname = (ExcelFile() & key).fetch1['filename', 'sheetname']
        df = pd.read_excel(filename, sheetname=sheetname)
        row = df.iloc[(Study() & key).fetch1['row_idx']]
        if not pd.isnull(row.species):
            row.species = row.species.replace('and', ',')
            for sp in row.species.split(','):
                key['species'] = sp.strip().lower()
                self.insert1(key)




@schema
class CellsPatched(dj.Imported):
    definition = """
    -> Species
    ---
    no_cells             : int       # number of cells patched
    guess=0              : boolean   # if it is guesswork or not
    upper_bound=0        : boolean   # if it is an upper bound or not
    """

    @property
    def populated_from(self):
        return Study()

    def _make_tuples(self, key):
        filename, sheetname = (ExcelFile() & key).fetch1['filename', 'sheetname']
        df = pd.read_excel(filename, sheetname=sheetname)
        row = df.iloc[(Study() & key).fetch1['row_idx']]
        cells = row['number of cells patched']
        pat =  re.compile("""\s*(?P<guess>~)?(?P<upper_bound>\<|\<\=)?(?P<no_cells>[0-9]+)\s*(?P<species>[a-zA-Z]+)?\s*""")
        if not pd.isnull(cells):
            if isinstance(cells, int):
                key['species'] = (Species() & key).fetch1['species']
                key['no_cells'] = cells
                key['guess'] = False
                key['upper_bound'] = False
                self.insert1(key)
            else:
                for no_cells in cells.split(','):
                    gr = {k:int(v) if k == 'no_cells' else v for k,v, in pat.match(no_cells).groupdict().items() if v is not None}
                    if 'guess' in gr:
                        gr['guess'] = True
                    if 'upper_bound' in gr:
                        gr['upper_bound'] = True
                    if 'species' not in gr:
                        gr['species'] = (Species() & key).fetch1['species']
                    key.update(gr)
                    self.insert1(key)


    def plot(self, vs, exclude=None):
        cm = plt.cm.get_cmap('viridis')
        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        if exclude is None:
            exclude = []
        # df = pd.DataFrame(self.fetch.order_by('no_cells')())
        n = (self - (vs + exclude)).fetch['no_cells']
        xlabel= ['field']
        cm = sns.color_palette("Set2", len(n))
        for j, nn in enumerate(np.cumsum(n)):
            ax.bar(0, nn, align='center',lw=0, zorder=-nn, color=cm[j])

        n = ((self & "species='mouse'") - (vs + exclude)).fetch['no_cells']
        n = n[np.argsort(-n)]
        cm = sns.color_palette("Set2", len(n))
        xlabel.append('mouse field')
        for j, nn in enumerate(np.cumsum(n)):
            ax.bar(1, nn, align='center',lw=0, zorder=-nn, color=cm[j])

        for i, key in enumerate(vs):
            print(key)
            nx = ((self - exclude) & key).fetch['no_cells']
            ax.bar(i+2, nx,color='steelblue', align='center',lw=0)
            xlabel.append(key['citation_key'])
        ax.set_xticks(np.arange(len(vs)+2))
        ax.set_xticklabels(xlabel)
        ax.xaxis.grid(False)
        ax.set_ylabel('Number of cells patched')
        print('Approx studies:', len(self))
        plt.show()



